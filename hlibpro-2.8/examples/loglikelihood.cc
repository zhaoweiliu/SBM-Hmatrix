#include <iostream>
#include <fstream>
#include <string>
#include <array>

#include <gsl/gsl_multimin.h>

#include "hlib.hh"

using namespace HLIB;

enum {
    IDX_SIGMA  = 0,
    IDX_LENGTH = 1,
    IDX_NU     = 2
};
    
void
read_data ( const std::string &       datafile,
            std::vector< T2Point > &  vertices,
            BLAS::Vector< double > &  Z_data )
{
    std::ifstream  in( datafile );
    
    if ( ! in ) // error
        exit( 1 );

    size_t  N_vtx = 0;
    
    in >> N_vtx;

    std::cout << "reading " << N_vtx << " datapoints" << std::endl;
    
    vertices.resize( N_vtx );
    Z_data = BLAS::Vector< double >( N_vtx );
        
    for ( idx_t  i = 0; i < idx_t(N_vtx); ++i )
    {
        int     index = i;
        double  x, y;
        double  v     = 0.0;

        in >> index >> x >> y >> v;

        vertices[ index ] = T2Point( x, y );
        Z_data( index )   = v;
    }// for

    //
    // visualization data with grid
    //

    using  triangle_t = std::array< idx_t, 3 >;
    
    size_t                     N_tri;
    std::vector< triangle_t >  triangles;

    in >> N_tri;

    triangles.resize( N_tri );
    
    for ( idx_t  i = 0; i < idx_t(N_tri); ++i )
    {
        int     index = i;
        int     v0, v1, v2;

        in >> index >> v0 >> v1 >> v2;

        triangles[ index ] = triangle_t{ v0, v1, v2 };
    }// for

    std::ofstream  out( "dataset.vtk" );

    out << "# vtk DataFile Version 2.0" << std::endl
        << "HLIBpro grid" << std::endl
        << "ASCII" << std::endl
        << "DATASET UNSTRUCTURED_GRID" << std::endl;
    
    out << "POINTS " << N_vtx << " FLOAT" << std::endl;
                                             
    for ( uint  i = 0; i < N_vtx; ++i )
        out << vertices[i].x() << " " << vertices[i].y() << " 0" << std::endl;

    out << "CELLS " << N_tri << " " << 4 * N_tri << std::endl;

    for ( uint i = 0; i < N_tri; i++ )
        out << "3 " << triangles[i][0] << " " << triangles[i][1] << " " << triangles[i][2] << std::endl;

    out << "CELL_TYPES " << N_tri << std::endl;
        
    for ( uint  i = 0; i < N_tri; ++i )
        out << "5 ";
    out << std::endl;
    
    out << "POINT_DATA " << N_vtx << std::endl
        << "SCALARS Z float 1" << std::endl
        << "LOOKUP_TABLE default " << std::endl;

    for ( uint  i = 0; i < N_vtx; ++i )
        out << Z_data(i) << " ";
    out << std::endl;
}

struct LogLikeliHoodProblem
{
    std::vector< T2Point >                vertices;
    std::unique_ptr< TCoordinate >        coord;
    std::unique_ptr< TClusterTree >       ct;
    std::unique_ptr< TBlockClusterTree >  bct;
    std::unique_ptr< TVector >            Z;
    double                                eps;

    LogLikeliHoodProblem ( const std::string &  coord_filename,
                           const double         accuracy )
            : eps( accuracy )
    {
        init( coord_filename );
    }

    void
    init ( const std::string &  datafile )
    {
        BLAS::Vector< double >  Z_data;

        read_data( datafile, vertices, Z_data );

        coord = std::make_unique< TCoordinate >( vertices );

        TAutoBSPPartStrat  part_strat;
        TBSPCTBuilder      ct_builder( & part_strat );
    
        ct = ct_builder.build( coord.get() );
    
        TStdGeomAdmCond    adm_cond( 2.0, use_min_diam );
        TBCBuilder         bct_builder;
    
        bct = bct_builder.build( ct.get(), ct.get(), & adm_cond );

        Z   = std::make_unique< TScalarVector >( *ct->root(), std::move( Z_data ) );

        ct->perm_e2i()->permute( Z.get() );
    }
    
    double
    eval ( const double  sigma,
           const double  length,
           const double  nu )
    {
        TMaternCovCoeffFn< T2Point >  matern_coefffn( sigma, length, nu, vertices );
        TPermCoeffFn< double >        coefffn( & matern_coefffn, ct->perm_i2e(), ct->perm_i2e() );

        TACAPlus< double >            aca( & coefffn );
        auto                          acc = fixed_prec( eps );
        TDenseMatBuilder< double >    h_builder( & coefffn, & aca );
    
        auto                          C     = h_builder.build( bct.get(), acc );
        auto                          C_fac = C->copy();

        chol( C_fac.get(), acc );
    
        auto                          C_inv = std::make_unique< TLLInvMatrix >( C_fac.get(), symmetric );
        const size_t                  N     = vertices.size();
        double                        log_det_C = 0.0;
    
        for ( idx_t  i = 0; i < idx_t(N); ++i )
            log_det_C += 2*std::log( C_fac->entry( i, i ) ); // two factors L!
        
        TStopCriterion                sstop( 250, 1e-16, 0.0 );
        TCG                           solver( sstop );
        auto                          sol = C->row_vector();
    
        solver.solve( C.get(), sol.get(), Z.get(), C_inv.get() );
        
        auto                          ZdotCZ = re( Z->dot( sol.get() ) );
        const double                  log2pi = std::log( 2.0 * Math::pi<double>() );
        auto                          LL     = -0.5 * ( N * log2pi + log_det_C + ZdotCZ );

        return LL;
    }
};

double
eval_logli ( const gsl_vector *  param,
             void *              data )
{
    double sigma  = gsl_vector_get( param, IDX_SIGMA );
    double length = gsl_vector_get( param, IDX_LENGTH );
    double nu     = gsl_vector_get( param, IDX_NU );

    LogLikeliHoodProblem *  problem = static_cast< LogLikeliHoodProblem * >( data );

    return - problem->eval( sigma, length, nu );
}

double
maximize_likelihood ( double &                sigma,
                      double &                length,
                      double &                nu,
                      LogLikeliHoodProblem &  problem )
{
    int        status   = 0;
    const int  max_iter = 200;

    const gsl_multimin_fminimizer_type *  T = gsl_multimin_fminimizer_nmsimplex2;
    gsl_multimin_fminimizer *             s = NULL;
    gsl_vector *                          ss;
    gsl_vector *                          x;
    gsl_multimin_function                 minex_func;
    double                                size;

    x = gsl_vector_alloc( 3 );        // start value
    gsl_vector_set( x, IDX_SIGMA,  sigma );
    gsl_vector_set( x, IDX_LENGTH, length );
    gsl_vector_set( x, IDX_NU,     nu );

    ss = gsl_vector_alloc( 3 );       // step sizes
    gsl_vector_set( ss, IDX_SIGMA,  0.001 );
    gsl_vector_set( ss, IDX_LENGTH, 0.001 );
    gsl_vector_set( ss, IDX_NU,     0.001 );

    // Initialize method and iterate
    minex_func.n      = 3;
    minex_func.f      = & eval_logli;
    minex_func.params = & problem;

    s = gsl_multimin_fminimizer_alloc( T, 3 );
    gsl_multimin_fminimizer_set(s, & minex_func, x, ss );

    int     iter = 0;
    double  LL = 0;

    do
    {
        iter++;
        status = gsl_multimin_fminimizer_iterate( s );

        if ( status != 0 )
            break;

        size   = gsl_multimin_fminimizer_size( s );    // return eps for stopping criteria
        status = gsl_multimin_test_size( size, 1e-3 ); // This function tests the minimizer specific characteristic size 

        if ( status == GSL_SUCCESS )
            std::cout << "converged to minimum" << std::endl;

        sigma  = gsl_vector_get( s->x, IDX_SIGMA );
        length = gsl_vector_get( s->x, IDX_LENGTH );
        nu     = gsl_vector_get( s->x, IDX_NU );
        LL     = -s->fval;

        std::cout << "  logliklihood at"
                  << "  σ = " << sigma
                  << ", ℓ = " << length
                  << ", ν = " << nu
                  << "   is " << LL
                  << std::endl;
    } while (( status == GSL_CONTINUE ) && ( iter < max_iter ));

    gsl_multimin_fminimizer_free( s );
    gsl_vector_free( ss );
    gsl_vector_free( x );

    return LL;
}

int
main ( int      argc,
       char **  argv )
{
    INIT();

    double  sigma  = 1.0;
    double  length = 1.29; 
    double  nu     = 0.325;
    
    LogLikeliHoodProblem  problem( "datafile.txt", 1e-7 );

    auto  LL = maximize_likelihood( sigma, length, nu, problem );

    DONE();
}
