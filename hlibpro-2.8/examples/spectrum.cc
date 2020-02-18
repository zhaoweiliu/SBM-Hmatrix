//
// Project     : HLib
// File        : spectrum.cc
// Description : computes spektrum of graph Laplacian
// Author      : Ronald Kriemann
// Copyright   : Max Planck Institute MIS 2004-2014. All Rights Reserved.
//

#include <iostream>

#include <boost/format.hpp>

#include <tbb/parallel_for.h>

#include "hlib.hh"

using namespace std;
using namespace HLIB;

using boost::format;

double
calc_abs_det ( const TMatrix *  L,
               const TMatrix *  Id,
               const double     lambda_0,
               const double     mu,
               const double     eps )
{
    auto           progress = unique_ptr< TProgressBar >( verbose(3) ? new TConsoleProgressBar : nullptr );
    auto           L_lmu_Id = L->copy();
    fac_options_t  opts( progress.get() );

    // ensure "real" LU instead of block-wise LU (with inversion at diagonal leaves)
    opts.eval = point_wise;

    add( HLIB::complex( -lambda_0, -mu ), Id, HLIB::complex( 1, 0 ), L_lmu_Id.get(), acc_exact );
    factorise( L_lmu_Id.get(), fixed_prec( eps ), opts );
    
    auto    diag = diagonal( L_lmu_Id.get() );
    double  sum  = 0.0;

    for ( auto  i : diag->is() )
    {
        const auto  u_ii     = diag->centry( i );
        const auto  norm_uii = std::sqrt( norm( u_ii ) );

        if ( norm_uii < 1e-12 )
        {
            std::cout << "detected near zero entry" << std::endl;
            return 0;
        }// if
        else
            sum += std::log( norm_uii );
    }// for

    return sum;
}

std::pair< BLAS::Vector< double >,
           BLAS::Vector< double > >
calc_y ( const TMatrix *  L,
         const TMatrix *  Id,
         const double     mu,
         const double     h,
         const uint       nsteps,
         const double     eps )
{
    auto                    progress = unique_ptr< TProgressBar >( verbose(1) ? new TConsoleProgressBar : nullptr );
    BLAS::Vector< double >  x( nsteps+1 );
    BLAS::Vector< double >  y( nsteps+1 );

    if ( verbose(1) )
        progress->init( 0, nsteps+1, 0 );
    
    // spectrum is in [0,2]
    const double  stepwidth = 2.0 / double(nsteps);

    tbb::parallel_for( uint(0), nsteps+1,
                       [=,&x,&y,&progress] ( const uint  i )
                       {
                           const auto  lambda_0     = i * stepwidth;
                           const auto  f_mu_plus_h  = calc_abs_det( L, Id, lambda_0, mu+h, eps );
                           const auto  f_mu_minus_h = calc_abs_det( L, Id, lambda_0, mu-h, eps );

                           x(i) = lambda_0;
                           y(i) = (f_mu_plus_h - f_mu_minus_h) / (2 * h);
                           
                           if ( verbose( 2 ) )
                               std::cout << "λ₀=" << format( "%.2e" ) % lambda_0
                                         << " : " << format( "%.4e" ) % y(i) << std::endl;

                           if ( verbose( 1 ) )
                               progress->advance( 1 );
                       } );

    return std::pair< BLAS::Vector< double >,
                      BLAS::Vector< double > >( std::move( x ),
                                                std::move( y ) );
}

int
main ( int      argc,
       char **  argv )
{
    std::string  graph_file = "graph.mat";
    double       mu         = 1e-2;
    double       h          = 1e-3;
    uint         nsteps     = 200;
    double       eps        = 1e-4;
    int          verbosity  = 1;
        
    if ( argc > 1 )
        graph_file = argv[1];
    
    try
    {
        HLIB::INIT();

        HLIB::CFG::set_verbosity( verbosity );
        
        auto  progress = unique_ptr< TProgressBar >( verbose(2) ? new TConsoleProgressBar : nullptr );
        auto  M = read_matrix( graph_file.c_str() );

        if ( ! IS_TYPE( M, TSparseMatrix ) )
        {
            std::cout << "input matrix is not sparse" << std::endl;
            return 1;
        }// if

        auto  S = ptrcast( M.get(), TSparseMatrix );

        S->set_unsymmetric();
        
        TBFSAlgPartStrat      part_strat;
        TAlgCTBuilder         base_ct_builder( & part_strat, 100 );
        TAlgNDCTBuilder       nd_ct_builder( & base_ct_builder );
        auto                  ct = nd_ct_builder.build( S );

        TWeakAlgAdmCond       adm_cond( S, ct->perm_i2e() );
        TBCBuilder            bct_builder;
        auto                  bct = bct_builder.build( ct.get(), ct.get(), & adm_cond );
        
        TSparseMBuilder       h_builder( S, ct->perm_i2e(), ct->perm_e2i() );
        auto                  L  = h_builder.build( bct->root(), acc_exact, progress.get() );
        auto                  Id = identity( bct.get() ); 

        // ensure both are complex valued
        L->set_complex( true );
        Id->set_complex( true );

        std::cout << "computing spectrum for " << graph_file << " with" << std::endl
                  << "  μ = " << format( "%.2e" ) % mu << ", "
                  << "h = " << format( "%.2e" ) % h << " and "
                  << nsteps << " steps" << std::endl;

        auto  starttime = Time::Wall::now();
        auto  spectrum  = calc_y( L.get(), Id.get(), mu, h, nsteps, eps );
        auto  runtime   = Time::Wall::since( starttime );

        std::cout << "  done in " << runtime << std::endl;

        TMatlabVectorIO  vio;

        vio.write( spectrum.first,  "x.mat", "x" );
        vio.write( spectrum.second, "y.mat", "y" );

        HLIB::DONE();
    }// try
    catch ( std::exception & e )
    {
        std::cout << e.what() << std::endl;
    }// catch
}
