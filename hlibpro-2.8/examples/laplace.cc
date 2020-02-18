//
// Project     : HLib
// File        : sparsealg.cc
// Description : example for sparse matrices with algebraic clustering
// Author      : Ronald Kriemann
// Copyright   : Max Planck Institute MIS 2004-2014. All Rights Reserved.
//

#include <iostream>

#include <boost/format.hpp>

#include "hlib.hh"

using namespace std;
using boost::format;
using namespace HLIB;
using HLIB::uint;

using real_t = HLIB::real;

namespace
{

//
// function for right-hand side
//
class TMyRHS : public TBEMFunction< real_t >
{
public:
    // actual RHS function
    virtual real_t
    eval  ( const T3Point &  pos,
            const T3Point & ) const
    {
        return sin( 0.5 * pos.x() * cos( 0.25 * pos.z() ) ) * cos( pos.y() );
    }
};

}// namespace anonymous

//
// main function
//
int
main ( int argc, char ** argv )
{
    //
    // read command line
    //
    
    string  gridfile  = "sphere-5";
    real_t  eps       = real_t(1e-4);

    if ( argc > 1 )
        gridfile = argv[1];
    
    try
    {
        //
        // init HLIBpro
        //
        
        INIT();

        CFG::set_verbosity( 3 );

        //
        // read grid
        //

        cout << "## constructing grid (" << gridfile << " )" << endl;
        
        auto  grid = make_grid( gridfile );
        
        TVTKGridVis  gvis;

        if ( verbose( 3 ) )
            gvis.print( grid.get(), "laplace_grid" );
        
        //
        // build ansatz and test space
        //

        cout << endl << "## building ansatz/test spaces" << endl;
        
        using  ansatzsp_t = TLinearFnSpace;
        using  testsp_t   = TConstFnSpace;
        using  dlpbf_t    = TLaplaceDLPBF< ansatzsp_t, testsp_t >;
        using  coeff_t    = TBFCoeffFn< dlpbf_t >;
        
        testsp_t    test_fnspace(   grid.get() );
        ansatzsp_t  ansatz_fnspace( grid.get() );
        
        cout << "    no. of ansatz dof = " << ansatz_fnspace.n_indices() << endl;
        cout << "    no. of test dof   = " << test_fnspace.n_indices() << endl;
        
        //
        // and the cluster trees
        //
        
        auto               coord1 = ansatz_fnspace.build_coord();
        auto               coord2 = test_fnspace.build_coord();

        TAutoBSPPartStrat  part_strat;
        TBSPCTBuilder      ct_builder( & part_strat );

        auto               ct1 = ct_builder.build( coord1.get() );
        auto               ct2 = ct_builder.build( coord2.get() );
        
        TStdGeomAdmCond    adm_cond;
        TBCBuilder         bct_builder;
        auto               bct = bct_builder.build( ct1.get(), ct2.get(), & adm_cond );

        if ( verbose( 3 ) )
        {
            TPSBlockClusterVis  bc_vis;

            bc_vis.print( bct->root(), "laplace_bct" );
        }// if

        //
        // build RHS
        //

        TMyRHS                           rhsfn;
        TQuadBEMRHS< testsp_t, real_t >  rhs_build( 4 );
        auto                             rhs = rhs_build.build( & test_fnspace, & rhsfn );

        if ( verbose( 3 ) )
            gvis.print( grid.get(), & test_fnspace, rhs.get(), "laplace_rhs" );

        // bring RHS into H-ordering
        ct2->perm_e2i()->permute( rhs.get() );
        
        //
        // build matrix
        //

        cout << endl << "## building Laplace DLP H-matrix K ( ε = " << format( "%.2e" ) % eps << " )" << endl;

        TConsoleProgressBar       progress;
        TWallTimer                timer;
        TTruncAcc                 acc( eps, 0.0 );
        dlpbf_t                   bf( & ansatz_fnspace, & test_fnspace );
        coeff_t                   bf_coeff_fn( & bf );
        TPermCoeffFn< real_t >    coeff_fn( & bf_coeff_fn, bct->row_ct()->perm_i2e(), bct->col_ct()->perm_i2e() );
        TACAPlus< real_t >        aca( & coeff_fn );
        TDenseMBuilder< real_t >  h_builder( & coeff_fn, & aca );
        TPSMatrixVis              mvis;

        // enable coarsening during construction
        // h_builder.set_coarsening( TTruncAcc( eps * 0.4, 0.0 ) );

        timer.start();

        auto  A = h_builder.build( bct.get(), acc, & progress );

        timer.pause();

        cout << "    done in " << timer << endl;
        cout << "    size of H-matrix = " << Mem::to_string( A->byte_size() ) << endl;
        cout << "    |A|₂             = " << format( "%.6e" ) % norm_2( A.get() ) << endl;
        cout << "    |A|_F            = " << format( "%.6e" ) % norm_F( A.get() ) << endl;

        if( verbose( 3 ) )
        {
            mvis.svd( true ).print( A.get(), "laplace_A" );
        }// if

        //
        // build mass matrix
        //
        
        using  massbf_t    = TMassBF< ansatzsp_t, testsp_t >;
        using  masscoeff_t = TBFCoeffFn< massbf_t >;
            
        massbf_t                  massbf( & ansatz_fnspace, & test_fnspace );
        masscoeff_t               bf_mcoeff_fn( & massbf );
        TPermCoeffFn< real_t >    mcoeff_fn( & bf_mcoeff_fn, bct->row_ct()->perm_i2e(), bct->col_ct()->perm_i2e() );
        TACAPlus< real_t >        maca( & mcoeff_fn );
        TDenseMBuilder< real_t >  mh_builder( & mcoeff_fn, & maca );

        cout << endl << "## building Mass H-matrix M " << endl;
        timer.start();

        auto  M = mh_builder.build( bct.get(), acc, & progress );

        timer.pause();
        cout << "    done in " << timer << endl;

        // test DLP
        auto  x = A->col_vector();
        auto  y = A->row_vector();

        x->fill( 1.0 );
        y->fill( 0.0 );

        M->mul_vec(  0.5, x.get(), 0.0, y.get(), apply_normal );
        A->mul_vec( -1.0, x.get(), 1.0, y.get(), apply_normal );

        cout << endl << "## approximation error " << endl
             << "    |(½M-K)·𝟏|/(|K|·|𝟏|) = "
             << format( "%.4e" ) % ( y->norm2() / ( norm_2( A.get() ) * x->norm2() ) )
             << endl;
            
        DONE();
    }// try
    catch ( Error & e )
    {
        cout << e.to_string() << endl;
    }// catch

    return 0;
}
