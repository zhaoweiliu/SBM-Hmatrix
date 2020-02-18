//
// Project     : HLib
// File        : sparsealg.cc
// Description : example for sparse matrices with algebraic clustering
// Author      : Ronald Kriemann
// Copyright   : Max Planck Institute MIS 2004-2014. All Rights Reserved.
//

#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>

#include <boost/format.hpp>

#include "hlib.hh"

using namespace std;
using boost::format;
using namespace HLIB;
using HLIB::uint;
using HLIB::real;

using  real_t    = HLIB::real;
using  complex_t = HLIB::complex;

namespace
{

//
// function for right-hand side
//
class TMyRHS : public TBEMFunction< complex_t >
{
public:
    virtual complex_t
    eval ( const T3Point &  pos,
           const T3Point &  ) const
    {
        return std::sin( pos.x() ) * std::cos( pos.y() );
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
    
    string     gridfile  = "sphere-5";
    real_t     eps       = real_t(1e-4);
    complex_t  kappa     = complex_t( 2, 1 );
    
    if ( argc > 1 )
        gridfile = argv[1];

    try
    {
        //
        // init HLIBpro
        //
        
        INIT();

        CFG::set_verbosity( 1 );
        
        //
        // read grid
        //

        cout << "## constructing grid (" << gridfile << " )" << endl;
        
        auto  grid = make_grid( gridfile );
        
        TVTKGridVis  gvis;
            
        if ( verbose( 3 ) )
            gvis.print( grid.get(), "helmholtz_grid" );
        
        //
        // build ansatz and test space
        //

        cout << endl << "## building ansatz/test spaces" << endl;
        
        using  ansatzsp_t = TLinearFnSpace;
        using  testsp_t   = TLinearFnSpace;
        using  slpbf_t    = THelmholtzSLPBF< ansatzsp_t, testsp_t >;
        using  coeff_t    = TBFCoeffFn< slpbf_t >;
        
        ansatzsp_t   ansatz_fnspace( grid.get() );
        testsp_t     test_fnspace(   grid.get() );

        cout << "    no. of ansatz dof = " << ansatz_fnspace.n_indices() << endl;
        cout << "    no. of test dof   = " << test_fnspace.n_indices() << endl;

        //
        // and the cluster trees
        //

        auto                  coord1 = ansatz_fnspace.build_coord();
        auto                  coord2 = test_fnspace.build_coord();

        TAutoBSPPartStrat     part_strat;
        TBSPCTBuilder         ct_builder( & part_strat );

        auto                  ct1 = ct_builder.build( coord1.get() );
        auto                  ct2 = ct_builder.build( coord2.get() );
        
        THiLoFreqGeomAdmCond  adm_cond( kappa, 10 );
        TBCBuilder            bct_builder;
        auto                  bct = bct_builder.build( ct1.get(), ct2.get(), & adm_cond );

        if ( verbose( 3 ) )
        {
            TPSBlockClusterVis  bc_vis;

            bc_vis.print( bct->root(), "helmholtz_bct" );
        }// if
        
        //
        // build RHS
        //

        TMyRHS                              rhsfn;
        TQuadBEMRHS< testsp_t, complex_t >  rhs_build( 4 );
        auto                                rhs = rhs_build.build( & test_fnspace, & rhsfn );

        if ( verbose( 3 ) )
        {
            auto  rhsre = rhs->restrict_re();
            
            gvis.print( grid.get(), & test_fnspace, rhsre.get(), "helmholtz_rhs" );
        }// if
        
        // bring RHS into H-ordering
        ct2->perm_e2i()->permute( rhs.get() );
        
        //
        // build matrix
        //
        
        cout << endl << "## building Helmholtz SLP H-matrix ("
             << " κ = " << kappa
             << ", ε = " << format( "%.2e" ) % eps
             << " )" << endl;

        TConsoleProgressBar          progress;
        TTimer                       timer( WALL_TIME );
        TTruncAcc                    acc( eps, 0.0 );
        slpbf_t                      bf( kappa, & ansatz_fnspace, & test_fnspace );
        coeff_t                      bf_coeff_fn( & bf );
        TPermCoeffFn< complex_t >    coeff_fn( & bf_coeff_fn, bct->row_ct()->perm_i2e(), bct->col_ct()->perm_i2e() );
        TACAPlus< complex_t >        aca( & coeff_fn );
        TDenseMBuilder< complex_t >  h_builder( & coeff_fn, & aca );
        TPSMatrixVis                 mvis;
        
        // enable coarsening during construction
        // h_builder.set_coarsening( TTruncAcc( eps * 0.4, 0.0 ) );

        timer.start();

        auto  A = unique_ptr< TMatrix >( h_builder.build( bct.get(), acc, & progress ) );

        timer.pause();
        cout << "    done in " << timer << endl;
        cout << "    size of H-matrix = " << Mem::to_string( A->byte_size() ) << endl;
        cout << "    |A|₂             = " << format( "%.8e" ) % norm_2( A.get() ) << endl;
        cout << "    |A|_F            = " << format( "%.8e" ) % norm_F( A.get() ) << endl;

        if( verbose( 3 ) )
        {
            mvis.svd( true ).print( A.get(), "helmholtz_A" );
        }// if

        //
        // LU decomposition
        //

        auto  B = A->copy();

        cout << endl << "## LU factorisation ( ε = " << format( "%.2e" ) % eps << " )" << endl;

        timer.start();

        auto  A_inv = factorise_inv( B.get(), acc, & progress );
    
        timer.pause();
        cout << "    done in " << timer << endl;

        if( verbose( 3 ) )
            mvis.print( B.get(), "helmholtz_LU" );
    
        cout << "    size of LU factor = " << Mem::to_string( B->byte_size() ) << endl;
        cout << "    inversion error   = " << format( "%.6e" ) % inv_approx_2( A.get(), A_inv.get() ) << endl;
        cout << "    |LU|              = " << format( "%.6e" ) % norm_2( B.get() ) << endl;
        cout << "    condest(LU)       = " << format( "%.6e" ) % condest( A_inv.get() ) << endl;

        //
        // build RHS and solve
        //

        auto         x = A->col_vector();
        TSolverInfo  solve_info( false, verbose( 4 ) );

        cout << endl << "## solving with H-LU preconditioner" << endl;
        
        solve( A.get(), x.get(), rhs.get(), A_inv.get(), & solve_info );
        
        cout << "    " << solve_info.to_string() << endl;

        if ( verbose( 3 ) )
        {
            auto  solre = x->restrict_re();
            
            ct1->perm_i2e()->permute( solre.get() );
            gvis.print( grid.get(), & ansatz_fnspace, solre.get(), "helmholtz_sol_pre" );
        }// if
        
        {
            auto  b2 = rhs->copy();

            A->mul_vec( 1.0, x.get(), -1.0, b2.get() );
            cout << "    |Ax-b| = " << b2->norm2() << endl;
        }
    
        DONE();
    }// try
    catch ( Error & e )
    {
        cout << e.to_string() << endl;
    }// catch
    
    return 0;
}
