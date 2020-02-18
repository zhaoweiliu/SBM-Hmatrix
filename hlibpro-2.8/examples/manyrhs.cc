//
// Project     : HLib
// File        : manyrhs.cc
// Description : example for solving many RHS
// Author      : Ronald Kriemann
// Copyright   : Max Planck Institute MIS 2004-2014. All Rights Reserved.
//

#include <iostream>
#include <string>

#include <boost/format.hpp>

#include <tbb/parallel_for.h>

#include "hlib.hh"

using namespace std;
using boost::format;
using namespace HLIB;

using  real_t    = HLIB::real;
using  fnspace_t = TConstFnSpace;
using  slpbf_t   = TLaplaceSLPBF< fnspace_t, fnspace_t >;

namespace
{

//
// function for right-hand side
//
class TMyRHS : public TBEMFunction< real_t >
{
private:
    const double  _alpha, _beta;
    
public:
    TMyRHS ( const double  alpha,
             const double  beta )
            : _alpha( alpha )
            , _beta( beta )
    {}
    
    // actual RHS function
    virtual real_t
    eval  ( const T3Point &  pos,
            const T3Point & ) const
    {
        return real_t( sin( 0.5 * pos.x() * cos( _alpha * 0.25 * pos.z() ) ) * cos( _beta * pos.y() ) );
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
    
    string     gridfile   = "sphere-5";
    double     eps        = 1e-4;
    double     fac_eps    = 1e-6;
    matform_t  mat_form   = unsymmetric;
    
    
    if ( argc > 1 )
        gridfile = argv[1];

    try
    {
        INIT();

        TConsoleProgressBar  progress;
        TWallTimer           timer;

        //
        // read grid and build function spaces, coordinates and cluster trees
        //

        TAutoBSPPartStrat  part_strat;
        TBSPCTBuilder      ct_builder( & part_strat );
        TStdGeomAdmCond    adm_cond( 2.0 );
        TBCBuilder         bct_builder;

        auto               grid    = make_grid( gridfile );
        auto               fnspace = make_unique< fnspace_t >( grid.get() );
        auto               coord   = fnspace->build_coord();
        auto               ct      = ct_builder.build( coord.get() );
        auto               bct     = bct_builder.build( ct.get(), ct.get(), & adm_cond );
    
        //
        // build matrix
        //

        cout << endl << "━━ building H-matrix" << endl;

        slpbf_t                   bf( fnspace.get(), fnspace.get() );
        TBFCoeffFn< slpbf_t >     bf_coeff_fn( & bf );
        TPermCoeffFn< real_t >    coeff_fn( & bf_coeff_fn, bct->row_ct()->perm_i2e(), bct->col_ct()->perm_i2e() );
        TACAPlus< real_t >        lrapx( & coeff_fn );
        TTruncAcc                 acc = fixed_prec( eps );
        TDenseMBuilder< real_t >  h_builder( & coeff_fn, & lrapx );

        auto  A = h_builder.build( bct.get(), mat_form, acc, & progress );

        //
        // build RHSs
        //

        cout << endl << "━━ building RHSs " << std::endl;

        const size_t                      nrhs = 500;
        TQuadBEMRHS< fnspace_t, real_t >  rhs_build( 4 );
        TDenseMatrix                      RHS( fnspace->n_indices(), nrhs );

        timer.start();
        
        tbb::parallel_for( size_t(0), nrhs,
                           [&RHS,&rhs_build,&fnspace] ( const size_t  i )
                           {
                               TMyRHS  rhsfn( 100 * drand48() - 50,
                                              100 * drand48() - 50 );
                               auto    b     = rhs_build.build( fnspace.get(), & rhsfn );
                               auto    rhs_i = RHS.blas_rmat().column( i );
                               
                               BLAS::copy( cptrcast( b.get(), TScalarVector )->blas_rvec(), rhs_i );
                           } );
    
        // bring RHS into H-ordering (only rows, columns correspond to number of rhs)
        RHS.permute( ct->perm_e2i(), nullptr );
        
        timer.pause();
        cout << "    done in " << timer << endl;

        //
        // LU decomposition
        //

        cout << endl << "━━ H-LU factorisation" << endl;

        auto       A_copy = A->copy();
        TTruncAcc  fac_acc( fac_eps, 0.0 );
        auto       A_inv  = factorise_inv( A_copy.get(), fac_acc, & progress );
    
        cout << "    inversion error = " << format( "%.6e" ) % inv_approx_2( A.get(), A_inv.get() ) << endl;

        ////////////////////////////////////////////////
        //
        // solve
        //
        ////////////////////////////////////////////////
    
        {
            cout << endl << "━━ direct matrix solve" << endl;
            
            TDenseMatrix  SOL( A->cols(), nrhs );
    
            timer.start();

            solve_option_t  opt_LL( block_wise, unit_diag );
            solve_option_t  opt_UR( block_wise, general_diag );
            
            RHS.copy_to( & SOL );
            
            if ( mat_form == symmetric )
            {
                solve_lower_left( apply_normal,     A_copy.get(), & SOL, fac_acc, opt_LL );
                solve_diag_left(  apply_normal,     A_copy.get(), & SOL, fac_acc, opt_UR );
                solve_lower_left( apply_transposed, A_copy.get(), & SOL, fac_acc, opt_LL );
            }// if
            else
            {
                solve_lower_left( apply_normal, A_copy.get(), & SOL, fac_acc, opt_LL );
                solve_upper_left( apply_normal, A_copy.get(), & SOL, fac_acc, opt_UR );
            }// else

            timer.pause();
            cout << "    done in " << timer << endl;

            auto  X = RHS.copy();

            multiply( real_t(-1), apply_normal, A.get(), apply_normal, & SOL, real_t(1), X.get(), acc_exact );
            cout << "    |AX-B| = " << norm_F( X.get() ) << endl;
        }
        
        {
            cout << endl << "━━ direct vector solve" << endl;
            
            TDenseMatrix  SOL( A->cols(), nrhs );
    
            timer.start();

            tbb::parallel_for( size_t(0), nrhs,
                               [&RHS,&SOL,&A,&A_inv] ( const size_t  i )
                               {
                                   BLAS::Vector< real_t >  rhs_i( RHS.blas_rmat().column( i ) );
                                   BLAS::Vector< real_t >  sol_i( SOL.blas_rmat().column( i ) );
                                   TScalarVector           b( A->row_is(), rhs_i );
                                   TScalarVector           x( A->col_is(), sol_i );
                                   
                                   A_inv->apply( & b, & x );
                               } );

            timer.pause();
            cout << "    done in " << timer << endl;

            auto  X = RHS.copy();

            multiply( real_t(-1), apply_normal, A.get(), apply_normal, & SOL, real_t(1), X.get(), acc_exact );
            cout << "    |AX-B| = " << norm_F( X.get() ) << endl;
        }

        {
            cout << endl << "━━ iterative vector solve" << endl;
            
            TDenseMatrix  SOL( A->cols(), nrhs );
    
            timer.start();

            tbb::parallel_for( size_t(0), nrhs,
                               [&RHS,&SOL,&A,&A_inv,fac_eps] ( const size_t  i )
                               {
                                   BLAS::Vector< real_t >  rhs_i( RHS.blas_rmat().column( i ) );
                                   BLAS::Vector< real_t >  sol_i( SOL.blas_rmat().column( i ) );
                                   TScalarVector           b( A->row_is(), rhs_i );
                                   TScalarVector           x( A->col_is(), sol_i );
                                   
                                   solve( A.get(), & x, & b, A_inv.get() );
                               } );

            timer.pause();
            cout << "    done in " << timer << endl;

            auto  X = RHS.copy();

            multiply( real_t(-1), apply_normal, A.get(), apply_normal, & SOL, real_t(1), X.get(), acc_exact );
            cout << "    |AX-B| = " << norm_F( X.get() ) << endl;
        }

        DONE();
    }// try
    catch ( Error & e )
    {
        cout << "HLIB::error : " << e.what() << endl;
    }// catch
    catch ( std::exception & e )
    {
        cout << "std::exception : " << e.what() << endl;
    }// catch

    return 0;
}
