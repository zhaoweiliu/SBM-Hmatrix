//
// Project     : HLib
// File        : sparsebsp.cc
// Description : example for sparse matrices with geometrical clustering
// Author      : Ronald Kriemann
// Copyright   : Max Planck Institute MIS 2004-2014. All Rights Reserved.
//
//

#include <iostream>
#include <iomanip>
#include <string>

#include <boost/format.hpp>

#include "hlib.hh"

using namespace std;
using namespace HLIB;
using HLIB::uint;
using boost::format;

using  real_t = HLIB::real;

//
// main function
//
int
main ( int argc, char ** argv )
{
    string   matrix_file = "matrix.hm";
    string   coord_file  = "coord.hm";
    string   rhs_file    = "";
    double   epslu       = 1e-8;  // <- corresponds to (almost) direct solve
    
    if ( argc > 1 )
        matrix_file = argv[1];
    if ( argc > 2 )
        coord_file  = argv[2];
    if ( argc > 3 )
        rhs_file    = argv[3];

    try
    {
        //
        // init HLIBpro
        //
    
        INIT();

        CFG::set_verbosity( 1 );

        //
        // load sparse matrix
        //
    
        cout << "## importing sparse matrix from " << matrix_file << endl;

        TPSMatrixVis  mvis;
        auto          M = read_matrix( matrix_file.c_str() );
    
        if ( M == nullptr )
        {
            cout << "error while reading matrix" << endl;
            exit( 1 );
        }// if
            
        if ( ! IS_TYPE( M, TSparseMatrix ) )
        {
            cout << "given matrix is not sparse (" << M->typestr() << ")" << endl;
            exit( 1 );
        }// if

        auto  S = ptrcast( M.get(), TSparseMatrix );
    
        cout << "  matrix has dimension " << S->rows() << " x " << S->cols() << endl;
        cout << "    no of non-zeroes    = " << S->n_non_zero() << " ("
             << format( "%.2f" ) % ( 1000.0 * double(S->n_non_zero()) / Math::square(double(S->rows())) ) << "‰"
             << ", avg=" << S->avg_entries_per_row()
             << ", max=" << S->max_entries_per_row() 
             << ")" << endl;
        cout << "    matrix is             " << ( S->is_complex() ? "complex" : "real" ) << " valued" << endl;
        cout << "    format              = ";
        if      ( S->is_nonsym()    ) cout << "non symmetric" << endl;
        else if ( S->is_symmetric() ) cout << "symmetric" << endl;
        else if ( S->is_hermitian() ) cout << "hermitian" << endl;
        cout << "  size of sparse matrix = " << Mem::to_string( S->byte_size() ) << endl;
        cout << "  |S|_2                 = " << norm_2( S ) << endl;

        if( verbose( 3 ) )
        {
            mvis.pattern( true ).print( S, "sparsebsp_S.ps" );
        }// if

        //
        // load coordinates
        //
    
        cout << endl << "## importing coordinates from " << coord_file << endl;
        
        auto  coord = read_coord( coord_file.c_str() );

        //
        // convert sparse to H
        //

        TTimer  timer( WALL_TIME );

        cout << endl << "## converting to H" << endl;

        timer.start();
                                       
        const uint         nmin = 80;
        TAutoBSPPartStrat  part_strat;
        TBSPNDCTBuilder    ct_builder( S, & part_strat, nmin );
        auto               ct = ct_builder.build( coord.get(), S );
        
        if( verbose( 3 ) )
        {
            TPSClusterVis  cvis;
            
            cvis.print( ct->root(), "sparsebsp_ct.ps" );
        }// if
        
        TStdGeomAdmCond    adm_cond;
        TBCBuilder         bct_builder;
        auto               bct = bct_builder.build( ct.get(), ct.get(), & adm_cond );
        
        if( verbose( 2 ) )
        {
            TPSBlockClusterVis  bcvis;
            
            bcvis.print( bct->root(), "sparsebsp_bct.ps" );
        }// if
        
        TConsoleProgressBar  progress;
        TSparseMBuilder      h_builder( S, ct->perm_i2e(), ct->perm_e2i() );
        TTruncAcc            acc( real_t(0.0) );
        auto                 A = h_builder.build( bct.get(), acc, & progress );

        timer.pause();
        cout << "    done in " << timer << endl;
        cout << "    size of H-matrix  = " << Mem::to_string( A->byte_size() ) << endl;
        cout << "    sparsity constant = " << bct->compute_c_sp() << endl;

        if ( verbose( 3 ) )
        {
            mvis.pattern( false ).svd( true ).print( A.get(), "sparsebsp_A.ps" );
        }// if
        
        //
        // LU decomposition
        //

        cout << endl << "## H-LU factorisation ( ε = " << format( "%.2e" ) % epslu << " )" << endl;
        
        const TTruncAcc  fac_acc( epslu );
        
        timer.start();
        
        auto  A_inv = factorise_inv( A.get(), real_t(epslu), & progress );
    
        timer.pause();
        cout << "    done in " << timer << endl;

        if ( verbose( 3 ) )
            mvis.print( A.get(), "sparsebsp_LU.ps" );

        auto  PA_inv = make_unique< TPermMatrix >( ct->perm_i2e(), A_inv.get(), ct->perm_e2i() );
        
        cout << "    size of LU factor = " << Mem::to_string( A->byte_size() ) << endl;
        cout << "    inv. error wrt S  = " << format( "%.6e" ) % inv_approx_2( S, PA_inv.get() ) << endl;
        cout << "    condest(LU)       = " << format( "%.6e" ) % condest( PA_inv.get() ) << endl;
        
        //
        // solve with LU decomposition
        //

        cout << endl << "## solving system" << endl;

        TAutoSolver            solver( 1000 );
        TSolverInfo            solve_info;
        unique_ptr< TVector >  b, x, sol;

        if ( rhs_file != "" )
            b = read_vector( rhs_file.c_str() );
        else
        {
            //
            // use solution identical to 1
            //
            
            b   = A->row_vector();
            sol = S->col_vector();
            
            sol->fill( 1.0 );
            S->mul_vec( 1.0, sol.get(), 0.0, b.get(), apply_normal );
        }// if
    
        x = S->col_vector();

        timer.start();
    
        solver.solve( S, x.get(), b.get(), PA_inv.get(), & solve_info );
    
        if ( solve_info.has_converged() )
            cout << "  converged in " << timer << " and "
                      << solve_info.n_iter() << " steps with rate " << solve_info.conv_rate()
                      << ", |r| = " << solve_info.res_norm() << endl;
        else
            cout << "  not converged in " << timer << " and "
                      << solve_info.n_iter() << " steps " << endl;

        {
            auto  b2 = b->copy();

            S->mul_vec( 1.0, x.get(), -1.0, b2.get() );
            cout << "  |Sx-b| = " << b2->norm2() << endl;
        }
    
        DONE();
    }// try
    catch ( Error & e )
    {
        cout << e.to_string() << endl;
    }// catch
    
    return 0;
}
