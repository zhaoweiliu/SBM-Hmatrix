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
using namespace HLIB;
using HLIB::uint;
using boost::format;

using real_t = HLIB::real;

//
// main function
//
int
main ( int argc, char ** argv )
{
    //
    // program options
    //
    
    string   matrix_file = "matrix.hm";
    string   rhs_file    = "";
    double   epslu       = 1e-8;  // <- corresponds to (almost) direct solve
    
    if ( argc > 1 )
        matrix_file = argv[1];
    if ( argc > 2 )
        rhs_file    = argv[2];
    
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

        // for easier access: convert to sparse
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
        cout << "  |S|_F                 = " << format( "%.6e" ) % norm_F( S ) << endl;

        if ( verbose( 3 ) )
        {
            mvis.pattern( true );
            mvis.print( S, "sparsealg_S" );
            mvis.pattern( false );
        }// if

        //
        // convert sparse to H
        //

        TTimer  timer( WALL_TIME );

        // set reasonable minimal cluster size
        size_t  nmin = max( size_t(80), min( S->avg_entries_per_row(), size_t(200) ) );

        cout << endl << "## converting to H ( nmin = " << nmin << " )" << endl;

        timer.start();

        TBFSAlgPartStrat    part_strat;    // <- good standard choice
        // TMLAlgPartStrat     part_strat;  // <- internal multi-level algorithm (might be better)  
        // TMETISAlgPartStrat  part_strat;  // <- prefer if available
        TAlgCTBuilder       ct_builder( & part_strat, nmin );
        TAlgNDCTBuilder     nd_ct_builder( & ct_builder, nmin );  // use nested dissection
        auto                ct = nd_ct_builder.build( S );
        
        if ( verbose( 1 ) )
            cout << "  built cluster tree (" << timer << ")" << endl;
        
        if( verbose( 3 ) )
        {
            TPSClusterVis  cvis;
            
            cvis.print( ct->root(), "sparsealg_ct" );
        }// if

        TWeakAlgAdmCond     adm_cond( S, ct->perm_i2e() );
        TBCBuilder          bct_builder;
        auto                bct( bct_builder.build( ct.get(), ct.get(), & adm_cond ) );
        
        if ( verbose( 1 ) )
            cout << "  built block cluster tree (" << timer << ")" << endl;
        cout << "    sparsity constant = " << bct->compute_c_sp() << endl;
        
        if( verbose( 3 ) )
        {
            TPSBlockClusterVis  bcvis;
            
            bcvis.print( bct->root(), "sparsealg_bct" );
        }// if

        TConsoleProgressBar  progress;
        TSparseMBuilder      h_builder( S, ct->perm_i2e(), ct->perm_e2i() );
        TTruncAcc            acc( real_t(0.0) );
        auto                 A = h_builder.build( bct.get(), acc, & progress );

        timer.pause();

        TPermMatrix  PA( ct->perm_i2e(), A.get(), ct->perm_e2i() );
        
        cout << "  built H-matrix in     " << timer << endl;
        cout << "    size of H-matrix  = " << Mem::to_string( A->byte_size() ) << endl;
        cout << "    |A|_F             = " << format( "%.6e" ) % norm_F( A.get() ) << endl;
        cout << "    |S-A|_2           = " << format( "%.6e" ) % diff_norm_2( S, & PA ) << endl;

        if( verbose( 3 ) )
            mvis.svd( true ).print( A.get(), "sparsealg_A" );

        //
        // optional: to fix most common error (singular diagonal block), try to
        //           add scaled identity to matrix before LU factorisation
        //           (this only changes preconditioner, not the system matrix!)
        //

        if ( false )
        {
            const real_t  factor = norm_F( A.get() ) * 1e-8;
            
            cout << endl
                 << "## fixing diagonal" << endl
                 << "    adding λ·I (λ = " << format( "%.6e" ) % factor << ")" << endl;

            add_identity( A.get(), factor );
        }// if
        
        //
        // LU decomposition
        //

        const TTruncAcc  fac_acc( epslu );

        cout << endl << "## H factorisation ( ε = " << format( "%.2e" ) % epslu << " )" << endl;

        timer.start();
        
        auto  A_inv = factorise_inv( A.get(), fac_acc, & progress );
    
        timer.pause();
        cout << "    done in " << timer << endl;
        
        if( verbose( 3 ) )
            mvis.print( A.get(), "sparsealg_LU" );

        auto  PA_inv = make_unique< TPermMatrix >( ct->perm_i2e(), A_inv.get(), ct->perm_e2i() );
        
        cout << "    size of LU factor = " << Mem::to_string( A->byte_size() ) << endl;
        cout << "    inv. error wrt S  = " << format( "%.6e" ) % inv_approx_2( S, PA_inv.get() ) << endl;
        cout << "    condest(LU)       = " << format( "%.6e" ) % condest( PA_inv.get() ) << endl;

        //
        // solve with LU decomposition
        //

        cout << endl << "## solving system" << endl;

        TAutoSolver            solver( 100 );
        TSolverInfo            solve_info( true );
        unique_ptr< TVector >  b, x, sol;

        x = S->col_vector();

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
    
        timer.start();
    
        solver.solve( S, x.get(), b.get(), PA_inv.get(), & solve_info );

        timer.pause();
        
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
            
            cout << "  |Sx-b|₂  = " << format( "%.6e" ) % b2->norm2() << endl;
        }
    
        DONE();
    }// try
    catch ( Error & e )
    {
        cout << e.to_string() << endl;
    }// catch
    
    return 0;
}
