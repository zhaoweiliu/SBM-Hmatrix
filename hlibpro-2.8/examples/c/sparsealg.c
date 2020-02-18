/*
 * Project     : HLib
 * File        : sparsealg.c
 * Description : example for solving sparse systems
 * Author      : Ronald Kriemann
 * Copyright   : Max Planck Institute MIS 2004-2014. All Rights Reserved.
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <hlib-c.h>

#define CHECK_INFO  { if ( info != HLIB_NO_ERROR ) \
                      { char buf[1024]; hlib_error_desc( buf, 1024 ); \
                        printf( "\n%s\n\n", buf ); exit(1); } }

/*
 * main function
 */

int
main ( int argc, char ** argv )
{
    hlib_matrix_t      S = NULL;
    hlib_vector_t      b = NULL, x = NULL, sol = NULL;
    int                info;
    int                n, out_n = 50000;
    char               fmatrix[256] = "\0";
    char               frhs[256] = "\0";
    char               fsol[256] = "\0";
    hlib_acc_t         acc       = hlib_acc_fixed_eps( 1e-4 );
    double             start;
    hlib_solver_t      solver;
    hlib_solve_info_t  solve_info;
    int                is_sym = 0;

    hlib_init( & info );                                                          CHECK_INFO;

    if ( argc > 3 )
        strcpy( fsol, argv[3] );
    if ( argc > 2 )
        strcpy( frhs, argv[2] );
    if ( argc > 1 )
        strcpy( fmatrix, argv[1] );
    else
    {
        printf( "usage: sparsealg <matrix> [<rhs> <sol>] \n" );
        exit( 1 );
    }

    hlib_set_verbosity( 2 );
    hlib_set_abs_eps( 0.0 );                                    
    
    printf( "importing sparse matrix from %s\n", fmatrix );
    S = hlib_load_matrix( fmatrix, & info );                                      CHECK_INFO;
    is_sym = ( hlib_matrix_is_sym( S, & info ) || hlib_matrix_is_herm( S, & info ) ? 1 : 0 );
        
    printf( "  matrix has dimension %lu x %lu\n", hlib_matrix_rows( S, & info ),
            hlib_matrix_cols( S, & info ) );
    printf( "  size of sparse matrix = %.2f MB\n",
            ((double) hlib_matrix_bytesize( S, & info )) / (1024.0 *1024.0) );
    printf( "  symmetric             = %d\n", is_sym );

    n = hlib_matrix_rows( S, & info );                                            CHECK_INFO;
    if ( n < out_n ) hlib_matrix_print_ps( S, "sparsealg_S", HLIB_MATIO_PATTERN, & info ); CHECK_INFO;

    if ( 0 )
    {
        /*
         * test S
         */

        hlib_vector_t  xone, xrand;
        hlib_vector_t  yone, yrand;

        xone  = hlib_matrix_col_vector( S, & info );
        xrand = hlib_matrix_col_vector( S, & info );
        yone  = hlib_matrix_row_vector( S, & info );
        yrand = hlib_matrix_row_vector( S, & info );

        hlib_vector_fill( xone, 1.0, & info );
        hlib_vector_scale( xone, 1.0 / hlib_vector_norm2( xone, & info ), & info );
        hlib_vector_fill_rand( xrand, & info );
        hlib_vector_scale( xrand, 1.0 / hlib_vector_norm2( xrand, & info ), & info );

        hlib_matrix_mulvec( 1.0, S, xone, 0.0, yone, HLIB_MATOP_NORM, & info );
        hlib_matrix_mulvec( 1.0, S, xrand, 0.0, yrand, HLIB_MATOP_NORM, & info );

        printf( "|S * 1| = %.4e\n", hlib_vector_norm2( yone, & info ) );
        printf( "|S * ?| = %.4e\n", hlib_vector_norm2( yrand, & info ) );
        printf( "|S|_2   = %.4e\n", hlib_matrix_norm_spectral( S, & info ) );
    }

    if ( frhs[0] != '\0' )
        b = hlib_load_vector( frhs, & info );                                     CHECK_INFO;
    
    if ( fsol[0] != '\0' )
        sol = hlib_load_vector( fsol, & info );                                   CHECK_INFO;
    
    x = hlib_matrix_col_vector( S, & info );                                      CHECK_INFO;

    if ( b == NULL )
    {
        b = hlib_matrix_row_vector( S, & info );                                  CHECK_INFO;
        hlib_vector_fill_rand( x, & info );                                       CHECK_INFO;
        hlib_matrix_mulvec( 1.0, S, x, 0.0, b, HLIB_MATOP_NORM, & info );         CHECK_INFO;
    }// if
    
    if ( sol != NULL )
    {
        hlib_matrix_mulvec( 1.0, S, sol, 0.0, x, HLIB_MATOP_NORM, & info );       CHECK_INFO;
        hlib_vector_axpy( -1.0, b, x, & info );
        printf( "  |S * sol - rhs| = %.2e\n", hlib_vector_norm2( x, & info ) );
    }

    solver = hlib_solver_auto( & info );                                          CHECK_INFO;

    /* turn start value initialisation on/off */
    hlib_solver_initialise_start_value( solver, 1, & info );                      CHECK_INFO;
    
    if ( 0 )
    {
        printf( "solving sparse system ... \n" );
        start = hlib_walltime();
        hlib_solver_solve( solver, (hlib_linearoperator_t) S, x, b, NULL, & solve_info, & info );         CHECK_INFO;
        if ( solve_info.converged )
            printf( "  converged in %.1f seconds and %u steps with rate %.2e, |r| = %.2e\n",
                    hlib_walltime() - start, solve_info.steps, solve_info.conv_rate,
                    solve_info.res_norm );
        else
            printf( "  not converged in %.1f seconds and %u steps\n",
                    hlib_walltime() - start, solve_info.steps );
    
        if ( sol != NULL )
        {
            hlib_vector_axpy( -1.0, sol, x, & info );
            printf( "  |x-sol|/|sol| = %.2e\n",
                    hlib_vector_norm2( x, & info ) / hlib_vector_norm2( sol, & info ) );
        }
    }

    /*
     * solve with LU decomposition and nested dissection
     */

    if ( 1 )
    {
        hlib_clustertree_t       ct;
        hlib_admcond_t           adm;
        hlib_blockclustertree_t  bct;
        hlib_matrix_t            A;
        hlib_linearoperator_t    LU, PLU;
        
        printf( "\n## LU decomposition with algebraic nested dissection\n" );

        printf( "  converting sparse matrix to H-matrix\n" );
        start = hlib_walltime();
        
        ct  = hlib_clt_build_alg_nd( S, HLIB_ALG_AUTO, 20, & info );              CHECK_INFO;
        adm = hlib_admcond_alg( HLIB_ADM_AUTO, 2.0, S,
                                hlib_clt_perm_i2e( ct, & info ),
                                hlib_clt_perm_i2e( ct, & info ),
                                & info );                                         CHECK_INFO;
        bct = hlib_bct_build( ct, ct, adm, & info );                              CHECK_INFO;
        A   = hlib_matrix_build_sparse( bct, S, acc, & info );                    CHECK_INFO;
        
        printf( "    done in %.1f seconds\n", hlib_walltime() - start );
        printf( "    size of H-matrix = %.2f MB\n",
                ((double) hlib_matrix_bytesize( A, & info )) / (1024.0 * 1024.0) );

        if ( n < out_n ) hlib_clt_print_ps( ct, "sparsealg_ct_nd", & info );      CHECK_INFO;
        if ( n < out_n ) hlib_bct_print_ps( bct, "sparsealg_bct_nd", & info );    CHECK_INFO;
        if ( n < out_n ) hlib_matrix_print_ps( A, "sparsealg_A_nd", HLIB_MATIO_SVD, & info ); CHECK_INFO;

        printf( "  LU factorising H-matrix ... \n" );
        start = hlib_walltime();
        LU = hlib_matrix_factorise_inv( A, acc, & info );                               CHECK_INFO;
        printf( "    done in %.1f seconds\n", hlib_walltime() - start );
        if ( n < out_n ) hlib_matrix_print_ps( A, "sparsealg_LU_nd", HLIB_MATIO_SVD, & info ); CHECK_INFO;

        printf( "  size of LU factor = %.2f MB\n",
                ((double) hlib_matrix_bytesize( A, & info )) / (1024.0 * 1024.0) );

        /* apply permutations to compare with S */
        PLU = hlib_perm_linearoperator( hlib_clt_perm_i2e( ct, & info ),
                                        LU,
                                        hlib_clt_perm_e2i( ct, & info ),
                                        & info );
        CHECK_INFO;
        
        printf( "  inversion error = %.4e\n",
                hlib_linearoperator_norm_inv_approx( (hlib_linearoperator_t) S, PLU, & info ) ); CHECK_INFO;

        printf( "  solving preconditioned system ...\n" );
        start = hlib_walltime();
        hlib_solver_solve( solver, (hlib_linearoperator_t) S, x, b, PLU, & solve_info, & info ); CHECK_INFO;
        if ( solve_info.converged )
            printf( "  converged in %.1f seconds and %u steps with rate %.2e, |r| = %.2e\n",
                    hlib_walltime() - start, solve_info.steps, solve_info.conv_rate,
                    solve_info.res_norm );
        else
            printf( "  not converged in %.1f seconds and %u steps\n",
                    hlib_walltime() - start, solve_info.steps );

        {
            hlib_vector_t  b2 = hlib_vector_copy( b, & info );

            hlib_matrix_mulvec( 1.0, S, x, -1.0, b2, HLIB_MATOP_NORM, & info );
            printf( "  |Sx-b| = %.4e\n", hlib_vector_norm2( b2, & info ) );

            hlib_vector_free( b2, & info );
        }
        
        if ( sol != NULL )
        {
            hlib_vector_axpy( -1.0, sol, x, & info );
            printf( "  |x-sol|/|sol| = %.2e\n",
                    hlib_vector_norm2( x, & info ) / hlib_vector_norm2( sol, & info ) );
        }

        hlib_linearoperator_free( PLU, & info );                                  CHECK_INFO;
        hlib_linearoperator_free( LU, & info );                                   CHECK_INFO;
        hlib_matrix_free( A, & info );                                            CHECK_INFO;
        hlib_bct_free( bct, & info );                                             CHECK_INFO;
        hlib_admcond_free( adm, & info );                                         CHECK_INFO;
        hlib_clt_free( ct, & info );                                              CHECK_INFO;
    }

    /*
     * solve with matrix inversion
     */
    
    if ( 0 )
    {
        hlib_clustertree_t       ct;
        hlib_admcond_t           adm;
        hlib_blockclustertree_t  bct;
        hlib_matrix_t            A;
        hlib_linearoperator_t    PA;

        printf( "\n## Inversion with algebraic bisection\n" );

        printf( "  converting sparse matrix to H-matrix ... \n" );
        start = hlib_walltime();
        
        ct  = hlib_clt_build_alg( S, HLIB_ALG_AUTO, 20, & info );                 CHECK_INFO;
        adm = hlib_admcond_alg( HLIB_ADM_AUTO, 2.0, S,
                                hlib_clt_perm_i2e( ct, & info ),
                                hlib_clt_perm_i2e( ct, & info ),
                                & info );                                         CHECK_INFO;
        bct = hlib_bct_build( ct, ct, adm, & info );                              CHECK_INFO;
        A   = hlib_matrix_build_sparse( bct, S, acc, & info );                    CHECK_INFO;

        printf( "    done in %.1f seconds\n", hlib_walltime() - start );
        printf( "    size of H-matrix = %.2f MB\n",
                ((double) hlib_matrix_bytesize( A, & info )) / (1024.0 * 1024.0) );
        
        if ( n < out_n ) hlib_clt_print_ps( ct, "sparsealg_ct", & info );         CHECK_INFO; 
        if ( n < out_n ) hlib_bct_print_ps( bct, "sparsealg_bct", & info );       CHECK_INFO;
        if ( n < out_n ) hlib_matrix_print_ps( A, "sparsealg_A", HLIB_MATIO_SVD, & info ); CHECK_INFO;

        printf( "  inverting H-matrix ... \n" );
        start = hlib_walltime();
        hlib_matrix_inv( A, acc, & info );                                        CHECK_INFO;
        printf( "    done in %.1f seconds\n", hlib_walltime() - start );
        if ( n < out_n ) hlib_matrix_print_ps( A, "sparsealg_I", HLIB_MATIO_SVD, & info );CHECK_INFO;

        printf( "  size of Inverse = %.2f MB\n",
                ((double) hlib_matrix_bytesize( A, & info )) / (1024.0 * 1024.0) );

        /* apply permutations to compare with S */
        PA = hlib_perm_matrix( hlib_clt_perm_i2e( ct, & info ),
                               A,
                               hlib_clt_perm_e2i( ct, & info ),
                               & info );

        CHECK_INFO;

        printf( "  inversion error = %.4e\n",
                hlib_linearoperator_norm_inv_approx( (hlib_linearoperator_t) S, 
                                                     PA, & info ) );              CHECK_INFO;
        
        printf( "  solving preconditioned system ...\n" );
        start = hlib_walltime();
        hlib_solver_solve( solver, (hlib_linearoperator_t) S, x, b,
                           (hlib_linearoperator_t) PA, & solve_info, & info );    CHECK_INFO;
        if ( solve_info.converged )
            printf( "  converged in %.1f seconds and %u steps with rate %.2e, |r| = %.2e\n",
                    hlib_walltime() - start, solve_info.steps, solve_info.conv_rate,
                    solve_info.res_norm );
        else
            printf( "  not converged in %.1f seconds and %u steps\n",
                    hlib_walltime() - start, solve_info.steps );

        {
            hlib_vector_t  b2 = hlib_vector_copy( b, & info );

            hlib_matrix_mulvec( 1.0, S, x, -1.0, b2, HLIB_MATOP_NORM, & info );
            printf( "  |Sx-b| = %.4e\n", hlib_vector_norm2( b2, & info ) );

            hlib_vector_free( b2, & info );
        }
        
        if ( sol != NULL )
        {
            hlib_vector_axpy( -1.0, sol, x, & info );
            printf( "  |x-sol|/|sol| = %.2e\n",
                    hlib_vector_norm2( x, & info ) / hlib_vector_norm2( sol, & info ) );
        }

        hlib_linearoperator_free( PA, & info );                                   CHECK_INFO;
        hlib_matrix_free( A, & info );                                            CHECK_INFO;
        hlib_bct_free( bct, & info );                                             CHECK_INFO;
        hlib_admcond_free( adm, & info );                                         CHECK_INFO;
        hlib_clt_free( ct, & info );                                              CHECK_INFO;
    }

    hlib_solver_free( solver, & info );                                           CHECK_INFO;
    hlib_vector_free( x, & info );                                                CHECK_INFO;
    hlib_vector_free( b, & info );                                                CHECK_INFO;
    if ( sol != NULL )
        hlib_vector_free( sol, & info );                                          CHECK_INFO;
    hlib_matrix_free( S, & info );                                                CHECK_INFO;

    hlib_done( & info );                                                          CHECK_INFO;
    
    return 0;
}
