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
    unsigned           n;
    unsigned           out_n = 50000;
    hlib_coord_t       coord;
    hlib_vector_t      b = NULL, x = NULL, sol = NULL;
    hlib_matrix_t      S;
    int                info;
    char               fmatrix[256] = "\0";
    char               fcoord[256] = "\0";
    char               frhs[256] = "\0";
    char               fsol[256] = "\0";
    hlib_acc_t         acc       = hlib_acc_fixed_eps( 1e-4 );
    double             start;
    hlib_solver_t      solver;
    hlib_solve_info_t  solve_info;
    int                is_sym = 0;

    hlib_init( & info );                                                           CHECK_INFO;

    if ( argc > 4 )
        strcpy( fsol, argv[3] );
    if ( argc > 3 )
        strcpy( frhs, argv[2] );
    if ( argc > 2 )
    {
        strcpy( fmatrix, argv[1] );
        strcpy( fcoord, argv[2] );
    }
    else
    {
        printf( "usage: sparsebsp <matrix> <coord> [rhs sol]\n" );
        exit( 1 );
    }
    
    hlib_set_verbosity( 2 );

    printf( "reading sparse system\n" );
    coord = hlib_load_coord( fcoord, & info );                                     CHECK_INFO;
    S = hlib_load_matrix( fmatrix, & info );                                       CHECK_INFO;

    n      = hlib_matrix_rows( S, & info );                                        CHECK_INFO;
    is_sym = ( hlib_matrix_is_sym( S, & info ) || hlib_matrix_is_herm( S, & info ) ? 1 : 0 );
    
    printf( "  matrix has dimension %lu x %lu\n", hlib_matrix_rows( S, & info ),
            hlib_matrix_cols( S, & info ) );
    if ( n < out_n ) hlib_matrix_print_ps( S, "sparsebsp_S", HLIB_MATIO_PATTERN, & info ); CHECK_INFO;
    if ( is_sym )
        printf( "  matrix is symmetric/hermitian\n" );

    if ( frhs[0] != '\0' )
        b = hlib_load_vector( frhs, & info );                                      CHECK_INFO;
    
    if ( fsol[0] != '\0' )
        sol = hlib_load_vector( fsol, & info );                                    CHECK_INFO;
    
    x = hlib_matrix_col_vector( S, & info );                                       CHECK_INFO;

    if ( b == NULL )
    {
        b = hlib_matrix_row_vector( S, & info );                                   CHECK_INFO;
        hlib_vector_fill_rand( x, & info );                                        CHECK_INFO;
        hlib_matrix_mulvec( 1.0, S, x, 0.0, b, HLIB_MATOP_NORM, & info );          CHECK_INFO;
    }// if
    
    if ( sol != NULL )
    {
        hlib_matrix_mulvec( 1.0, S, sol, 0.0, x, HLIB_MATOP_NORM, & info );        CHECK_INFO;
        hlib_vector_axpy( -1.0, b, x, & info );
        printf( "  |S * sol - rhs| = %.2e\n", hlib_vector_norm2( x, & info ) );
    }

    solver = hlib_solver_auto( & info );                                           CHECK_INFO;

    printf( "solving sparse system ...\n" );
    start = hlib_walltime();
    hlib_solver_solve( solver, (hlib_linearoperator_t) S, x, b, NULL, & solve_info, & info ); CHECK_INFO;

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
        
        printf( "\n## LU decomposition with nested dissection\n" );

        printf( "  building cluster tree\n" );
        ct = hlib_clt_build_bsp_nd( coord, S, HLIB_BSP_AUTO, 20, & info );           CHECK_INFO;
        if ( n < out_n ) hlib_clt_print_ps( ct, "sparsebsp_ct_nd", & info );      CHECK_INFO;
        
        printf( "  building block cluster tree\n" );
        adm = hlib_admcond_geom( HLIB_ADM_AUTO, 2.0, & info );                       CHECK_INFO;
        bct = hlib_bct_build( ct, ct, adm, & info );                                 CHECK_INFO;
        if ( n < out_n ) hlib_bct_print_ps( bct, "sparsebsp_bct_nd", & info );    CHECK_INFO;

        printf( "  converting sparse matrix to H-matrix\n" );
        A = hlib_matrix_build_sparse( bct, S, acc, & info );                         CHECK_INFO;
        printf( "    size of H-matrix = %.2f MB\n",
                ((double) hlib_matrix_bytesize( A, & info )) / (1024.0 * 1024.0) );
        if ( n < out_n ) hlib_matrix_print_ps( A, "sparsebsp_A_nd", HLIB_MATIO_PATTERN, & info ); CHECK_INFO;

        printf( "  LU factorising H-matrix\n" );
        start = hlib_walltime();
        LU = hlib_matrix_factorise_inv( A, acc, & info );                            CHECK_INFO;
        printf( "    done in %.1f seconds\n", hlib_walltime() - start );
        if ( n < out_n ) hlib_matrix_print_ps( A, "sparsebsp_LU_nd", HLIB_MATIO_SVD, & info ); CHECK_INFO;

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
     * solve by matrix inversion
     */
    
    if ( 1 )
    {
        hlib_clustertree_t       ct;
        hlib_admcond_t           adm;
        hlib_blockclustertree_t  bct;
        hlib_matrix_t            A;
        hlib_linearoperator_t    PA;
        
        printf( "\n## Inversion with Geometrical Bissection\n" );

        printf( "  building cluster tree\n" );
        ct = hlib_clt_build_bsp( coord, HLIB_BSP_AUTO, 20, & info );                 CHECK_INFO;
        if ( n < out_n ) hlib_clt_print_ps( ct, "sparsebsp_ct", & info );         CHECK_INFO;

        printf( "  building block cluster tree\n" );
        adm = hlib_admcond_geom( HLIB_ADM_AUTO, 2.0, & info );                       CHECK_INFO;
        bct = hlib_bct_build( ct, ct, adm, & info );                                 CHECK_INFO;
        if ( n < out_n ) hlib_bct_print_ps( bct, "sparsebsp_bct", & info );       CHECK_INFO;

        printf( "  converting sparse matrix to H-matrix\n" );
        A = hlib_matrix_build_sparse( bct, S, acc, & info );                         CHECK_INFO;
        printf( "    size of H-matrix = %.2f MB\n",
                ((double) hlib_matrix_bytesize( A, & info )) / (1024.0 * 1024.0) );
        if ( n < out_n ) hlib_matrix_print_ps( A, "sparsebsp_A", HLIB_MATIO_PATTERN, & info ); CHECK_INFO;

        printf( "  inverting H-matrix\n" );
        start = hlib_walltime();
        hlib_matrix_inv( A, acc, & info );                                           CHECK_INFO;
        printf( "    done in %.1f seconds\n", hlib_walltime() - start );
        if ( n < out_n ) hlib_matrix_print_ps( A, "sparsebsp_I", HLIB_MATIO_SVD, & info );CHECK_INFO;

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
                                                     PA, & info ) );               CHECK_INFO;
        
        printf( "  solving preconditioned system ...\n" );
        start = hlib_walltime();
        hlib_solver_solve( solver, (hlib_linearoperator_t) S, x, b,
                           PA, & solve_info, & info );                             CHECK_INFO;
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

        hlib_linearoperator_free( PA, & info );                                    CHECK_INFO;
        hlib_matrix_free( A, & info );                                             CHECK_INFO;
        hlib_bct_free( bct, & info );                                              CHECK_INFO;
        hlib_admcond_free( adm, & info );                                          CHECK_INFO;
        hlib_clt_free( ct, & info );                                               CHECK_INFO;
    }
    
    hlib_solver_free( solver, & info );                                            CHECK_INFO;
    hlib_vector_free( x, & info );                                                 CHECK_INFO;
    if ( sol != NULL )
        hlib_vector_free( sol, & info );                                           CHECK_INFO;
    hlib_vector_free( b, & info );                                                 CHECK_INFO;
    hlib_matrix_free( S, & info );                                                 CHECK_INFO;
    hlib_coord_free( coord, & info );                                              CHECK_INFO;

    hlib_done( & info );                                                           CHECK_INFO;
    
    return 0;
}
