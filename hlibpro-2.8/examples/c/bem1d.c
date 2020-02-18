/*
 * Project     : HLib
 * File        : bem1d.c
 * Description : example for 1d BEM problem
 * Author      : Ronald Kriemann
 * Copyright   : Max Planck Institute MIS 2004-2014. All Rights Reserved.
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <hlib-c.h>

#define CHECK_INFO  { if ( info != HLIB_NO_ERROR ) \
                      { char buf[1024]; hlib_error_desc( buf, 1024 ); \
                        printf( "\n%s\n\n", buf ); exit(1); } }

/*
 * function for evaluating kernel: log|x-y| in [0,1]
 */

void
kernel ( const size_t   n,
         const int *    rowidx,
         const size_t   m,
         const int *    colidx,
         hlib_real_t *  matrix,
         void *         arg )
{
    size_t  rowi, colj;
    double  h = *((double*) arg);

    for ( colj = 0; colj < m; colj++ )
    {
        const int idx1 = colidx[colj];
        
        for ( rowi = 0; rowi < n; rowi++ )
        {
            const int idx0 = rowidx[rowi];
            double    value;

            if ( idx0 == idx1 ) 
                value = -1.5*h*h + h*h*log(h);
            else
            {
                const double dist = h * ( fabs( (double) (idx0-idx1) ) - 1.0 );
                const double t1   = dist+1.0*h;
                const double t2   = dist+2.0*h;
                
                value = ( - 1.5*h*h + 0.5*t2*t2*log(t2) - t1*t1*log(t1) );
                
                if ( fabs(dist) > 1e-8 )
                    value += 0.5*dist*dist*log(dist);
            }
            
            matrix[(colj*n) + rowi] = -value;
        }
    }
}

/*
 * function for evaluating rhs (for solution u = -1)
 */
double
rhs ( int i, int n )
{
    const double  a     = ((double)i)     / ((double) n);
    const double  b     = ((double)i+1.0) / ((double) n);
    double        value = -1.5 * (b - a);
    
    if ( fabs( b )       > 1e-8 ) value += 0.5*b*b*log(b);
    if ( fabs( a )       > 1e-8 ) value -= 0.5*a*a*log(a);
    if ( fabs( 1.0 - b ) > 1e-8 ) value -= 0.5*(1.0-b)*(1.0-b)*log(1.0-b);
    if ( fabs( 1.0 - a ) > 1e-8 ) value += 0.5*(1.0-a)*(1.0-a)*log(1.0-a); 
    
    return value;
}

/*
 * main function
 */
int
main ( int argc, char ** argv )
{
    int                      n;
    double                   h;
    double **                vertices;
    int                      i, bytesize;
    int                      info;
    hlib_real_t *            b_arr;
    hlib_coord_t             coord;
    hlib_vector_t            x, b, one;
    hlib_clustertree_t       ct;
    hlib_admcond_t           adm;
    hlib_blockclustertree_t  bct;
    hlib_matrix_t            A;
    double                   error, start;
    hlib_solver_t            solver;
    hlib_solve_info_t        solve_info;
    hlib_acc_t               acc   = hlib_acc_fixed_eps( 1e-4 );
    
    if ( argc > 1 ) n = atoi( argv[1] );
    else            n = 8;
    if ( n <= 0 )   n = 8;  /* sanity check */

    h = 1.0 / ((double) n);
    
    hlib_init( & info );                                    CHECK_INFO;
    hlib_set_verbosity( 2 );
    
    /*
     * build coordinates of 1d example
     */

    vertices = (double**) malloc( n * sizeof(double*) );

    for ( i = 0; i < n; i++ )
    {
        vertices[i]    = (double*) malloc( sizeof(double) );
        vertices[i][0] = h * ((double) i) + (h / 2.0);
    }
    
    printf( "building cluster tree\n" );
    coord = hlib_coord_import( n, 1, vertices, NULL, & info );              CHECK_INFO;
    ct    = hlib_clt_build_bsp( coord, HLIB_BSP_AUTO, 20, & info );         CHECK_INFO;
    hlib_clt_print_ps( ct, "bem1d_ct", & info );                         CHECK_INFO;
    
    printf( "building block cluster tree\n" );
    adm = hlib_admcond_geom( HLIB_ADM_STD_MIN, 2.0, & info );               CHECK_INFO;
    bct = hlib_bct_build( ct, ct, adm, & info );                            CHECK_INFO;
    hlib_bct_print_ps( bct, "bem1d_bct", & info );                       CHECK_INFO;

    printf( "building BEM matrix using ACA+ ...\n" );
    start = hlib_walltime();

    A = hlib_matrix_build_coeff( bct, kernel, & h, HLIB_LRAPX_ACAPLUS,
                                 acc, 1, & info );                          CHECK_INFO;
    
    printf( "  done in %.1f seconds\n", hlib_walltime() - start );
    hlib_matrix_print_ps( A, "bem1d_A", HLIB_MATIO_SVD, & info );        CHECK_INFO;

    bytesize = hlib_matrix_bytesize( A, & info );
    
    printf( "  compression ratio = %.2f%% (%.2f MB compared to %.2f MB)\n",
            100.0 * ((double) bytesize) /
            (((double) n) * ((double) n) * ((double) sizeof(double))),
            ((double) bytesize) / (1024.0 * 1024.0),
            (((double) n) * ((double) n) * ((double) sizeof(double))) / (1024.0 * 1024.0) );

    printf( "  |A|_F = %.4e\n", hlib_matrix_norm_frobenius( A, & info ) );
    printf( "  |A|_2 = %.4e\n", hlib_matrix_norm_spectral(  A, & info ) );
        
    if ( n < 20 )
    {
        hlib_real_t *  D   = (hlib_real_t *) malloc( n * n * sizeof(hlib_real_t) );
        int *          idx = (int *)         malloc( n * sizeof(int) );
        int            j;
        hlib_matrix_t  A_D;
        
        for ( j = 0; j < n; j++ )
            idx[j] = j;

        printf( "building exact BEM matrix ...\n" );
        
        start = hlib_walltime();
        kernel( n, idx, n, idx, D, & h );
        
        A_D = hlib_matrix_import_dense( n, n, D, 0, & info );              CHECK_INFO;

        printf( "  done in %.1f seconds\n", hlib_walltime() - start );
        printf( "  |D - A| = %.4e\n", hlib_matrix_norm_spectral_diff( A_D, A, & info ) );

        hlib_matrix_free( A_D, & info );
        free( idx );
        free( D );
    }
    
    if ( n < 20 )
    {
        hlib_matrix_t  B;
        hlib_acc_t     acc2 = acc;

        acc2.fixed_eps.eps = 1e-12;

        printf( "building BEM matrix using SVD ...\n" );
        start = hlib_walltime();
        
        B = hlib_matrix_build_coeff( bct, kernel, & h, HLIB_LRAPX_SVD,
                                     acc2, 1, & info );                         CHECK_INFO;
        
        printf( "  done in %.1f seconds\n", hlib_walltime() - start );
        hlib_matrix_print_ps( B, "bem1d_A_svd", HLIB_MATIO_SVD, & info );    CHECK_INFO;

        printf( "  |A-A~|_F/|A|_F = %.4e\n", hlib_matrix_norm_spectral_diff( A, B, & info ) );

        hlib_matrix_free( B, & info );
    }
    
    x = hlib_vector_build( n, & info );                                         CHECK_INFO;
    b = hlib_vector_build( n, & info );                                         CHECK_INFO;

    /*
     * set up rhs
     */

    b_arr = (hlib_real_t *) malloc( n * sizeof(hlib_real_t) );

    for ( i = 0; i < n; i++ )
        b_arr[i] = rhs( i, n );

    hlib_vector_import( b, b_arr, & info );                                CHECK_INFO;

    /* bring into H-ordering */
    hlib_vector_permute( b, hlib_clt_perm_e2i( ct, & info ), & info );     CHECK_INFO;
    
    printf( "  solving system ...\n" );

    solver = hlib_solver_auto( & info );                                   CHECK_INFO;
    
    start = hlib_walltime();
    hlib_solver_solve( solver, (hlib_linearoperator_t) A, x, b, NULL, & solve_info, & info );      CHECK_INFO;
    if ( solve_info.converged )
        printf( "  converged in %.1f seconds and %u steps with rate %.2e, |r| = %.2e\n",
                hlib_walltime() - start, solve_info.steps, solve_info.conv_rate,
                solve_info.res_norm );
    else
        printf( "  not converged in %.1f seconds and %u steps\n",
                hlib_walltime() - start, solve_info.steps );
    
    /* bring into external ordering to compare with exact solution */
    hlib_vector_permute( x, hlib_clt_perm_i2e( ct, & info ), & info );      CHECK_INFO;

    /* compute |x-1|_2 */
    one = hlib_vector_copy( x, & info );                                    CHECK_INFO;
    hlib_vector_fill( one, 1.0, & info );                                   CHECK_INFO;

    hlib_vector_axpy( 1.0, one, x, & info );                                CHECK_INFO;
    error = hlib_vector_norm_inf( x, & info );                              CHECK_INFO;
    
    printf( "  error of solution = %.4e\n", error );
        
    /*
     * use LU decomposition to solve system
     */

    if ( 1 )
    {
        hlib_matrix_t          B;
        hlib_linearoperator_t  LU;
        
        printf( "\n## LU decomposition\n" );

        printf( "  LU factorising H-matrix ... \n" );
        start = hlib_walltime();
        
        B  = hlib_matrix_copy( A, & info );                                 CHECK_INFO;
        LU = hlib_matrix_factorise_inv( B, acc, & info );                   CHECK_INFO;

        printf( "    done in %.1f seconds\n", hlib_walltime() - start );

        hlib_matrix_print_ps( B, "bem1d_LU", HLIB_MATIO_SVD, & info );  CHECK_INFO;

        printf( "  size of LU factor = %.2f MB\n",
                ((double) hlib_matrix_bytesize( B, & info )) / (1024.0 *1024.0) );

        printf( "  inversion error = %.4e\n",
                hlib_linearoperator_norm_inv_approx( (hlib_linearoperator_t) A, LU, & info ) );     CHECK_INFO;
        
        printf( "  solving preconditioned system\n" );
        start = hlib_walltime();
        hlib_solver_solve( solver, (hlib_linearoperator_t) A, x, b, LU, & solve_info, & info );     CHECK_INFO;
        if ( solve_info.converged )
            printf( "  converged in %.1f seconds and %u steps with rate %.2e, |r| = %.2e\n",
                    hlib_walltime() - start, solve_info.steps, solve_info.conv_rate,
                    solve_info.res_norm );
        else
            printf( "  not converged in %.1f seconds and %u steps\n",
                    hlib_walltime() - start, solve_info.steps );

        /* bring into external ordering to compare with exact solution */
        hlib_vector_permute( x, hlib_clt_perm_i2e( ct, & info ), & info );  CHECK_INFO;
        
        hlib_vector_axpy( 1.0, one, x, & info );                            CHECK_INFO;
        error = hlib_vector_norm_inf( x, & info );                          CHECK_INFO;
        
        printf( "  error of solution = %.4e\n", error );
        
        hlib_linearoperator_free( LU, & info );                             CHECK_INFO;
        hlib_matrix_free( B, & info );                                      CHECK_INFO;
    }

    if ( 1 )
    {
        hlib_matrix_t          B;
        hlib_linearoperator_t  WZ;
        
        printf( "\n## WAZ decomposition\n" );

        printf( "  WAZ factorising H-matrix ... \n" );
        start = hlib_walltime();
        
        B  = hlib_matrix_copy( A, & info );                                 CHECK_INFO;
        WZ = hlib_matrix_waz( B, acc, & info );                             CHECK_INFO;

        printf( "    done in %.1f seconds\n", hlib_walltime() - start );

        printf( "  inversion error = %.4e\n",
                hlib_linearoperator_norm_inv_approx( (hlib_linearoperator_t) A, WZ, & info ) );     CHECK_INFO;
        
        printf( "  solving preconditioned system\n" );
        start = hlib_walltime();
        hlib_solver_solve( solver, (hlib_linearoperator_t) A, x, b, WZ, & solve_info, & info );     CHECK_INFO;
        if ( solve_info.converged )
            printf( "  converged in %.1f seconds and %u steps with rate %.2e, |r| = %.2e\n",
                    hlib_walltime() - start, solve_info.steps, solve_info.conv_rate,
                    solve_info.res_norm );
        else
            printf( "  not converged in %.1f seconds and %u steps\n",
                    hlib_walltime() - start, solve_info.steps );

        /* bring into external ordering to compare with exact solution */
        hlib_vector_permute( x, hlib_clt_perm_i2e( ct, & info ), & info );  CHECK_INFO;
        
        hlib_vector_axpy( 1.0, one, x, & info );                            CHECK_INFO;
        error = hlib_vector_norm_inf( x, & info );                          CHECK_INFO;
        
        printf( "  error of solution = %.4e\n", error );
        
        hlib_linearoperator_free( WZ, & info );                             CHECK_INFO;
    }

    /*
     * solve system with inversion
     */
    
    if ( 1 )
    {
        hlib_matrix_t   IA;
        
        printf( "\n## Inversion\n" );

        printf( "  inverting H-matrix ... \n" );
        start = hlib_walltime();

        IA = hlib_matrix_copy( A, & info );                                 CHECK_INFO;
        hlib_matrix_inv( IA, acc, & info );                                 CHECK_INFO;
        
        printf( "    done in %.1f seconds\n", hlib_walltime() - start );
        
        hlib_matrix_print_ps( IA, "bem1d_I", HLIB_MATIO_SVD, & info );   CHECK_INFO;

        printf( "  size of Inverse = %.2f MB\n",
                ((double) hlib_matrix_bytesize( IA, & info )) / (1024.0 *1024.0) );
        
        printf( "  inversion error = %.4e\n",
                hlib_matrix_norm_inv_approx( A, IA, & info ) );             CHECK_INFO;
        
        printf( "  solving preconditioned system\n" );
        start = hlib_walltime();
        hlib_solver_solve( solver, (hlib_linearoperator_t) A, x, b,
                           (hlib_linearoperator_t) IA, & solve_info, & info ); CHECK_INFO;
        if ( solve_info.converged )
            printf( "  converged in %.1f seconds and %u steps with rate %.2e, |r| = %.2e\n",
                    hlib_walltime() - start, solve_info.steps, solve_info.conv_rate,
                    solve_info.res_norm );
        else
            printf( "  not converged in %.1f seconds and %u steps\n",
                    hlib_walltime() - start, solve_info.steps );

        /* bring into external ordering to compare with exact solution */
        hlib_vector_permute( x, hlib_clt_perm_i2e( ct, & info ), & info );  CHECK_INFO;
        
        hlib_vector_axpy( 1.0, one, x, & info );                            CHECK_INFO;
        error = hlib_vector_norm_inf( x, & info );                          CHECK_INFO;
        
        printf( "  error of solution = %.4e\n", error );
        
        hlib_matrix_free( IA, & info );                                     CHECK_INFO;
    }
    
    hlib_vector_free( x, & info );                                          CHECK_INFO;
    hlib_vector_free( b, & info );                                          CHECK_INFO;
    hlib_vector_free( one, & info );                                        CHECK_INFO;
    hlib_solver_free( solver, & info );                                     CHECK_INFO;
    hlib_matrix_free( A, & info );                                          CHECK_INFO;
    hlib_bct_free( bct, & info );                                           CHECK_INFO;
    hlib_admcond_free( adm, & info );                                       CHECK_INFO;
    hlib_clt_free( ct, & info );                                            CHECK_INFO;
    hlib_coord_free( coord, & info );                                       CHECK_INFO;

    free( b_arr );
    
    for ( i = 0; i < n; i++ )
        free( vertices[i] );
    free( vertices );

    hlib_done( & info );                                                    CHECK_INFO;
    
    return 0;
}
