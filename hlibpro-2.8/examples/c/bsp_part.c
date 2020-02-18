/*
 * Project     : HLib
 * File        : bsp_part.c
 * Description : example for 1d BEM problem
 * Author      : Ronald Kriemann
 * Copyright   : Max Planck Institute MIS 2004-2014. All Rights Reserved.
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <hlib-c.h>

/*
 * error handler function
 */
void
errorfn ( const int errcode, const char * errmsg )
{
    printf( errmsg );
    exit( errcode );
}

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
    const double  a     = ((double) i)       / ((double) n);
    const double  b     = ((double) i + 1.0) / ((double) n);
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
    hlib_real_t            * b_arr;
    hlib_coord_t             coord;
    hlib_vector_t            x, b, one;
    unsigned int *           partition;
    hlib_clustertree_t       ct;
    hlib_admcond_t           adm;
    hlib_blockclustertree_t  bct;
    hlib_matrix_t            A, B;
    hlib_linearoperator_t    LU;
    double                   error, start;
    hlib_solver_t            solver;
    hlib_solve_info_t        solve_info;
    hlib_acc_t               acc;
    hlib_acc_t               subacc[4];
    
    if ( argc > 1 ) n = atoi( argv[1] );
    else            n = 8;
    if ( n <= 0 )   n = 8;  /* sanity check */

    h = 1.0 / ((double) n);
    
    hlib_init( & info );

    hlib_set_error_fn( errorfn );
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

    /*
     * split indexset into two subsets 
     */
    
    partition = (unsigned int *) malloc( sizeof(int) * n );

    for ( i = 0; i < 7*n/8; i++ ) partition[i] = 0;
    for (      ; i < n;     i++ ) partition[i] = 1;

    /*
     * build cluster trees accoridng to coordinates and
     * predefined splitting
     */
    
    printf( "building cluster tree\n" );
    coord = hlib_coord_import( n, 1, vertices, NULL, NULL );
    ct    = hlib_clt_build_bsp_part( coord, partition, 1, HLIB_BSP_AUTO, 20, NULL );
    hlib_clt_print_ps( ct, "bsp_part_ct", NULL );
    
    printf( "building block cluster tree\n" );
    adm = hlib_admcond_geom( HLIB_ADM_STD_MAX, 2.0, NULL );
    bct = hlib_bct_build( ct, ct, adm, NULL );
    hlib_bct_print_ps( bct, "bsp_part_bct", NULL );

    /*
     * build matrix with block wise accuracy
     */

    subacc[0] = hlib_acc_fixed_eps( 1e-6 );
    subacc[1] = hlib_acc_fixed_eps( 1e-4 );
    subacc[2] = hlib_acc_fixed_eps( 1e-0 );
    subacc[3] = hlib_acc_fixed_eps( 1e-6 );
    acc       = hlib_acc_blocked( subacc );
    
    printf( "building BEM matrix using ACA+ ...\n" );
    start = hlib_walltime();

    A = hlib_matrix_build_coeff( bct, kernel, & h, HLIB_LRAPX_ACAPLUS,
                                 acc, 1, NULL );
    
    printf( "  done in %.1f seconds\n", hlib_walltime() - start );
    hlib_matrix_print_ps( A, "bsp_part_A", HLIB_MATIO_SVD, NULL );

    bytesize = hlib_matrix_bytesize( A, NULL );
    
    printf( "  memory per dof    = %.2f kb\n", ((double) bytesize) / ( 1024.0 * n ) );
    printf( "  compression ratio = %.2f%%\n",
            100.0 * ((double) bytesize) / (((double) n) * ((double) n) * ((double) sizeof(double))) );

    printf( "  |A|_F = %.4e\n", hlib_matrix_norm_frobenius( A, NULL ) );
    printf( "  |A|_2 = %.4e\n", hlib_matrix_norm_spectral(  A, NULL ) );
        
    /*
     * use LU decomposition to solve system
     */

    printf( "  LU factorising H-matrix ... \n" );
    start = hlib_cputime();
        
    B = hlib_matrix_copy_acc( A, acc, 0, NULL );
    hlib_matrix_print_ps( B, "bsp_part_B", HLIB_MATIO_SVD, NULL );
    
    LU = hlib_matrix_factorise_inv( B, acc, NULL );

    printf( "    done in %.1f seconds\n", hlib_cputime() - start );

    hlib_matrix_print_ps( B, "bsp_part_LU", HLIB_MATIO_SVD, NULL );

    printf( "    memory per dof  = %.2f kb\n",
            ((double) hlib_matrix_bytesize( B, NULL )) / ( 1024.0 * n ) );

    printf( "    inversion error = %.4e\n",
            hlib_linearoperator_norm_inv_approx( (hlib_linearoperator_t) A, LU, NULL ) );
        
    x = hlib_vector_build( n, NULL );
    b = hlib_vector_build( n, NULL );

    /*
     * set up rhs
     */

    b_arr = (hlib_real_t *) malloc( n * sizeof(hlib_real_t) );
    
    for ( i = 0; i < n; i++ )
        b_arr[i] = rhs( i, n );
    
    hlib_vector_import( b, b_arr, NULL );

    /* bring into H-ordering */
    hlib_vector_permute( b, hlib_clt_perm_e2i( ct, NULL ), NULL );
    
    printf( "  solving preconditioned system\n" );
    solver = hlib_solver_auto( NULL );
    start = hlib_cputime();
    hlib_solver_solve( solver, (hlib_linearoperator_t) A, x, b, LU, & solve_info, NULL );
    if ( solve_info.converged )
        printf( "    converged in %.1f seconds and %u steps with rate %.2e, |r| = %.2e\n",
                hlib_cputime() - start, solve_info.steps, solve_info.conv_rate,
                solve_info.res_norm );
    else
        printf( "    not converged in %.1f seconds and %u steps\n",
                hlib_cputime() - start, solve_info.steps );

    /* compute |x-1|_2 */
    one = hlib_vector_copy( x, NULL );
    hlib_vector_fill( one, 1.0, NULL );

    hlib_vector_axpy( 1.0, one, x, NULL );
    error = hlib_vector_norm_inf( x, NULL );
        
    printf( "  error of solution = %.4e\n", error );
        
    hlib_linearoperator_free( LU, NULL );
    hlib_matrix_free( B, NULL );
    
    hlib_solver_free( solver, NULL );
    
    hlib_vector_free( x, NULL );
    hlib_vector_free( b, NULL );
    hlib_vector_free( one, NULL );
    hlib_matrix_free( A, NULL );
    hlib_bct_free( bct, NULL );
    hlib_admcond_free( adm, NULL );
    hlib_clt_free( ct, NULL );
    hlib_coord_free( coord, NULL );

    free( b_arr );
    
    for ( i = 0; i < n; i++ )
        free( vertices[i] );
    free( vertices );

    hlib_done( NULL );
    
    return 0;
}
