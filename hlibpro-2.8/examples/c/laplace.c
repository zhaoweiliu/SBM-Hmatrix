/*
 * Project     : HLib
 * File        : laplace.c
 * Description : example for a 3d BEM problem
 * Author      : Ronald Kriemann
 * Copyright   : Max Planck Institute MIS 2004-2014. All Rights Reserved.
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include <hlib-c-bem.h>

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
 * main function
 */
int
main ( int argc, char ** argv )
{
    int                      n;
    char                     gridfile[1024];
    hlib_grid_t              grid;
    hlib_fnspace_t           fnspace;
    hlib_coord_t             coord;
    hlib_clustertree_t       ct;
    hlib_admcond_t           adm;
    hlib_blockclustertree_t  bct;
    hlib_matrix_t            M = NULL;
    double                   start, error, norm;
    const double             view[3] = { 0.25, -0.35, 1 };
    unsigned long            bytesize;
    hlib_acc_t               acc = hlib_acc_fixed_eps( 1e-4 );

    if ( argc > 2 )
        acc = hlib_acc_fixed_eps( atof( argv[2] ) );
    if ( argc > 1 )
        strcpy( gridfile, argv[1] );
    else
    {
        printf( "usage: laplace <gridname> [<accuracy>]\n" );
        exit( 1 );
    }
    
    /*
     * init HLIBpro
     */
    
    hlib_init( NULL );

    hlib_set_error_fn( errorfn );
    hlib_set_verbosity( 2 );
    hlib_set_coarsening( 1, 0 );
    
    /*
     * read grid and build functions space
     */

    printf( "reading grid\n" );
    
    grid = hlib_load_grid( gridfile, NULL );

    hlib_grid_print_ps( grid, view, "laplace_grid", NULL );
    hlib_grid_print_vrml( grid, "laplace_grid", NULL );

    printf( "  grid uses %.2f MB of memory\n",
            ((double) hlib_grid_bytesize( grid, NULL )) / (1024.0*1024.0) );

    fnspace = hlib_fnspace_build_linear( grid, NULL );
    n       = hlib_fnspace_dim( fnspace, NULL );
    
    printf( "building cluster tree\n" );
    coord = hlib_fnspace_coord( fnspace, NULL );
    ct    = hlib_clt_build_bsp( coord, HLIB_BSP_AUTO, 20, NULL );
    hlib_clt_print_ps( ct, "laplace_ct", NULL );
    
    printf( "building block cluster tree\n" );
    adm = hlib_admcond_geom( HLIB_ADM_STD_MIN, 2.0, NULL );
    bct = hlib_bct_build( ct, ct, adm, NULL );
    hlib_bct_print_ps( bct, "laplace_bct", NULL );

    /*
     * build mass matrix
     */

    {
        hlib_bemrbf_t  massbf = hlib_bembf_mass( fnspace, fnspace, NULL );
        
        printf( "building mass matrix\n" );
    
        start = hlib_walltime();

        M = hlib_matrix_build_coeff( bct, hlib_eval_bemrbf_cb, massbf, HLIB_LRAPX_ACAPLUS,
                                     acc, 1, NULL );

        printf( "  done in %.2f sec\n", hlib_walltime() - start );
    
        printf( "  matrix has dimension %lu x %lu\n", hlib_matrix_rows( M, NULL ),
                hlib_matrix_cols( M, NULL ) );
    
        hlib_matrix_print_ps( M, "laplace_M", HLIB_MATIO_SVD, NULL );
    
        bytesize = hlib_matrix_bytesize( M, NULL );
    
        printf( "  compression ratio = %.2f%% (%.2f MB compared to %.2f MB)\n",
                100.0 * ((double) bytesize) /
                (((double) n) * ((double) n) * ((double) sizeof(double))),
                ((double) bytesize) / (1024.0 * 1024.0),
                (((double) n) * ((double) n) * ((double) sizeof(double))) / (1024.0 * 1024.0) );

        hlib_bemrbf_free( massbf, NULL );
    }

    /*
     * build DLP matrix
     */

    {
        hlib_bemrbf_t  dlpbf = hlib_bembf_laplace_dlp( fnspace, fnspace, NULL );
        hlib_matrix_t  K;
        
        printf( "building Laplace DLP matrix\n" );
    
        start = hlib_walltime();
        
        K = hlib_matrix_build_coeff( bct, hlib_eval_bemrbf_cb, dlpbf, HLIB_LRAPX_ACAPLUS,
                                     acc, 0, NULL );
        
        printf( "  done in %.2f sec\n", hlib_walltime() - start );
        
        hlib_matrix_print_ps( K, "laplace_K", HLIB_MATIO_SVD, NULL );
        
        bytesize = hlib_matrix_bytesize( K, NULL );
        
        printf( "  compression ratio = %.2f%% (%.2f MB compared to %.2f MB)\n",
                100.0 * ((double) bytesize) /
                (((double) n) * ((double) n) * ((double) sizeof(double))),
                ((double) bytesize) / (1024.0 * 1024.0),
                (((double) n) * ((double) n) * ((double) sizeof(double))) / (1024.0 * 1024.0) );
        printf( "  memory per dof    = %.2f kB\n", ((double) bytesize) / ( 1024.0 * ((double) n) ) );
        
        norm  = hlib_matrix_norm_spectral( K, NULL );
    
        printf( "  |K| = %.4e\n", norm );
        
        /*
         * check approximation by evaluating (K+½M) · 𝟏
         */

        {
            hlib_vector_t  x, y;
            
            x = hlib_matrix_col_vector( K, NULL );
            y = hlib_matrix_row_vector( K, NULL );
            
            hlib_vector_fill( x, 1.0, NULL );

            hlib_matrix_mulvec(  0.5, M, x, 0.0, y, HLIB_MATOP_NORM, NULL );
            hlib_matrix_mulvec( -1.0, K, x, 1.0, y, HLIB_MATOP_NORM, NULL );
            
            error = hlib_vector_norm2( y, NULL );
            norm  = hlib_matrix_norm_spectral( K, NULL );
            
            printf( "approximation error:\n" );
            printf( "  |(½M - K)*1|             = %.4e\n", error );
            printf( "  |(½M - K)*1| / (|K| |1|) = %.4e\n", error / (norm * sqrt(n)) );
            
            hlib_vector_free( x, NULL );
            hlib_vector_free( y, NULL );
        }
        
        hlib_matrix_free( K, NULL );
        hlib_bemrbf_free( dlpbf, NULL );
    }
    
    /*
     * finish HLIBpro
     */

    hlib_matrix_free( M, NULL );
    hlib_bct_free( bct, NULL );
    hlib_admcond_free( adm, NULL );
    hlib_clt_free( ct, NULL );
    hlib_coord_free( coord, NULL );
    hlib_fnspace_free( fnspace, NULL );
    hlib_grid_free( grid, NULL );
    hlib_done( NULL );
    
    return 0;
}
