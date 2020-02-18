//
// Project     : HLib
// File        : sparsealg.cc
// Description : Acoustic Scattering example
// Author      : Ronald Kriemann
// Copyright   : Max Planck Institute MIS 2004-2014. All Rights Reserved.
//

#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <string>

#include <boost/format.hpp>

#include "hlib.hh"

using namespace std;
using namespace HLIB;

using boost::format;

using  real_t    = HLIB::real;
using  complex_t = HLIB::complex;
using  fnspace_t = HLIB::TConstFnSpace;

/////////////////////////////////////////////////////////////////////////
//
// Acoustic Scattering
//
// Two spheres are placed next to each other with positive distance.
// An incident wave is reflected from first to second sphere and back.
// If the reflected wave is not changing anymore, the iteration stops.
//
// Limitations: assumes grid with sphere centered at (0,0,0) and radius 1
//
/////////////////////////////////////////////////////////////////////////

namespace
{

//
// function for right-hand side
//
class TMyRHS : public TBEMFunction< complex_t >
{
private:
    const T3Point    _inc_dir;  // incident direction
    const complex_t  _ikappa;
    
public:
    TMyRHS ( const T3Point &  inc_dir,
             const complex_t  kappa )
            : _inc_dir( inc_dir )
            , _ikappa( complex_t( 0, 1 ) * kappa )
    {}
    
    virtual complex_t
    eval ( const T3Point &  pos,
           const T3Point &  normal ) const
    {
        // 2 i k exp( i k <d,x> ) <d,n>
        return  real_t(2) * _ikappa * exp( _ikappa * dot( _inc_dir, pos ) ) * dot( _inc_dir, normal );
    }
};

//
// assuming sphere with center at (0,0,0) trivially add vertex normals
//
void
add_vtx_normal ( TGrid *  grid )
{
    size_t                  nvtx = grid->n_vertices();
    std::vector< T3Point >  normals( nvtx );

    for ( size_t i = 0; i < nvtx; ++i )
    {
        T3Point  n( grid->vertex( i ) );

        n.normalise2();
        normals[i] = n;
    }// for

    grid->add_vtx_normal( normals );
}

//
// build acoustic scattering matrix
//
unique_ptr< TMatrix >
build_acoustic ( const complex_t            kappa,
                 const TBlockClusterTree *  bct,
                 const fnspace_t *          ansatzsp,
                 const fnspace_t *          testsp,
                 const TTruncAcc &          acc,
                 TProgressBar *             progress )
{
    using  bf_t = TAcousticScatterBF< fnspace_t, fnspace_t >;

    bf_t                       bf( kappa, ansatzsp, testsp );
    TBFCoeffFn< bf_t >         bf_coeff_fn( & bf );
    TPermCoeffFn< complex_t >  coeff_fn( & bf_coeff_fn, bct->row_ct()->perm_i2e(), bct->col_ct()->perm_i2e() );
    
    // choose low-rank approximation
    TACAPlus< complex_t >        lrapx( & coeff_fn );
    TDenseMBuilder< complex_t >  h_builder( & coeff_fn, & lrapx );
    TWallTimer                   timer;
            
    timer.start();

    auto  K = h_builder.build( bct, acc, progress );
    
    timer.pause();
    cout << "    done in " << timer << endl
         << "    size of H-matrix = " << Mem::to_string( K->byte_size() ) << endl
         << "    |A|₂ / |A|_F     = "
         << format( "%.4e" ) % norm_2( K.get() )
         << " / "
         << format( "%.4e" ) % norm_F( K.get() )
         << endl;

    return K;
}

//
// build Mass matrix
//
unique_ptr< TMatrix >
build_mass ( const TBlockClusterTree *  bct,
             const fnspace_t *          ansatzsp,
             const fnspace_t *          testsp,
             const TTruncAcc &          acc,
             TProgressBar *             progress )
{
    using  bf_t = TMassBF< fnspace_t, fnspace_t >;

    bf_t                      bf( ansatzsp, testsp );
    TBFCoeffFn< bf_t >        bf_coeff_fn( & bf );
    TPermCoeffFn< real_t >    coeff_fn( & bf_coeff_fn, bct->row_ct()->perm_i2e(), bct->col_ct()->perm_i2e() );
    
    // choose low-rank approximation
    TACAPlus< real_t >        lrapx( & coeff_fn );
    TDenseMBuilder< real_t >  h_builder( & coeff_fn, & lrapx );
    TWallTimer                timer;
            
    timer.start();

    auto  M = h_builder.build( bct, acc, progress );
    
    M->set_complex( true );
    
    timer.pause();
    cout << "    done in " << timer << endl
         << "    size of H-matrix = " << Mem::to_string( M->byte_size() ) << endl
         << "    |A|₂ / |A|_F     = "
         << format( "%.4e" ) % norm_2( M.get() )
         << " / "
         << format( "%.4e" ) % norm_F( M.get() )
         << endl;

    return M;
}

//
// compare old and new solution (x/x_old) entrywise
//
real_t
compare ( const TVector * x,
          const TVector * x_old )
{
    real_t  vmax = 0.0;
    real_t  vmin = 1e100;
    size_t  n   = x->size();

    for ( idx_t  j = 0; j < idx_t(n); ++j )
    {
        const complex_t  v     = x->centry( j );
        const complex_t  v_old = x_old->centry( j );
        const real_t     f     = abs( v / v_old );
        
        if ( f > vmax ) vmax = f;
        if ( f < vmin ) vmin = f;
    }// for

    if ( verbose(4) )
        cout << "    x/x_old = "
             << format( "%.4f" ) % vmin
             << " / "
             << format( "%.4f" ) % vmax
             << endl;

    return vmax;
}

//
// print complex valued grid function
//
void
print_vec ( const TVector *   vec,
            const TGrid *     grid,
            const TFnSpace *  fnspace,
            const char *      vecname,
            const int         it = -1 )
{
    TVTKGridVis  gvis;

    gvis.set_colourmap( "coolwarm" );

    if ( false )
    {
        auto                vecre( vec->restrict_re() );
        std::ostringstream  filename;

        filename << "acoustic_" << vecname << "_re";
        if ( it >= 0 ) filename << "_" << format( "%02d" ) % it;

        gvis.print( grid, fnspace, vecre.get(), filename.str() );
    }

    if ( false )
    {
        auto                vecim( vec->restrict_im() );
        std::ostringstream  filename;

        filename << "acoustic_" << vecname << "_im";
        if ( it >= 0 ) filename << "_" << format( "%02d" ) % it;

        gvis.print( grid, fnspace, vecim.get(), filename.str() );
    }

    if ( true )
    {
        auto                vecabs = make_unique< TScalarVector >( vec->is() );
        std::ostringstream  filename;

        for ( idx_t  i = 0; i < idx_t(vecabs->size()); ++i )
            vecabs->set_entry( i, Math::abs( vec->centry( i ) ) );
        
        filename << "acoustic_" << vecname << "_abs";
        if ( it >= 0 ) filename << "_" << format( "%02d" ) % it;

        gvis.print( grid, fnspace, vecabs.get(), filename.str() );
    }
}

}// namespace

//
// main function
//
int
main ( int argc, char ** argv )
{
    //
    // read command line
    //

    const char *  gridname  = "sphere-5";
    double        eps       = 1e-4;
    double        fac_eps   = 1e-4;
    complex_t     kappa     = complex_t( 4, 0 );
    
    try
    {
        INIT();

        TConsoleProgressBar  progress;
        TWallTimer           timer;
        
        //
        // read grids
        // - only works for spheres, which are assumed to be centered at (0,0,0) with radius 1
        //
        
        cout << "━━ building grid " << endl;
        
        auto  grid1 = make_grid( gridname );

        grid1->scale( T3Point( 4, 4, 1 ) );
        
        if ( ! grid1->has_vtx_normal() )
            add_vtx_normal( grid1.get() );
        
        cout << "    no. of vertices   = " << grid1->n_vertices() << endl;
        cout << "    no. of triangles  = " << grid1->n_triangles() << endl;

        if ( verbose(3) )
        {
            TVTKGridVis  gvis;

            gvis.print( grid1.get(), "acoustic_grid1" );
        }// if
        
        // read second sphere and translate it: assuming sphere has center at (0,0,0) with radius 1
        auto  grid2 = make_grid( gridname );

        if ( ! grid2->has_vtx_normal() )
            add_vtx_normal( grid2.get() );

        // set up second sphere next to first one (distance: 0.5)
        grid2->scale( T3Point( 4, 4, 1 ) );
        grid2->translate( T3Point( 0, 0, 2.5 ) );
        
        if ( verbose(3) )
        {
            TVTKGridVis  gvis;

            gvis.print( grid2.get(), "acoustic_grid2" );
        }// if
        
        //
        // build ansatz and test space (both equal)
        //

        cout << endl << "━━ building ansatz/test spaces" << endl;
        
        fnspace_t                        fnspace1( grid1.get() );
        fnspace_t                        fnspace2( grid2.get() );

        cout << "    no. of dof = " << fnspace1.n_indices() << endl;

        //
        // and the cluster trees
        //

        auto               coord1 = fnspace1.build_coord();
        auto               coord2 = fnspace2.build_coord();

        TAutoBSPPartStrat  part_strat( adaptive_split_axis );
        TBSPCTBuilder      ct_builder( & part_strat );
        auto               ct1 = ct_builder.build( coord1.get() );
        auto               ct2 = ct_builder.build( coord2.get() );
        
        TStdGeomAdmCond    adm_cond( 2.0 );
        TBCBuilder         bct_builder;
        auto               bct11 = bct_builder.build( ct1.get(), ct1.get(), & adm_cond );
        auto               bct12 = bct_builder.build( ct1.get(), ct2.get(), & adm_cond );
        auto               bct21 = bct_builder.build( ct2.get(), ct1.get(), & adm_cond );

        if ( verbose( 3 ) )
        {
            TPSBlockClusterVis  bc_vis;

            bc_vis.print( bct11->root(), "acoustic_bct11" );
            bc_vis.print( bct12->root(), "acoustic_bct12" );
            bc_vis.print( bct21->root(), "acoustic_bct21" );
        }// if

        //
        // build matrix
        //
        
        cout << endl << "━━ building Acoustic H-matrix ("
             << " κ = " << kappa
             << ", ε = " << format( "%.2e" ) % eps
             << " )" << endl;

        auto  acc = fixed_prec( eps );
            
        cout << endl << "  ∙ building K11 " << endl;
        auto  K11 = build_acoustic( kappa, bct11.get(), & fnspace1, & fnspace1, acc, & progress );

        cout << endl << "  ∙ building M11 " << endl;
        auto  M11 = build_mass( bct11.get(), & fnspace1, & fnspace1, acc, & progress );

        cout << endl << "  ∙ computing M11 - K11 " << endl;
        add( 1.0, M11.get(), -1.0, K11.get(), acc );

        cout << endl << "  ∙ building K12 " << endl;
        auto  K12 = build_acoustic( kappa, bct12.get(), & fnspace1, & fnspace2, acc, & progress );

        cout << endl << "  ∙ building K21 " << endl;
        auto  K21 = build_acoustic( kappa, bct21.get(), & fnspace2, & fnspace1, acc, & progress );

        if ( verbose( 3 ) )
        {
            TPSMatrixVis  mvis;

            mvis.print( K11.get(), "acoustic_K11" );
            mvis.print( K12.get(), "acoustic_K12" );
            mvis.print( K21.get(), "acoustic_K21" );
        }// if

        //
        // construct preconditioner
        //

        cout << endl << "━━ LU factorisation ( ε = " << format( "%.2e" ) % fac_eps << " )" << endl;

        auto  acc_lu = fixed_prec( fac_eps );

        timer.start();
        
        auto  K11_copy = K11->copy();
        auto  K11_pre  = factorise_inv( K11_copy.get(), acc_lu, & progress );

        timer.pause();
        cout << "    done in " << timer << endl;
        
        cout << "    size of LU factor = " << Mem::to_string( K11_copy->byte_size() ) << endl;
        cout << "    inversion error   = " << format( "%.6e" ) % inv_approx_2( K11.get(), K11_pre.get() ) << endl;

        if( verbose( 3 ) )
        {
            TPSMatrixVis  mvis;
            
            mvis.print( K11_copy.get(), "acoustic_LU" );
        }// if

        //
        // build RHS
        //

        cout << endl << "━━ building RHS " << endl;
            
        TMyRHS                               rhsfn( T3Point( 1, 1, 0 ), kappa );
        TQuadBEMRHS< fnspace_t, complex_t >  rhs_build( 4 );
        auto                                 rhs( rhs_build.build( & fnspace1, & rhsfn ) );
        
        if ( verbose( 2 ) )
            print_vec( rhs.get(), grid1.get(), & fnspace1, "rhs" );

        // bring RHS into H-ordering
        ct2->perm_e2i()->permute( rhs.get() );
        
        //
        // acoustic scattering
        //

        auto         x         = K11->col_vector();
        auto         x_old     = K11->col_vector();
        TSolverInfo  solve_info( false, verbose( 4 ) );
        real_t       old_ratio = real_t(1);
        
        cout << endl << "━━ Acoustic Scattering Iteration " << endl;

        timer.start();
        
        for ( int it = 0; it < 100; ++it )
        {
            cout << "  ∙ step " << format( "%2d" ) % it;

            //
            // solve on first sphere
            //
            
            solve( K11.get(), x.get(), rhs.get(), K11_pre.get(), & solve_info );

            {
                auto  t = x->copy();

                ct1->perm_i2e()->permute( t.get() );
                print_vec( t.get(), grid1.get(), & fnspace1, "sol1", it );
            }

            //
            // compare old and new solution and stop if converged or zero
            //
            
            if ( it > 0 )
            {
                const real_t  ratio = compare( x.get(), x_old.get() );
                const real_t  diff  = ( it > 1 ? std::abs( ratio - old_ratio ) / std::abs( old_ratio ) : 0.0 );

                old_ratio = ratio;

                cout << ", |x_i/x_i-1| = " << format( "%.4e" ) % ratio;
                if ( it > 1 ) cout << ", diff = " << format( "%.4e" ) % diff;
                cout << endl;

                // stop condition
                if (( it > 1 ) && ( diff < 1e-4 ))
                    break;
            }// if
            else
                cout << endl;

            const real_t  normx = x->norm2();

            if ( normx < 1e-20 )
            {
                cout << "    |x| = 0"<< endl;
                break;
            }// if
            
            x->scale( real_t(1) / normx );
            x_old->assign( real_t(1), x.get() );

            //
            // compute right-hand side for second sphere: rhs = - K12 x
            //

            mul_vec( real_t(-1), K12.get(), x.get(), real_t(0), rhs.get(), apply_normal );
        
            //
            // solve on second sphere (identical with first sphere)
            //

            solve( K11.get(), x.get(), rhs.get(), K11_pre.get(), & solve_info );

            {
                auto  t = x->copy();

                ct2->perm_i2e()->permute( t.get() );
                print_vec( t.get(), grid2.get(), & fnspace2, "sol2", it );
            }

            //
            // compute right-hand side for first sphere: rhs = - K21 x
            //

            mul_vec( real_t(-1), K21.get(), x.get(), real_t(0), rhs.get(), apply_normal );
        }// for

        timer.pause();
        cout << "  done in " << timer << endl;
        
        //
        // finish
        //
        
        DONE();
    }// try
    catch ( Error & e )
    {
        cout << e.to_string() << endl;
    }// catch
    
    return 0;
}
