
#include <iostream>

#include <boost/format.hpp>
#include <boost/program_options.hpp>

#include "hlib.hh"

using namespace HLIB;
using namespace boost::program_options;

using std::string;

using real_t    = HLIB::real;
using complex_t = HLIB::complex;

int
main ( int     argc,
       char ** argv )
{
    //
    // define command line options
    //

    real_t   eps       = 1e-4;
    string   method    = "lu";
    string   matrix    = "A.hm";
    int      nthreads  = 0;
    int      verbosity = 1;
    
    options_description             all_opts;
    options_description             vis_opts( "usage: factorize [options] [matrix]\n  where options include" );
    options_description             hid_opts( "Hidden options" );
    positional_options_description  pos_opts;
    variables_map                   vm;

    // standard options
    vis_opts.add_options()
        ( "eps,e",       value<double>(), ": set H-algebra accuracy ε (default = 1e-4)" )
        ( "help,h",                       ": print this help text" )
        ( "method,m",    value<string>(), ": inversion method:\n"
                                          "      lu    : \tLU factorisation\n"
                                          "      waz   : \tLU+triangular inversion (default)\n"
                                          "      inv   : \tLU+tri.inv.+U·L\n"
                                          "      gauss : \tGaussian elimination"
          )
        ( "threads,t",   value<int>(),    ": number of parallel threads" )
        ( "verbosity,v", value<int>(),    ": verbosity level" )
        ;
        
    hid_opts.add_options()
        ( "matrix",      value<string>()->default_value("A.hm"), ": filename containing H-matrix (default = A.hm)" )
        ;

    // options for command line parsing
    all_opts.add( vis_opts ).add( hid_opts );

    // all "non-option" arguments should be "--matrix" arguments
    pos_opts.add( "matrix", -1 );

    //
    // parse command line options
    //

    try
    {
        store( command_line_parser( argc, argv ).options( all_opts ).positional( pos_opts ).run(), vm );
        notify( vm );
    }// try
    catch ( required_option &  e )
    {
        std::cout << e.get_option_name() << " requires an argument, try \"-h\"" << std::endl;
        exit( 1 );
    }// catch
    catch ( unknown_option &  e )
    {
        std::cout << e.what() << ", try \"-h\"" << std::endl;
        exit( 1 );
    }// catch

    //
    // eval command line options
    //

    if ( vm.count( "help") )
    {
        std::cout << vis_opts << std::endl;
        exit( 1 );
    }// if

    if ( vm.count( "eps"       ) ) eps       = vm["eps"].as<double>();
    if ( vm.count( "threads"   ) ) nthreads  = vm["threads"].as<int>();
    if ( vm.count( "verbosity" ) ) verbosity = vm["verbosity"].as<int>();
    if ( vm.count( "method"    ) ) method    = vm["method"].as<string>();

    if ( vm.count( "matrix" ) )
        matrix = vm["matrix"].as<string>();
    else
    {
        std::cout << "usage: factorize [options] [matrix]" << std::endl;
        exit( 1 );
    }// if

    //
    // factorisation
    //
    
    INIT();

    if ( nthreads != 0 )
        CFG::set_nthreads( nthreads );

    CFG::set_verbosity( verbosity );

    // activating accumulator arithmetic and RRQR approximation/truncation
    CFG::Arith::use_accu     = true;
    CFG::BLAS::approx_method = use_rrqr;
    CFG::BLAS::trunc_method  = use_rrqr;
    
    std::cout << "┎────────────────────────────────────────────────────────" << std::endl
              << "┃ reading matrix from " << matrix                          << std::endl
              << "┖────────────────────────────────────────────────────────" << std::endl;

    auto          A   = read_matrix( matrix.c_str() );
    auto          acc = fixed_prec( eps );
    TPSMatrixVis  mvis;

    std::cout << "  #rows / #cols  = " << A->rows() << " / " << A->cols() << std::endl
              << "  symmetric      = " << ( A->is_symmetric() || A->is_hermitian() ? "yes" : "no"  ) << std::endl
              << "  complex valued = " << ( A->is_complex() ? "yes" : "no" ) << std::endl
              << std::endl;
    
    mvis.svd( false );
    mvis.structure( true );
    
    if ( verbose( 3 ) )
        mvis.print( A.get(), "A" );
    
    if ( method == "lu" )    
    {
        std::cout << "┎────────────────────────────────────────────────────────" << std::endl
                  << "┃ LU factorisation"                                        << std::endl
                  << "┖────────────────────────────────────────────────────────" << std::endl;
        
        auto  A_copy = A->copy();
        auto  tic    = Time::Wall::now();

        if ( A->is_unsymmetric() )
            lu( A_copy.get(), acc );
        else
            ldl( A_copy.get(), acc );

        auto  toc = Time::Wall::since( tic );

        std::unique_ptr< TLinearOperator >  A_inv;

        if ( A->is_unsymmetric() )
            A_inv.reset( new TLUInvMatrix( A_copy.get() ) );
        else
            A_inv.reset( new TLDLInvMatrix( A_copy.get(), A->form() ) );

        std::cout << "lu()" << std::endl
                  << "  done in     " << toc << std::endl;

        std::cout 
                  << "  inv error = " << boost::format( "%.4e" ) % inv_approx_2( A.get(), A_inv.get() ) << std::endl
                  << "  memory    = " << Mem::to_string( A_copy->byte_size() ) << std::endl;

        if ( verbose( 3 ) )
            mvis.print( A_copy.get(), "LU" );

        auto  x = A_inv->domain_vector();
        auto  y = A_inv->range_vector();

        x->fill( 1 );
        y->fill( 0 );

        tic = Time::Wall::now();

        for ( uint i = 0; i < 100; ++i )
        {
            A_inv->apply( x.get(), y.get(), apply_normal );
        }// for

        toc = Time::Wall::since( tic );
        std::cout << "  100x MVM in " << toc << std::endl;
    }// if
    else if ( method == "waz" )
    {
        std::cout << "┎────────────────────────────────────────────────────────" << std::endl
                  << "┃ WAZ"                                                     << std::endl
                  << "┖────────────────────────────────────────────────────────" << std::endl;
        
        auto  A_copy = A->copy();
        auto  tic    = Time::Wall::now();
        auto  WZ     = WAZ::factorise( A_copy.release(), acc );
        auto  toc    = Time::Wall::since( tic );
        auto  W      = std::unique_ptr< TMatrix >( WZ.first );
        auto  Z      = std::unique_ptr< TMatrix >( WZ.second );
        auto  A_inv  = std::unique_ptr< TLinearOperator >( A->is_unsymmetric()
                                                           ? ptrcast( matrix_product( Z.get(), W.get() ).release(),
                                                                      TLinearOperator )
                                                           : ptrcast( matrix_product( 1.0, apply_adjoint, W.get(),
                                                                                      1.0, apply_normal,  Z.get(),
                                                                                      1.0, apply_normal,  W.get() ).release(),
                                                                      TLinearOperator ) );

        std::cout << "factored_inverse()" << std::endl
                  << "  done in     " << toc << std::endl
                  << "  inv error = " << boost::format( "%.4e" ) % inv_approx_2( A.get(), A_inv.get() ) << std::endl
                  << "  memory    = " << Mem::to_string( W->byte_size() + Z->byte_size() ) << std::endl;

        if ( verbose( 3 ) )
        {
            mvis.print( W.get(), "W" );
            mvis.print( Z.get(), "Z" );
        }// if

        auto  x = A_inv->domain_vector();
        auto  y = A_inv->range_vector();

        x->fill( 1 );
        y->fill( 0 );
        
        tic = Time::Wall::now();

        for ( uint i = 0; i < 100; ++i )
        {
            A_inv->apply( x.get(), y.get(), apply_normal );
        }// for
        
        toc = Time::Wall::since( tic );
        std::cout << "  100x MVM in " << toc << std::endl;
    }// if
    else if ( method == "inv" )
    {
        std::cout << "┎────────────────────────────────────────────────────────" << std::endl
                  << "┃ Inversion"                                               << std::endl
                  << "┖────────────────────────────────────────────────────────" << std::endl;
        
        auto  A_inv = A->copy();
        auto  tic   = Time::Wall::now();

        invert( A_inv.get(), acc );

        auto  toc = Time::Wall::since( tic );
    
        std::cout << "invert()" << std::endl
                  << "  done in     " << toc << std::endl
                  << "  inv error = " << boost::format( "%.4e" ) % inv_approx_2( A.get(), A_inv.get() ) << std::endl
                  << "  memory    = " << Mem::to_string( A_inv->byte_size() ) << std::endl;

        if ( verbose( 3 ) )
            mvis.print( A_inv.get(), "A_inv" );

        auto  x = A_inv->domain_vector();
        auto  y = A_inv->range_vector();

        x->fill( 1 );
        y->fill( 0 );
        
        tic = Time::Wall::now();

        for ( uint i = 0; i < 100; ++i )
        {
            A_inv->apply( x.get(), y.get(), apply_normal );
        }// for
        
        toc = Time::Wall::since( tic );
        std::cout << "  100x MVM in " << toc << std::endl;
    }// if
    else if ( method == "gauss" )
    {
        std::cout << "┎────────────────────────────────────────────────────────" << std::endl
                  << "┃ Gaussian Elimination"                                    << std::endl
                  << "┖────────────────────────────────────────────────────────" << std::endl;
        
        auto  A_copy = A->copy();
        auto  A_inv  = A->copy();
        auto  tic    = Time::Wall::now();

        gauss_elim( A_copy.get(), A_inv.get(), acc );

        auto  toc = Time::Wall::since( tic );
    
        std::cout << "invert()" << std::endl
                  << "  done in     " << toc << std::endl
                  << "  inv error = " << boost::format( "%.4e" ) % inv_approx_2( A.get(), A_inv.get() ) << std::endl
                  << "  memory    = " << Mem::to_string( A_inv->byte_size() ) << std::endl;

        if ( verbose( 3 ) )
            mvis.print( A_inv.get(), "A_inv" );

        auto  x = A_inv->domain_vector();
        auto  y = A_inv->range_vector();

        x->fill( 1 );
        y->fill( 0 );
        
        tic = Time::Wall::now();

        for ( uint i = 0; i < 100; ++i )
        {
            A_inv->apply( x.get(), y.get(), apply_normal );
        }// for
        
        toc = Time::Wall::since( tic );
        std::cout << "  100x MVM in " << toc << std::endl;
    }// if
    else
        std::cout << "unknown method, try -h" << std::endl;
    
    DONE();
}
