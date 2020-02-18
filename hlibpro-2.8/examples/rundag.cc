
#include <iostream>

#include <boost/format.hpp>
#include <boost/program_options.hpp>

#include "hlib.hh"

using namespace HLIB;
using namespace HLIB::Time::Wall;
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
    options_description             vis_opts( "usage: rundag [options] [matrix]\n  where options include" );
    options_description             hid_opts( "Hidden options" );
    positional_options_description  pos_opts;
    variables_map                   vm;

    // standard options
    vis_opts.add_options()
        ( "eps,e",       value<double>(), ": set H-algebra accuracy ε (default = 1e-4)" )
        ( "help,h",                       ": print this help text" )
        ( "method,m",    value<string>(), ": inversion method:\n"
                                          "      lu    : \tLU factorisation\n"
                                          "      ldl   : \tLDL factorisation\n"
                                          "      waz   : \tLU+triangular inversion (default)\n"
                                          "      invll : \tlower left triangular inversion.\n"
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
        std::cout << "usage: rundag [options] [matrix]" << std::endl;
        exit( 1 );
    }// if

    //
    // factorisation
    //
    
    INIT();

    if ( nthreads != 0 )
        CFG::set_nthreads( nthreads );

    CFG::set_verbosity( verbosity );

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
    
    // init timer
    auto  tic = now();
    auto  toc = since( tic );
    
    if ( method == "lu" )    
    {
        std::cout << "┎────────────────────────────────────────────────────────" << std::endl
                  << "┃ LU factorisation"                                        << std::endl
                  << "┖────────────────────────────────────────────────────────" << std::endl;
        
        auto  A_copy = A->copy();

        tic    = Time::Wall::now();

        auto  dag = LU::gen_dag( A_copy.get(), fac_options_t() );

        toc = since( tic );
        
        std::cout << "  dag in      " << toc << std::endl;
        std::cout << "    #nodes  = " << dag.nnodes() << std::endl;
        std::cout << "    #edges  = " << dag.nedges() << std::endl;

        dag.print_dot( "lu.dot" );
        
        tic = now();

        dag.run( acc );

        toc = since( tic );
        
        auto  A_inv = std::make_unique< TLUInvMatrix >( A_copy.get() );

        std::cout << "  done in     " << toc << std::endl;
        std::cout << "  inv error = " << boost::format( "%.4e" ) % inv_approx_2( A.get(), A_inv.get() ) << std::endl
                  << "  memory    = " << Mem::to_string( A_copy->byte_size() ) << std::endl;
    }// if
    else if ( method == "ldl" )
    {
        std::cout << "┎────────────────────────────────────────────────────────" << std::endl
                  << "┃ LU factorisation"                                        << std::endl
                  << "┖────────────────────────────────────────────────────────" << std::endl;
        
        auto  A_copy = A->copy();
        auto  D      = copy_diag( A.get() );

        tic = now();

        auto  dag = LDL::gen_dag( A_copy.get(), D.get(), fac_options_t() );

        toc = since( tic );
        
        std::cout << "  dag in      " << toc << std::endl;
        std::cout << "    #nodes  = " << dag.nnodes() << std::endl;
        std::cout << "    #edges  = " << dag.nedges() << std::endl;

        dag.print_dot( "ldl.dot" );
        
        tic = now();

        dag.run( acc );

        toc = since( tic );
        
        auto  A_inv = std::make_unique< TLDLInvMatrix >( A_copy.get(), A->form() );

        std::cout << "  done in     " << toc << std::endl;
        std::cout << "  inv error = " << boost::format( "%.4e" ) % inv_approx_2( A.get(), A_inv.get() ) << std::endl
                  << "  memory    = " << Mem::to_string( A_copy->byte_size() ) << std::endl;
    }// if
    else if ( method == "waz" )
    {
        std::cout << "┎────────────────────────────────────────────────────────" << std::endl
                  << "┃ WAZ"                                                     << std::endl
                  << "┖────────────────────────────────────────────────────────" << std::endl;
        

        auto  A_copy = A->copy();
        
        tic = now();

        auto  dag = WAZ::gen_dag( A_copy.get(), fac_options_t() );

        toc = since( tic );

        std::cout << "  dag in      " << toc << std::endl;
        std::cout << "    #nodes  = " << dag.nnodes() << std::endl;
        std::cout << "    #edges  = " << dag.nedges() << std::endl;

        dag.print_dot( "waz.dot" );

        tic = now();

        dag.run( acc );
        
        toc = since( tic );
        
        auto  WZ     = WAZ::split( A_copy.get() );
        auto  W      = WZ.first;
        auto  Z      = WZ.second;
        auto  A_inv  = matrix_product( Z, W );

        std::cout << "  done in     " << toc << std::endl;
        std::cout << "  inv error = " << boost::format( "%.4e" ) % inv_approx_2( A.get(), A_inv.get() ) << std::endl
                  << "  memory    = " << Mem::to_string( W->byte_size() + Z->byte_size() ) << std::endl;
    }// if
    else if ( method == "invll" )
    {
        std::cout << "┎────────────────────────────────────────────────────────" << std::endl
                  << "┃ lower left triangular inversion"                         << std::endl
                  << "┖────────────────────────────────────────────────────────" << std::endl;

        // just to be sure
        A->set_unsymmetric();
        
        auto  A_inv = A->copy();

        tic = now();

        auto  dag = gen_dag_invert_ll( A_inv.get(), inv_options_t() );

        toc = since( tic );
        
        std::cout << "  dag in      " << toc << std::endl;
        std::cout << "    #nodes  = " << dag.nnodes() << std::endl;
        std::cout << "    #edges  = " << dag.nedges() << std::endl;

        dag.print_dot( "invll.dot" );
        
        tic = now();

        dag.run( acc );

        toc = since( tic );
    
        std::cout << "  done in     " << toc << std::endl;
        std::cout << "  inv error = " << boost::format( "%.4e" ) % inv_approx_2( A.get(), A_inv.get() ) << std::endl
                  << "  memory    = " << Mem::to_string( A_inv->byte_size() ) << std::endl;
    }// if
    else
        std::cout << "unknown method, try -h" << std::endl;
    
    DONE();
}
