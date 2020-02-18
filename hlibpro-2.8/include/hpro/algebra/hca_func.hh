#ifndef __HLIB_HCA_FUNC_HH
#define __HLIB_HCA_FUNC_HH
//
// Project     : HLib
// File        : hca_func.hh
// Description : special HCA functions for SLP, DLP and generic function interface
// Author      : Ronald Kriemann
// Copyright   : Max Planck Institute MIS 2004-2018. All Rights Reserved.
//

#include "hpro/base/TPoint.hh"

#include "hpro/cluster/TPermutation.hh"

#include "hpro/algebra/TLowRankApx.hh"

#include "hpro/bem/TFnSpace.hh"

namespace HLIB
{

#if 0

////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////
//
// generic call-back based function interface for HCA
//
////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////

//
// function to evaluate kernel at given coordinates
//
class THCACBGeneratorFn : public THCA::TGeneratorFn
{
public:
    // function type
    using  func_t  = real    (*) ( const double * x, const double * y, void * arg );
    using  cfunc_t = complex (*) ( const double * x, const double * y, void * arg );
    
protected:
    // function to call
    func_t   _func;
    cfunc_t  _cfunc;
    
    // additional argument
    void  * _arg;
    
public:
    //
    // constructor
    //

    THCACBGeneratorFn ( func_t  fn,
                        void *  arg )
    {
        _func  = fn;
        _cfunc = NULL;
        _arg   = arg;
        
        if ( fn == NULL )
            HERROR( ERR_ARG, "(THCACBGeneratorFn)", "kernel function is NULL" );
    }
    
    THCACBGeneratorFn ( cfunc_t  fn,
                        void *   arg )
    {
        _cfunc = fn;
        _func  = NULL;
        _arg   = arg;
        
        if ( fn == NULL )
            HERROR( ERR_ARG, "(THCACBGeneratorFn)", "kernel function is NULL" );
    }

    // return info about real or complex valued function
    bool is_complex () const { return _cfunc == NULL ? false : true; }
    
    //
    // evaluate function
    //
    
    real eval ( const T3Point & x, const T3Point & y ) const
    {
        if ( _func != NULL )
            return _func( x.vector(), y.vector(), _arg );
        else
            HERROR( ERR_COMPLEX, "(THCACBGeneratorFn) eval", "" );
    }

    complex ceval ( const T3Point & x, const T3Point & y ) const
    {
        if ( _cfunc != NULL )
            return _cfunc( x.vector(), y.vector(), _arg );
        else
            HERROR( ERR_REAL, "(THCACBGeneratorFn) ceval", "" );
    }
};

//
// function to return integral over kernel and basisfunction
//
class THCACBIntegralFn : public THCA::TIntegralFn
{
public:
    // function type
    using  func_t  = real    (*) ( const uint i, const real * x, void * arg );
    using  cfunc_t = complex (*) ( const uint i, const real * x, void * arg );
    
protected:
    // function to call
    func_t               _func;
    cfunc_t              _cfunc;
    
    // additional argument
    void *               _arg;
    
    // mapping from internal to external numering
    const TPermutation * _perm;
    
public:
    //
    // constructor
    //
    
    THCACBIntegralFn ( func_t                fn,
                       void *                arg,
                       const TPermutation *  perm )
    {
        _func  = fn;
        _cfunc = NULL;
        _arg   = arg;
        _perm  = perm;
        
        if ( fn == NULL )
            HERROR( ERR_ARG, "(THCACBIntegralFn)", "integral function is NULL" );
    }
    
    THCACBIntegralFn ( cfunc_t               fn,
                       void *                arg,
                       const TPermutation *  perm )
    {
        _cfunc = fn;
        _func  = NULL;
        _arg   = arg;
        _perm  = perm;
        
        if ( fn == NULL )
            HERROR( ERR_ARG, "(THCACBIntegralFn)", "integral function is NULL" );
    }
    
    // return info about real or complex valued function
    bool is_complex () const { return _cfunc == NULL ? false : true; }
    
    //
    // evaluate function
    //
    
    void eval ( const TIndexSet &,
                const std::vector< T3Point > &,
                BLAS::Matrix< real > &,
                const uint ) const
    {
        /*
        if ( _perm == NULL )
            return _func( i, x.vector(), _arg );
        else
            return _func( _perm->permute(i), x.vector(), _arg );
        */
    }

    void ceval ( const TIndexSet &,
                 const std::vector< T3Point > &,
                 BLAS::Matrix< complex > &,
                 const uint ) const
    {
        /*
        if ( _perm == NULL )
            return _cfunc( i, x.vector(), _arg );
        else
            return _cfunc( _perm->permute(i), x.vector(), _arg );
        */
    }
};


#endif

}// namespace

#endif  // __HLIB_HCA_FUNC_HH
