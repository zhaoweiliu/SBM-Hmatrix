#ifndef __HLIB_BEM_TBFCOEFFFN_H
#define __HLIB_BEM_TBFCOEFFFN_H
//
// Project     : HLib
// File        : TBFCoeffFn.hh
// Description : matrix coefficient functions for bilinear forms
// Author      : Ronald Kriemann
// Copyright   : Max Planck Institute MIS 2004-2018. All Rights Reserved.
//

#include <tbb/parallel_for.h>
#include <tbb/blocked_range2d.h>

#include "hpro/matrix/TCoeffFn.hh"
#include "hpro/bem/TBEMBF.hh"

namespace HLIB
{

////////////////////////////////////////////////////////////////////
//!
//! \class  TBFCoeffFn
//! \brief  Provide matrix coefficients defined by bilinear forms
//!
template < typename  T_bf > 
class TBFCoeffFn : public TCoeffFn< typename T_bf::value_t >
{
public:
    //
    // template types as internal types
    //
    using  bf_t       = T_bf;
    using  ansatzsp_t = typename bf_t::ansatzsp_t;
    using  testsp_t   = typename bf_t::testsp_t;
    using  value_t    = typename bf_t::value_t;

protected:
    //! bilinear form to evaluate
    const bf_t *  _bf;
    
public:
    //////////////////////////////////////
    //
    // constructor and destructor
    //

    //! construct coefficient function for bilinear form \a bf
    TBFCoeffFn ( const bf_t *  bf )
            : _bf( bf )
    {
        if ( _bf == NULL )
            HERROR( ERR_ARG, "(TBFCoeffFn)", "invalid bilinear form supplied (NULL)" );
    }

    //! destructor
    virtual ~TBFCoeffFn () {}

    //! return true if function is complex valued
    bool       is_complex     () const { return _bf->is_complex(); }
        
    //! return format of matrix, e.g. symmetric or hermitian
    matform_t  matrix_format  () const { return _bf->format(); }
    
    //////////////////////////////////////
    //
    // evaluation of coefficient
    //

    //! evaluate matrix coefficients in \a rowis × \a colis and store values
    //! in \a matrix (real valued)
    virtual void eval ( const std::vector< idx_t > &  rowidxs,
                        const std::vector< idx_t > &  colidxs,
                        value_t *                     matrix ) const
    {
        //
        // call bilinear form
        //
        
        const size_t  n = rowidxs.size();
        const size_t  m = colidxs.size();

        if ( std::min( n, m ) >= 1024 )
        {
            //
            // parallel mode
            //

            const size_t  block_size = 256;

            tbb::parallel_for( tbb::blocked_range2d< size_t >( 0, n, block_size,
                                                               0, m, block_size ),
                               [&] ( const tbb::blocked_range2d< size_t > &  r )
                               {
                                   const size_t  n_start = r.rows().begin();
                                   const size_t  n_end   = r.rows().end();
                                   const size_t  n_sub   = n_end - n_start;
                                   const size_t  m_start = r.cols().begin();
                                   const size_t  m_end   = r.cols().end();
                                   const size_t  m_sub   = m_end - m_start;

                                   std::vector< idx_t >     sub_rowidxs( n_sub );
                                   std::vector< idx_t >     sub_colidxs( m_sub );
                                   BLAS::Matrix< value_t >  sub_M( n_sub, m_sub );

                                   for ( size_t  i = 0; i < n_sub; ++i ) sub_rowidxs[i] = rowidxs[ i + n_start ];
                                   for ( size_t  i = 0; i < m_sub; ++i ) sub_colidxs[i] = colidxs[ i + m_start ];
                                   
                                   _bf->eval( sub_rowidxs, sub_colidxs, sub_M );

                                   for ( idx_t  j = 0; j < idx_t(m_sub); ++j )
                                       for ( idx_t  i = 0; i < idx_t(n_sub); ++i )
                                           matrix[ ( j + m_start ) * n + ( i + n_start ) ] = sub_M( i, j );
                               } );
        }// if
        else
        {
            //
            // sequential mode
            //
            
            BLAS::Matrix< value_t >  M( n, m );
        
            _bf->eval( rowidxs, colidxs, M );
            
            for ( idx_t  i = 0; i < idx_t(n*m); ++i )
                matrix[i] = M.data()[i];
        }// else
    }

    using TCoeffFn< value_t >::eval;

    DISABLE_COPY_OP( TBFCoeffFn );
};

}// namespace HLIB

#endif  // __HLIB_BEM_TBFCOEFFFN_H
