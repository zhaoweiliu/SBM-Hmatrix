#ifndef __HLIB_TCOARSEN_HH
#define __HLIB_TCOARSEN_HH
//
// Project     : HLib
// File        : TCoarsen.hh
// Description : class for coarsening H-matrices
// Author      : Ronald Kriemann
// Copyright   : Max Planck Institute MIS 2004-2018. All Rights Reserved.
//

#include "hpro/matrix/TMatrix.hh"
#include "hpro/misc/TProgressBar.hh"

namespace HLIB
{

//!
//! \class  TCoarsen
//! \brief  Implements coarsening, e.g. agglomeration of blocked matrices
//!         into dense or low-rank matrices if the later use less memory.
//!
class TCoarsen
{
protected:
    //! @cond
    
    //! maximal size of block-matrices, still worth to try coarsening
    size_t  _max_block_size;

    //! maximal size of dense matrices, still worth to try coarsening
    size_t  _max_dense_size;

    //! allow conversion to dense, e.g. low-rank to dense or block to dense
    bool    _allow_dense;
    
    //! recompress low-rank blocks if true
    bool    _recompress_lr;
    
    //! do not convert diagonal blocks to low-rank if true
    bool    _only_offdiag_lr;

    //! @endcond
    
public:
    //
    // constructor and destructor
    //

    //! constructs coarsening object
    TCoarsen ( const bool    aallow_dense    = true,
               const bool    arecompress_lr  = false,
               const size_t  amax_block_size = 50,
               const size_t  amax_dense_size = 500 );

    //! destructor
    virtual ~TCoarsen () {}

    //
    // options
    //

    //! set dense conversion mode, e.g. low-rank to dense or block to dense
    TCoarsen &  allow_dense      ( const bool    b );
    
    //! set recompression of low-rank blocks
    TCoarsen &  recompress_lr    ( const bool    b );
    
    //! turn on/off low-rank compression of diagonal blocks
    TCoarsen &  only_offdiag_lr  ( const bool    b );
    
    //! set maximal size of block-matrices, still worth to try coarsening
    TCoarsen &  max_block_size   ( const size_t  n );

    //! set maximal size of dense matrices, still worth to try coarsening
    TCoarsen &  max_dense_size   ( const size_t  n );
    
    //
    // coarsen matrix up to given precision
    //

    //! coarsen matrix \a A with recursion, i.e. first try to coarsen subblocks.
    //! \param  A          matrix to coarsen
    //! \param  acc        accuracy to maintain during coarsening
    virtual TMatrix * rec_coarsen ( TMatrix *          A,
                                    const TTruncAcc &  acc ) const;

    //! coarsen matrix \a A without recursion, i.e. \b no coarsening of subblocks.
    //! \param  A    matrix to coarsen
    //! \param  acc  accuracy to maintain during coarsening
    virtual TMatrix * coarsen     ( TMatrix *          A,
                                    const TTruncAcc &  acc ) const;

protected:
    
    //! return "corrected" bytesizes, e.g. in symmetric case, of matrices for internal use
    size_t byte_size ( const TMatrix * A ) const;

    DISABLE_COPY_OP( TCoarsen );
};

//
// functional version
//
inline
TMatrix *
rec_coarsen ( TMatrix *          A,
              const TTruncAcc &  acc )
{
    TCoarsen  mcoarsen;

    return mcoarsen.rec_coarsen( A, acc );
}

//
// try to replace dense by lowrank without changing structure
//
TMatrix *
compress ( TMatrix *          A,
           const TTruncAcc &  acc );

}// namespace HLIB

#endif  // __HLIB_TCOARSEN_HH
