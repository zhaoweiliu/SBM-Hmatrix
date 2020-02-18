#ifndef __HLIB_TMATRIXVIS_HH
#define __HLIB_TMATRIXVIS_HH
//
// Project     : HLib
// File        : TMatrixVis.hh
// Description : matrix visualisation classes
// Author      : Ronald Kriemann
// Copyright   : Max Planck Institute MIS 2004-2018. All Rights Reserved.
//

#include "hpro/matrix/TMatrix.hh"

namespace HLIB
{
    
///////////////////////////////////////////////////////////////
//!
//! \ingroup  IO_Module
//! \class    TMatrixVis
//! \brief    Base class for matrix visualisation.
//!
class TMatrixVis
{
public:
    //////////////////////////////////////
    //
    // constructor and destructor
    //

    TMatrixVis () {}

    virtual ~TMatrixVis () {}

    //////////////////////////////////////
    //
    // visualise matrix
    //

    //! print matrix \a A to file \a filename
    virtual void print ( const TMatrix *      A,
                         const std::string &  filename ) const = 0;
};

///////////////////////////////////////////////////////////////
//!
//! \ingroup  IO_Module
//! \class    T2DMatrixVis
//! \brief    Implements 2D based matrix visualisation.
//!

// forward decl.
class T2DPrinter;
struct option_t;

class T2DMatrixVis : public TMatrixVis
{
private:
    //! @cond

    // internal options
    option_t *  _opt;
    
    //! @endcond
    
public:
    ///////////////////////////////////////////////
    //
    // constructor and destructor
    //

    //! construct matrix visualisation object
    T2DMatrixVis ();

    virtual ~T2DMatrixVis ();

    ///////////////////////////////////////////////
    //
    // option management
    //
    
    //! turn on/off printing of matrix block border (default: on)
    T2DMatrixVis &  border          ( const bool  b );

    //! turn on/off color usage for visualization (default: on)
    T2DMatrixVis &  color           ( const bool  b );

    //! turn on/off printing of structure of matrix (default: on)
    T2DMatrixVis &  structure       ( const bool  b );

    //! turn on/off printing of indexsets (default: off)
    T2DMatrixVis &  indexset        ( const bool  b );

    //! turn on/off printing of id (default: off)
    T2DMatrixVis &  id              ( const bool  b );

    //! turn on/off printing of only nonempty blocks (default: off)
    T2DMatrixVis &  nonempty        ( const bool  b );

    //! turn on/off printing of local/non-local matrix blocks (default: off)
    T2DMatrixVis &  only_local      ( const bool  b );

    //! turn on/off printing of neighbourhood relation (default: off)
    T2DMatrixVis &  neighbours      ( const bool  b );

    //! turn on/off printing of matrix coefficients (default: off)
    T2DMatrixVis &  entries         ( const bool  b );

    //! turn on/off printing of sparsity pattern (default: off)
    T2DMatrixVis &  pattern         ( const bool  b );

    //! turn on/off printing of sparsity pattern of sparse matrices only (default: off)
    T2DMatrixVis &  sparse_pattern  ( const bool  b );

    //! turn on/off printing of SVD of matrix-blocks (default: off)
    T2DMatrixVis &  svd             ( const bool  b );

    //! turn on/off colouring matrix blocks according to rank compared with \a max_rank (default: off)
    //! if \a max_rank is 0, it is compared with block local maximal rank
    T2DMatrixVis &  rank_col        ( const bool           b,
                                      const size_t         max_rank = 0,
                                      const std::string &  cmap     = "default" );

    //! turn on/off colouring matrix blocks according to memory consumption (default: off)
    T2DMatrixVis &  mem_col         ( const bool  b,
                                      const std::string &  cmap = "default" );

    //! set maximal allowed ratio of block size compared to largest block size
    //! (default: 1000); if ratio is exceeded, block will not be printed to 
    //! limit file size for very large matrices
    T2DMatrixVis &  max_size_ratio  ( const double  r );

    //! set maximal level in matrix (no blocks on levels below will be printed) (default: 0 [off])
    T2DMatrixVis &  max_level       ( const size_t  l );

    //! set value of largest/smallest singular value to compare with for each block;
    //! - if \a ref_max/\a ref_min ≤ 0, the largest/smallest local singular values per block
    //!   are used
    //! - good choice for \a ref_max is ∥A∥₂
    T2DMatrixVis &  svd_ref         ( const double  ref_max,
                                      const double  ref_min = Limits::epsilon<double>() );

    ///////////////////////////////////////////////
    //
    // interface for visualisation
    //

    //! print matrix \a A to file \a filename
    virtual void print ( const TMatrix *      A,
                         const std::string &  filename ) const;

protected:
    // return printer object for given format
    virtual T2DPrinter *  get_printer ( const double         width,
                                        const double         height,
                                        const std::string &  filename ) const = 0;
};

//
//! \ingroup  IO_Module
//! \class    TPSMatrixVis
//! \brief    class for matrix visualisation in PostScript format
//
class TPSMatrixVis : public T2DMatrixVis
{
protected:
    // return printer object for given format
    virtual T2DPrinter *  get_printer ( const double         width,
                                        const double         height,
                                        const std::string &  filename ) const;
};

//
//! \ingroup  IO_Module
//! \fn       print_ps
//! \brief    functional version of TPSMatrixVis
//
void
print_ps ( const TMatrix *      cl,
           const std::string &  filename );

//
//! \ingroup  IO_Module
//! \class    TPDFMatrixVis
//! \brief    class for matrix visualisation in PDF format
//
class TPDFMatrixVis : public T2DMatrixVis
{
protected:
    // return printer object for given format
    virtual T2DPrinter *  get_printer ( const double         width,
                                        const double         height,
                                        const std::string &  filename ) const;
};

//
//! \ingroup  IO_Module
//! \fn       print_pdf
//! \brief    functional version of TPDFMatrixVis
//
void
print_pdf ( const TMatrix *      cl,
            const std::string &  filename );

}// namespace HLIB

#endif  // __HLIB_TMATRIXVIS_HH
