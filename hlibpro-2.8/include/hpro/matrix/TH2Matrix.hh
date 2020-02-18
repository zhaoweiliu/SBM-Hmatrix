#ifndef __HLIB_TH2MATRIX_HH
#define __HLIB_TH2MATRIX_HH
//
// Project     : HLib
// File        : TH2Matrix.hh
// Description : class for H²-matrices
// Author      : Ronald Kriemann
// Copyright   : Max Planck Institute MIS 2004-2018. All Rights Reserved.
//

#include "hpro/cluster/TClusterBasis.hh"
#include "hpro/cluster/TPermutation.hh"

#include "hpro/matrix/TBlockMatrix.hh"

namespace HLIB
{

// local matrix type
DECLARE_TYPE( TH2Matrix );

//!
//! \ingroup  Matrix_Module
//! \class    TH2Matrix
//! \brief    Class for an H²-matrix, which extends block matrices
//!           with additional functionality, e.g. permutations and
//!           uniform vectors.
//!
class TH2Matrix : public TBlockMatrix
{
private:
    // cluster bases for rows and columns
    const TClusterBasis< real > *     _rrow_cb;
    const TClusterBasis< real > *     _rcol_cb;
    const TClusterBasis< complex > *  _crow_cb;
    const TClusterBasis< complex > *  _ccol_cb;
    
    // mappings of indices from external to
    // internal numbering (for mul_vec)
    TPermutation                      _row_perm_e2i;
    TPermutation                      _col_perm_e2i;
    TPermutation                      _row_perm_i2e;
    TPermutation                      _col_perm_i2e;
    
public:
    ///////////////////////////////////////////////
    //
    // constructor and destructor
    //

    //! default ctor without cluster bases
    TH2Matrix ();
    
    //! construct H²-matrix over block index set \a bis
    //! with cluster bases \a row_cb and \a col_cb
    TH2Matrix ( const TBlockIndexSet &            bis,
                const TClusterBasis< real > *     row_cb,
                const TClusterBasis< real > *     col_cb );

    //! construct H²-matrix over block index set \a bis
    //! with cluster bases \a row_cb and \a col_cb
    TH2Matrix ( const TBlockIndexSet &            bis,
                const TClusterBasis< complex > *  row_cb,
                const TClusterBasis< complex > *  col_cb );

    //! construct H²-matrix over block index set \a bis
    //! with cluster bases \a row_cb and \a col_cb and
    //! permutations \a row_perm_e2i and \a col_perm_e2i
    TH2Matrix ( const TBlockIndexSet &            bis,
                const TClusterBasis< real > *     row_cb,
                const TClusterBasis< real > *     col_cb,
                const TPermutation &              row_perm_e2i,
                const TPermutation &              col_perm_e2i );

    //! construct H²-matrix over block index set \a bis
    //! with cluster bases \a row_cb and \a col_cb and
    //! permutations \a row_perm_e2i and \a col_perm_e2i
    TH2Matrix ( const TBlockIndexSet &            bis,
                const TClusterBasis< complex > *  row_cb,
                const TClusterBasis< complex > *  col_cb,
                const TPermutation &              row_perm_e2i,
                const TPermutation &              col_perm_e2i );

    virtual ~TH2Matrix () {}

    /////////////////////////////////////////////////
    //
    // access internal data
    //

    //! access single matrix coefficient
    virtual real           entry         ( const idx_t i, const idx_t j ) const;
    virtual const complex  centry        ( const idx_t i, const idx_t j ) const;

    
    //! access row cluster basis
    const TClusterBasis< real > *     rrow_cb () const { return _rrow_cb; }
    const TClusterBasis< complex > *  crow_cb () const { return _crow_cb; }
    
    //! access column cluster basis
    const TClusterBasis< real > *     rcol_cb () const { return _rcol_cb; }
    const TClusterBasis< complex > *  ccol_cb () const { return _ccol_cb; }

    //! assign cluster bases to local matrix and to all uniform sub blocks
    //! - rank and value type of bases must be identical to corresponding
    //!   dimension and value type of coefficient matrices (also in sub blocks!)
    void  assign_cb    ( const TClusterBasis< real > *     row_cb,
                         const TClusterBasis< real > *     col_cb );
    void  assign_cb    ( const TClusterBasis< complex > *  row_cb,
                         const TClusterBasis< complex > *  col_cb );

    //! set row permutations
    void set_row_perm ( const TPermutation & perm_e2i,
                        const TPermutation & perm_i2e )
    {
        _row_perm_e2i = perm_e2i;
        _row_perm_i2e = perm_i2e;
    }

    //! set column permutations
    void set_col_perm ( const TPermutation & perm_e2i,
                        const TPermutation & perm_i2e )
    {
        _col_perm_e2i = perm_e2i;
        _col_perm_i2e = perm_i2e;
    }

    //! access row permutation (extern to intern)
    const TPermutation &   row_perm_e2i  () const { return _row_perm_e2i; }

    //! access row permutation (intern to extern)
    const TPermutation &   row_perm_i2e  () const { return _row_perm_i2e; }

    //! access column permutation (extern to intern)
    const TPermutation &   col_perm_e2i  () const { return _col_perm_e2i; }

    //! access column permutation (intern to extern)
    const TPermutation &   col_perm_i2e  () const { return _col_perm_i2e; }

    //! return true if \b all mappings are present and of correct size
    bool                   has_perm      () const
    {
        return (( _row_perm_e2i.size() == rows() ) &&
                ( _row_perm_i2e.size() == rows() ) &&
                ( _col_perm_e2i.size() == cols() ) &&
                ( _col_perm_i2e.size() == cols() ));
    }
    
    /////////////////////////////////////////////////
    //
    // BLAS-routines (real valued)
    //

    //! compute y ≔ β·y + α·op(M)·x, with M = this
    virtual void mul_vec ( const real      alpha,
                           const TVector * x,
                           const real      beta,
                           TVector       * y,
                           const matop_t   op = MATOP_NORM ) const;

    /////////////////////////////////////////////////
    //
    // BLAS-routines (complex valued)
    //

    //! compute y ≔ β·y + α·op(M)·x, with M = this
    virtual void cmul_vec ( const complex   alpha,
                            const TVector * x,
                            const complex   beta,
                            TVector       * y,
                            const matop_t   op = MATOP_NORM ) const;
    
    /////////////////////////////////////////////////
    //
    // misc.
    //

    //
    // virtual constructor
    //

    //! return matrix of same class (but no content)
    virtual auto  create       () const -> std::unique_ptr< TMatrix >
    {
        return std::make_unique< TH2Matrix >();
    }
    
    //! return copy of matrix
    virtual auto  copy         () const -> std::unique_ptr< TMatrix >;
    
    //! copy matrix wrt. accuracy \a acc and optional coarsening
    virtual auto  copy         ( const TTruncAcc &  acc,
                                 const bool         coarsen = false ) const -> std::unique_ptr< TMatrix >;

    //! return structural copy of matrix
    virtual auto  copy_struct  () const -> std::unique_ptr< TMatrix >;
    
    //! copy matrix into \a A
    virtual void  copy_to      ( TMatrix *          A ) const;

    //! copy matrix into \a A with accuracy \a acc and optional coarsening
    virtual void  copy_to      ( TMatrix *          A,
                                 const TTruncAcc &  acc,
                                 const bool         coarsen = false ) const;

    //! copy complete structural information from given matrix
    virtual void  copy_struct_from  ( const TMatrix *    A );
    
    //! return appropriate row vector type for matrix
    virtual auto  row_vector   () const -> std::unique_ptr< TVector >;

    //! return appropriate column vector type for matrix
    virtual auto  col_vector   () const -> std::unique_ptr< TVector >;
    
    //
    // size of object
    //
    
    //! return size in bytes used by this object
    virtual size_t byte_size () const;

    //
    // RTTI
    //

    HLIB_RTTI_DERIVED( TH2Matrix, TBlockMatrix )
    
    //
    // serialisation
    //

    //! read data from stream \a s and copy to matrix
    virtual void read  ( TByteStream &  s );

    //! use data from stream \a s to build matrix
    virtual void build ( TByteStream &  s );

    //! write data to stream \a s
    virtual void write ( TByteStream &  s ) const;

    //! returns size of object in bytestream
    virtual size_t bs_size () const;
};

}// namespace HLIB

#endif  // __HLIB_TH2MATRIX_HH
