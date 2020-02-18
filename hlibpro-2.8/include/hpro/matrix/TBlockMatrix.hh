#ifndef __HLIB_TBLOCKMATRIX_HH
#define __HLIB_TBLOCKMATRIX_HH
//
// Project     : HLib
// File        : TBlockMatrix.hh
// Description : class for a matrix consisting of submatrices
// Author      : Ronald Kriemann
// Copyright   : Max Planck Institute MIS 2004-2018. All Rights Reserved.
//

#include <list>

#include "hpro/matrix/TMatrix.hh"
#include "hpro/matrix/TRkMatrix.hh"

namespace HLIB
{

// local matrix type
DECLARE_TYPE( TBlockMatrix );

//!
//! \ingroup  Matrix_Module
//! \class    TBlockMatrix
//! \brief    Class for a n×m block matrix of TMatrix sub matrices.
//!
class TBlockMatrix : public TMatrix
{
private:
    //! @cond

    //! number of rows of matrix
    size_t      _rows;

    //! number of columns of matrix
    size_t      _cols;
    
    //! number of block rows of matrix
    uint        _block_rows;

    //! number of block columns of matrix
    uint        _block_cols;
    
    //! array of sub-blocks of this matrix (column-wise storage !!!)
    TMatrix  ** _blocks;

    //! @endcond
    
public:
    ///////////////////////////////////////////
    //
    // constructor and destructor
    //

    //! construct block matrix with size and block structure defined by \a bct
    TBlockMatrix ( const TBlockCluster * bct = nullptr )
            : TMatrix( bct )
            , _block_rows(0)
            , _block_cols(0)
            , _blocks(nullptr)
    {
        if ( bct != nullptr ) set_cluster( bct );
        else                  set_block_struct( 1, 1 );
    }

    //! construct block matrix with over block index set \a row_is × \a col_is
    TBlockMatrix ( const TIndexSet &  row_is,
                   const TIndexSet &  col_is )
            : TMatrix()
            , _block_rows(0)
            , _block_cols(0)
            , _blocks(nullptr)
    {
        set_ofs(  row_is.first(), col_is.first() );
        set_size( row_is.size(),  col_is.size() );
    }

    //! construct block matrix with over block index set \a bis
    TBlockMatrix ( const TBlockIndexSet &  bis )
            : TMatrix()
            , _block_rows(0)
            , _block_cols(0)
            , _blocks(nullptr)
    {
        set_ofs(  bis.row_is().first(), bis.col_is().first() );
        set_size( bis.row_is().size(),  bis.col_is().size() );
    }

    //! dtor
    virtual ~TBlockMatrix ();

    //! set cluster of matrix
    virtual void          set_cluster ( const TBlockCluster * c );
    
    ///////////////////////////////////////////
    //
    // block-handling
    //

    //! return number of rows
    virtual size_t        rows () const { return _rows; }

    //! return number of columns
    virtual size_t        cols () const { return _cols; }

    //! set size of matrix
    void                  set_size ( const size_t  r, const size_t  c ) { _rows = r; _cols = c; }
    
    //! return number of block-rows
    uint                  block_rows   ()                    const { return _block_rows; }
    uint                  block_rows   ( const matop_t  op ) const { return ( op == apply_normal
                                                                              ? block_rows()
                                                                              : block_cols() ); }
    uint                  nblock_rows  ()                    const { return block_rows(); }
    uint                  nblock_rows  ( const matop_t  op ) const { return block_rows( op ); }
    
    //! return number of block-columns
    uint                  block_cols   ()                    const { return _block_cols; }
    uint                  block_cols   ( const matop_t  op ) const { return ( op == apply_normal
                                                                              ? block_cols()
                                                                              : block_rows() ); }
    uint                  nblock_cols  ()                    const { return block_cols(); }
    uint                  nblock_cols  ( const matop_t  op ) const { return block_cols( op ); }

    //! return number of sub blocks of matrix
    uint                  no_of_blocks () const { return _block_rows*_block_cols; }

    //! set block structure
    void                  set_block_struct ( const uint n, const uint m );

    //! set matrix format
    virtual void          set_form         ( const matform_t  f );
    
    //! set matrix format
    virtual void          set_form         ( const matform_t         f,
                                             const recursion_type_t  rec_type );
    
    //! return matrix coefficient (\a i, \a j)
    virtual real          entry  ( const idx_t i, const idx_t j ) const;

    //! return matrix coefficient (\a i, \a j)
    virtual const complex centry ( const idx_t i, const idx_t j ) const;

    //! return block at block index (\a i, \a j)
    TMatrix *             block ( const uint i, const uint j )       { return _blocks[(j*_block_rows)+i]; }

    //! return block at block index (\a i, \a j)
    const TMatrix *       block ( const uint i, const uint j ) const { return _blocks[(j*_block_rows)+i]; }

    //! return block at block index (\a i, \a j)
    TMatrix *             block ( const uint     i,
                                  const uint     j,
                                  const matop_t  op )                { return ( op == apply_normal
                                                                                ? _blocks[(j*_block_rows)+i]
                                                                                : _blocks[(i*_block_rows)+j] ); }

    //! return block at block index (\a i, \a j)
    const TMatrix *       block ( const uint     i,
                                  const uint     j,
                                  const matop_t  op ) const          { return ( op == apply_normal
                                                                                ? _blocks[(j*_block_rows)+i]
                                                                                : _blocks[(i*_block_rows)+j] ); }

    //! set matrix block at block index (\a i,\a j) to matrix \a A
    void                  set_block ( const uint  i,
                                      const uint  j,
                                      TMatrix *   A )
    {
        _blocks[(j*_block_rows)+i] = A;

        if ( A != nullptr )
            A->set_parent( this );
    }

    //! replace matrix block \a A by matrix \a B (don't delete A !)
    void                  replace_block ( TMatrix * A, TMatrix * B );
    
    //! delete block (i,j)
    void                  delete_block ( const uint i, const uint j )
    {
        delete _blocks[(j*_block_rows)+i];
        _blocks[(j*_block_rows)+i] = nullptr;
    }
    
    //! return subblock of matrix corresponding to block cluster \a t
    TMatrix *             bc_block ( const TBlockCluster * t ) const;

    //! return subblock of matrix corresponding to block cluster (\a tau, \a sigma)
    TMatrix *             bc_block ( const TCluster * tau, const TCluster * sigma ) const;

    //! clear pointers to all sub blocks
    void                  clear_blocks ();
    
    //! convert data to real valued representation (if possible)
    virtual void          to_real    ();

    //! convert data to complex valued representation (if possible)
    virtual void          to_complex ();

    //! make value type of this and all sub blocks consistent
    virtual void          adjust_value_type ();
    
    //! return true, if matrix is blocked
    virtual bool          is_blocked () const { return true; }

    //! set processor set of matrix
    virtual void          set_procs  ( const TProcSet &        ps,
                                       const recursion_type_t  rec_type = nonrecursive );

    ///////////////////////////////////////////
    //
    // management of update accumulator
    //

    //! apply stored updates U to local matrix M, e.g., M = M + U,
    //! with accuracy \a acc
    virtual void  apply_updates ( const TTruncAcc &       acc,
                                  const recursion_type_t  recursion );
    
    //! return true, if matrix has updates not yet applied
    virtual bool  has_updates   ( const recursion_type_t  recursion ) const;

    /////////////////////////////////////////////////
    //
    // BLAS-routines (real valued)
    //

    //! compute this ≔ α·this
    virtual void scale ( const real alpha );
    
    //! compute this ≔ this + α · matrix
    virtual void add ( const real alpha, const TMatrix * matrix );
    
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

    //! compute this ≔ α·this
    virtual void cscale ( const complex alpha );
    
    //! compute this ≔ this + α · matrix
    virtual void cadd ( const complex alpha, const TMatrix * matrix );

    //! compute y ≔ β·y + α·op(M)·x, with M = this
    virtual void cmul_vec ( const complex   alpha,
                            const TVector * x,
                            const complex   beta,
                            TVector       * y,
                            const matop_t   op = MATOP_NORM ) const;
    
    ///////////////////////////////////////////////////////////
    //
    // linear operator mapping
    //

    //! same as above but only the dimension of the vector spaces is tested,
    //! not the corresponding index sets
    virtual void  apply_add   ( const real                       alpha,
                                const BLAS::Vector< real > &     x,
                                BLAS::Vector< real > &           y,
                                const matop_t                    op = apply_normal ) const;
    virtual void  apply_add   ( const complex                    alpha,
                                const BLAS::Vector< complex > &  x,
                                BLAS::Vector< complex > &        y,
                                const matop_t                    op = apply_normal ) const;

    using TMatrix::apply_add;

    /////////////////////////////////////////////////
    //
    // misc.
    //

    //! collect matrix blocks corresponding to leaves in list \a leaf_list
    template < typename T_list >
    void collect_leaves ( T_list &  leaf_list ) const
    {
        for ( uint i = 0; i < block_rows(); i++ )
            for ( uint j = 0; j < block_cols(); j++ )
            {
                const TMatrix  * A_ij = block(i,j);
            
                if ( A_ij == nullptr )
                    continue;
            
                if ( ! IS_TYPE( A_ij, TBlockMatrix ) )
                    leaf_list.push_back( const_cast< TMatrix * >( A_ij ) );
                else
                    cptrcast( A_ij, TBlockMatrix )->collect_leaves( leaf_list );
            }// for
    }

    //! transpose matrix
    virtual void transpose ();
    
    //! conjugate matrix coefficients
    virtual void conjugate ();
    
    //! truncate all rank-blocks in matrix to accuracy \a acc
    virtual void truncate ( const TTruncAcc & acc );
    
    //! print matrix to stdout
    virtual void print ( const uint ofs = 0 ) const;
    
    //! return matrix of same class (but no content)
    virtual auto create () const -> std::unique_ptr< TMatrix >
    {
        return std::make_unique< TBlockMatrix >();
    }

    //! return copy of matrix
    virtual auto copy         () const -> std::unique_ptr< TMatrix >;

    //! copy matrix wrt. accuracy \a acc and optional coarsening
    virtual auto copy         ( const TTruncAcc & acc,
                                const bool        coarsen = false ) const -> std::unique_ptr< TMatrix >;

    //! return structural copy of matrix
    virtual auto copy_struct  () const -> std::unique_ptr< TMatrix >;

    //! copy matrix into \a A
    virtual void copy_to      ( TMatrix * A ) const;

    //! copy matrix into \a A with accuracy \a acc and optional coarsening
    virtual void copy_to      ( TMatrix         * A,
                                const TTruncAcc & acc,
                                const bool        coarsen = false ) const;

    //! return size in bytes used by this object
    virtual size_t byte_size () const;

    //! copy complete structural information from given matrix
    virtual void copy_struct_from ( const TMatrix * M );
    
    //
    // RTTI
    //

    HLIB_RTTI_DERIVED( TBlockMatrix, TMatrix )

    //
    // serialisation
    //

    //! read data from stream \a s and copy to matrix
    virtual void read  ( TByteStream & s );

    //! use data from stream \a s to build matrix
    virtual void build ( TByteStream & s );

    //! write data to stream \a s
    virtual void write ( TByteStream & s ) const;

    //! returns size of object in bytestream
    virtual size_t bs_size () const;

    //
    // tests
    //

    //! test data for invalid values, e.g. INF and NAN
    virtual void check_data () const;
};

}// namespace

#endif  // __HLIB_TBLOCKMATRIX_HH
