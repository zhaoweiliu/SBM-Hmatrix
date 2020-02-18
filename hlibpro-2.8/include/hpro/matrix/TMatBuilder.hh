#ifndef __HLIB_TMATBUILDER_HH
#define __HLIB_TMATBUILDER_HH
//
// Project     : HLib
// File        : TMatBuilder.hh
// Description : class for building H- or other matrices
// Author      : Ronald Kriemann
// Copyright   : Max Planck Institute MIS 2004-2018. All Rights Reserved.
//

#include <map>
#include <vector>
#include <unordered_map>

#include "hpro/cluster/TBlockClusterTree.hh"
#include "hpro/cluster/TClusterBasis.hh"
#include "hpro/cluster/TBCBuilder.hh"

#include "hpro/matrix/TMatrix.hh"
#include "hpro/matrix/TSparseMatrix.hh"
#include "hpro/matrix/TBlockMatrix.hh"
#include "hpro/matrix/THMatrix.hh"
#include "hpro/matrix/TCoeffFn.hh"

#include "hpro/misc/TProgressBar.hh"

#include "hpro/algebra/TLowRankApx.hh"

namespace HLIB
{

////////////////////////////////////////////////////////////////
//!
//! \class  TMatBuilder
//! \brief  Base class for building matrices implementing basic
//!         management and parallel construction.
//!
class TMatBuilder
{
public:

    //
    // statistic information of matrix building
    //
    struct stat_t {
        // time for nearfield/farfield blocks and coarsening
        double  time_nearfield;
        double  time_farfield;
        double  time_coarsen;

        // detailed timing for blocks addressed by their id
        std::unordered_map< int, double >  time_block;
        
        // number of nearfield/farfield blocks and number of coarsening
        size_t  count_nearfield;
        size_t  count_farfield;
        size_t  count_coarsen;

        stat_t ();

        friend std::ostream &  operator << ( std::ostream & os, const stat_t &  stat );
    };
    
protected:
    //! @cond
    
    // if true, coarsening is applied during construction
    bool       _coarsening;
    
    // if true, the accuracy during coarsening is identical to the
    // accuracy during block construction (standard accuracy)
    bool       _use_construct_acc;

    // if true, ghost matrices will be built for all non-local
    // matrices
    bool       _build_ghosts;
    
    // defines coarsening accuracy
    TTruncAcc  _coarse_acc;

    // block cluster tree builder for adaptive refinement
    const TBCBuilder *  _bc_builder;
    
    // statistics
    mutable stat_t  _statistics;
    
    //! @endcond
    
public:
    //////////////////////////////////////////////
    //
    // constructor and destructor
    //

    //! construct matrix construction object with standard coarsening,
    //! e.g. same precision as in block construction
    TMatBuilder ( const bool  coarsening = CFG::Build::coarsen );

    virtual ~TMatBuilder () {}

    //////////////////////////////////////////////
    //
    // change build behaviour
    //

    //! (de-) activate coarsening with standard accuracy (default: off)
    void set_coarsening  ( const bool  b )
    {
        _coarsening        = b;
        _use_construct_acc = true;
    }
    
    //! activate coarsening with accuracy \a acc
    void set_coarsening  ( const TTruncAcc &  acc )
    {
        _coarsening        = true;
        _use_construct_acc = false;
        _coarse_acc        = acc;
    }

    //! activate/deactivate construction of ghost matrices for
    //! non-local matrix blocks (default: off)
    void set_build_ghosts ( const bool  b )
    {
        _build_ghosts = b;
    }

    //! set block cluster tree builder object (default: nullptr)
    void set_bc_builder ( const TBCBuilder *  bc_builder )
    {
        _bc_builder = bc_builder;
    }
    
    //////////////////////////////////////////////
    //
    // build H-matrix
    //

    //! build the H-matrix with block-wise accuracy defined by \a acc
    virtual std::unique_ptr< TMatrix >  build ( const TBlockClusterTree *  bct,
                                                const TTruncAcc &          acc,
                                                TProgressBar *             progress = nullptr ) const;

    //! same as standard \see build, but build block matrices for given 
    //! block cluster without permutations etc.
    virtual std::unique_ptr< TMatrix >  build ( const TBlockCluster *      bc,
                                                const TTruncAcc &          acc,
                                                TProgressBar *             progress = nullptr ) const;

    
    //! same as \see build but with user defined matrix format \a matformat
    virtual std::unique_ptr< TMatrix >  build ( const TBlockClusterTree *  cluster,
                                                const matform_t            matformat,
                                                const TTruncAcc &          acc,
                                                TProgressBar *             progress = nullptr ) const;

    //! same as \see build but with user defined matrix format \a matformat
    virtual std::unique_ptr< TMatrix >  build ( const TBlockCluster *      cluster,
                                                const matform_t            matformat,
                                                const TTruncAcc &          acc,
                                                TProgressBar *             progress = nullptr ) const;

protected:
    
    //! threaded building process
    virtual std::unique_ptr< TMatrix >       thr_build     ( const TBlockCluster *  bc,
                                                             const matform_t        matformat,
                                                             const TTruncAcc &      acc,
                                                             TProgressBar *         progress ) const;

public:
    
    //! build blocked matrix
    virtual std::unique_ptr< TBlockMatrix >  build_blocked ( const TBlockCluster *  bc ) const;
    
    //! build matrix corresponding to leaves in the block cluster tree
    virtual std::unique_ptr< TMatrix >       build_leaf    ( const TBlockCluster *  bc,
                                                             const matform_t        matformat,
                                                             const TTruncAcc &      acc ) const = 0;

    //! build placeholder matrix for remote blocks
    virtual std::unique_ptr< TMatrix >       build_ghost   ( const TBlockCluster *  bc ) const = 0;

    //! return matrix format
    virtual matform_t      matrix_format () const = 0;

    //! return statistics information
    const stat_t &         statistics    () const { return _statistics; }
};

////////////////////////////////////////////////////////////////
//!
//! \class  TSparseMatBuilder
//! \brief  Creates H-matrices out of sparse matrices.
//!
class TSparseMatBuilder : public TMatBuilder
{
protected:
    //! sparse matrix holding coefficients
    const TSparseMatrix *  _sparse_mat;
    
    //! mapping of of internal numbering (in row cluster tree) to external numbering (in sparse matrix)
    const TPermutation *   _row_perm_i2e;

    //! mapping of of external numbering (in sparse matrix) to internal numbering (in column cluster tree)
    const TPermutation *   _col_perm_e2i;

    //! if true, build sparse matrices in leaves (default: false)
    bool                   _build_sparse;
    
    //! if true, use TZeroMatrix for off-diagonal domain-domain blocks (default: true)
    bool                   _use_zero_mat;

    //! value type of matrices to create (default: defined by sparse matrix)
    value_type_t           _value_type;
    
public:
    //////////////////////////////////////////////
    //
    // constructor and destructor
    //

    //! construct TSparseMatBuilder object with given sparse matrix and permutations
    TSparseMatBuilder ( const TSparseMatrix *  sparse_mat,
                        const TPermutation *   row_perm_i2e,
                        const TPermutation *   col_perm_e2i );

    //! construct TSparseMatBuilder object with given sparse matrix and permutations
    //! and user-defined value type
    TSparseMatBuilder ( const TSparseMatrix *  sparse_mat,
                        const TPermutation *   row_perm_i2e,
                        const TPermutation *   col_perm_e2i,
                        const value_type_t     value_type );

    virtual ~TSparseMatBuilder ();

    //! (de-) activate construction of sparse leaves
    void set_sparse_mode ( const bool  b )
    {
        _build_sparse = b;
    }
    
    //! (de-) activate construction of sparse leaves
    void set_use_zero_mat ( const bool  b )
    {
        _use_zero_mat = b;
    }
    
    //////////////////////////////////////////////
    //
    // build subblocks for given h-matrix
    //

protected:
    
    //! threaded building process
    virtual std::unique_ptr< TMatrix >  thr_build   ( const TBlockCluster * bc,
                                                      const matform_t       matformat,
                                                      const TTruncAcc &     acc,
                                                      TProgressBar *        progress ) const;

public:
    
    //! construct matrices for leaves in block cluster tree
    virtual std::unique_ptr< TMatrix >  build_leaf  ( const TBlockCluster * bc,
                                                      const matform_t       matformat,
                                                      const TTruncAcc &     acc ) const;

    //! build placeholder matrix for remote blocks
    virtual std::unique_ptr< TMatrix >  build_ghost ( const TBlockCluster * bc ) const;
    
    //! return matrix format
    virtual matform_t  matrix_format () const
    {
        return _sparse_mat->form();
    }

private:
    //! build dense matrix for given block cluster
    std::unique_ptr< TMatrix >  build_dense  ( const TBlockCluster * bc,
                                               const matform_t       matformat,
                                               const TTruncAcc &     acc ) const;
    
    //! build low-rank matrix for given block cluster
    std::unique_ptr< TMatrix >  build_rank   ( const TBlockCluster * bc,
                                               const matform_t       matformat,
                                               const TTruncAcc &     acc ) const;

    //! build sparse matrices for given block cluster
    std::unique_ptr< TMatrix >  build_sparse ( const TBlockCluster * bc,
                                               const matform_t       matformat,
                                               const TTruncAcc &     acc ) const;

    //! return true, if sparse matrix as non-zero entry in index set
    bool                        has_entry    ( const TIndexSet &     rowis,
                                               const TIndexSet &     colis ) const;
};

// define old name
using TSparseMBuilder = TSparseMatBuilder;

////////////////////////////////////////////////////////////////
//!
//! \class  TDenseMatBuilder
//! \brief  Creates matrices by computing low rank approximations
//!         of dense sub matrices, e.g. via ACA or SVD.
//!
template < typename T >
class TDenseMatBuilder : public TMatBuilder
{
public:
    //
    // template arguments as internal types
    //

    using  value_t    = T;
    using  coeff_fn_t = TCoeffFn< value_t >;

    //!
    //! \struct  stat_t
    //! \brief   statistical data of computations
    //!
    struct stat_t
    {
        //! number of requested matrix coefficients
        size_t  ncoeff;

        //! constructor
        stat_t ()
                : ncoeff(0)
        {}
    };
    
protected:
    //! function for computing matrix entries
    const coeff_fn_t *   _coeff_fn;
    
    //! low-rank approximation algorithm
    const TLowRankApx *  _lr_apx;
    
    //! if true, use recompression for low-rank blocks
    bool                 _recompress;

    //! if true, convert low-rank to dense if rank too high
    const bool           _convert_to_dense;

    //! statistics object
    stat_t *             _stat;
    
public:
    //////////////////////////////////////////////
    //
    // constructor and destructor
    //

    //! ctor for construction of H-matrix out of given coefficient function
    //! and low-rank approximation method.
    TDenseMatBuilder ( const coeff_fn_t *   coeff_fn,
                       const TLowRankApx *  lr_apx,
                       const bool           recompress   = CFG::Build::recompress,
                       const bool           cvt_to_dense = CFG::Build::to_dense,
                       stat_t *             stat         = nullptr );

    virtual ~TDenseMatBuilder ();
    
    //////////////////////////////////////////////
    //
    // change build behaviour
    //

    //! (de-) activate recompression of low-rank blocks
    void set_recompression ( const bool  b )
    {
        _recompress = b;
    }
    
    //! assign statistics object
    void set_stat ( stat_t *  stat )
    {
        _stat = stat;
    }
    
    //////////////////////////////////////////////
    //
    // build subblocks for given H-matrix
    //

public:
    
    //! construct matrices for leaves in the block cluster tree
    virtual std::unique_ptr< TMatrix > build_leaf  ( const TBlockCluster * bc,
                                                     const matform_t       matformat,
                                                     const TTruncAcc &     acc ) const;

    //! build placeholder matrix for remote blocks
    virtual std::unique_ptr< TMatrix > build_ghost ( const TBlockCluster * bc ) const;

protected:
    
    //! build dense matrix for given block cluster by directly evaluating
    //! the matrix coefficient via the coefficient function
    std::unique_ptr< TMatrix >         build_dense ( const TBlockCluster * bc,
                                                     const matform_t       matformat,
                                                     const TTruncAcc &     acc ) const;
    
    //! build low-rank matrix for given block cluster by using
    //! internal low-rank approximation object
    std::unique_ptr< TMatrix >         build_rank  ( const TBlockCluster * bc,
                                                     const TTruncAcc &     acc ) const;

public:
    
    //! return matrix format
    virtual matform_t                  matrix_format () const { return _coeff_fn->matrix_format(); }
};

// define old name
template < typename T > using TDenseMBuilder = TDenseMatBuilder< T >;

////////////////////////////////////////////////////////////////
//!
//! \class  TH2MatBuilder
//! \brief  Base class for H² matrix builders providing
//!         leaf bulding function with corresponding
//!         cluster bases
//!
template < typename T >
class TH2MatBuilder : public TMatBuilder
{
public:
    //
    // provide template types
    //

    using  value_t    = T;
    using  cl_basis_t = TClusterBasis< value_t >;
    using  cb_map_t   = std::map< TIndexSet, const cl_basis_t *, TIndexSet::map_cmp_t >;
    
protected:

    //! @cond
    
    // row and column cluster bases
    const cl_basis_t *  _row_cb;
    const cl_basis_t *  _col_cb;

    // mappings from index sets to cluster bases
    mutable cb_map_t    _rowcb_map, _colcb_map;
    
    //! @endcond
    
public:
    //////////////////////////////////////////////
    //
    // ctor
    //

    //! construct H2 matrix builder with supplied cluster bases
    TH2MatBuilder ( const cl_basis_t *  row_cb,
                    const cl_basis_t *  col_cb );
    
    //////////////////////////////////////////////
    //
    // build subblocks for given H²-matrix
    //

    //! construct matrices for leaves in the block cluster tree
    virtual std::unique_ptr< TMatrix >      build_leaf         ( const TBlockCluster *  bc,
                                                                 const matform_t        matformat,
                                                                 const TTruncAcc &      acc ) const;

    //! build blocked matrix
    virtual std::unique_ptr< TBlockMatrix > build_blocked      ( const TBlockCluster *  bc ) const;
    
    //! build placeholder matrix for remote blocks
    virtual std::unique_ptr< TMatrix >      build_ghost        ( const TBlockCluster *  bc ) const;

protected:
    
    //! construct uniform matrices for leaves in the block cluster tree
    virtual std::unique_ptr< TMatrix >      build_uniform_leaf ( const TBlockCluster *  bc,
                                                                 const cl_basis_t *     row_cb,
                                                                 const cl_basis_t *     col_cb,
                                                                 const matform_t        matformat,
                                                                 const TTruncAcc &      acc ) const = 0;
};

////////////////////////////////////////////////////////////////
//!
//! \class  TIdMatBuilder
//! \brief  Construct identity matrix for given block cluster trees
//!
class TIdMatBuilder : public TMatBuilder
{
private:
    //! value type of matrices to create (default: real_valued)
    value_type_t  _value_type;
    
public:
    //////////////////////////////////////////////
    //
    // ctor and dtor
    //

    //! construct identity matrix builder with coarsening disabled
    TIdMatBuilder ( const value_type_t  avalue_type = real_valued );

    //////////////////////////////////////////////
    //
    // build subblocks for given H-matrix
    //

    //! build matrix corresponding to leaves in the block cluster tree
    virtual std::unique_ptr< TMatrix >  build_leaf  ( const TBlockCluster *  bc,
                                                      const matform_t        matformat,
                                                      const TTruncAcc &      acc ) const;

    //! build placeholder matrix for remote blocks
    virtual std::unique_ptr< TMatrix >  build_ghost ( const TBlockCluster *  bc ) const;

    //! return matrix format
    virtual matform_t                   matrix_format () const { return MATFORM_HERM; }
};

//!
//! identity construction in functional form
//!
std::unique_ptr< TMatrix >
identity ( const TBlockClusterTree *  cluster,
           const value_type_t         avalue_type = real_valued );

////////////////////////////////////////////////////////////////
//!
//! \class  TZeroMatBuilder
//! \brief  Construct empty matrix for given block cluster trees
//!
class TZeroMatBuilder : public TMatBuilder
{
private:
    //! value type of matrices to create (default: real_valued)
    value_type_t  _value_type;
    
public:
    //////////////////////////////////////////////
    //
    // ctor and dtor
    //

    //! construct zero matrix builder with coarsening disabled
    TZeroMatBuilder ( const value_type_t  avalue_type = real_valued );

    //////////////////////////////////////////////
    //
    // build subblocks for given H-matrix
    //

    //! build matrix corresponding to leaves in the block cluster tree
    virtual std::unique_ptr< TMatrix >  build_leaf  ( const TBlockCluster *  bc,
                                                      const matform_t        matformat,
                                                      const TTruncAcc &      acc ) const;

    //! build placeholder matrix for remote blocks
    virtual std::unique_ptr< TMatrix >  build_ghost ( const TBlockCluster *  bc ) const;

    //! return matrix format (assuming general situation)
    virtual matform_t                   matrix_format () const { return unsymmetric; }
};

//!
//! construct empty matrix for given block cluster tree
//!
std::unique_ptr< TMatrix >
build_zero_mat ( const TBlockClusterTree *  cluster,
                 const value_type_t         avalue_type = real_valued );

////////////////////////////////////////////////////////////////
//!
//! \fn     assemble_block
//! \brief  Construct block matrix out of given set of matrices
//!
//!         This function creates a \a block_rows × \a block_cols
//!         block matrix with sub matrices defined (in column
//!         wise ordering) by \a submatrices. If \a copy is true
//!         the matrices will be copied. Otherwise, the given
//!         matrices itself will be used. Furthermore, all index sets
//!         and (possibly) permutations will be adjusted. Only if
//!         \b all matrices are of type THMatrix, an H-matrix is returned.
//!         It is assumed, that all cluster trees are compatible, e.g.
//!         equal in each block row and column.
//!
std::unique_ptr< TMatrix >
assemble_block ( const size_t  block_rows,
                 const size_t  block_cols,
                 TMatrix **    submatrices,
                 const bool    copy = false );

}// namespace

#endif  // __HLIB_TMATBUILDER_HH
