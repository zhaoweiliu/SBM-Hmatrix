#ifndef __HLIB_MAT_INV_HH
#define __HLIB_MAT_INV_HH
//
// Project     : HLib
// File        : mat_inv.hh
// Description : matrix inversion
// Author      : Ronald Kriemann
// Copyright   : Max Planck Institute MIS 2004-2018. All Rights Reserved.
//

#include "hpro/matrix/TMatrix.hh"
#include "hpro/matrix/TDenseMatrix.hh"

#include "hpro/misc/TProgressBar.hh"

#include "hpro/algebra/solve_types.hh"

#include "hpro/parallel/dag.hh"

namespace HLIB
{

/////////////////////////////////////////////////////////////////////////////////////
//
// option type for inversion
//
/////////////////////////////////////////////////////////////////////////////////////

//!
//! \ingroup Algebra_Module
//! \struct  inv_options_t
//! \brief   options for matrix inversion
//!
struct inv_options_t
{
    //
    // inversion properties
    //

    //! is diagonal unit or not
    diag_type_t     diag;

    //! defines storage format for diagonal blocks (e.g., in invert_ll/ur)
    storage_type_t  storage;
    
    //! if true, coarsening is applied during factorisation
    bool            do_coarsen;
        
    //! progress bar for inversion
    TProgressBar *  progressbar;

    //
    // ctor
    //
        
    //! default constructor
    inv_options_t ( const diag_type_t     adiag     = general_diag,
                    const storage_type_t  astorage  = store_normal,
                    const bool            acoarsen  = CFG::Arith::coarsen,
                    TProgressBar *        aprogress = nullptr )
            : diag( adiag )
            , storage( astorage )
            , do_coarsen( acoarsen )
            , progressbar( aprogress )
    {}

    //! default constructor (with simplification of progress bar usage)
    inv_options_t ( TProgressBar *  progress )
    {
        diag        = general_diag;
        storage     = store_normal;
        do_coarsen  = CFG::Arith::coarsen;
        progressbar = progress;
    }

};

/////////////////////////////////////////////////////////////////////////////////////
//
// matrix inversion functions
//
/////////////////////////////////////////////////////////////////////////////////////

//!
//! \ingroup Algebra_Module
//! \brief   compute \f$ A^{-1} \f$
//!
//!          Compute inverse of \f$ A \f$ with block-wise accuracy \a acc.
//!          \a A is overwritten by \f$ A^{-1} \f$ during inversion.
//!
void
invert ( TMatrix *              A,
         const TTruncAcc &      acc,
         const inv_options_t &  opts = inv_options_t() );

//!
//! \ingroup Algebra_Module
//! \brief   generate DAG for computing \f$ A^{-1} \f$
//!
DAG::Graph
gen_dag_invert ( TMatrix *              A,
                 const inv_options_t &  opts = inv_options_t() );

//!
//! \ingroup Algebra_Module
//! \brief   compute \f$ A^{-1} \f$
//!
//!          Compute inverse of lower-left triangular matrix \f$ A \f$ with 
//!          block-wise accuracy \a acc. \a A is overwritten by \f$ A^{-1} \f$
//!          during inversion.
//!
void
invert_ll     ( TMatrix *             A,
                const TTruncAcc &     acc,
                const inv_options_t & opts = inv_options_t() );

//!
//! same as above but enforce recursive algorithm
//!
void
invert_ll_rec ( TMatrix *             A,
                const TTruncAcc &     acc,
                const inv_options_t & opts = inv_options_t() );

//!
//! same as above but enforce DAG algorithm
//!
void
invert_ll_dag ( TMatrix *             A,
                const TTruncAcc &     acc,
                const inv_options_t & opts = inv_options_t() );

//!
//! \ingroup Algebra_Module
//! \brief   generate DAG for computing \f$ A^{-1} \f$ for lower-left part of \f$A\f$
//!
DAG::Graph
gen_dag_invert_ll ( TMatrix *              A,
                    const inv_options_t &  opts = inv_options_t() );

//!
//! \ingroup Algebra_Module
//! \brief   compute \f$ X = L^{-1} \f$
//!
//!          Compute inverse of lower-left triangular matrix \f$ L \f$ with 
//!          block-wise accuracy \a acc. \a L may be modified during inversion.
//!
void
invert_ll ( TMatrix *             L,
            TMatrix *             X,
            const TTruncAcc &     acc,
            const inv_options_t & opts = inv_options_t() );

//!
//! \fn      invert_ll_steps
//! \brief   report number of steps in \see invert_ll for progress meter initialisation
//!
size_t
invert_ll_steps ( const TMatrix * A );
    
//!
//! \ingroup Algebra_Module
//! \brief   compute \f$ A^{-1} \f$
//!
//!          Compute inverse of upper-right triangular matrix \f$ A \f$ with 
//!          block-wise accuracy \a acc. \a A is overwritten by \f$ A^{-1} \f$
//!          during inversion.
//!
void
invert_ur     ( TMatrix *             A,
                const TTruncAcc &     acc,
                const inv_options_t & opts = inv_options_t() );

//!
//! same as above but enforce recursive algorithm
//!
void
invert_ur_rec ( TMatrix *             A,
                const TTruncAcc &     acc,
                const inv_options_t & opts = inv_options_t() );

//!
//! same as above but enforce DAG algorithm
//!
void
invert_ur_dag ( TMatrix *             A,
                const TTruncAcc &     acc,
                const inv_options_t & opts = inv_options_t() );

//!
//! \ingroup Algebra_Module
//! \brief   generate DAG for computing \f$ A^{-1} \f$ for upper-right part of \f$A\f$
//!
DAG::Graph
gen_dag_invert_ur ( TMatrix *              A,
                    const inv_options_t &  opts = inv_options_t() );

//!
//! \fn      invert_ur_steps
//! \brief   report number of steps in \see invert_ur for progress meter initialisation
//!
size_t
invert_ur_steps ( const TMatrix * A );
    
//!
//! \ingroup Algebra_Module
//! \brief   compute \f$ A^{-1} \f$
//!
//!          Compute inverse of diagonal matrix \f$ A \f$ with block-wise
//!          accuracy \a acc. \a A is overwritten by \f$ A^{-1} \f$
//!          during inversion.
//!
void
invert_diag ( TMatrix *             A,
              const TTruncAcc &     acc,
              const inv_options_t & opts = inv_options_t() );

//!
//! \fn      invert_diag_steps
//! \brief   report number of steps in \see invert_diag for progress meter initialisation
//!
size_t
invert_diag_steps ( const TMatrix * A );
    
//!
//! \ingroup Algebra_Module
//! \brief   Compute diagonal of \f$ A^{-1} \f$
//!
//!          Compute only the diagonal of \a A and return the resulting
//!          vector. \a A is modified during computation.
//!
TVector *
inverse_diag ( TMatrix *              A,
               const TTruncAcc &      acc,
               const inv_options_t &  opts = inv_options_t() );

//!
//! \ingroup Algebra_Module
//! \brief   Compute \f$A^{-1} \f$ via Gaussian elimination
//!
//!          \a A is overwritten by \f$A^{-1} \f$. \a C is optional
//!          and may be used during computation as temporary space.
//!          If present, it should have same format as \a A.
//!
void
gauss_elim   ( TMatrix *              A,
               TMatrix *              C,
               const TTruncAcc &      acc,
               const inv_options_t &  opts = inv_options_t() );

//!
//! \ingroup Algebra_Module
//! \brief   Compute \f$A^{-1} \f$ via Gaussian elimination
//!
//!          \a A is overwritten by \f$A^{-1} \f$. The elimination
//!          is performed in-place, although local copy operations
//!          are needed. (ATTENTION: experimental)
//!
void
gauss_elim   ( TMatrix *              A,
               const TTruncAcc &      acc,
               const inv_options_t &  opts = inv_options_t() );

//!
//! \ingroup Algebra_Module
//! \brief   compute \f$ A^{-1} \f$
//!
//!          Compute inverse of \f$ A \f$ (dense matrix version).
//!          \a A is overwritten by \f$ A^{-1} \f$ during inversion.
//!
void
invert       ( TDenseMatrix *         A );

//!
//! \ingroup Algebra_Module
//! \brief   compute and return \f$ A^{-1} \f$
//!
//!          Compute inverse of \f$ A \f$ (dense matrix version)
//!          but do not modify \f$ A \f$.
//!
std::unique_ptr< TDenseMatrix >
inverse      ( const TDenseMatrix *   A );

}// namespace HLIB

#endif  // __HLIB_MAT_INV_HH
