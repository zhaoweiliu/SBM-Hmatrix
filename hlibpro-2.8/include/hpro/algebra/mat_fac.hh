#ifndef __HLIB_MAT_FAC_HH
#define __HLIB_MAT_FAC_HH
//
// Project     : HLib
// File        : mat_fac.hh
// Description : matrix factorisation functions
// Author      : Ronald Kriemann
// Copyright   : Max Planck Institute MIS 2004-2018. All Rights Reserved.
//

#include <tuple>

#include "hpro/matrix/TMatrix.hh"
#include "hpro/matrix/TFacMatrix.hh"
#include "hpro/matrix/TFacInvMatrix.hh"

#include "hpro/misc/TProgressBar.hh"

#include "hpro/algebra/solve_types.hh"

#include "hpro/parallel/dag.hh"

namespace HLIB
{

//!
//! \{
//! \name Matrix Factorisation
//!       Functions related to matrix factorisation, e.g. LU or Cholesky factorisation.
//!

/////////////////////////////////////////////////////////////////////////////////////
//
// option type for factorisation algorithms
//
/////////////////////////////////////////////////////////////////////////////////////

//!
//! \ingroup Algebra_Module
//! \struct  fac_options_t
//! \brief   options for matrix factorisations
//!
struct fac_options_t
{
    //
    // factorisation type
    //
        
    //! factorise block (default) or point wise
    eval_type_t     eval;

    //! store inverse diagonal blocks or leave unchanged
    storage_type_t  storage;
    
    //! if true, coarsening is applied during factorisation (default: off)
    bool            do_coarsen;
        
    //! progress bar for LU
    TProgressBar *  progressbar;

    //
    // ctor
    //
        
    //! default constructor (with simplification of progress bar usage)
    fac_options_t ( TProgressBar *  progress = nullptr )
    {
        eval           = CFG::Arith::eval_type;
        storage        = CFG::Arith::storage_type;
        do_coarsen     = CFG::Arith::coarsen;
        progressbar    = progress;
    }

    //! default constructor (with simplification of progress bar usage)
    fac_options_t ( const eval_type_t     aeval,
                    const storage_type_t  astorage,
                    const bool            ado_coarsen,
                    TProgressBar *        aprogress = nullptr )
            : eval( aeval )
            , storage( astorage )
            , do_coarsen( ado_coarsen )
            , progressbar( aprogress )
    {}
};

/////////////////////////////////////////////////////////////////////////////////////
//
// LU factorisation
//
/////////////////////////////////////////////////////////////////////////////////////

namespace LU
{
//!
//! \ingroup Algebra_Module
//! \class   TLU
//! \brief   Computes LU factorisation \f$ A = LU \f$
//! 
//!          This class computes the LU factorisation \f$A = LU\f$ of a matrix \f$A\f$
//!          with lower, unit triangular matrix \f$L\f$ and upper triangular matrix \f$U\f$.
//!
//!          The factorisation may be either point wise, i.e. a real LU factorisation,
//!          or block wise in which case, dense diagonal matrix blocks are inverted.
//!
//!          Support for multiple threads is available, although the expectable speedup is
//!          limited, e.g. best suited for at most 4 threads.
//!

//! compute LU factorisation of given matrix
//! \param  A          on input matrix to factorise; on output factors L and U
//! \param  acc        accuracy of factorisation
void
factorise    ( TMatrix *            A,
               const TTruncAcc &    acc,
               const fac_options_t  opts = fac_options_t() );

//! compute LU factorisation of given matrix using DAG algorithm
//! \param  A          on input matrix to factorise; on output factors L and U
//! \param  acc        accuracy of factorisation
void
factorise_dag ( TMatrix *            A,
                const TTruncAcc &    acc,
                const fac_options_t  opts = fac_options_t(),
                TMatrix *            D    = nullptr );

//! compute LU factorisation of given matrix using recursive algorithm
//! \param  A          on input matrix to factorise; on output factors L and U
//! \param  acc        accuracy of factorisation
void
factorise_rec ( TMatrix *            A,
                const TTruncAcc &    acc,
                const fac_options_t  opts = fac_options_t() );

//! return suitable representation for evaluating factors L and U
//! \param  A          LU factors to be represented
std::unique_ptr< TFacMatrix >
eval_matrix  ( TMatrix *            A,
               const fac_options_t  opts = fac_options_t() ) ;
    
//! return suitable inverse representation of factors L and U
//! \param  A          LU factors to be represented
std::unique_ptr< TFacInvMatrix >
inv_matrix   ( TMatrix *            A,
               const fac_options_t  opts = fac_options_t() );
    
//! split given matrix A into individual factors L and U
//! \param  A          joined LU factors
//! \param  L          matrix pointer to store factor L
//! \param  U          matrix pointer to store factor U
//! \param  eval_type  determines block-wise or point-wise evaluation
void
split        ( const TMatrix *      A,
               TMatrix * &          L,
               TMatrix * &          U,
               const fac_options_t  opts = fac_options_t() );

//! split matrix \a A into individual factors L and U; A is destroyed
//! and L/U will contain the former submatrices of A
std::pair< TMatrix *, TMatrix * >
split        ( TMatrix *            A,
               const fac_options_t  opts = fac_options_t() );

//! return number of factorisation steps for A for progress meter 
size_t
pm_steps     ( const TMatrix *      A,
               const fac_options_t  opts = fac_options_t() );


//
// compute DAG for LU(A)
//
DAG::Graph
gen_dag      ( TMatrix *              A,
               const fac_options_t &  options  = fac_options_t() );

}// namespace LU

/////////////////////////////////////////////////////////////////////////////////////
//
// LDU factorisation
//
/////////////////////////////////////////////////////////////////////////////////////

//!
//! \ingroup Algebra_Module
//! \class  TLDU
//! \brief  Computes LDU factorisation \f$ A = LDU \f$
//! 
//!          This class computes the LU factorisation \f$A = LDU\f$ of a matrix \f$A\f$
//!          with lower, unit triangular matrix \f$L\f$, upper triangular matrix \f$U\f$
//!          and diagonal matrix \f$D\f$.
//!
//!          The factorisation may be either point wise, i.e. a real LDU factorisation,
//!          or block wise in which case, dense diagonal matrix blocks are inverted.
//!
//!          Support for multiple threads is not available.
//!
class TLDU
{
private:
    //! local options for factorisation
    fac_options_t  _options;
    
public:
    //////////////////////////////////////////////
    //
    // constructor and destructor
    //

    //! standard constructor with optional LDU settings
    TLDU  ( const fac_options_t  opts = fac_options_t() )
            : _options( opts )
    {}
    
    //! destructor
    ~TLDU () {}

    //! set options for factorisation
    void set_options ( const fac_options_t  opts )
    {
        _options = opts;
    }
    
    //
    // compute LDU factorisation of given matrix
    // - A is overwritten with L,U
    //

    //! compute LDU factorisation of given matrix
    //! \param  A         on input matrix to factorise; on output factors L and U 
    //! \param  acc       accuracy of factorisation
    void factorise ( TMatrix *          A,
                     const TTruncAcc &  acc ) const;

    //
    // misc.
    //

    //! return suitable representation for evaluating factors L, D and U
    //! \param  A        LDU factors to be represented
    auto  eval_matrix  ( TMatrix *  A ) const -> std::unique_ptr< TFacMatrix >;
    
    //! return suitable inverse representation of factors L, D and U
    //! \param  A        LDU factor to be represented
    auto  inv_matrix   ( TMatrix *   A ) const -> std::unique_ptr< TFacInvMatrix >;
    
    //! split given matrix A into individual factors L, D and U
    //! \param  A          joined LDU factors
    //! \param  L          matrix pointer to store factor L
    //! \param  D          matrix pointer to store factor D
    //! \param  U          matrix pointer to store factor U
    //! \param  eval_type  determines block-wise or point-wise evaluation
    void       split       ( const TMatrix *    A,
                             TMatrix * &        L,
                             TMatrix * &        D,
                             TMatrix * &        U ) const;
    
    //! return number of factorisation steps for A for progress meter 
    size_t     pm_steps    ( const TMatrix * A ) const;

    DISABLE_COPY_OP( TLDU );
};

/////////////////////////////////////////////////////////////////////////////////////
//
// Cholesky factorisation
//
/////////////////////////////////////////////////////////////////////////////////////

namespace LL
{
//!
//! \ingroup Algebra_Module
//! \class  TLL
//! \brief  computes Cholesky factorisation \f$ A = LL^T \f$ or \f$ A=LL^H \f$ 
//! 
//!          This class computes the Cholesky factorisation \f$A = LL^T\f$ (\f$A = LL^H\f$)
//!          of a symmetric (hermitian) matrix \f$A\f$ with lower triangular matrix \f$L\f$.
//!
//!          Support for multiple threads is available, although the expectable speedup is
//!          limited, e.g. best suited for at most 4 threads.
//!

//! compute LL factorisation of given matrix
//! \param  A          on input matrix to factorise; on output factor L
//! \param  acc        accuracy of factorisation
void
factorise    ( TMatrix *            A,
               const TTruncAcc &    acc,
               const fac_options_t  opts = fac_options_t() );

//! compute LL factorisation of given matrix using DAG algorithm
//! \param  A          on input matrix to factorise; on output factor L
//! \param  acc        accuracy of factorisation
void
factorise_dag ( TMatrix *            A,
                const TTruncAcc &    acc,
                const fac_options_t  opts = fac_options_t() );

//! compute LL factorisation of given matrix using recursive algorithm
//! \param  A          on input matrix to factorise; on output factor L
//! \param  acc        accuracy of factorisation
void
factorise_rec ( TMatrix *            A,
                const TTruncAcc &    acc,
                const fac_options_t  opts = fac_options_t() );

//! return suitable representation for evaluating factor L
//! \param  A          LL factors to be represented
//! \param  matform    symmetric or hermitian
std::unique_ptr< TFacMatrix >
eval_matrix   ( TMatrix *            A,
                const matform_t      matform ) ;
    
//! return suitable inverse representation of factor L
//! \param  A          LL factors to be represented
//! \param  matform    symmetric or hermitian
std::unique_ptr< TFacInvMatrix >
inv_matrix    ( TMatrix *            A,
                const matform_t      matform );
    
//! return number of factorisation steps
size_t
pm_steps      ( const TMatrix *      A );

//
// compute DAG for LL(A)
//
DAG::Graph
gen_dag       ( TMatrix *              A,
                const fac_options_t &  options  = fac_options_t() );

}// namespace LL

/////////////////////////////////////////////////////////////////////////////////////
//
// LDL factorisation
//
/////////////////////////////////////////////////////////////////////////////////////

namespace LDL
{
//!
//! \ingroup Algebra_Module
//! \class   TLDL
//! \brief   computes LDL factorisation \f$ A = LDL^T \f$ or \f$ A = LDL^H \f$
//! 
//!          This class computes the LDL factorisation \f$A = LDL^T\f$ (\f$A = LDL^H\f$)
//!          of a symmetric (hermitian) matrix \f$A\f$ with lower, unit triangular matrix \f$L\f$
//!          and diagonal matrix \f$D\f$.
//!
//!          The factorisation may be either point wise, i.e. a real LU factorisation,
//!          or block wise in which case, dense diagonal matrix blocks are inverted.
//!
//!          Support for multiple threads is available, although the expectable speedup is
//!          limited, e.g. best suited for at most 4 threads.
//!

//! compute LDL factorisation of given matrix
//! \param  A          on input matrix to factorise; on output factors L and D
//! \param  acc        accuracy of factorisation
void
factorise     ( TMatrix *            A,
                const TTruncAcc &    acc,
                const fac_options_t  opts = fac_options_t() );

//! compute LDL factorisation of given matrix using DAG algorithm
//! \param  A          on input matrix to factorise; on output factors L and D
//! \param  acc        accuracy of factorisation
void
factorise_dag ( TMatrix *            A,
                const TTruncAcc &    acc,
                const fac_options_t  opts = fac_options_t(),
                TMatrix *            D    = nullptr,
                TMatrix *            DI   = nullptr );

//! compute LDL factorisation of given matrix using recursive algorithm
//! \param  A          on input matrix to factorise; on output factors L and D
//! \param  acc        accuracy of factorisation
void
factorise_rec ( TMatrix *            A,
                const TTruncAcc &    acc,
                const fac_options_t  opts = fac_options_t(),
                TMatrix *            D    = nullptr,
                TMatrix *            DI   = nullptr );

//! return suitable representation for evaluating factors L and D
//! \param  A          LDL factors to be represented
//! \param  matform    symmetric or hermitian
std::unique_ptr< TFacMatrix >
eval_matrix  ( TMatrix *            A,
               const matform_t      matform,
               const fac_options_t  opts = fac_options_t() ) ;
    
//! return suitable inverse representation of factors L and D
//! \param  A          LDL factors to be represented
//! \param  matform    symmetric or hermitian
std::unique_ptr< TFacInvMatrix >
inv_matrix   ( TMatrix *            A,
               const matform_t      matform,
               const fac_options_t  opts = fac_options_t() );
    
//! split given matrix A into individual factors L and D
//! \param  A          joined LDL factors
//! \param  L          matrix pointer to store factor L
//! \param  D          matrix pointer to store factor D
//! \param  eval_type  determines block-wise or point-wise evaluation
void
split        ( const TMatrix *      A,
               TMatrix * &          L,
               TMatrix * &          D,
               const fac_options_t  opts = fac_options_t() );

//! split matrix \a A into individual factors L and D; A is destroyed
//! and L/D will contain the former submatrices of A
std::pair< TMatrix *, TMatrix * >
split        ( TMatrix *            A,
               const fac_options_t  opts = fac_options_t() );

//! return number of factorisation steps for A for progress meter 
size_t
pm_steps     ( const TMatrix *      A,
               const fac_options_t  opts = fac_options_t() );

//
// compute DAG for LDL(A)
//
DAG::Graph
gen_dag      ( TMatrix *              A,
               TMatrix *              D,
               const fac_options_t &  options  = fac_options_t() );

}// namespace LDL

/////////////////////////////////////////////////////////////////////////////////////
//
// WAZ factorisation
//
/////////////////////////////////////////////////////////////////////////////////////

namespace WAZ
{

//! compute WAZ factorisation of given matrix
//! \param  A          on input matrix to factorise; on output destroyed
//! \param  acc        accuracy of factorisation
std::pair< TMatrix *, TMatrix * >
factorise      ( TMatrix *            A,
                 const TTruncAcc &    acc,
                 const fac_options_t  opts = fac_options_t() );

//! compute WAZ factorisation of given matrix using recursive algorithm
std::pair< TMatrix *, TMatrix * >
factorise_rec  ( TMatrix *            A,
                 const TTruncAcc &    acc,
                 const fac_options_t  opts = fac_options_t() );

//! compute WAZ factorisation of given matrix using DAG algorithm
std::pair< TMatrix *, TMatrix * >
factorise_dag  ( TMatrix *            A,
                 const TTruncAcc &    acc,
                 const fac_options_t  opts = fac_options_t() );

//! split matrix \a A into individual factors W and Z; A is destroyed
//! and L/U will contain the former submatrices of A
std::pair< TMatrix *, TMatrix * >
split          ( TMatrix *            A );

//! return number of factorisation steps for A for progress meter 
size_t
pm_steps       ( const TMatrix *      A );

//
// compute DAG for WAZ(A)
//
DAG::Graph
gen_dag        ( TMatrix *              A,
                 const fac_options_t &  options  = fac_options_t() );

}// namespace WAZ

/////////////////////////////////////////////////////////////////////////////////////
//
// functional versions of factorisation 
//
/////////////////////////////////////////////////////////////////////////////////////

//!
//! \ingroup Algebra_Module
//! \brief   compute factorisation of \a A
//!
//!          Compute triangular factorisation of \a A while choosing
//!          appropriate factorisation method depending on the format
//!          of \a A, e.g. if unsymmetric, symmetric or hermitian.
//!
//!          The return value is an operator object representing the
//!          factorised form and suitable for evaluation, e.g.
//!          matrix-vector multiplication (see TFacMatrix).
//!          \a A will be overwritten with the actual factorisation data.
//!
std::unique_ptr< TFacMatrix >
factorise     ( TMatrix *              A,
                const TTruncAcc &      acc,
                const fac_options_t &  options = fac_options_t() );

//!
//! \ingroup Algebra_Module
//! \brief   compute factorisation of \a A and return inverse operator
//!
//!          Compute triangular factorisation of \a A as in factorise
//!          but instead of an operator for evaluation of the factorised
//!          \a A, an operator for evaluation of the inverse of \a A is
//!          returned (see TFacInvMatrix).
//!
std::unique_ptr< TFacInvMatrix >
factorise_inv ( TMatrix *              A,
                const TTruncAcc &      acc,
                const fac_options_t &  options = fac_options_t() );

//!
//! \ingroup Algebra_Module
//! \brief   compute LU factorisation
//!
void
lu ( TMatrix *              A,
     const TTruncAcc &      acc,
     const fac_options_t &  options = fac_options_t() );

//!
//! \ingroup Algebra_Module
//! \brief   compute LU factorisation and return inverse operator
//!
std::unique_ptr< TFacInvMatrix >
lu_inv ( TMatrix *              A,
         const TTruncAcc &      acc,
         const fac_options_t &  options = fac_options_t() );

//!
//! \ingroup Algebra_Module
//! \brief   compute LDU factorisation
//!
void
ldu ( TMatrix *              A,
      const TTruncAcc &      acc,
      const fac_options_t &  options = fac_options_t() );

//!
//! \ingroup Algebra_Module
//! \brief   compute LDU factorisation and return inverse operator
//!
std::unique_ptr< TFacInvMatrix >
ldu_inv ( TMatrix *              A,
          const TTruncAcc &      acc,
          const fac_options_t &  options = fac_options_t() );

//!
//! \ingroup Algebra_Module
//! \brief   compute Cholesky factorisation
//!
void
ll ( TMatrix *              A,
     const TTruncAcc &      acc,
     const fac_options_t &  options = fac_options_t() );

inline
void
chol ( TMatrix *              A,
       const TTruncAcc &      acc,
       const fac_options_t &  options = fac_options_t() )
{
    ll( A, acc, options );
}

//!
//! \ingroup Algebra_Module
//! \brief   compute Cholesky factorisation and return inverse operator
//!
std::unique_ptr< TFacInvMatrix >
ll_inv ( TMatrix *              A,
         const TTruncAcc &      acc,
         const fac_options_t &  options = fac_options_t() );

inline
std::unique_ptr< TFacInvMatrix >
chol_inv ( TMatrix *              A,
           const TTruncAcc &      acc,
           const fac_options_t &  options = fac_options_t() )
{
    return ll_inv( A, acc, options );
}

//!
//! \ingroup Algebra_Module
//! \brief   compute LDL^H factorisation
//!
void
ldl ( TMatrix *              A,
      const TTruncAcc &      acc,
      const fac_options_t &  options = fac_options_t() );

//!
//! \ingroup Algebra_Module
//! \brief   compute LDL factorisation and return inverse operator
//!
std::unique_ptr< TFacInvMatrix >
ldl_inv ( TMatrix *              A,
          const TTruncAcc &      acc,
          const fac_options_t &  options = fac_options_t() );

//!
//! \ingroup Algebra_Module
//! \brief   test, if given matrix B_ij is lowrank with too large rank
//!          and convert to dense
//!
void
test_convert_dense ( TBlockMatrix *  B,
                     const uint      i,
                     const uint      j );

//! \}

}// namespace HLIB

#endif  // __HLIB_MAT_FAC_HH
