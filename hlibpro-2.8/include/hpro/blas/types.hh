#ifndef __HLIB_BLAS_TYPES_HH
#define __HLIB_BLAS_TYPES_HH
//
// Project     : HLib
// File        : types.hh
// Description : provide basic types for BLAS operations
// Author      : Ronald Kriemann
// Copyright   : Max Planck Institute MIS 2004-2018. All Rights Reserved.
//

namespace HLIB
{

//! \enum  matop_t
//! \brief matrix transformations, e.g. transposed or conjugate transpose
enum matop_t
{
    //! do not change matrix
    MATOP_NORM       = 'N',
    apply_normal     = 'N',

    //! use transposed matrix
    MATOP_TRANS      = 'T',
    apply_trans      = 'T',
    apply_transposed = 'T',

    //! use adjoint, e.g. conjugate transposed matrix
    MATOP_ADJ        = 'C',        
    MATOP_CONJTRANS  = 'C',        
    apply_adj        = 'C',
    apply_adjoint    = 'C',
    apply_conjtrans  = 'C'        
};

//! \enum   approx_t
//! \brief  different types of low-rank approximation methods
enum approx_t : int
{
    use_svd      = 0,
    use_rrqr     = 1,
    use_rand     = 2,
    use_randsvd  = 2,
    use_randlr   = 3,
    use_aca      = 4,
    use_svd_pair = 5
};

namespace BLAS
{

//! \enum   diag_type_t
//! \brief  defines whether diagonal is unit or non-unit
enum diag_type_t
{
    unit_diag        = 'U',
    general_diag     = 'N'
};

//! \enum   tri_type_t
//! \brief  defines whether matrix is upper or lower triangular
enum tri_type_t
{
    lower_triangular = 'L',
    upper_triangular = 'U'
};

//! \enum   eval_side_t
//! \brief  defines whether triangular system is multiplied from left or right
enum eval_side_t
{
    from_left        = 'L',
    from_right       = 'R'
};

}// namespace BLAS

}// namespace HLIB

#endif  // __HLIB_BLAS_TYPES_HH
