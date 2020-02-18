#ifndef __HLIB_TMATRIXSTAT_HH
#define __HLIB_TMATRIXSTAT_HH
//
// Project     : HLib
// File        : TMatrixStat.hh
// Description : class for obtaining statistical information about a matrix
// Author      : Ronald Kriemann
// Copyright   : Max Planck Institute MIS 2004-2018. All Rights Reserved.
//

#include <vector>

#include "hpro/matrix/TMatrix.hh"
#include "hpro/matrix/TBlockMatrix.hh"
#include "hpro/matrix/TRkMatrix.hh"

namespace HLIB
{

//
// computes various statistical information about the matrix
//
class TMatrixStat
{
public:
    ////////////////////////////////////////////////////
    //
    // constructor and destructor
    //

    TMatrixStat  () {}
    ~TMatrixStat () {}
    
    ////////////////////////////////////////////////////
    //
    // compute rank distribution of low-rank
    // blocks in matrix
    //

    // count how many submatrices have a certain rank
    // - rank[i] contains the number of blocks with rank i
    void    rank_distr  ( const TMatrix *        A,
                          std::vector< uint > &  ranks ) const;

    // return number of coefficients in matrix
    size_t  ncoeff      ( const TMatrix *        A ) const;
    
private:
    ////////////////////////////////////////////////////
    //
    // special methods for rank distribution
    //

    void rank_distr_gen   ( const TMatrix *        A,
                            std::vector< uint > &  ranks ) const;
    
    void rank_distr_block ( const TBlockMatrix *   A,
                            std::vector< uint > &  ranks ) const;

    void rank_distr_rank  ( const TRkMatrix *      A,
                            std::vector< uint > &  ranks ) const;
};

}// namespace

#endif // __TMATRIXSTAT_HH
