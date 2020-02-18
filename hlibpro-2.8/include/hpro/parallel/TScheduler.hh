#ifndef __HLIB_TSCHEDULER_HH
#define __HLIB_TSCHEDULER_HH
//
// Project     : HLib
// File        : TScheduler.hh
// Description : provides base class for scheduling algorithms
// Author      : Ronald Kriemann
// Copyright   : Max Planck Institute MIS 2004-2018. All Rights Reserved.
//

#include <vector>

#include "hpro/base/types.hh"

namespace HLIB
{

// types of scheduling
enum sched_t
{
    SCHED_LIST,
    SCHED_LPT,
    SCHED_BSSEQ,
    SCHED_SEQ,
    SCHED_MFIT
};

//
// baseclass for all scheduling algorithms
//
class TScheduler
{
public:
    ///////////////////////////////////////////////
    //
    // constructor
    //

    virtual ~TScheduler () {}
    
    ///////////////////////////////////////////////
    //
    // scheduling algorithms
    //

    //! Schedule entries in \a costs onto \a p processors and store
    //! assignment to processors in \a sched. Return true if scheduling
    //! was successful and false otherwise.
    virtual bool schedule ( const uint                     p,
                            std::vector< int > &           sched,
                            const std::vector< double > &  costs ) const = 0;
};

}// namespace

#endif  // __TSCHEDULER_HH
