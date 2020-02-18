#ifndef __HLIB_TIMER_HH
#define __HLIB_TIMER_HH
//
// Project     : HLib
// File        : TTimer.hh
// Description : class for measuring time
// Author      : Ronald Kriemann
// Copyright   : Max Planck Institute MIS 2004-2018. All Rights Reserved.
//

#include <ostream>
#include <string>

namespace HLIB
{

// type for different timer-measurment
enum timetype_t
{
    WALL_TIME,        // real clock
    CPU_TIME,         // cpu time of process
    CPU_TIME_THREAD   // cpu time of thread
};

//!
//! \class  TTimer
//! \brief  Timer class
//!
//!         Provides class for timing either using wall time or cpu time, i.e.
//!         actual time doing something without waiting for free resources.
//!
class TTimer 
{
protected:
    // start and overall time
    double      _start, _sum;

    // what kind of time we should stop
    timetype_t  _type;
    
public:
    /////////////////////////////////////////////////
    //
    // constructor
    //

    //! construct timer of type \a type
    TTimer ( const timetype_t type = CPU_TIME )
            : _start(0), _sum(0), _type(type)
    {}
    
    /////////////////////////////////////////////////
    //
    // compute times
    //
    
    //! start timer, resets time-sum
    TTimer &  start    ()       { _start = current_time(); _sum = 0.0; return *this; }

    //! pause timer, add elapsed time to sum
    TTimer &  pause    ()       { _sum  += current_time() - _start; _start = 0.0; return *this; }

    //! continue timer, set start time to current time
    TTimer &  cont     ()       { _start = current_time(); return *this; }

    //! returns elapsed time without stopping timer
    double    elapsed  () const { return (_start > 0 ? _sum + current_time() - _start : _sum); }

    /////////////////////////////////////////////////
    //
    // output
    //

    //! set type of time to wall clock time
    TTimer &  wall       () { _type = WALL_TIME;       return *this; }

    //! set type of time to wall cpu time
    TTimer &  cpu        () { _type = CPU_TIME;        return *this; }
    
    //! set type of time to wall cpu time
    TTimer &  cpu_thread () { _type = CPU_TIME_THREAD; return *this; }
    
    //! convert to string
    std::string  to_string () const;

    //! stream output
    friend std::ostream & operator << ( std::ostream & os, const TTimer & timer )
    {
        return os << timer.to_string();
    }

private:
    
    /////////////////////////////////////////////////
    //
    // misc functions
    //

    // get current time of system (cputime or wall clock)
    double  current_time  () const;
};

//!
//! \class  TCPUTimer
//! \brief  Timer class measuring CPU time
//!
class TCPUTimer : public TTimer
{
public:
    // ctor
    TCPUTimer ()
            : TTimer( CPU_TIME )
    {}
};

//!
//! \class  TThreadCPUTimer
//! \brief  Timer class measuring CPU time of current thread
//!
class TThreadCPUTimer : public TTimer
{
public:
    // ctor
    TThreadCPUTimer ()
            : TTimer( CPU_TIME_THREAD )
    {}
};

//!
//! \class  TWallTimer
//! \brief  Timer class measuring wall clock time
//!
class TWallTimer : public TTimer
{
public:
    // ctor
    TWallTimer ()
            : TTimer( WALL_TIME )
    {}
};

}// namespace HLIB

#endif  // __HLIB_TIMER_HH
