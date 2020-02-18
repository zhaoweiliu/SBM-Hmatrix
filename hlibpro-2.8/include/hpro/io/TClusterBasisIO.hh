#ifndef __HLIB_TCLUSTERBASISIO_HH
#define __HLIB_TCLUSTERBASISIO_HH
//
// Project     : HLib
// File        : TClusterBasisIO.hh
// Description : classes for visualisation of cluster basis
// Author      : Ronald Kriemann
// Copyright   : Max Planck Institute MIS 2004-2018. All Rights Reserved.
//

#include "hpro/cluster/TClusterBasis.hh"

namespace HLIB
{

//!
//! \ingroup  IO_Module
//! \class    TClusterBasisIO
//! \brief    base class for cluster basis input/output
//!
template <typename T>
class TClusterBasisIO
{
public:
    ///////////////////////////////////////////////
    //
    // visualisation function
    //

    //! write cluster basis \a cb to file \a filename
    virtual void
    write ( const TClusterBasis< T > *  cb,
            const std::string &         filename ) const = 0;

    //! read cluster basis from file \a filename
    virtual std::unique_ptr< TClusterBasis< T > >
    read  ( const std::string &         filename ) const = 0;
};

/////////////////////////////////////////////////////////////////
//
// HLIBpro file format
//
/////////////////////////////////////////////////////////////////

//
//! \ingroup  IO_Module
//! \class    THLibClusterBasisIO
//! \brief    cluster basis I/O in HLIBpro file format
//
template <typename T>
class THLibClusterBasisIO : public TClusterBasisIO< T >
{
public:
    //! write cluster basis \a cb to file \a filename
    virtual void
    write ( const TClusterBasis< T > *  cb,
            const std::string &         filename ) const;

    //! read cluster basis from file \a filename
    virtual std::unique_ptr< TClusterBasis< T > >
    read  ( const std::string &         filename ) const;
};

}// namespace HLIB

#endif  // __TCLUSTERBASISVIS_HH
