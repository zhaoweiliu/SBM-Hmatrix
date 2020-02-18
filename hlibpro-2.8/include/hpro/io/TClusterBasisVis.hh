#ifndef __HLIB_TCLUSTERBASISVIS_HH
#define __HLIB_TCLUSTERBASISVIS_HH
//
// Project     : HLib
// File        : TClusterBasisVis.hh
// Description : classes for visualisation of cluster basis
// Author      : Ronald Kriemann
// Copyright   : Max Planck Institute MIS 2004-2018. All Rights Reserved.
//

#include "hpro/cluster/TClusterBasis.hh"

namespace HLIB
{

//!
//! \ingroup  IO_Module
//! \class    TClusterBasisVis
//! \brief    base class for cluster basis visualisation
//!
template <typename T>
class TClusterBasisVis
{
public:
    ///////////////////////////////////////////////
    //
    // visualisation function
    //

    //!
    //! visualise cluster basis \a cl
    //!
    virtual void visualise ( const TClusterBasis< T > *  cb ) const = 0;
};

/////////////////////////////////////////////////////////////////
//
// PostScript visualisation
//
/////////////////////////////////////////////////////////////////

//
//! \ingroup  IO_Module
//! \class    TPSClusterBasisVis
//! \brief    cluster basis visualisation in PostScript format
//
template <typename T>
class TPSClusterBasisVis : public TClusterBasisVis< T >
{
public:
    //! @cond
    
    //
    // visualisation options (see below)
    //
    struct option_t
    {
        bool  basis_rank;
        bool  colourise;
        bool  legend;

        option_t ();
    };
    
private:
    
    // output stream
    std::ostream *  _output;

    // options for visualisation
    option_t        _opt;
    
    //!@endcond
    
public:
    ///////////////////////////////////////////////
    //
    // constructor and destructor
    //

    //!
    //! construct visualisation
    //!
    TPSClusterBasisVis ();
    
    //!
    //! construct visualisation to stream \a output
    //!
    TPSClusterBasisVis ( std::ostream &  output );


    ///////////////////////////////////////////////
    //
    // option managment
    //
    
    //! turn on/off printing of rank of cluster bases
    TPSClusterBasisVis< T > &  basis_rank  ( const bool  b );
    
    //! turn on/off colourising tree according to rank compared to full rank
    TPSClusterBasisVis< T > &  colourise   ( const bool  b );
    
    //! turn on/off showing of legend (with colourmap)
    TPSClusterBasisVis< T > &  legend      ( const bool  b );

    ///////////////////////////////////////////////
    //
    // visualisation functions
    //

    //!
    //! visualise cluster basis \a cl
    //!
    virtual void visualise ( const TClusterBasis< T > *  cb ) const;

    //!
    //! visualise cluster basis \a cl to stream \a output
    //!
    virtual void visualise ( const TClusterBasis< T > *  cb,
                             std::ostream &              output ) const;
};

}// namespace HLIB

#endif  // __TCLUSTERBASISVIS_HH
