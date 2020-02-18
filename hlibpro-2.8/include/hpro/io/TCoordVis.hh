#ifndef __HLIB_TCOORDVIS_HH
#define __HLIB_TCOORDVIS_HH
//
// Project     : HLib
// File        : TCoordVis.hh
// Description : coordinate visualisation classes
// Author      : Ronald Kriemann
// Copyright   : Max Planck Institute MIS 2004-2018. All Rights Reserved.
//

#include <vector>

#include "hpro/cluster/TCoordinate.hh"
#include "hpro/cluster/TPermutation.hh"
#include "hpro/cluster/TCluster.hh"

#include "hpro/matrix/TSparseMatrix.hh"

#include "hpro/vector/TScalarVector.hh"

namespace HLIB
{
    
///////////////////////////////////////////////////////////////
//!
//! \ingroup  IO_Module
//! \class    TCoordVis
//! \brief    Base class for coordinate visualisation.
//!
class TCoordVis
{
public:
    //////////////////////////////////////
    //
    // constructor and destructor
    //

    TCoordVis () {}

    virtual ~TCoordVis () {}

    //////////////////////////////////////
    //
    // print coord
    //

    //! print \a coord to file \a filename
    virtual void print ( const TCoordinate &  coord,
                         const std::string &  filename ) const = 0;
};

///////////////////////////////////////////////////////////////
//!
//! \ingroup  IO_Module
//! \class    TPSCoordVis
//! \brief    Coordinate visualisation in PostScript format
//!
class TPSCoordVis
{
private:
    // viewing direction
    T3Point  _view;
    
public:
    //////////////////////////////////////
    //
    // constructor and destructor
    //

    TPSCoordVis ();
    TPSCoordVis ( const T3Point view_dir );

    virtual ~TPSCoordVis () {}

    //////////////////////////////////////
    //
    // print coord
    //
    
    //! print \a coord to file \a filename
    virtual void print ( const TCoordinate *   coord,
                         const std::string &   filename ) const;

    //! print coordinates \a coord and choose coordinate color based on label
    virtual void print ( const TCoordinate *          coord,
                         const std::vector< uint > &  label,
                         const std::string &          filename ) const;

    //! print only coordinates \a coord with filter(i)=true
    virtual void print ( const TCoordinate *          coord,
                         const std::vector< uint > &  label,
                         const std::vector< bool > &  filter,
                         const std::string &          filename ) const;
};

//
//! \ingroup  IO_Module
//! \fn       print_ps
//! \brief    functional version of TPSCoordVis
//
void
print_ps ( const TCoordinate *  coord,
           const std::string &  filename );

///////////////////////////////////////////////////////////////
//!
//! \ingroup  IO_Module
//! \class    TVRMLCoordVis
//! \brief    Coordinate visualisation in VRML format
//!
class TVRMLCoordVis
{
public:
    //////////////////////////////////////
    //
    // constructor and destructor
    //

    TVRMLCoordVis () {}

    virtual ~TVRMLCoordVis () {}

    //////////////////////////////////////
    //
    // print coord
    //
    
    //! print \a coord to file \a filename
    virtual void print ( const TCoordinate *   coord,
                         const std::string &   filename ) const;

    //! print coordinates \a coo and mark those in cluster \a cl
    //! - conversion of indices from coordinate numbering to
    //!   cluster numbering is performed by \a perm
    virtual void print ( const TCoordinate *          coord,
                         const std::vector< uint > &  label,
                         const std::string &          filename ) const;

    //! same as \see print, but also print connectivity between indices
    //! defined by sparse matrix \a S
    //! - ordering in matrix MUST be the same as in coordinates 
    virtual void print ( const TCoordinate *   coord,
                         const TBlockCluster * cluster,
                         const TSparseMatrix * S,
                         const TPermutation *  perm,
                         const std::string &   filename ) const;
};

//
//! \ingroup  IO_Module
//! \fn       print_vrml
//! \brief    functional version of TVRMLCoordVis
//
void
print_vrml ( const TCoordinate *  coord,
             const std::string &  filename );

///////////////////////////////////////////////////////////////
//!
//! \ingroup  IO_Module
//! \class    TVTKCoordVis
//! \brief    Coordinate visualisation in VTK format
//!
class TVTKCoordVis
{
public:
    //////////////////////////////////////
    //
    // colouring type for connectivity
    //

    //! colourisation options for drawing connectivity
    enum  coltype_t
    {
        //! do not colour connections
        NO_COLOUR,

        //! use absolute value of matrix coefficients
        //! to determine colour
        VAL_COLOUR,

        //! use logarithm of absolute value instead
        LOG_COLOUR   
    };
    
public:
    //////////////////////////////////////
    //
    // constructor and destructor
    //

    TVTKCoordVis () {}

    virtual ~TVTKCoordVis () {}

    //////////////////////////////////////
    //
    // print coord
    //
    
    //! print \a coord to file \a filename
    virtual void print ( const TCoordinate *    coord,
                         const std::string &    filename ) const;

    //! print coordinates \a coord and colour coordinates with same
    //! label with same colour; the labels are defined by \a label and
    //! must begin with 0
    virtual void print ( const TCoordinate *          coord,
                         const std::vector< uint > &  label,
                         const std::string &          filename ) const;

    //! print coordinates \a coord and the connections between
    //! them as defined by the coefficients of the sparse matrix \a S
    //! - the ordering of \a S is assumed to be equal to \a coord
    //! - if \a log_scale is true, coefficient scaling is logarithmic
    virtual void print ( const TCoordinate *    coord,
                         const TSparseMatrix *  S,
                         const std::string &    filename,
                         const coltype_t        col_type = TVTKCoordVis::NO_COLOUR ) const;
};

//
//! \ingroup  IO_Module
//! \fn       print_vtk
//! \brief    functional version of TVTKCoordVis
//
void
print_vtk ( const TCoordinate *  coord,
            const std::string &  filename );

void
print_vtk ( const TCoordinate *          coord,
            const std::vector< uint > &  label,
            const std::string &          filename );

}// namespace

#endif  // __HLIB_TCOORDVIS_HH
