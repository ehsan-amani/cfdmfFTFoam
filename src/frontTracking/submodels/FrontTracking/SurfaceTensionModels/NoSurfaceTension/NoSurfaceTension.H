/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2013-2020 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Class
    Foam::SurfaceTensionModels::NoSurfaceTension

Description

SourceFiles
    NoSurfaceTension.C

\*---------------------------------------------------------------------------*/

#ifndef NoSurfaceTension_H
#define NoSurfaceTension_H

#include "SurfaceTensionModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace SurfaceTensionModels
{

/*---------------------------------------------------------------------------*\
                         Class NoSurfaceTension Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType>
class NoSurfaceTension
:
    public SurfaceTensionModel<CloudType>
{
public:

    //- Runtime type information
    TypeName("none");

    // Constructors

        //- Construct from components
        NoSurfaceTension(const dictionary& dict, CloudType& owner);

        //- Construct copy
        NoSurfaceTension(const NoSurfaceTension<CloudType>& cm);

        //- Construct and return a clone
        virtual autoPtr<SurfaceTensionModel<CloudType>> clone() const
        {
            return autoPtr<SurfaceTensionModel<CloudType>>
            (
                new NoSurfaceTension<CloudType>(*this)
            );
        }


    //- Destructor
    virtual ~NoSurfaceTension();


    //- Member Functions

        //- Calculate velocities
        virtual void calculate();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace SurfaceTensionModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "NoSurfaceTension.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
