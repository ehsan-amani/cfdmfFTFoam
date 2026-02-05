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
    Foam::SurfaceTensionModels::DirectElementBased

Description
    DirectElementBased return-to-surfaceTension model.

SourceFiles
    DirectElementBased.C

\*---------------------------------------------------------------------------*/

#ifndef DirectElementBased_H
#define DirectElementBased_H

#include "SurfaceTensionModel.H"
#include "bubbleData.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace SurfaceTensionModels
{

/*---------------------------------------------------------------------------*\
                         Class DirectElementBased Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType>
class DirectElementBased
:
    public SurfaceTensionModel<CloudType>
{
    // Private Member Functions


public:

    //- Runtime type information
    TypeName("directElementBased");

    // Constructors

        //- Construct from components
        DirectElementBased(const dictionary& dict, CloudType& owner);

        //- Construct copy
        DirectElementBased(const DirectElementBased<CloudType>& cm);

        //- Construct and return a clone
        virtual autoPtr<SurfaceTensionModel<CloudType>> clone() const
        {
            return autoPtr<SurfaceTensionModel<CloudType>>
            (
                new DirectElementBased<CloudType>(*this)
            );
        }


    //- Destructor
    virtual ~DirectElementBased();


    //- Member Functions

        //- Calculate velocities
        virtual void calculate();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace SurfaceTensionModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "DirectElementBased.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
