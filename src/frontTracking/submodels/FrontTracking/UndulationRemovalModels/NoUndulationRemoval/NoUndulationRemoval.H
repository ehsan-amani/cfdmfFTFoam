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
    Foam::UndulationRemovalModels::NoUndulationRemoval

Description

SourceFiles
    NoUndulationRemoval.C

\*---------------------------------------------------------------------------*/

#ifndef NoUndulationRemoval_H
#define NoUndulationRemoval_H

#include "UndulationRemovalModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace UndulationRemovalModels
{

/*---------------------------------------------------------------------------*\
                         Class NoUndulationRemoval Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType>
class NoUndulationRemoval
:
    public UndulationRemovalModel<CloudType>
{
public:

    //- Runtime type information
    TypeName("none");

    // Constructors

        //- Construct from components
        NoUndulationRemoval(const dictionary& dict, CloudType& owner);

        //- Construct copy
        NoUndulationRemoval(const NoUndulationRemoval<CloudType>& cm);

        //- Construct and return a clone
        virtual autoPtr<UndulationRemovalModel<CloudType>> clone() const
        {
            return autoPtr<UndulationRemovalModel<CloudType>>
            (
                new NoUndulationRemoval<CloudType>(*this)
            );
        }


    //- Destructor
    virtual ~NoUndulationRemoval();


    //- Member Functions

        //- Calculate velocities
        virtual void calculate();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace UndulationRemovalModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "NoUndulationRemoval.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
