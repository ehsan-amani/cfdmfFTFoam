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
    Foam::UndulationRemovalModels::TSUR3D

Description
    TSUR3D return-to-undulationRemoval model.
	// note: function to remove all front undulation based on TSUR method :
	// Trapezoidal Sub-grid Undulations Removal. A front-tracking/front-capturing
	// method for the simulation of 3D multi-fluid flows with free surfaces, F.S. de Sousa et al,
	// Journal of Computational Physics 198 (2004) 469–499

SourceFiles
    TSUR3D.C

\*---------------------------------------------------------------------------*/

#ifndef TSUR3D_H
#define TSUR3D_H

#include "UndulationRemovalModel.H"
#include "bubbleData.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace UndulationRemovalModels
{

/*---------------------------------------------------------------------------*\
                         Class TSUR3D Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType>
class TSUR3D
:
    public UndulationRemovalModel<CloudType>
{
    // Private Member Functions


public:

    //- Runtime type information
    TypeName("tsur3d");


    // Constructors

        //- Construct from components
        TSUR3D(const dictionary& dict, CloudType& owner);

        //- Construct copy
        TSUR3D(const TSUR3D<CloudType>& cm);

        //- Construct and return a clone
        virtual autoPtr<UndulationRemovalModel<CloudType>> clone() const
        {
            return autoPtr<UndulationRemovalModel<CloudType>>
            (
                new TSUR3D<CloudType>(*this)
            );
        }


    //- Destructor
    virtual ~TSUR3D();


    //- Member Functions

        //- Calculate velocities
        virtual void calculate();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace UndulationRemovalModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "TSUR3D.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
