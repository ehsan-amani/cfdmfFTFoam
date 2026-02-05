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
    Foam::UndulationRemovalModels::VCSIV

Description
    VCSIV return-to-undulationRemoval model.
	// note: function to smooth the front mesh, algorithm IV.
	// Volume Conserving Smoothing for Piecewise Linear Curves, Surfaces, and Triple Lines
	//A. Kuprat et al, Journal of Computational Physics 172, 99–118 (2001)

SourceFiles
    VCSIV.C

\*---------------------------------------------------------------------------*/

#ifndef VCSIV_H
#define VCSIV_H

#include "UndulationRemovalModel.H"
#include "bubbleData.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace UndulationRemovalModels
{

/*---------------------------------------------------------------------------*\
                         Class VCSIV Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType>
class VCSIV
:
    public UndulationRemovalModel<CloudType>
{
    // Private Member Functions


public:

    //- Runtime type information
    TypeName("vcsiv");

	Switch combined_;


    // Constructors

        //- Construct from components
        VCSIV(const dictionary& dict, CloudType& owner);

        //- Construct copy
        VCSIV(const VCSIV<CloudType>& cm);

        //- Construct and return a clone
        virtual autoPtr<UndulationRemovalModel<CloudType>> clone() const
        {
            return autoPtr<UndulationRemovalModel<CloudType>>
            (
                new VCSIV<CloudType>(*this)
            );
        }


    //- Destructor
    virtual ~VCSIV();


    //- Member Functions

        //- Calculate velocities
        virtual void calculate();

		//inline vector smoothTo(pointData<CloudType> currentPt);
		inline vector smoothTo(pointData currentPt);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace UndulationRemovalModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "VCSIV.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
