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
    Foam::FrontToFieldModels::Poisson

Description
    Poisson return-to-frontToField model.

SourceFiles
    Poisson.C

\*---------------------------------------------------------------------------*/

#ifndef Poisson_H
#define Poisson_H

#include "FrontToFieldModel.H"
#include "BaseF2FModel.H"
#include "bubbleData.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace FrontToFieldModels
{

/*---------------------------------------------------------------------------*\
                         Class Poisson Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType>
class Poisson
:
    public BaseF2FModel<CloudType>
{
    // Private Member Functions


public:

    //- Runtime type information
    TypeName("poisson");

	autoPtr<volVectorField> IndicatorGrad_;

	autoPtr<volVectorField> densityIGrad_;

	autoPtr<volVectorField> viscosityIGrad_;

    // Constructors

        //- Construct from components
        Poisson(const dictionary& dict, CloudType& owner);

        //- Construct copy
        Poisson(const Poisson<CloudType>& cm);

        //- Construct and return a clone
        virtual autoPtr<FrontToFieldModel<CloudType>> clone() const 
        {
            return autoPtr<FrontToFieldModel<CloudType>>
            (
                new Poisson<CloudType>(*this)
            );
        }


    //- Destructor
    virtual ~Poisson();


    //- Member Functions

        //- Calculate velocities
		virtual void IndicatorsConstruction();

		virtual void collectFields(label cellI, label bDI, scalar wBar, vector elAreaVec); 

		virtual void initFields();

/*
		//EA add- return the density indicator from the frontDynamic class.
        inline volScalarField& IndicatorGradRef();

        //EA add- return the density indicator from the frontDynamic class.
        inline tmp<volScalarField> IndicatorGrad() const;

		//EA add- return the density indicator from the frontDynamic class.
        inline volScalarField& densityIGradRef();

        //EA add- return the density indicator from the frontDynamic class.
        inline tmp<volScalarField> densityIGrad() const;

		//EA add- return the density indicator from the frontDynamic class.
        inline volScalarField& viscosityIGradRef();

        //EA add- return the density indicator from the frontDynamic class.
        inline tmp<volScalarField> viscosityIGrad() const;
*/
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace FrontToFieldModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "Poisson.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
