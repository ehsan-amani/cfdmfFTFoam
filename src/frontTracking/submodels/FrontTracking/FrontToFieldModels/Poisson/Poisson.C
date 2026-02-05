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

\*---------------------------------------------------------------------------*/

#include "Poisson.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::FrontToFieldModels::Poisson<CloudType>::Poisson
(
    const dictionary& dict,
    CloudType& owner
)
:
    BaseF2FModel<CloudType>(dict, owner, typeName), //EA d2
	IndicatorGrad_(nullptr),
	densityIGrad_(nullptr),
	viscosityIGrad_(nullptr)
{
	if (this->twoFluidFlow_)
	{
		IndicatorGrad_.reset
		(
		    new volVectorField
		    (
		        IOobject
				(
				    "IndicatorGrad",
				    this->owner().db().time().timeName(),
				    this->owner().mesh(),
					IOobject::NO_READ, //rem IOobject::MUST_READ, 
					IOobject::AUTO_WRITE //rem IOobject::NO_WRITE
				),
				this->owner().mesh(),
				dimensionedVector
				(
				    "IndicatorGrad",
				    dimless/dimLength, //EA d2
				    vector::zero
				),
				zeroGradientFvPatchVectorField::typeName //rem
		    )
		);
	}
	else
	{
		densityIGrad_.reset
		(
		    new volVectorField
		    (
		        IOobject
				(
				    "densityIGrad",
				    this->owner().db().time().timeName(),
				    this->owner().mesh(),
					IOobject::NO_READ, //rem IOobject::MUST_READ, 
					IOobject::AUTO_WRITE //rem IOobject::NO_WRITE
				),
				this->owner().mesh(),
				dimensionedVector
				(
				    "densityIGrad",
				    this->owner().rho().dimensions()/dimLength,
				    vector::zero
				),
				zeroGradientFvPatchVectorField::typeName //rem
		    )
		);
		viscosityIGrad_.reset
		(
		    new volVectorField
		    (
		        IOobject
				(
				    "viscosityIGrad",
				    this->owner().db().time().timeName(),
				    this->owner().mesh(),
					IOobject::NO_READ, //rem IOobject::MUST_READ, 
					IOobject::AUTO_WRITE //rem IOobject::NO_WRITE
				),
				this->owner().mesh(),
				dimensionedVector
				(
				    "viscosityIGrad",
				    this->owner().mu().dimensions()/dimLength,
				    vector::zero
				),
				zeroGradientFvPatchVectorField::typeName //rem
		    )
		);
	}
}


template<class CloudType>
Foam::FrontToFieldModels::Poisson<CloudType>::Poisson
(
    const Poisson<CloudType>& cm
)
:
    BaseF2FModel<CloudType>(cm),
	IndicatorGrad_(cm.IndicatorGrad_),
	densityIGrad_(cm.densityIGrad_),
	viscosityIGrad_(cm.viscosityIGrad_)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class CloudType>
Foam::FrontToFieldModels::Poisson<CloudType>::~Poisson()
{}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType> 
void Foam::FrontToFieldModels::Poisson<CloudType>::IndicatorsConstruction()
{
    //const fvMesh& mesh = this->owner().mesh(); //EA d2
    //const pointField& ctrs = mesh.cellCentres(); //EA d2
	//DynamicList<bubbleData<CloudType>>& bDataL = this->owner().bDataL(); //EA d2

    if (this->twoFluidFlow_)
    {
		volScalarField& Indicator = this->IndicatorRef(); //EA d2
		volVectorField& IndicatorGrad = IndicatorGrad_(); //EA d2

        IndicatorGrad.correctBoundaryConditions(); //EA d2

        fvScalarMatrix IndicatorEqn
        (
            	fvm::laplacian(Indicator) //EA d2
      		==
           		fvc::div(IndicatorGrad) //EA d2
        );
        IndicatorEqn.solve();
/*
		//EA d2
       	forAll(Indicator, cI) 
        {
			scalar& Ind = Indicator[cI];
			if( Ind > 0.998 ) //EA4
			{
			    Ind = 1.0;
			}
			else if( Ind < 0.002 ) //EA4
			{
			    Ind = 0.0;
			}
        }
*/
        Indicator = max(Indicator, 0.0);
        Indicator = min(Indicator, 1.0); 
    }
    else //EA d2
    {
		volScalarField& densityIndicator = this->densityIndicatorRef(); //EA d2
		volScalarField& viscosityIndicator = this->viscosityIndicatorRef(); //EA d2
		volVectorField& densityIGrad = densityIGrad_(); //EA d2
		volVectorField& viscosityIGrad = viscosityIGrad_(); //EA d2

        densityIGrad.correctBoundaryConditions(); //EA d2
        viscosityIGrad.correctBoundaryConditions(); //EA d2

        fvScalarMatrix densityIEqn
        (
            	fvm::laplacian(densityIndicator) //EA d2
      		==
            	fvc::div(densityIGrad) //EA d2
        );
        densityIEqn.solve();

        fvScalarMatrix viscosityIEqn
        (
            	fvm::laplacian(viscosityIndicator) //EA d2
      		==
            	fvc::div(viscosityIGrad) //EA d2
        );
        viscosityIEqn.solve();

		densityIndicator = max(densityIndicator, 0.0); //EA d2
        viscosityIndicator = max(viscosityIndicator, 0.0); //EA d2
    }
}


template<class CloudType>
void Foam::FrontToFieldModels::Poisson<CloudType>::initFields()
{
	if (this->twoFluidFlow_)
	{
		IndicatorGrad_().primitiveFieldRef() = vector::zero;
	}
	else
	{
		densityIGrad_().primitiveFieldRef() = vector::zero;
		viscosityIGrad_().primitiveFieldRef() = vector::zero;
	}
}


template<class CloudType>
void Foam::FrontToFieldModels::Poisson<CloudType>::collectFields
(
	label cellI, 
	label bDI, 
	scalar wBar,
	vector elAreaVec
)
{
	//DynamicList<bubbleData<CloudType>>& bDataL_= this->owner().bDataL();
	DynamicList<bubbleData>& bDataL_= this->owner().bDataL();

	if (this->twoFluidFlow_)
	{
		IndicatorGrad_().primitiveFieldRef()[cellI] += (1.0-0.0) * wBar * elAreaVec; // /dVolume; //EA
	}
	else
	{
		densityIGrad_().primitiveFieldRef()[cellI] += (bDataL_[bDI].outerFluidDensity-bDataL_[bDI].density)
		                                       * wBar * elAreaVec; // /dVolume; //EA
		viscosityIGrad_().primitiveFieldRef()[cellI] += (bDataL_[bDI].outerFluidViscosity-bDataL_[bDI].viscosity)
		                                       * wBar * elAreaVec; // /dVolume; //EA
	}
}


// ************************************************************************* //
