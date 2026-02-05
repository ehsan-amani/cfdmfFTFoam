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

#include "FrontToFieldModel.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::FrontToFieldModel<CloudType>::FrontToFieldModel(CloudType& owner)
:
    CloudSubModelBase<CloudType>(owner),
	//EA add
	twoFluidFlow_(true), //EA d2
	h_(0),
	hFactor_(2.0),
	cellsToRefine_
	(
        IOobject
        (
            "cellsToRefine_",
            this->owner().db().time().timeName(),
            this->owner().mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        this->owner().mesh(),
        dimensionedScalar("cellsToRefine", dimless, 0.0)
    )
{}


template<class CloudType>
Foam::FrontToFieldModel<CloudType>::FrontToFieldModel
(
    const dictionary& dict,
    CloudType& owner,
    const word& type
)
:
    CloudSubModelBase<CloudType>(owner, dict, typeName, type),
	//EA add 
	twoFluidFlow_(this->coeffDict().lookup("twoFluidFlow")),
	h_(0),
	hFactor_(2.0),
	cellsToRefine_
	(
        IOobject
        (
            "cellsToRefine_",
            this->owner().db().time().timeName(),
            this->owner().mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        this->owner().mesh(),
        dimensionedScalar("cellsToRefine", dimless, 0.0)
    ),
	Indicator_(nullptr),
	densityIndicator_(nullptr),
	viscosityIndicator_(nullptr)
{
	if (this->twoFluidFlow_) //EA d2
	{
		Indicator_.reset
		(
		    new volScalarField 
		    (
		        IOobject
				(
				    "indicator",
				    this->owner().db().time().timeName(),
				    this->owner().mesh(),
					IOobject::MUST_READ, 
					IOobject::AUTO_WRITE 
				),
				this->owner().mesh()/*,
				dimensionedScalar 
				(
				    "Indicator",
				    dimless, //EA d2
				    Zero
				),
				zeroGradientFvPatchScalarField::typeName */
		    )
		);
	}
	else
	{
		densityIndicator_.reset
		(
		    new volScalarField
		    (
		        IOobject
				(
				    "densityIndicator",
				    this->owner().db().time().timeName(),
				    this->owner().mesh(),
					IOobject::MUST_READ, 
					IOobject::AUTO_WRITE 
				),
				this->owner().mesh()/*,
				dimensionedScalar 
				(
				    "densityIndicator",
				    this->owner().rho().dimensions(),
				    Zero
				),
				zeroGradientFvPatchScalarField::typeName */
		    )
		);
		viscosityIndicator_.reset
		(
		    new volScalarField
		    (
		        IOobject
				(
				    "viscosityIndicator",
				    this->owner().db().time().timeName(),
				    this->owner().mesh(),
					IOobject::MUST_READ, 
					IOobject::AUTO_WRITE 
				),
				this->owner().mesh()/*,
				dimensionedScalar
				(
				    "viscosityIndicator",
				    this->owner().mu().dimensions(),
				    Zero
				),
				zeroGradientFvPatchScalarField::typeName */
		    )
		);
	}
}


template<class CloudType>
Foam::FrontToFieldModel<CloudType>::FrontToFieldModel
(
    const FrontToFieldModel<CloudType>& cm
)
:
    CloudSubModelBase<CloudType>(cm),
	twoFluidFlow_(cm.twoFluidFlow_), 
	h_(cm.h_),
	hFactor_(cm.hFactor_),
	cellsToRefine_(cm.cellsToRefine_),
	Indicator_(cm.Indicator_),
	densityIndicator_(cm.densityIndicator_),
	viscosityIndicator_(cm.viscosityIndicator_)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class CloudType>
Foam::FrontToFieldModel<CloudType>::~FrontToFieldModel()
{}


// * * * * * * * * * * * * * * * *  Selector * * * * * * * * * * * * * * * * //

template<class CloudType>
Foam::autoPtr<Foam::FrontToFieldModel<CloudType>>
Foam::FrontToFieldModel<CloudType>::New
(
    const dictionary& dict,
    CloudType& owner
)
{
    word modelType(dict.lookup(typeName));

    Info<< "Selecting frontToField model " << modelType << endl;

    typename dictionaryConstructorTable::iterator cstrIter =
        dictionaryConstructorTablePtr_->find(modelType);

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        FatalErrorInFunction
            << "Unknown frontToField model type " << modelType
            << ", constructor not in hash table" << nl << nl
            << "    Valid frontToField model types are:" << nl
            << dictionaryConstructorTablePtr_->sortedToc()
            << exit(FatalError);
    }

    return
        autoPtr<FrontToFieldModel<CloudType>>
        (
            cstrIter()(dict, owner)
        );
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType> //EA d2
inline Foam::scalar Foam::FrontToFieldModel<CloudType>::bubbleHl(label bDI)
{
	return this->owner().lengthScaleNearTheFront()[bDI];
}

template<class CloudType>
inline Foam::scalar Foam::FrontToFieldModel<CloudType>::hFactor() const
{
    return hFactor_;
}

template<class CloudType>
inline Foam::tmp<Foam::volScalarField>
Foam::FrontToFieldModel<CloudType>::Indicator() const
{
    return Indicator_();
}

template<class CloudType>
inline Foam::volScalarField&
Foam::FrontToFieldModel<CloudType>::IndicatorRef()
{
    return Indicator_();
}

template<class CloudType>
inline Foam::tmp<Foam::volScalarField>
Foam::FrontToFieldModel<CloudType>::densityIndicator() const
{
    return densityIndicator_();
}

template<class CloudType>
inline Foam::volScalarField&
Foam::FrontToFieldModel<CloudType>::densityIndicatorRef()
{
    return densityIndicator_();
}

template<class CloudType>
inline Foam::tmp<Foam::volScalarField>
Foam::FrontToFieldModel<CloudType>::viscosityIndicator() const
{
    return viscosityIndicator_();
}

template<class CloudType>
inline Foam::volScalarField&
Foam::FrontToFieldModel<CloudType>::viscosityIndicatorRef()
{
    return viscosityIndicator_();
}

template<class CloudType> //EA d2
inline bool
Foam::FrontToFieldModel<CloudType>::twoFluidFlow() const
{
    return twoFluidFlow_;
}

template<class CloudType>
inline const Foam::volScalarField&
Foam::FrontToFieldModel<CloudType>::cellsToRefine() const
{
    return cellsToRefine_;
}

// ************************************************************************* //
