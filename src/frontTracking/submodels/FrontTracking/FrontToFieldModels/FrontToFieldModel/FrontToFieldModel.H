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
    Foam::FrontToFieldModel

Description
    Base class for collisional return-to-frontToField models.

SourceFiles
    FrontToFieldModel.C

\*---------------------------------------------------------------------------*/

#ifndef FrontToFieldModel_H
#define FrontToFieldModel_H

#include "IOdictionary.H"
#include "autoPtr.H"
#include "runTimeSelectionTables.H"
#include "CloudSubModelBase.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class FrontToFieldModel Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType>
class FrontToFieldModel
:
    public CloudSubModelBase<CloudType>
{
public:

    //- Runtime type information
    TypeName("frontToFieldModel");

	// operations for phisical margins to density and viscosity
    bool  twoFluidFlow_;

	//base length  scale 
	scalar h_;

	//length  scale factor
	scalar hFactor_; 

    volScalarField cellsToRefine_;

	autoPtr<volScalarField> Indicator_;

	autoPtr<volScalarField> densityIndicator_;

	autoPtr<volScalarField> viscosityIndicator_;

    //- Declare runtime constructor selection table
    declareRunTimeSelectionTable
    (
        autoPtr,
        FrontToFieldModel,
        dictionary,
        (
            const dictionary& dict,
            CloudType& owner
        ),
        (dict, owner)
    );


    // Constructors

        //- Construct null from owner
        FrontToFieldModel(CloudType& owner);

        //- Construct from components
        FrontToFieldModel
        (
            const dictionary& dict,
            CloudType& owner,
            const word& type
        );

        //- Construct a copy
        FrontToFieldModel(const FrontToFieldModel<CloudType>& cm);

        //- Construct and return a clone
        virtual autoPtr<FrontToFieldModel<CloudType>> clone() const = 0;


    //- Destructor
    virtual ~FrontToFieldModel();


    //- Selector
    static autoPtr<FrontToFieldModel<CloudType>> New
    (
        const dictionary& dict,
        CloudType& owner
    );


    //- Member Functions

        //- Calculate velocities
        virtual void calculate() = 0;

	//acess
		inline scalar bubbleHl(label bDI); //EA d2

		//- Return const-access to hFActor_
        inline scalar hFactor() const;

		//EA add- return the density indicator from the frontDynamic class.
        inline volScalarField& IndicatorRef();

        //EA add- return the density indicator from the frontDynamic class.
        inline tmp<volScalarField> Indicator() const;

		//EA add- return the density indicator from the frontDynamic class.
        inline volScalarField& densityIndicatorRef();

        //EA add- return the density indicator from the frontDynamic class.
        inline tmp<volScalarField> densityIndicator() const;

		//EA add- return the density indicator from the frontDynamic class.
        inline volScalarField& viscosityIndicatorRef();

        //EA add- return the density indicator from the frontDynamic class.
        inline tmp<volScalarField> viscosityIndicator() const;

		//EA add- return the density indicator from the frontDynamic class.
        inline bool twoFluidFlow() const; //EA d2

		//EA add- return cellsToRefine_.
        inline const volScalarField& cellsToRefine() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define makeFrontToFieldModel(CloudType)                                           \
                                                                               \
    typedef Foam::CloudType::FrontTrackingCloudType FrontTrackingCloudType;                    \
    defineNamedTemplateTypeNameAndDebug                                        \
    (                                                                          \
        Foam::FrontToFieldModel<FrontTrackingCloudType>,                                   \
        0                                                                      \
    );                                                                         \
                                                                               \
    namespace Foam                                                             \
    {                                                                          \
        defineTemplateRunTimeSelectionTable                                    \
        (                                                                      \
            FrontToFieldModel<FrontTrackingCloudType>,                                     \
            dictionary                                                         \
        );                                                                     \
    }


#define makeFrontToFieldModelType(SS, CloudType)                                   \
                                                                               \
    typedef Foam::CloudType::FrontTrackingCloudType FrontTrackingCloudType;                    \
    defineNamedTemplateTypeNameAndDebug                                        \
        (Foam::FrontToFieldModels::SS<FrontTrackingCloudType>, 0);                         \
                                                                               \
    Foam::FrontToFieldModel<FrontTrackingCloudType>::                                      \
        adddictionaryConstructorToTable                                        \
        <Foam::FrontToFieldModels::SS<FrontTrackingCloudType>>                             \
            add##SS##CloudType##FrontTrackingCloudType##ConstructorToTable_;


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "FrontToFieldModel.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
