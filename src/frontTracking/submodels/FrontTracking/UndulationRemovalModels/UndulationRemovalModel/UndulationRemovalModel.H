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
    Foam::UndulationRemovalModel

Description
    Base class for collisional return-to-undulationRemoval models.

SourceFiles
    UndulationRemovalModel.C
    UndulationRemovalModelNew.C

\*---------------------------------------------------------------------------*/

#ifndef UndulationRemovalModel_H
#define UndulationRemovalModel_H

#include "IOdictionary.H"
#include "autoPtr.H"
#include "runTimeSelectionTables.H"
#include "CloudSubModelBase.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class UndulationRemovalModel Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType>
class UndulationRemovalModel
:
    public CloudSubModelBase<CloudType>
{
public:

    //- Runtime type information
    TypeName("undulationRemovalModel");

	label undulationRemovalInterval_;

	label URIntervalCounter_;

	label URRepeatNum_;

    //- Declare runtime constructor selection table
    declareRunTimeSelectionTable
    (
        autoPtr,
        UndulationRemovalModel,
        dictionary,
        (
            const dictionary& dict,
            CloudType& owner
        ),
        (dict, owner)
    );


    // Constructors

        //- Construct null from owner
        UndulationRemovalModel(CloudType& owner);

        //- Construct from components
        UndulationRemovalModel
        (
            const dictionary& dict,
            CloudType& owner,
            const word& type
        );

        //- Construct a copy
        UndulationRemovalModel(const UndulationRemovalModel<CloudType>& cm);

        //- Construct and return a clone
        virtual autoPtr<UndulationRemovalModel<CloudType>> clone() const = 0;


    //- Destructor
    virtual ~UndulationRemovalModel();


    //- Selector
    static autoPtr<UndulationRemovalModel<CloudType>> New
    (
        const dictionary& dict,
        CloudType& owner
    );


    //- Member Functions

        //- Calculate velocities
        virtual void calculate() = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define makeUndulationRemovalModel(CloudType)                                           \
                                                                               \
    typedef Foam::CloudType::FrontTrackingCloudType FrontTrackingCloudType;                    \
    defineNamedTemplateTypeNameAndDebug                                        \
    (                                                                          \
        Foam::UndulationRemovalModel<FrontTrackingCloudType>,                                   \
        0                                                                      \
    );                                                                         \
                                                                               \
    namespace Foam                                                             \
    {                                                                          \
        defineTemplateRunTimeSelectionTable                                    \
        (                                                                      \
            UndulationRemovalModel<FrontTrackingCloudType>,                                     \
            dictionary                                                         \
        );                                                                     \
    }


#define makeUndulationRemovalModelType(SS, CloudType)                                   \
                                                                               \
    typedef Foam::CloudType::FrontTrackingCloudType FrontTrackingCloudType;                    \
    defineNamedTemplateTypeNameAndDebug                                        \
        (Foam::UndulationRemovalModels::SS<FrontTrackingCloudType>, 0);                         \
                                                                               \
    Foam::UndulationRemovalModel<FrontTrackingCloudType>::                                      \
        adddictionaryConstructorToTable                                        \
        <Foam::UndulationRemovalModels::SS<FrontTrackingCloudType>>                             \
            add##SS##CloudType##FrontTrackingCloudType##ConstructorToTable_;


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "UndulationRemovalModel.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
