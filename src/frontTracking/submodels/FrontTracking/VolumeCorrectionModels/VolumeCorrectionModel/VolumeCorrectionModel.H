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
    Foam::VolumeCorrectionModel

Description
    Base class for collisional return-to-volumeCorrection models.

SourceFiles
    VolumeCorrectionModel.C
    VolumeCorrectionModelNew.C

\*---------------------------------------------------------------------------*/

#ifndef VolumeCorrectionModel_H
#define VolumeCorrectionModel_H

#include "IOdictionary.H"
#include "autoPtr.H"
#include "runTimeSelectionTables.H"
#include "CloudSubModelBase.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class VolumeCorrectionModel Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType>
class VolumeCorrectionModel
:
    public CloudSubModelBase<CloudType>
{
public:

    //- Runtime type information
    TypeName("volumeCorrectionModel");

    //- Declare runtime constructor selection table
    declareRunTimeSelectionTable
    (
        autoPtr,
        VolumeCorrectionModel,
        dictionary,
        (
            const dictionary& dict,
            CloudType& owner
        ),
        (dict, owner)
    );


    // Constructors

        //- Construct null from owner
        VolumeCorrectionModel(CloudType& owner);

        //- Construct from components
        VolumeCorrectionModel
        (
            const dictionary& dict,
            CloudType& owner,
            const word& type
        );

        //- Construct a copy
        VolumeCorrectionModel(const VolumeCorrectionModel<CloudType>& cm);

        //- Construct and return a clone
        virtual autoPtr<VolumeCorrectionModel<CloudType>> clone() const = 0;


    //- Destructor
    virtual ~VolumeCorrectionModel();


    //- Selector
    static autoPtr<VolumeCorrectionModel<CloudType>> New
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

#define makeVolumeCorrectionModel(CloudType)                                           \
                                                                               \
    typedef Foam::CloudType::FrontTrackingCloudType FrontTrackingCloudType;                    \
    defineNamedTemplateTypeNameAndDebug                                        \
    (                                                                          \
        Foam::VolumeCorrectionModel<FrontTrackingCloudType>,                                   \
        0                                                                      \
    );                                                                         \
                                                                               \
    namespace Foam                                                             \
    {                                                                          \
        defineTemplateRunTimeSelectionTable                                    \
        (                                                                      \
            VolumeCorrectionModel<FrontTrackingCloudType>,                                     \
            dictionary                                                         \
        );                                                                     \
    }


#define makeVolumeCorrectionModelType(SS, CloudType)                                   \
                                                                               \
    typedef Foam::CloudType::FrontTrackingCloudType FrontTrackingCloudType;                    \
    defineNamedTemplateTypeNameAndDebug                                        \
        (Foam::VolumeCorrectionModels::SS<FrontTrackingCloudType>, 0);                         \
                                                                               \
    Foam::VolumeCorrectionModel<FrontTrackingCloudType>::                                      \
        adddictionaryConstructorToTable                                        \
        <Foam::VolumeCorrectionModels::SS<FrontTrackingCloudType>>                             \
            add##SS##CloudType##FrontTrackingCloudType##ConstructorToTable_;


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "VolumeCorrectionModel.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
