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
    Foam::SurfaceTensionModel

Description
    Base class for collisional return-to-surfaceTension models.

SourceFiles
    SurfaceTensionModel.C
    SurfaceTensionModelNew.C

\*---------------------------------------------------------------------------*/

#ifndef SurfaceTensionModel_H
#define SurfaceTensionModel_H

#include "IOdictionary.H"
#include "autoPtr.H"
#include "runTimeSelectionTables.H"
#include "CloudSubModelBase.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class SurfaceTensionModel Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType>
class SurfaceTensionModel
:
    public CloudSubModelBase<CloudType>
{
public:

    //- Runtime type information
    TypeName("surfaceTensionModel");

    //- Declare runtime constructor selection table
    declareRunTimeSelectionTable
    (
        autoPtr,
        SurfaceTensionModel,
        dictionary,
        (
            const dictionary& dict,
            CloudType& owner
        ),
        (dict, owner)
    );


    // Constructors

        //- Construct null from owner
        SurfaceTensionModel(CloudType& owner);

        //- Construct from components
        SurfaceTensionModel
        (
            const dictionary& dict,
            CloudType& owner,
            const word& type
        );

        //- Construct a copy
        SurfaceTensionModel(const SurfaceTensionModel<CloudType>& cm);

        //- Construct and return a clone
        virtual autoPtr<SurfaceTensionModel<CloudType>> clone() const = 0;


    //- Destructor
    virtual ~SurfaceTensionModel();


    //- Selector
    static autoPtr<SurfaceTensionModel<CloudType>> New
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

#define makeSurfaceTensionModel(CloudType)                                           \
                                                                               \
    typedef Foam::CloudType::FrontTrackingCloudType FrontTrackingCloudType;                    \
    defineNamedTemplateTypeNameAndDebug                                        \
    (                                                                          \
        Foam::SurfaceTensionModel<FrontTrackingCloudType>,                                   \
        0                                                                      \
    );                                                                         \
                                                                               \
    namespace Foam                                                             \
    {                                                                          \
        defineTemplateRunTimeSelectionTable                                    \
        (                                                                      \
            SurfaceTensionModel<FrontTrackingCloudType>,                                     \
            dictionary                                                         \
        );                                                                     \
    }


#define makeSurfaceTensionModelType(SS, CloudType)                                   \
                                                                               \
    typedef Foam::CloudType::FrontTrackingCloudType FrontTrackingCloudType;                    \
    defineNamedTemplateTypeNameAndDebug                                        \
        (Foam::SurfaceTensionModels::SS<FrontTrackingCloudType>, 0);                         \
                                                                               \
    Foam::SurfaceTensionModel<FrontTrackingCloudType>::                                      \
        adddictionaryConstructorToTable                                        \
        <Foam::SurfaceTensionModels::SS<FrontTrackingCloudType>>                             \
            add##SS##CloudType##FrontTrackingCloudType##ConstructorToTable_;


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "SurfaceTensionModel.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
