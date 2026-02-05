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

#include "UndulationRemovalModel.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::UndulationRemovalModel<CloudType>::UndulationRemovalModel(CloudType& owner)
:
    CloudSubModelBase<CloudType>(owner),
	undulationRemovalInterval_(0),
	URIntervalCounter_(0),
	URRepeatNum_(1)
{}


template<class CloudType>
Foam::UndulationRemovalModel<CloudType>::UndulationRemovalModel
(
    const dictionary& dict,
    CloudType& owner,
    const word& type
)
:
    CloudSubModelBase<CloudType>(owner, dict, typeName, type),
	undulationRemovalInterval_(0),
	URIntervalCounter_(0),
	URRepeatNum_(1)
{}


template<class CloudType>
Foam::UndulationRemovalModel<CloudType>::UndulationRemovalModel
(
    const UndulationRemovalModel<CloudType>& cm
)
:
    CloudSubModelBase<CloudType>(cm),
	undulationRemovalInterval_(cm.undulationRemovalInterval_),
	URIntervalCounter_(cm.URIntervalCounter_),
	URRepeatNum_(cm.URRepeatNum_)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class CloudType>
Foam::UndulationRemovalModel<CloudType>::~UndulationRemovalModel()
{}


// * * * * * * * * * * * * * * * *  Selector * * * * * * * * * * * * * * * * //

template<class CloudType>
Foam::autoPtr<Foam::UndulationRemovalModel<CloudType>>
Foam::UndulationRemovalModel<CloudType>::New
(
    const dictionary& dict,
    CloudType& owner
)
{
    word modelType(dict.lookup(typeName));

    Info<< "Selecting undulationRemoval model " << modelType << endl;

    typename dictionaryConstructorTable::iterator cstrIter =
        dictionaryConstructorTablePtr_->find(modelType);

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        FatalErrorInFunction
            << "Unknown undulationRemoval model type " << modelType
            << ", constructor not in hash table" << nl << nl
            << "    Valid undulationRemoval model types are:" << nl
            << dictionaryConstructorTablePtr_->sortedToc()
            << exit(FatalError);
    }

    return
        autoPtr<UndulationRemovalModel<CloudType>>
        (
            cstrIter()(dict, owner)
        );
}


// ************************************************************************* //
