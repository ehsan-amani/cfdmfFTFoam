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
    Foam::VolumeCorrectionModels::NoVolumeCorrection

Description

SourceFiles
    NoVolumeCorrection.C

\*---------------------------------------------------------------------------*/

#ifndef NoVolumeCorrection_H
#define NoVolumeCorrection_H

#include "VolumeCorrectionModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace VolumeCorrectionModels
{

/*---------------------------------------------------------------------------*\
                         Class NoVolumeCorrection Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType>
class NoVolumeCorrection
:
    public VolumeCorrectionModel<CloudType>
{
public:

    //- Runtime type information
    TypeName("none");

    // Constructors

        //- Construct from components
        NoVolumeCorrection(const dictionary& dict, CloudType& owner);

        //- Construct copy
        NoVolumeCorrection(const NoVolumeCorrection<CloudType>& cm);

        //- Construct and return a clone
        virtual autoPtr<VolumeCorrectionModel<CloudType>> clone() const
        {
            return autoPtr<VolumeCorrectionModel<CloudType>>
            (
                new NoVolumeCorrection<CloudType>(*this)
            );
        }


    //- Destructor
    virtual ~NoVolumeCorrection();


    //- Member Functions

        //- Calculate velocities
        virtual void calculate();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace VolumeCorrectionModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "NoVolumeCorrection.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
