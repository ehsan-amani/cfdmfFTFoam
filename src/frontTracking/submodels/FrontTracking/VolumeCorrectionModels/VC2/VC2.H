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
    Foam::VolumeCorrectionModels::VC2

Description
    VC2 return-to-volumeCorrection model.

SourceFiles
    VC2.C

\*---------------------------------------------------------------------------*/

#ifndef VC2_H
#define VC2_H

#include "VolumeCorrectionModel.H"
#include "bubbleData.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace VolumeCorrectionModels
{

/*---------------------------------------------------------------------------*\
                         Class VC2 Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType>
class VC2
:
    public VolumeCorrectionModel<CloudType>
{
    // Private Member Functions


public:

    //- Runtime type information
    TypeName("vc2");

	scalar volumeCorrectionTolerance_;
	int VCMaxIter_;

    // Constructors

        //- Construct from components
        VC2(const dictionary& dict, CloudType& owner);

        //- Construct copy
        VC2(const VC2<CloudType>& cm);

        //- Construct and return a clone
        virtual autoPtr<VolumeCorrectionModel<CloudType>> clone() const
        {
            return autoPtr<VolumeCorrectionModel<CloudType>>
            (
                new VC2<CloudType>(*this)
            );
        }


    //- Destructor
    virtual ~VC2();


    //- Member Functions

        //- Calculate velocities
        virtual void calculate();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace VolumeCorrectionModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "VC2.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
