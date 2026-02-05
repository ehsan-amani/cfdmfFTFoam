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
    Foam::FrontToFieldModels::NoFrontToField

Description

SourceFiles
    NoFrontToField.C

\*---------------------------------------------------------------------------*/

#ifndef NoFrontToField_H
#define NoFrontToField_H

#include "FrontToFieldModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace FrontToFieldModels
{

/*---------------------------------------------------------------------------*\
                         Class NoFrontToField Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType>
class NoFrontToField
:
    public FrontToFieldModel<CloudType>
{
public:

    //- Runtime type information
    TypeName("none");

    // Constructors

        //- Construct from components
        NoFrontToField(const dictionary& dict, CloudType& owner);

        //- Construct copy
        NoFrontToField(const NoFrontToField<CloudType>& cm);

        //- Construct and return a clone
        virtual autoPtr<FrontToFieldModel<CloudType>> clone() const
        {
            return autoPtr<FrontToFieldModel<CloudType>>
            (
                new NoFrontToField<CloudType>(*this)
            );
        }


    //- Destructor
    virtual ~NoFrontToField();


    //- Member Functions

        //- Calculate velocities
        virtual void calculate();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace FrontToFieldModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "NoFrontToField.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
