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
    Foam::FrontToFieldModels::CPT

Description
    CPT return-to-frontToField model.

SourceFiles
    CPT.C

\*---------------------------------------------------------------------------*/

#ifndef CPT_H
#define CPT_H

#include "FrontToFieldModel.H"
#include "BaseF2FModel.H"
#include "bubbleData.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace FrontToFieldModels
{

/*---------------------------------------------------------------------------*\
                         Class CPT Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType>
class CPT
:
    public BaseF2FModel<CloudType>
{
    // Private Member Functions


public:

    //- Runtime type information
    TypeName("cpt");

	//EA add3
	class CPTData // closest point transform data
    {
       public:
              // list of density and viscosity induced by each bubble. 
              List<scalar> inducedDensity;
              List<scalar> inducedViscosity;

			  //constructor
              CPTData()
              {
              };
    };	

	List<CPTData> CPTList_;

    // Constructors

        //- Construct from components
        CPT(const dictionary& dict, CloudType& owner);

        //- Construct copy
        CPT(const CPT<CloudType>& cm);

        //- Construct and return a clone
        virtual autoPtr<FrontToFieldModel<CloudType>> clone() const 
        {
            return autoPtr<FrontToFieldModel<CloudType>>
            (
                new CPT<CloudType>(*this)
            );
        }


    //- Destructor
    virtual ~CPT();


    //- Member Functions

        //- Calculate velocities
		virtual void IndicatorsConstruction();

		virtual void collectFields(label cellI, label bDI, scalar wBar, vector elAreaVec); 

		virtual void initFields();

		inline void initialAllocatingTheCPTList();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace FrontToFieldModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "CPT.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
