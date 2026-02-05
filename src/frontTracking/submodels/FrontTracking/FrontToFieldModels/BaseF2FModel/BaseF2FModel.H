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
    Foam::FrontToFieldModels::BaseF2FModel

Description
    BaseF2FModel return-to-frontToField model.

SourceFiles
    BaseF2FModel.C

\*---------------------------------------------------------------------------*/

#ifndef BaseF2F_H
#define BaseF2F_H

#include "FrontToFieldModel.H"
#include "bubbleData.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace FrontToFieldModels
{

/*---------------------------------------------------------------------------*\
                         Class BaseF2FModel Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType>
class BaseF2FModel
:
    public FrontToFieldModel<CloudType>
{
    // Private Member Functions


public:

    //- Runtime type information
    //TypeName("baseF2F"); //EA d2

	// Enumerations
    enum spreadingOption
    {
        fixed,
        bubbleLocal,
        gridLocal
    };

	//Switch spreadingOption_;
	spreadingOption spreadingOption_;

	word front2field_; 

	List<label> p_; 
	List<label> q_; 
	List<label> r_;
	List<scalar> beta_;

	label numberOfIndicatorFiltering_ = 0;
	scalar filteringWeight_ = 0.0;

    // Constructors

        //- Construct from components //EA d2
        BaseF2FModel
		(
			const dictionary& dict, 
			CloudType& owner,
			const word& modelName //= typeName
		);

        //- Construct copy
        BaseF2FModel(const BaseF2FModel<CloudType>& cm);

		//EA d2
		/*        
		//- Construct and return a clone
        virtual autoPtr<FrontToFieldModel<CloudType>> clone() const
        {
            return autoPtr<FrontToFieldModel<CloudType>>
            (
                new BaseF2FModel<CloudType>(*this)
            );
        }
		*/
		virtual autoPtr<FrontToFieldModel<CloudType>> clone() const = 0;


    //- Destructor
    virtual ~BaseF2FModel();


    //- Member Functions

        //- Calculate velocities
        virtual void calculate();

		virtual void F2FCommunication();

		virtual void IndicatorsConstruction() = 0;

		virtual void IndicatorsFiltering();

		virtual void collectFields(label cellI, label bDI, scalar wBar, vector elAreaVec) = 0; 

		virtual void initFields() = 0;

		inline void adjustingTheImportantVarsForIndicators();

		inline scalar bubbleHl(label bDI); //EA d2

		inline scalar PeskinD(scalar r,scalar h); //EA d2

		inline scalar windowFunc(vector del, scalar h); //EA d2

		inline scalar windowFuncBar(vector del, scalar h, List<scalar> beta); //EA d2

		inline scalar mPQR(vector del, scalar h, label p, label q, label r); //EA d2

		inline scalar LiuHeviside(scalar outerFluidProperty, scalar property, vector minDistance, vector elementSurfaceArea, scalar gama); //EA d2
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace FrontToFieldModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "BaseF2FModel.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
