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

#include "CPT.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::FrontToFieldModels::CPT<CloudType>::CPT
(
    const dictionary& dict,
    CloudType& owner
)
: //EA d2
    BaseF2FModel<CloudType>(dict, owner, typeName) //EA d2
{
	this->twoFluidFlow_ = false; //EA d2
}


template<class CloudType>
Foam::FrontToFieldModels::CPT<CloudType>::CPT
(
    const CPT<CloudType>& cm
)
:
    BaseF2FModel<CloudType>(cm) //EA d2
{}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class CloudType>
Foam::FrontToFieldModels::CPT<CloudType>::~CPT()
{}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
/*
template<class CloudType>
void Foam::FrontToFieldModels::CPT<CloudType>::calculate()
{
	this->F2FCommunication();
	this->IndicatorsConstruction();
}
*/

//note: func to initialize the CPTList_.
template<class CloudType>
inline void FrontToFieldModels::CPT<CloudType>::initialAllocatingTheCPTList()
{
    const fvMesh& mesh = this->mesh();
    const pointField& ctrs = mesh.cellCentres();
	//DynamicList<bubbleData<CloudType>>& bDataL_= this->owner().bDataL();
	DynamicList<bubbleData>& bDataL_= this->owner().bDataL();

    CPTList_.clear();

    CPTList_.setSize(ctrs.size());

    forAll(CPTList_, i)
    {
        CPTList_[i].inducedDensity.setSize(bDataL_.size());
        CPTList_[i].inducedViscosity.setSize(bDataL_.size());
    }

    forAll(bDataL_,bDI)
    {//10
        forAll(ctrs, cI)
        {//2
            CPTData& CPT = CPTList_[cI]; 
/* seems not required
            if( (xMin <= ctrs[cI].x() && ctrs[cI].x() <= xMax) &&
                (yMin <= ctrs[cI].y() && ctrs[cI].y() <= yMax) &&
                (zMin <= ctrs[cI].z() && ctrs[cI].z() <= zMax)  )
            {
                CPT.inducedDensity[bDI] = bDataL_[bDI].density;
                CPT.inducedViscosity[bDI] = bDataL_[bDI].viscosity; 
            }
*/
            CPT.inducedDensity[bDI] = bDataL_[bDI].outerFluidDensity;
            CPT.inducedViscosity[bDI] = bDataL_[bDI].outerFluidViscosity; 
        }//2
    }//10
}

template<class CloudType>
void Foam::FrontToFieldModels::CPT<CloudType>::IndicatorsConstruction()
{
    const fvMesh& mesh = this->owner().mesh();
    const pointField& ctrs = mesh.cellCentres();
	//DynamicList<bubbleData<CloudType>>& bDataL_= this->owner().bDataL();
	DynamicList<bubbleData>& bDataL_= this->owner().bDataL();

    forAll(bDataL_, bDI)
    {//10
        DynamicList<elementInfo>& eL = bDataL_[bDI].eL;
        DynamicList<label> ctrsLabel = this->owner().allCellsInEeachMasket()[bDI]; //EA d2

        scalar hl = this->bubbleHl(bDI); //EA d2

        forAll(ctrsLabel,I)
        {//5
            const label& cI = ctrsLabel[I];

            CPTData& CPT = CPTList_[cI];

            vector minDistance;
            label targetElement = 0; //EA d2
            minDistance.x() = GREAT;
            minDistance.y() = GREAT;
            minDistance.z() = GREAT;

			forAll(eL, elementI)
			{//4
			    const elementInfo& elInfo = eL[elementI];
				const point& pos = elInfo.centrePosInDomain;
				vector distance = ctrs[cI] - pos;
	            if( mag(distance) <= mag(minDistance) )
	            {
	                targetElement = elementI;
	                minDistance = distance;
	            }
            }//4

            //calculating the density and viscosity induced by each bubble.
            CPT.inducedDensity[bDI]= this->LiuHeviside(bDataL_[bDI].outerFluidDensity, bDataL_[bDI].density, //EA d2
                                                  minDistance, eL[targetElement].elementSurfaceArea, this->hFactor_*hl); //EA2 //EA4
            CPT.inducedViscosity[bDI]= this->LiuHeviside(bDataL_[bDI].outerFluidViscosity, bDataL_[bDI].viscosity, //EA d2
                                                  minDistance, eL[targetElement].elementSurfaceArea, this->hFactor_*hl); //EA2 //EA4
        }//5
    }//10

    // constructing final indicators
    forAll(ctrs, cI)
    {//5
		CPTData& CPT = CPTList_[cI];

		scalar finalDensity = 0.0;
		scalar finalViscosity = 0.0;
		scalar densityDiffMax = 0.0;
		scalar viscosityDiffMax = 0.0;

		for(int bDI = 0; bDI < bDataL_.size() ; bDI++)
		{//0
			scalar densityDiff = mag(CPT.inducedDensity[bDI] -bDataL_[bDI].outerFluidDensity);
			scalar viscosityDiff = mag(CPT.inducedViscosity[bDI] -bDataL_[bDI].outerFluidViscosity);

			if( densityDiff >= densityDiffMax )
			{
				densityDiffMax = densityDiff;
			    finalDensity = CPT.inducedDensity[bDI];
			}
			if( viscosityDiff >= viscosityDiffMax )
			{
				viscosityDiffMax = viscosityDiff;
			    finalViscosity = CPT.inducedViscosity[bDI];
			}
		}//0
		this->densityIndicatorRef().primitiveFieldRef()[cI] = finalDensity; //EA d2
		this->viscosityIndicatorRef().primitiveFieldRef()[cI] = finalViscosity;  //EA d2
    }//5    
}


template<class CloudType>
void Foam::FrontToFieldModels::CPT<CloudType>::initFields()
{}


template<class CloudType>
void Foam::FrontToFieldModels::CPT<CloudType>::collectFields
(
	label cellI, 
	label bDI, 
	scalar wBar,
	vector elAreaVec
)
{}


// ************************************************************************* //
