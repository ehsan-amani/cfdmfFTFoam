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

#include "DirectElementBased.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::SurfaceTensionModels::DirectElementBased<CloudType>::DirectElementBased
(
    const dictionary& dict,
    CloudType& owner
)
:
    SurfaceTensionModel<CloudType>(dict, owner, typeName)
{}


template<class CloudType>
Foam::SurfaceTensionModels::DirectElementBased<CloudType>::DirectElementBased
(
    const DirectElementBased<CloudType>& cm
)
:
    SurfaceTensionModel<CloudType>(cm)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class CloudType>
Foam::SurfaceTensionModels::DirectElementBased<CloudType>::~DirectElementBased()
{}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType>
void Foam::SurfaceTensionModels::DirectElementBased<CloudType>::calculate()
{
	//Inputs: pL.currentPoint, eL.elementIndex, el.pointIndex, el.elementSurfaceArea
	//Outputs (updated): eL.averageSurfaceTension, bData.pressureJumpAtTheInterface

    //const fvMesh& mesh = this->owner().mesh();
    //const scalar deltaT(this->owner().db().time().deltaTValue());

	//DynamicList<bubbleData<CloudType>>& bDataL_= this->owner().bDataL();
	DynamicList<bubbleData>& bDataL_= this->owner().bDataL();

	forAll(bDataL_,bDI)
    {//10
		//necessary updates
		if (!bDataL_[bDI].elSurfaceFlag) bDataL_[bDI].updateBubbleSurface();

        vector sumSurfaceTension = vector::zero;
		DynamicList<elementInfo>& eL = bDataL_[bDI].eL;

		//forAll(eL, elementI)
		//{//4
	    //	elementInfo& elInfo = eL[elementI];
        //    elInfo.elementSurfaceArea = elInfo.calcElementSurfaceArea();
        //}//4

        //Info << "\nST coeff\n\n " << bDataL_[bDI].surfaceTensionCoeff;
		forAll(eL, elementI)
        {
            elementInfo& elInfo = eL[elementI];
            label p1I = elInfo.pointIndex[0];
            label p2I = elInfo.pointIndex[1];
            label p3I = elInfo.pointIndex[2];

            point p1 = bDataL_[bDI].pL[p1I].currentPoint;
            point p2 = bDataL_[bDI].pL[p2I].currentPoint;
            point p3 = bDataL_[bDI].pL[p3I].currentPoint;

            point t1 = p2 - p1;
            point t2 = p3 - p2;
            point t3 = p1 - p3;
            
            label e1I = elInfo.elementIndex[0];
            label e2I = elInfo.elementIndex[1];
            label e3I = elInfo.elementIndex[2];

            point ne1 = eL[e1I].elementSurfaceArea/mag(eL[e1I].elementSurfaceArea);
            point ne2 = eL[e2I].elementSurfaceArea/mag(eL[e2I].elementSurfaceArea);
            point ne3 = eL[e3I].elementSurfaceArea/mag(eL[e3I].elementSurfaceArea);

            elInfo.averageSurfaceTension = 0.5 * bDataL_[bDI].surfaceTensionCoeff * ( (t1 ^ ne1) + (t2 ^ ne2) + (t3 ^ ne3) );
            
            sumSurfaceTension = sumSurfaceTension + elInfo.averageSurfaceTension;
        }
    	//Info << "\nSum surface Tension\n\n " << mag(sumSurfaceTension);	

        bDataL_[bDI].pressureJumpAtTheInterface = sumSurfaceTension/bDataL_[bDI].surfaceArea;
    }//10
}


// ************************************************************************* //
