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

#include "VCSIII.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::UndulationRemovalModels::VCSIII<CloudType>::VCSIII
(
    const dictionary& dict,
    CloudType& owner
)
:
    UndulationRemovalModel<CloudType>(dict, owner, typeName)
{
	this->coeffDict().lookup("undulationRemovalInterval") >> this->undulationRemovalInterval_;
	this->coeffDict().lookup("URRepeatNum") >> this->URRepeatNum_;
	Info << "---------> The undulation removal interval is " << this->undulationRemovalInterval_ << endl;
	Info << "---------> The number of undulation removal loops over fronts is " << this->URRepeatNum_ << endl;
}


template<class CloudType>
Foam::UndulationRemovalModels::VCSIII<CloudType>::VCSIII
(
    const VCSIII<CloudType>& cm
)
:
    UndulationRemovalModel<CloudType>(cm)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class CloudType>
Foam::UndulationRemovalModels::VCSIII<CloudType>::~VCSIII()
{}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType>
void Foam::UndulationRemovalModels::VCSIII<CloudType>::calculate()
{
	//Inputs: pL.currentIndex, pL.currentPoint, pL.connectedPoints, pL.connextedPointsIndex
	//Outputs (updated): pL.currentPoint
	//Changed (not updated): pL.PosInDomain, pL.connectedPoints, el.Points, el.centerPosIndomain, el.elementSurfaceArea, bData.volume

	if( this->URIntervalCounter_ == this->undulationRemovalInterval_)
    {//0
        Info << "\n---------> Performing the undulation removal ... " << endl;

		//DynamicList<bubbleData<CloudType>>& bDataL_= this->owner().bDataL();
		DynamicList<bubbleData>& bDataL_= this->owner().bDataL();

        forAll(bDataL_,bDI)
        {//1
			//necessary updates
			if (!bDataL_[bDI].plConnectionFlag) bDataL_[bDI].ptConnectedPoints();

        	//DynamicList<pointData<CloudType>>& pL = bDataL_[bDI].pL;
        	DynamicList<pointData>& pL = bDataL_[bDI].pL;
			//DynamicList<elementInfo>& eL = bDataL_[bDI].eL;

			for( int repeat =0; repeat <= this->URRepeatNum_-1; repeat++)
			{//repeat
				//if (repeat>0) bDataL_[bDI].ptConnectedPoints();
				// Vertex balance procedure
				forAll(pL, ptI)
				{//2
			   	    //pointData<CloudType>& currentPt = pL[ptI];
			   	    pointData& currentPt = pL[ptI];
					point& x = currentPt.currentPoint;
					DynamicList<point>& cPoints = currentPt.connectedPoints;

					// calculation of the
					point xs = vector::zero;
					forAll(cPoints, cPI)
					{
						xs = xs + cPoints[cPI];
					}
					xs = xs/cPoints.size();

					point xn = cPoints[cPoints.size()-1];
					point x0 = cPoints[0];
					vector normalSum = vector::zero;
					normalSum = normalSum + ( (xn - x) ^ (x0 - x) );
					for( int cPI = 0; cPI <= cPoints.size()-2; cPI++)
					{
						point xi = cPoints[cPI];
						point xii = cPoints[cPI+1];
						normalSum = normalSum + ( (xi - x) ^ (xii - x) );
					}
					vector n = normalSum/(mag(normalSum)+ROOTVSMALL);

					// calculation of the dxs
					vector dxs = xs - x;

					// calculation of the smoothed x
					x = x + dxs - ( dxs & n ) * n;
				 }//2
			}//repeat

			//bDataL_[bDI].pointToElementMapping();        	

        }//1
		this->URIntervalCounter_ = 0;
		this->owner().setMotionFlags();
    }//0
    this->URIntervalCounter_++;

}


// ************************************************************************* //
