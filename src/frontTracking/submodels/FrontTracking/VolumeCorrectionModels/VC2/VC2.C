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

#include "VC2.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::VolumeCorrectionModels::VC2<CloudType>::VC2
(
    const dictionary& dict,
    CloudType& owner
)
:
    VolumeCorrectionModel<CloudType>(dict, owner, typeName),
	volumeCorrectionTolerance_(this->coeffDict().template lookup<scalar>("volumeCorrectionTolerance")),
	VCMaxIter_(this->coeffDict().template lookup<int>("VCMaxIter"))
{}


template<class CloudType>
Foam::VolumeCorrectionModels::VC2<CloudType>::VC2
(
    const VC2<CloudType>& cm
)
:
    VolumeCorrectionModel<CloudType>(cm),
	volumeCorrectionTolerance_(cm.volumeCorrectionTolerance_),
	VCMaxIter_(cm.VCMaxIter_)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class CloudType>
Foam::VolumeCorrectionModels::VC2<CloudType>::~VC2()
{}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType>
void Foam::VolumeCorrectionModels::VC2<CloudType>::calculate()
{
	Info << "\n---------> Performing the volume correction ... " << endl;

	//Inputs: pL.currentPoint, pL.connectedPoints, el.pointIndex, bData.volume (initial is sufficient)
	//Outputs (updated): pL.currentPoint
	//Changed (not updated): pL.PosInDomain, pL.connectedPoints, el.Points, el.centerPosIndomain, el.elementSurfaceArea, bData.volume

    //const fvMesh& mesh = this->owner().mesh();
    //const scalar deltaT(this->owner().db().time().deltaTValue());

	//DynamicList<bubbleData<CloudType>>& bDataL_= this->owner().bDataL();
	DynamicList<bubbleData>& bDataL_= this->owner().bDataL();

	forAll(bDataL_,bDI)
    {//1
		//necessary updates
		if (!bDataL_[bDI].plConnectionFlag) bDataL_[bDI].ptConnectedPoints();
		if (!bDataL_[bDI].volumeFlag) bDataL_[bDI].updateBubbleVolume();

		//DynamicList<pointData<CloudType>>& pL = bDataL_[bDI].pL;
		DynamicList<pointData>& pL = bDataL_[bDI].pL;
		DynamicList<elementInfo>& eL = bDataL_[bDI].eL;

		scalar Vloss = bDataL_[bDI].volume0 - bDataL_[bDI].volume;
		if (mag(Vloss)/bDataL_[bDI].volume0 > volumeCorrectionTolerance_)
		{	
			// normal velocity calculation at points
				DynamicList<vector> dp; //could be static list

			/*
						if (Pstream::parRun())
						{//2
						 //have to be implemented
						}//2
						if (!Pstream::parRun())
					{//2
					forAllIter(typename CloudType, this->owner(), iter1)
					{//3
						parcelType& p1 = iter1();

						label bubbleDI = p1.bubbleIndex();
						label pointDI = p1.currentIndex();

						pointData& currentPt = pL[pointDI];

						point& x = currentPt.currentPoint;
						DynamicList<point>& cPoints = currentPt.connectedPoints;

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
						vector dn = normalSum/(mag(normalSum)+ROOTVSMALL);
						//if (normalVelocity_)
						//{
							//interpolate u at p
							vector up = p1.U(); //
						//vector up = dn;
						dn = (dn & up) * dn;
						//}
					}//3
					}//2
			*/

				forAll(pL, ptI)
				{//2
			   	    //pointData<CloudType>& currentPt = pL[ptI];
			   	    pointData& currentPt = pL[ptI];
					point& x = currentPt.currentPoint;
					vector dn = vector::zero;

					//if( x.z() > thresholdTh_ ) //EA9 uncomment condition //3 for threshold cases
					//{//3
					DynamicList<point>& cPoints = currentPt.connectedPoints;

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
					dn = normalSum/(mag(normalSum)+ROOTVSMALL);
					//if (normalVelocity_)
					//{
						//interpolate u at p
						//vector up = ;        
					//vector up = dn;
					//dn = (dn & up) * dn;
					//}
					//}//3         
					dp.append(dn);       
				 }//2

				//Calculating a, b, c
				scalar a = 0.0;
				scalar b = 0.0;
				scalar c = 0.0;

				forAll(eL, elementI)
				{//4
					elementInfo& elInfo = eL[elementI];
					point x1 = pL[elInfo.pointIndex[0]].currentPoint;
					point x2 = pL[elInfo.pointIndex[1]].currentPoint;
					point x3 = pL[elInfo.pointIndex[2]].currentPoint;
					vector n1 = dp[elInfo.pointIndex[0]];
					vector n2 = dp[elInfo.pointIndex[1]];
					vector n3 = dp[elInfo.pointIndex[2]];

					a += (n1 & (n2 ^ n3));
					b += (x1 & (n2 ^ n3)) + (x2 & (n3 ^ n1)) + (x3 & (n1 ^ n2));
					c += (n1 & (x2 ^ x3)) + (n2 & (x3 ^ x1)) + (n3 & (x1 ^ x2));
				}//4
				a /= 6;
				b /= 6;
				c /= 6;

				//Solving for eps in Vloss=a*eps^3+b*eps^2+c*eps
				scalar eps = 0; //initial guess
				scalar guess = eps;
				scalar error = 1.e9;

				for (label nIter = 0; nIter < VCMaxIter_; ++nIter)
				{
					scalar f = a*pow(eps,3)+b*eps*eps+c*eps-Vloss;
					scalar d = 3*a*eps*eps+2*b*eps+c;

					eps = eps - f/d;
					
					error = mag(eps - guess)/(mag(eps)+VSMALL);
					
					if (error <= volumeCorrectionTolerance_) break;
					
					guess = eps;
				}
				Info << "\n				 Volume correction relative error = " << error << endl;
				Info << "\n				 Vloss = " << Vloss << "   eps = " << eps << endl;

			   //correct vertices
				forAll(pL, ptI)
				{//2
					pL[ptI].currentPoint += eps*dp[ptI];
					//pL[ptI].posInDomain = pL[ptI].currentPoint;
				}//2

				//bDataL_[bDI].pointToElementMapping(); //EA9 for updating el.points
		}
    }//1

	this->owner().setMotionFlags();
}


// ************************************************************************* //
