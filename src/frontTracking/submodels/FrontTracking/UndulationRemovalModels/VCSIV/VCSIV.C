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

#include "VCSIV.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::UndulationRemovalModels::VCSIV<CloudType>::VCSIV
(
    const dictionary& dict,
    CloudType& owner
)
:
    UndulationRemovalModel<CloudType>(dict, owner, typeName),
	combined_(this->coeffDict().lookup("combined"))
{
	this->coeffDict().lookup("undulationRemovalInterval") >> this->undulationRemovalInterval_;
	this->coeffDict().lookup("URRepeatNum") >> this->URRepeatNum_;
	Info << "---------> The undulation removal interval is " << this->undulationRemovalInterval_ << endl;
	Info << "---------> The number of undulation removal loops over fronts is " << this->URRepeatNum_ << endl;
	Info << "---------> The VCSIV combined mode " << combined_ << endl;
}


template<class CloudType>
Foam::UndulationRemovalModels::VCSIV<CloudType>::VCSIV
(
    const VCSIV<CloudType>& cm
)
:
    UndulationRemovalModel<CloudType>(cm),
	combined_(cm.combined_)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class CloudType>
Foam::UndulationRemovalModels::VCSIV<CloudType>::~VCSIV()
{}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType>
void Foam::UndulationRemovalModels::VCSIV<CloudType>::calculate()
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


			//scalar volDiff = (calcBubbleVolume(eL) - initialBubbleVolume_)/eL.size();
			scalar volDiff = 0.0; //(calcBubbleVolume(eL) - initialBubbleVolume_)/pL.size();

			//EA9
			//Info << "Bubble clacVolume is: " << calcBubbleVolume(eL) << endl;
			//Info << "Bubble initialVolume is: " << initialBubbleVolume_ << endl;
			//Info << "Bubble volDiff is: " << volDiff << endl;

			scalar omega = 1.0;
			label f;

			for( int repeat =0; repeat <= this->URRepeatNum_-1; repeat++)
			{//repeat
				 //if (repeat>0) bDataL_[bDI].ptConnectedPoints();
				 // Undulation removal procedure
				 forAll(pL, ptI)
				 {//2
			  	     //pointData<CloudType>& curPt = pL[ptI];
			  	     pointData& curPt = pL[ptI];
					 DynamicList<point> cPoints = curPt.connectedPoints;
					 DynamicList<label> cPointsIndex = curPt.connectedPointsIndex;
				     point& x2 = curPt.currentPoint;

				     // calculation of the A2
				     point xn = cPoints[cPoints.size()-1];
				     point x0 = cPoints[0];
				     vector A2 = vector::zero;
				     A2 = A2 + ( (xn - x2) ^ (x0 - x2) );
				     for( int cPI = 0; cPI <= cPoints.size()-2; cPI++)
				     {
				         point xi = cPoints[cPI];
				         point xii = cPoints[cPI+1];
				         A2 = A2 + ( (xi - x2) ^ (xii - x2) );
				     }

				     forAll(cPoints, cPI)
				     {//4
				         //pointData<CloudType>& otherPt = pL[cPointsIndex[cPI]];
				         pointData& otherPt = pL[cPointsIndex[cPI]];
				         DynamicList<point> otherCPoints = otherPt.connectedPoints;
				         point& x1 = otherPt.currentPoint;

				         // calculation of the A1
				         xn = otherCPoints[otherCPoints.size()-1];
				         x0 = otherCPoints[0];
				         vector A1 = vector::zero;
				         A1 = A1 + ( (xn - x1) ^ (x0 - x1) );
				         for( int oCPI = 0; oCPI <= otherCPoints.size()-2; oCPI++)
				         {
				             point xi = otherCPoints[oCPI];
				             point xii = otherCPoints[oCPI+1];
				             A1 = A1 + ( (xi - x1) ^ (xii - x1) );
				         }

				         // finding e22 for v
				         label up, down;
				         for( f=0; f <= cPoints.size()-1; f++)
				         {
				             if( curPt.connectedPointsIndex[f] == otherPt.currentIndex )
				             {
				                 break;
				             }
				         }
				         if ( f == 0)
				         {
				             up = 1;
				             down = cPoints.size()-1;
				         }
				         else if ( f == (cPoints.size()-1) )
				         {
				             up = 0;
				             down = cPoints.size()-2;
				         }
				         else
				         {
				             up = f+1;
				             down = f-1;
				         }
				         vector e2n2 = cPoints[down] - x2;
				         vector e22 = cPoints[up] - x2;
				         vector v = e2n2 - e22;

				         // calculation of the dx1s and dx2s
				         vector dx1s = vector::zero;
				         vector dx2s = vector::zero;

				         if(!combined_)
				         {
				             // calculation of the sumX2
				             vector sumX2 = vector::zero;
				             for( int oCPI = 0; oCPI <= cPoints.size()-1; oCPI++)
				             {
				                 sumX2 = sumX2 + cPoints[oCPI];
				             }
				             sumX2 = sumX2 - x1;

				             // calculation of the sumX1
				             vector sumX1 = vector::zero;
				             for( int oCPI = 0; oCPI <= otherCPoints.size()-1; oCPI++)
				             {
				                 sumX1 = sumX1 + otherCPoints[oCPI];
				             }
				             sumX1 = sumX1 - x2;

				             // calculation of the x1s and x2s
				             scalar n1 = otherCPoints.size();
				             scalar n2 = cPoints.size();
				             vector x1s = 1./( n1*n2-1.) * (sumX2 + n2 * sumX1);
				             vector x2s = (1./n2) * (x1s + sumX2);

				             // calculation of the dx1s and dx2s
				             dx1s = omega * (x1s - x1);
				             dx2s = omega * (x2s - x2);
				         }
				         else
				         {
				             // calculation of the dx1s and dx2s
				             dx1s = omega * (smoothTo(otherPt) - x1);
				             dx2s = omega * (smoothTo(curPt) - x2);
				         }

				         // calculation of the A
				         vector A = A1 + A2 + (v ^ (dx1s - dx2s));

				         // calculation of the normal,n
				         vector n = A/(mag(A)+ROOTVSMALL);

				         // calculation of the h
				         scalar h = -( (dx1s & A1) + (dx2s & A2) + ( dx2s & (v ^ dx1s) ) )/( (n & A) + volDiff + ROOTVSMALL);

				         // calculation of the new x2 and x1
				         x2 = x2 + dx2s + h * n;
				         x1 = x1 + dx1s + h * n;
				     }//4
				 }//2
			}//repeat

			//bDataL_[bDI].pointToElementMapping();  
			      	
        }//1
		this->URIntervalCounter_ = 0;
		this->owner().setMotionFlags();
    }//0
    this->URIntervalCounter_++;

}


//--------------------------------------------------------------------------
//-----------------------------Undulation Uttility--------------------------
//--------------------------------------------------------------------------
//- note: func to smooth a pointData var.
template<class CloudType>
inline Foam::vector Foam::UndulationRemovalModels::VCSIV<CloudType>::smoothTo
(
	pointData currentPt //pointData<CloudType> currentPt
)
{
    point ver = currentPt.currentPoint;
	DynamicList<point> cPoints = currentPt.connectedPoints;

    // Normal vector calculation
    // centre of the base for current point
    point p = vector::zero;
    forAll(cPoints, cPI)
    {
        p = p + cPoints[cPI];
    }
    p = p/cPoints.size();

    point xn = cPoints[cPoints.size()-1];
    point x0 = cPoints[0];
    vector normalSum = vector::zero;
    normalSum = normalSum + ( (xn - p) ^ (x0 - p) );
    for( int cPI = 0; cPI <= cPoints.size()-2; cPI++)
    {
        point xi = cPoints[cPI];
        point xii = cPoints[cPI+1];
        normalSum = normalSum + ( (xi - p) ^ (xii - p) );
    }
    vector n = normalSum/(mag(normalSum)+ROOTVSMALL);

    // calculation the V1 ,volume of the unitary polyhedron
    // (p+n, x0,..., xn, p)
    point pPlusN = p + n;
    scalar V1 = 0;
    vector Sf = 0.5 * (xn - p) ^ (x0 - p);
    vector he = (p + xn + x0)/3.0 - pPlusN;
    V1 = V1 + ( he & Sf)/3.0;
    for( int cPI = 0; cPI <= cPoints.size()-2; cPI++)
    {
        // volume of tetrahedron of p+n, p, xi, xii
        // he is a vector from p+n to the average of p, xi, and xii.
        point xi = cPoints[cPI];
        point xii = cPoints[cPI+1];

        Sf = 0.5 * (xi - p) ^ (xii - p);
        he = (p + xi + xii)/3.0 - pPlusN;
        V1 = V1 + ( he & Sf)/3.0;
    }

    // calculation the V ,volume of the polyhedron
    // (ver, x0,..., xn, p)
    scalar V = 0;
    Sf = 0.5 * (xn - p) ^ (x0 - p);
    he = (p + xn + x0)/3.0 - ver;
    V = V + ( he & Sf)/3.0;
	for( int cPI = 0; cPI <= cPoints.size()-2; cPI++)
    {
        // volume of tetrahedron of p+n, p, xi, xii
        // he is a vector from p+n to the average of p, xi, and xii.
        point xi = cPoints[cPI];
        point xii = cPoints[cPI+1];

        Sf = 0.5 * (xi - p) ^ (xii - p);
        he = (p + xi + xii)/3.0 - ver;
        V = V + ( he & Sf)/3.0;
     }

     scalar h = V/(V1+ROOTVSMALL);
     ver = p + h * n;
     return ver;
}


// ************************************************************************* //
