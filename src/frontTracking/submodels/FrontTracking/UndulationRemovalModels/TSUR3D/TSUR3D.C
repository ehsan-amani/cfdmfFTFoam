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

#include "TSUR3D.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::UndulationRemovalModels::TSUR3D<CloudType>::TSUR3D
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
Foam::UndulationRemovalModels::TSUR3D<CloudType>::TSUR3D
(
    const TSUR3D<CloudType>& cm
)
:
    UndulationRemovalModel<CloudType>(cm)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class CloudType>
Foam::UndulationRemovalModels::TSUR3D<CloudType>::~TSUR3D()
{}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType>
void Foam::UndulationRemovalModels::TSUR3D<CloudType>::calculate()
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
				    point& ver = currentPt.currentPoint;
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
				 }//2

				 // Undulation removal procedure
				 forAll(pL, ptI)
				 {//3
			  	     //pointData<CloudType>& curPt = pL[ptI];
			  	     pointData& curPt = pL[ptI];
					 DynamicList<point> cPoints = curPt.connectedPoints;
					 DynamicList<label> cPointsIndex = curPt.connectedPointsIndex;

					 point& ver1 = curPt.currentPoint;

					 // calculation of  p1
			  	     point p1 = vector::zero;
					 forAll(cPoints, cPI)
					 {
						 p1 = p1 + cPoints[cPI];
					 }
					 p1 = p1 /cPoints.size();

					 forAll(cPoints, cPI)
					 {//4
					     //posInList = quickSelect(pL, cPointsIndex[cPI], 0, pL.size()-1);
					     //pointData<CloudType>& otherPt = pL[posInList];
					     //pointData<CloudType>& otherPt = pL[cPointsIndex[cPI]];
					     pointData& otherPt = pL[cPointsIndex[cPI]];
					     DynamicList<point> otherCPoints = otherPt.connectedPoints;

						 //-------------------------------------------------------------------
					 	 // calculation of V1(volume of the polyhedron), n1(normal vector)

					 	 // Normal vector calculation for curPt = ver1
						 point x1n = cPoints[cPoints.size() -1];
						 point x10 = cPoints[0];
						 vector normalSum1 = vector::zero;
						 normalSum1 = normalSum1 + ( (x1n - p1) ^ (x10 - p1) );
						 for( int cPI = 0; cPI <= cPoints.size()-2; cPI++)
						 {
							 point x1i = cPoints[cPI];
							 point x1ii = cPoints[cPI+1];
							 normalSum1 = normalSum1 + ( (x1i - p1) ^ (x1ii - p1) );
						 }
						 vector n1 = normalSum1/(mag(normalSum1)+ROOTVSMALL);

					     //----------------------------------------------------------------------
					     // calculation of V2(volume of the polyhedron), p2, n2(normal vector)
					     point& ver2 = otherPt.currentPoint;

					     point p2 = vector::zero;
					 	 forAll(otherCPoints, oCPI)
					 	 {
						 	p2 = p2 + otherCPoints[oCPI];
						 }
					 	 p2 = p2 /otherCPoints.size();

					 	 // Normal vector calculation for otherPt = ver2
					   	 point x2n = otherCPoints[otherCPoints.size() -1];
					   	 point x20 = otherCPoints[0];
						 vector normalSum2 = vector::zero;
						 normalSum2 = normalSum2 + ( (x2n - p2) ^ (x20 - p2) );
						 for( int oCPI = 0; oCPI <= otherCPoints.size()-2; oCPI++)
						 {
							 point x2i = otherCPoints[oCPI];
							 point x2ii = otherCPoints[oCPI+1];
							 normalSum2 = normalSum2 + ( (x2i - p2) ^ (x2ii - p2) );
						 }
						 vector n2 = normalSum2/(mag(normalSum2)+ROOTVSMALL);

					     //----------------------------------------------------------------------
					     //----------------------------------------------------------------------
						 vector n = (n1 + n2)/mag(n1 + n2);

						 vector m = 0.5 * (ver1 + ver2);

						 scalar h1 = mag( (ver1 - m) & n );
						 scalar h2 = mag( (ver2 - m) & n );

						 // new definitions for p1 and p2
						 p1 = ver1 - (h1 * n);
						 p2 = ver2 - (h2 * n);

						 //----------------------------------------------------------------------
				 		 // calculation the V1 ,volume of the polyhedron
						 // ver1, x10,..., x1n, p1)
						 scalar V1 = 0;
						 vector Sf1 = 0.5 * (x1n - p1) ^ (x10 - p1);
						 vector he1 = (p1 + x1n + x10)/3.0 - ver1;
						 V1 = V1 + mag( he1 & Sf1)/3.0;
					 	 for( int cPI = 0; cPI <= cPoints.size()-2; cPI++)
					 	 {
							 point x1i = cPoints[cPI];
							 point x1ii = cPoints[cPI+1];

							 //volume of tetrahedron of ver1, p1, x1i, x1ii
							 //he1 is a vector from ver1 to the average of p1, x1i, and x1ii.
							 Sf1 = 0.5 * (x1i - p1) ^ (x1ii - p1);
							 he1 = (p1 + x1i + x1ii)/3.0 - ver1;

							 V1 = V1 + mag( he1 & Sf1)/3.0;
					 	 }

						 //----------------------------------------------------------------------
						 // calculation the V2 ,volume of the polyhedron
						 // (ver, x20,..., x2n, p2)
						 scalar V2 = 0;
					     vector Sf2 = 0.5 * (x2n - p2) ^ (x20 - p2);
					     vector he2 = (p2 + x2n + x20)/3.0 - ver2;
					     V2 = V2 + mag( he2 & Sf2)/3.0;
						 for( int oCPI = 0; oCPI <= otherCPoints.size()-2; oCPI++)
						 {
							 point x2i = otherCPoints[oCPI];
							 point x2ii = otherCPoints[oCPI+1];

							 // volume of tetrahedron of ver2, p2, xi, xii
							 // h is a vector from p+n to the average of p, xi, and xii.
						     vector Sf2 = 0.5 * (x2i - p2) ^ (x2ii - p2);
							 vector he2 = (p2 + x2i + x2ii)/3.0 - ver2;

						     V2 = V2 + mag( he2 & Sf2)/3.0;
					 	 }

						 //----------------------------------------------------------------------
						 scalar h = (V1 + V2) / ( mag( V1/(h1 + ROOTVSMALL) ) + mag( V2/(h2 + ROOTVSMALL) ) );
						 ver1 = p1 + h * n;
			 	         ver2 = p2 + h * n;
				 	}//4
				 }//3

			}//repeat

			//bDataL_[bDI].pointToElementMapping();        	

        }//1
		this->URIntervalCounter_ = 0;
		this->owner().setMotionFlags();
    }//0
    this->URIntervalCounter_++;

}


// ************************************************************************* //
