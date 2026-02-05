/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2020 OpenFOAM Foundation
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
    Foam::bubbleData

Description
    Stores all relevant solution info for cloud

SourceFiles
    bubbleDataI.H
    bubbleData.C

\*---------------------------------------------------------------------------*/

#ifndef bubbleData_H
#define bubbleData_H

#include "fvMesh.H"
//#include "Switch.H"
//#include "Tuple2.H"
//EA add
#include "fileName.H"
#include "pointField.H"
#include "point.H"
#include "DLList.H"
#include "sphereToCell.H"
#include "cellSet.H"
#include "simpleMatrix.H"
#include "mathematicalConstants.H" //e,pi,twoPi,piByTwo
#include "fvCFD.H"

#include "IFstream.H"
#include "OFstream.H"
#include "IStringStream.H"
#include "OStringStream.H"

#include <fstream> 
#include <iostream> 
#include <string>
#include <cmath>
#include <iomanip>
using namespace Foam::constant;
using std::fstream;
using std::ifstream;
using std::ofstream;
using std::iostream;
using std::ios;
using std::string;
using std::setw;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class bubbleData Declaration
\*---------------------------------------------------------------------------*/

//template<class CloudType>
class pointData
{
protected:

    // Protected Data

        //- Convenience typedef to the cloud's parcel type
        //typedef typename CloudType::parcelType parcelType;


   public:
          // The index of the current point. 
          label currentIndex;

          // The list of connected points to the current point as a DynamicList<pointData>. 
          //DynamicList<pointData> connectedPoints;

          DynamicList<point> connectedPoints;
          DynamicList<label> connectedPointsIndex;

          // positions of the current point as a vector. 
          vector currentPoint;
          vector posInDomain;

          // this is used when an element must be eliminated from DynamicList of elementData.
          bool keepPoint;

          label elOwner;

		  //autoPtr<parcelType> parcel;
//		  parcelType* parcel; //just to refer not own

		  //constructor
          pointData()
          {
              currentIndex = -1;
              currentPoint = vector::zero;
              posInDomain = vector::zero;
              keepPoint = true;
              elOwner = -1;
//			  parcel = nullptr;
          };
};

class elementInfo
{
   public:
          // The index of the current element. 
          label currentIndex;

          // indexes of the elements which are neighbours of the current one in counter clockwise order. 
          labelList elementIndex;

          // indexes of the points which are ends of the current element in counter clockwise order. 
          labelList pointIndex;

          // position of the points which are ends of the current element in counter clockwise order. 
          List<point> points;

          // this is used when an element must be eliminated from DynamicList of elementData.
          bool keepElement;

          // required data 
          vector centrePosInDomain;
          vector elementSurfaceArea;
          vector velocity;
          scalar pressure;
          symmTensor shearStress;
          vector averageSurfaceTension; 

		  //constructor
          elementInfo()
          {
              elementIndex.setSize(3);
              pointIndex.setSize(3);
              points.setSize(3);
              currentIndex = -1;
              for( int i=0; i<=2; i++)
              {
                  elementIndex[i] = -1;
                  pointIndex[i] = -1;
                  points[i] = vector::zero;
              }
              keepElement = true;
              centrePosInDomain = vector::zero;
              elementSurfaceArea = vector::zero;
              velocity = vector::zero;
              pressure = 0.0;
              shearStress = symmTensor::zero;
              averageSurfaceTension = vector::zero;
          };
		  
		  //- Destructor
		  virtual ~elementInfo() {}

		  // Member Functions
		  inline void updateElementSurfaceArea();

		  inline void findPointsOrder( label f, label& cNI, label& cNII, label& cNIII );
};


//template<class CloudType>
class bubbleData
{
    // Private Data

protected:

    // Protected Data

        //- Convenience typedef to the cloud's parcel type
        //typedef typename CloudType::parcelType parcelType;


public: 

		// The pointData List for the current bubble. 
		//DynamicList<pointData<CloudType>> pL;
		DynamicList<pointData> pL;

		// The elementInfo List for the current bubble. 
		DynamicList<elementInfo> eL;

		// The pressure force exerted on the current bubble. 
		vector pressureForce;

		// The shear stress (friction force) force exerted on the current bubble. 
		vector shearStressForce;

		// The pressure drag force value exerted on the current bubble. 
		scalar pressureDragForce;

		// The friction drag force value exerted on the current bubble. 
		scalar frictionDragForce;

		// The pressure lift force value exerted on the right side of the current bubble. 
		scalar pressureLiftForceI;

		// The pressure lift force value exerted on the last side of the current bubble. 
		scalar pressureLiftForceII;

		// The friction lift force exerted on the right side of the current bubble. 
		scalar frictionLiftForceI;

		// The friction lift force exerted on the last side of the current bubble. 
		scalar frictionLiftForceII;

		// The total drag force exerted on the current bubble. 
		scalar totalDragForce;

		// The total lift force value exerted on the right side of the current bubble. 
		scalar totalLiftForceI;

		// The total lift force value exerted on the last side of the current bubble.  
		scalar totalLiftForceII;

		// The buoyant force exerted on the current bubble. 
		vector buoyantForce;

		// The gravity force exerted on the current bubble. 
		vector gravityForce;

		// The accelaration of the current bubble. 
		vector accelaration;

		// The velocity of the current bubble. 
		vector velocity;

		// The velocity of the current bubble at old time. 
		vector velocityOld;

		// The centre of the current bubble. 
		vector centre;

		// The centre of the current bubble at old time. 
		vector centreOld;

		// The volume of the current bubble. 
		scalar volume;

		// initial volume of the current bubble. //EA9
		scalar volume0;

		// The surface area of the current bubble. 
		scalar surfaceArea;

		// The current bubble mass. 
		scalar mass;

		// The bubble density. 
		scalar density;

		// The bubble density ratio. 
		scalar outerFluidDensity;

		// The bubble viscosity. 
		scalar viscosity;

		// The bubble viscosity ratio. 
		scalar outerFluidViscosity;

		// The bubble drag coefficient. 
		scalar dragCoeff;

		// The bubble lift coefficient of right side of motion direction. 
		scalar LiftCoeffI;

		// The bubble lift coefficient of last side direction. 
		scalar LiftCoeffII;

		// The index of the current bubble. 
		label currentIndex;

		// The surface tension Coeff of the current bubble. 
		scalar surfaceTensionCoeff;

		// The Morton number of the current bubble. 
		scalar MortonNumber;

		// The Eotvos number of the current bubble. 
		scalar EotvosNumber;

		// The equalized diameter of the current bubble when it has the same volume with a sphere. 
		scalar sphereBubbleDiameter;

		// The maximum diameter, each diameter component show the maximum value of bubble in that direction. 
		List<scalar> diameter;

		// The sphericity showing how the bubble shape diffes from a sphere with the same volume. 
		scalar sphericity;

		vector pressureJumpAtTheInterface;

		// this is used when an bubble must 
		// be eliminated from DynamicList of bubbleData.
		bool keepBubble; 

		//flag: updated el.elementIndex
		bool elNeighborFlag;

		//flag: updated el.points and pL.elOwner
		bool elPointsFlag;

		//flag: updated pL.connectedPoints and pL.connectedPointsIndex
		bool plConnectionFlag;

		//flag: updated el.elementSurfaceArea and bData.surfaceArea
		bool elSurfaceFlag;

		//flag: updated bData.volume
		bool volumeFlag;

		//flag: updated pl.posInDomain
		bool pLPosInDomainFlag;

		//flag: updated el.centerPosInDomain
		bool elCPosInDomainFlag;

		//constructor
		bubbleData()
		{
			// The pointData List for the current bubble. 
			pL.clearStorage();
			eL.clearStorage();
			pressureForce = vector::zero;
			shearStressForce = vector::zero;
			pressureDragForce = 0.0;
			frictionDragForce = 0.0;
			pressureLiftForceI = 0.0;
			pressureLiftForceII = 0.0;
			frictionLiftForceI = 0.0;
			frictionLiftForceII = 0.0;
			totalDragForce = 0.0;
			totalLiftForceI = 0.0;
			totalLiftForceII = 0.0;
			buoyantForce = vector::zero;
			gravityForce = vector::zero;
			accelaration = vector::zero;
			velocity = vector::zero;
			velocityOld = vector::zero;
			centre = vector::zero;
			centreOld = vector::zero;
			volume = 0.0;
			volume0 = 0.0; //EA9
			surfaceArea = 0.0;
			mass = 0.0;
			density = 0.0;
			outerFluidDensity = 0.0;
			viscosity = 0.0;
			outerFluidViscosity = 0.0;
			dragCoeff = 0.0;
			LiftCoeffI = 0.0;
			LiftCoeffII = 0.0;
			currentIndex = 0.0;
			surfaceTensionCoeff = 0.0;
			MortonNumber = 0.0;
			EotvosNumber = 0.0;
			sphereBubbleDiameter = 0.0;
			diameter.setSize(3);
			sphericity = 0.0;
			pressureJumpAtTheInterface = vector::zero;
			keepBubble = true;
			elNeighborFlag = false;
			elPointsFlag = false;
			plConnectionFlag = false;
			elSurfaceFlag = false;
			volumeFlag = false;
			pLPosInDomainFlag = false;
			elCPosInDomainFlag = false;
		};

        // Run-time options
/*
public:

    // Constructors

        //- Construct null from mesh reference
        bubbleData(const fvMesh& mesh);

        //- Construct from mesh and dictionary
        bubbleData(const fvMesh& mesh, const dictionary& dict);

        //- Construct copy
        bubbleData(const bubbleData& cs);

*/
    //- Destructor
    virtual ~bubbleData() {}


    // Member Functions

		void startPtIndexFromZero();

		void findNighboures();

		void pointToElementMapping();

		void ptConnectedPoints();

		labelList adjustingTheMeshDataAfterChanges();

		void printInitialFronts(const fileName& outfileName);

		inline void updateBubbleVolume();

		inline void updateBubbleSurface();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "bubbleDataI.H"
/*
#ifdef NoRepository //for template
#   include "bubbleData.C"
#endif
*/
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
