/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2013-2021 OpenFOAM Foundation
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
    Foam::FrontTrackingCloud

Description
    Adds FrontTracking modelling to clouds

SourceFiles
    FrontTrackingCloudI.H
    FrontTrackingCloud.C

\*---------------------------------------------------------------------------*/

#ifndef FrontTrackingCloud_H
#define FrontTrackingCloud_H

#include "volFieldsFwd.H"
#include "fvMatricesFwd.H"
#include "dimensionedTypes.H"
#include "fvMesh.H"
#include "fluidThermo.H"
#include "Cloud.H"

#include "bubbleData.H" //EA add

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

#include "UDictionary.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes

//EA rem
//template<class CloudType>
//class PackingModel;

//template<class CloudType>
//class DampingModel;

//template<class CloudType>
//class IsotropyModel;

//EA add
template<class CloudType>
class SurfaceTensionModel;

template<class CloudType>
class FrontToFieldModel;

template<class CloudType>
class VolumeCorrectionModel;

template<class CloudType>
class UndulationRemovalModel;


/*---------------------------------------------------------------------------*\
                       Class FrontTrackingCloudName Declaration
\*---------------------------------------------------------------------------*/

TemplateName(FrontTrackingCloud);


/*---------------------------------------------------------------------------*\
                         Class FrontTrackingCloud Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType>
class FrontTrackingCloud
:
    public CloudType,
    public FrontTrackingCloudName
{
public:

    // Public Typedefs

        //- Type of cloud this cloud was instantiated for
        typedef CloudType cloudType;

        //- Type of parcel the cloud was instantiated for
        typedef typename CloudType::parcelType parcelType;

        //- Convenience typedef for this cloud type
        typedef FrontTrackingCloud<CloudType> FrontTrackingCloudType;


private:

    // Private Data

        //- Cloud copy pointer
        autoPtr<FrontTrackingCloud<CloudType>> cloudCopyPtr_;


protected:

    // Protected data

		// References to the carrier gas fields

		// Sources

            //EA add- field for surface tension force.
            autoPtr<volVectorField::Internal> sTensionForceFromFT_;


            //EA add- field for pressure jump at interface to reduce spurious current near the interface.
            autoPtr<volVectorField::Internal> pressureJumpAtTheInterfaceFromFT_;

			word initialMakingFileOfTheBubbles_;

        	//scalar timeTobeSteady_ = 0;

			scalar timeOld_ = 0.0; // container of old time

			scalar magFlowVelocity_ = 0.0;

			// The bubbleData list of the bubbleData class.
            //DynamicList<bubbleData<FrontTrackingCloud<CloudType>>> bDataL_; 
            DynamicList<bubbleData> bDataL_; 

/*          
			// The periodical options for curved duct. note: At theta =0, the r( or polar) axis
			// is the very x-axis and y axis is alinged to flow direction

		    //- the centre of periodical axis ( z-axis can be rotational and/or translational periodical axis)
		    vector rCentre_ = vector::zero;

		    //- the periodical angle about rotational axis
		    scalar thetaP_ = 0.0;

		    //- the translational periodic length in the direction of z-axis
		    scalar tPLength_ = 0.0;

		    //- the periodic length of straight duct in the direction of z-axis
		    scalar sDPLength_ = 0.0;

		    //- periodical option
		    word periodicalOption_;// rotationalPeriodicalInCurvedDuct/translationalPeriodicalInCurvedDuct/
		                           // combinedRotationalPeriodicalInCurvedDuct/combinedTranslationalPeriodicalInCurvedDuct/none/straightDuctPeriodical

		    vector flowD_ = vector::zero; // flow direction.
		    vector normalToflowDI_ = vector::zero; // the first (radial direction) normal to flow direction.
		    vector normalToflowDII_ = vector::zero; // the second (z direction) normal to flow direction.
		    DynamicList<vector> flowDList_;// list of all flow direction at each cell.


		    //scalar EotvosNumber_ = 0.0;
		    //scalar MortonNumber_ = 0.0;
		    scalar ReDh_ = 0.0;
		    scalar Dh_ = 0.0; //length used in ReDh
		    scalar surfaceTensionCoeff_ = 0.0;
		    
		    scalar gradPAtflowDirection_ = 0.0;
*/
    	// The important input parameters of the front mesh.

			// Enumerations
			enum fCellScaleOption
			{
				fixed,
				bubbleAverage,
				frontAverage,
				gridLocal
			};

			//Switch fCellScaleOption_;
			fCellScaleOption fCellScaleOption_;

			//base length  scale 
			scalar h_;

			//front cell length scale calculation relaxation
			scalar fcsRelaxation_;

		    //- periodical option
		    word frontMeshInputOption_;// roughlyCalculated/finerCalculated/manually

		    scalar lengthScaleOfTheMesh_;

		    scalar minEdge_;

		    scalar maxEdge_;

		    scalar maxAspectRatio_;

		    scalar howFinner_; // in the range 0.05 (high finner)- 0.3 (low finner)

		    scalar factorForMinEdge_; // in the range 0.2 (high finner)- 0.3 (low finner)
		    scalar factorForMaxEdge_; // in the range 1.0 (high finner)- 2.0 (low finner)

		   //threshold options //EA
		   //scalar thresholdTh_;
		   List<DynamicList<label> > allCellsInEeachMasket_;
		   List<scalar> lengthScaleNearTheFront_;
		   List<scalar> minEdgeB_;
		   List<scalar> maxEdgeB_;

        // References to the cloud sub-models

			//EA rem
            //- Packing model
            //autoPtr<PackingModel<FrontTrackingCloud<CloudType>>>
            //    packingModel_;

            //- Damping model
            //autoPtr<DampingModel<FrontTrackingCloud<CloudType>>>
            //    dampingModel_;

            //- Exchange model
            //autoPtr<IsotropyModel<FrontTrackingCloud<CloudType>>>
            //    isotropyModel_;

			//EA add
			//- Surface tension model
			autoPtr<SurfaceTensionModel<FrontTrackingCloud<CloudType>>>
                surfaceTensionModel_;

			//- front to field model
			autoPtr<FrontToFieldModel<FrontTrackingCloud<CloudType>>>
                frontToFieldModel_;

			//- volume correction model
			autoPtr<VolumeCorrectionModel<FrontTrackingCloud<CloudType>>>
                volumeCorrectionModel_;

			//- undulation removal model
			autoPtr<UndulationRemovalModel<FrontTrackingCloud<CloudType>>>
                undulationRemovalModel_;
			

	// Protected Member Functions

		    //EA add
			//- reading initial bubbles
		    void readInitialBubbles();


        // Initialisation

            //- Set cloud sub-models
            void setModels();


public:

    // Constructors

        //- Construct given carrier fields
        FrontTrackingCloud
        (
            const word& cloudName,
            const volScalarField& rho,
            const volVectorField& U,
            const volScalarField& mu,
            const dimensionedVector& g,
            const bool readFields = true
        );

        //- Construct given carrier fields and thermo
        FrontTrackingCloud
        (
            const word& cloudName,
            const volScalarField& rho,
            const volVectorField& U,
            const dimensionedVector& g,
            const fluidThermo& carrierThermo,
            const bool readFields = true
        );

        //- Copy constructor with new name
        FrontTrackingCloud
        (
            FrontTrackingCloud<CloudType>& c,
            const word& name
        );

        //- Copy constructor with new name - creates bare cloud
        FrontTrackingCloud
        (
            const fvMesh& mesh,
            const word& name,
            const FrontTrackingCloud<CloudType>& c
        );

        //- Disallow default bitwise copy construction
        FrontTrackingCloud(const FrontTrackingCloud&) = delete;

        //- Construct and return clone based on (this) with new name
        virtual autoPtr<Cloud<parcelType>> clone(const word& name)
        {
            return autoPtr<Cloud<parcelType>>
            (
                new FrontTrackingCloud(*this, name)
            );
        }

        //- Construct and return bare clone based on (this) with new name
        virtual autoPtr<Cloud<parcelType>> cloneBare(const word& name) const
        {
            return autoPtr<Cloud<parcelType>>
            (
                new FrontTrackingCloud(this->mesh(), name, *this)
            );
        }


    //- Destructor
    virtual ~FrontTrackingCloud();


    // Member Functions

        // Access

            //- Return a reference to the cloud copy
            inline const FrontTrackingCloud& cloudCopy() const;

            //EA add- return the density indicator from the frontDynamic class.
            inline void updateDensityFromFT(volScalarField& rho); //EA d2

            //EA add- return the viscosity indicator from the frontDynamic class.
            //inline volScalarField::Internal& viscosityFromFTRef(); //EA d2

            //EA add- return the viscosity indicator from the frontDynamic class.
            inline void updateViscosityFromFT(volScalarField& mu); //EA d2

            //EA add- return the surface tension force from the frontDynamic class.
            inline volVectorField::Internal& sTensionForceFromFTRef();

            //EA add- return the surface tension force from the frontDynamic class.
            inline tmp<volVectorField::Internal> sTensionForceFromFT() const;

            //EA add- field for pressure jump at interface to reduce spurious current near the interface.
            inline volVectorField::Internal& pressureJumpAtTheInterfaceFromFTRef();

            //EA add- field for pressure jump at interface to reduce spurious current near the interface.
            inline tmp<volVectorField::Internal> pressureJumpAtTheInterfaceFromFT() const;

			//EA rem
            //- Return const access to the packing model
            //inline const PackingModel<FrontTrackingCloud<CloudType>>&
            //    packingModel() const;

            //- Return a reference to the packing model
            //inline PackingModel<FrontTrackingCloud<CloudType>>& packingModel();

            //- Return condt access to the damping model
            //inline const DampingModel<FrontTrackingCloud<CloudType>>&
            //    dampingModel() const;

            //- Return a reference to the damping model
            //inline DampingModel<FrontTrackingCloud<CloudType>>& dampingModel();

            //- Return condt access to the isotropy model
            //inline const IsotropyModel<FrontTrackingCloud<CloudType>>&
            //    isotropyModel() const;

            //- Return a reference to the isotropy model
            //inline IsotropyModel<FrontTrackingCloud<CloudType>>& isotropyModel();

			//EA add3
			//- Return condt access to the surface tension model
            inline const SurfaceTensionModel<FrontTrackingCloud<CloudType>>&
                surfaceTensionModel() const;

            //- Return a reference to the surface tension model
            inline SurfaceTensionModel<FrontTrackingCloud<CloudType>>& surfaceTensionModel();

			//- Return condt access to the surface tension model
            inline const FrontToFieldModel<FrontTrackingCloud<CloudType>>&
                frontToFieldModel() const;

            //- Return a reference to the surface tension model
            inline FrontToFieldModel<FrontTrackingCloud<CloudType>>& frontToFieldModel();

			//- Return condt access to the surface tension model
            inline const VolumeCorrectionModel<FrontTrackingCloud<CloudType>>&
                volumeCorrectionModel() const;

            //- Return a reference to the surface tension model
            inline VolumeCorrectionModel<FrontTrackingCloud<CloudType>>& volumeCorrectionModel();

			//- Return condt access to the surface tension model
            inline const UndulationRemovalModel<FrontTrackingCloud<CloudType>>&
                undulationRemovalModel() const;

            //- Return a reference to the surface tension model
            inline UndulationRemovalModel<FrontTrackingCloud<CloudType>>& undulationRemovalModel();

			//- Return const access to bubble Data
            //inline const DynamicList<bubbleData<FrontTrackingCloud<CloudType>>>&
            inline const DynamicList<bubbleData>&
                bDataL() const;

            //- Return reference to bubble Data
            //inline DynamicList<bubbleData>& bDataL();
            inline DynamicList<bubbleData>& bDataL();

			//- Return const-access to the length Scale Of The Mesh
            inline scalar lengthScaleOfTheMesh() const;

			inline List<DynamicList<label>> allCellsInEeachMasket() const;  

			inline List<scalar> lengthScaleNearTheFront() const; 

			inline scalar fCellL(label bDI); //EA d2

        // Cloud evolution functions

            //- Store the current cloud state
            void storeState();

            //- Reset the current cloud to the previously stored state
            void restoreState();

            //- Evolve the cloud
            void evolve();

            //- Particle motion
            template<class TrackCloudType>
            void motion
            (
                TrackCloudType& cloud,
                typename parcelType::trackingData& td
            );

		//EA add Front Member Functions

		    void cloudCleaning();

			void mapFrontToParcel();

			template<class TrackCloudType>
			void updateParcelsFromFront
			(
                TrackCloudType& cloud,
                typename parcelType::trackingData& td
            );

			template<class TrackCloudType>
			void frontCloudAdvection
			(
                TrackCloudType& cloud,
                typename parcelType::trackingData& td
            );

			//void mapPointToParcel(label bDI, pointData<FrontTrackingCloudType>& ptD);
			void mapPointToParcel(label bDI, pointData& ptD);
										//label& cellI, label& tetFaceI, label& tetPtI, label& posInList);

//			void mapPointToParcel(label bDI, pointData<FrontTrackingCloudType>& ptN, pointData<FrontTrackingCloudType>& ptC);
//			void mapPointToParcel(label bDI, pointData& ptN, pointData& ptC);

			void constructingInitialFront();

			void readDataOfTheFrontsAtTheCurrentTimeIntimeName();

			inline Foam::scalar calcLengthScaleOfTheMesh();

			inline void setLengthScalesNearTheFront();

			void communicationsOfTheFrontAndEulerianGrid();

			void searchAllCellsInEeachMasket();

//			void flowDirectionField();

			void processingThefrontMeshes();

			//void updatePlPosInDomain(bubbleData<FrontTrackingCloudType>& bData);
			void updatePlPosInDomain(bubbleData& bData);

			//void updateElCentrePosInDomain(bubbleData<FrontTrackingCloudType>& bData);
			void updateElCentrePosInDomain(bubbleData& bData);

			void adjustingTheCloudDataAfterChanges
			( 
				labelList& locInPointDataList,
				label bDI
			);

			//void frontRefining( label bDI, DynamicList<pointData<FrontTrackingCloudType>>& ptDataList, 
			void frontRefining( label bDI, DynamicList<pointData>& ptDataList, 
									   DynamicList<elementInfo>& elInfoList,
		                               scalar minEdge, scalar maxEdge, scalar maxAspectRatio);

		    inline bool isElementSmall(elementInfo el, scalar minEdge, scalar maxEdge,
						               scalar maxAspectRatio, scalar& s1, scalar& s2,
						               scalar& s3, label& n1, label& n2, label& n3);

			inline bool isElementRefined(elementInfo el, scalar minEdge, scalar maxEdge, scalar maxAspectRatio,
		                                 scalar& s1, scalar& s2, scalar& s3, label& n1, label& n2, label& n3);

			//void frontCoarsening(label bDI, DynamicList<pointData<FrontTrackingCloudType>>& ptDataList, 
			void frontCoarsening(label bDI, DynamicList<pointData>& ptDataList, 
										DynamicList<elementInfo>& elInfoList,
		                                scalar minEdge, scalar maxEdge, scalar maxAspectRatio);

			//void thresholdSet(); //EA

			//update pL.currentPoint from p.position(), p.currentIndex(), p.bubbleIndex() and pL.currentIndex(=p.currentIndex() ?)
			void mapParcelToFront();

			inline void setMotionFlags();

			inline void setCRFlags();

			//Foam::vector calcThePositionInDomain( point pos);

			static string getLineNoComment(IFstream& is)
			{
				string line;
				do
				{
					is.getLine(line);
				}
				while ((line.empty() || line[0] == '#') && is.good());
				return line;
			}


/*
		    inline void cacheFields();

		    inline Foam::tmp<Foam::volScalarField> p() const;   
*/

        //- I-O

            //- Print cloud information
            //void info(); //EA rem
			void printingAllImportantRunTimeData();

			void outputTimePostProcessing();

			void bubblePostProcessing();

			void writeBubblePositions();

			void writeBubbleVelocityAndAccelaration();

		    void writeBubbleSphericity();

		    //void writeBubbleDragCoefs();

		    //void writeBubbleLiftCoefs();

		    //void writeBubbleDiameterComponents();

		    inline void calcBubbleNonDimensionalNumbers();

			void printFrontsForAllTimes();

		    void writeDataOfTheFrontsAtTheOutputTimeIntimeName();


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const FrontTrackingCloud&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "FrontTrackingCloudI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "FrontTrackingCloud.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
