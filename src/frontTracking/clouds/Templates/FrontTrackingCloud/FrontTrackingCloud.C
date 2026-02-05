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

\*---------------------------------------------------------------------------*/

#include "FrontTrackingCloud.H"
//EA rem
//#include "NoPacking.H"
//#include "ParticleStressModel.H"
//#include "NoDamping.H"
//#include "NoIsotropy.H"
//#include "TimeScaleModel.H"
//EA add
#include "SurfaceTensionModel.H" 
#include "FrontToFieldModel.H"
#include "VolumeCorrectionModel.H"
#include "UndulationRemovalModel.H"

using namespace Foam::constant; //EA

// * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * //

template<class CloudType>
void Foam::FrontTrackingCloud<CloudType>::setModels()
{
    // calction of the important input parameters of the front mesh.
    lengthScaleOfTheMesh_ = calcLengthScaleOfTheMesh(); //EA

	const word fCellScaleOption = this->subModelProperties().lookup("fCellScaleOption"); // to be moved to coarsening/refining and source spreading 
    dictionary bubbleAverageCoeffsDict = this->subModelProperties().subOrEmptyDict("bubbleAverageCoeffs");
    dictionary fixedCoeffsDict = this->subModelProperties().subOrEmptyDict("fixedCoeffs");
	dictionary frontAverageCoeffsDict = this->subModelProperties().subOrEmptyDict("frontAverageCoeffs");

	h_ = lengthScaleOfTheMesh_;
    if (fCellScaleOption == "bubbleAverage")
    {
        fCellScaleOption_ = bubbleAverage;
		bubbleAverageCoeffsDict.lookup("factorForMinEdge") >> factorForMinEdge_; 
		bubbleAverageCoeffsDict.lookup("factorForMaxEdge") >> factorForMaxEdge_;
        bubbleAverageCoeffsDict.lookup("maxAspectRatio") >> maxAspectRatio_;
		fcsRelaxation_ = bubbleAverageCoeffsDict.lookupOrDefault<Foam::scalar>("fcsRelaxation", 1.0);
    }
    else if (fCellScaleOption == "frontAverage")
    {
        fCellScaleOption_ = frontAverage;
		frontAverageCoeffsDict.lookup("factorForMinEdge") >> factorForMinEdge_; 
		frontAverageCoeffsDict.lookup("factorForMaxEdge") >> factorForMaxEdge_;
		frontAverageCoeffsDict.lookup("maxAspectRatio") >> maxAspectRatio_;
		fcsRelaxation_ = frontAverageCoeffsDict.lookupOrDefault<Foam::scalar>("fcsRelaxation", 1.0);
    }
    else if (fCellScaleOption == "fixed")
    {
        fCellScaleOption_ = fixed;
		h_ = fixedCoeffsDict.lookupOrDefault<Foam::scalar>("h", lengthScaleOfTheMesh_);
		fixedCoeffsDict.lookup("factorForMinEdge") >> factorForMinEdge_; 
		fixedCoeffsDict.lookup("factorForMaxEdge") >> factorForMaxEdge_;
		fixedCoeffsDict.lookup("maxAspectRatio") >> maxAspectRatio_;
		fcsRelaxation_ = 1.0;
    }
    //else if (fCellScaleOption == "gridLocal")
    //{
    //    fCellScaleOption_ = gridLocal;
    //}
    else
    {
        FatalErrorInFunction
            << "fCellScaleOption must be either 'fixed', 'frontAverage', 'bubbleAverage', or 'gridLocal (future)'"
            << nl << exit(FatalError);
    }

	Info << "\nFront coarseining/refining parameters: " << endl;
	Info << "	Selecting front cell length scale " << fCellScaleOption << endl;
    //Info << "---------> The lengh scale of the mesh is " << this->owner().lengthScaleOfTheMesh() << endl; //EA d2

    Info << "---------> The current value of front cell length scale, h, " << h_ << endl;
	Info << "---------> The minEdge/h (factorForMinEdge_) " << factorForMinEdge_ << endl;
	Info << "---------> The maxEdge/h (factorForMaxEdge_) " << factorForMaxEdge_ << endl;
    Info << "---------> maxAspectRatio " << maxAspectRatio_ <<  endl; 
	Info << "---------> fcsRelaxation " << fcsRelaxation_ <<  endl; 

	//EA
//EA add
	this->subModelProperties().lookup("initialMakingFileOfTheBubbles") >> initialMakingFileOfTheBubbles_;
//	this->coeffDict().lookup("periodicalOption") >> periodicalOption_;
    this->subModelProperties().lookup("frontMeshInput") >> frontMeshInputOption_; // to be moved to coarsening/refining
    //this->subModelProperties().lookup("timeTobeSteady") >> timeTobeSteady_;

    //this->coeffDict().lookup("thresholdTh") >> thresholdTh_; //EA

    // Here, the dict is the same dictionary as subModels in (cloudName + "Properties") file.
//    dictionary periodicalOptionCoeffsDict = FrontDynamicCoeffsDict.subOrEmptyDict("periodicalOptionCoeffs");
    dictionary manuallyMeshInputCoeffsDict = this->subModelProperties().subOrEmptyDict("manuallyMeshInputCoeffs");
    dictionary finerCalculatedMeshInputCoeffsDict = this->subModelProperties().subOrEmptyDict("finerCalculatedMeshInputCoeffs");
	// to be moved to coarsening/refining
    if( frontMeshInputOption_ == "roughlyCalculated" || frontMeshInputOption_ == "finerCalculated" ) 
    {
        minEdge_ = factorForMinEdge_ * lengthScaleOfTheMesh_;
        maxEdge_ = factorForMaxEdge_ * lengthScaleOfTheMesh_;             

        Info << "\nFront mesh parameter calculation option " << frontMeshInputOption_ << endl;
        Info << "			 The calculated lengh scale of the mesh is " << lengthScaleOfTheMesh_ << endl;
        Info << "			" << " minEdge = " << minEdge_ << " factorForMinEdge = " << factorForMinEdge_ <<  endl; 
        Info << "			" << " maxEdge = " << maxEdge_ << " factorForMaxEdge = " << factorForMaxEdge_ <<  endl; 

        if( frontMeshInputOption_ == "finerCalculated" )
        {
            finerCalculatedMeshInputCoeffsDict.lookup("howFinner") >> howFinner_;
            Info << "			 howFinner = " << howFinner_ << endl;
        }
    }     
    else if( frontMeshInputOption_ == "manually")
    {
        manuallyMeshInputCoeffsDict.lookup("howFinner") >> howFinner_; 
        manuallyMeshInputCoeffsDict.lookup("minEdge") >> minEdge_;
        manuallyMeshInputCoeffsDict.lookup("maxEdge") >> maxEdge_; 

        Info << "\nFront mesh parameter calculation option " << frontMeshInputOption_ << endl;
        Info << "			 The calculated lengh scale of the mesh is " << lengthScaleOfTheMesh_ << endl;
        Info << "			" << " minEdge = " << minEdge_ << " maxEdge = " << maxEdge_ <<  endl; 
		Info << "			" << " howFinner = " << howFinner_ <<  endl;             
    }

    //Info << "---------> The lengh scale of the mesh is " << lengthScaleOfTheMesh_ << endl; //EA
    //EA rem
    /*
    packingModel_.reset
    (
        PackingModel<FrontTrackingCloud<CloudType>>::New
        (
            this->subModelProperties(),
            *this
        ).ptr()
    );
    dampingModel_.reset
    (
        DampingModel<FrontTrackingCloud<CloudType>>::New
        (
            this->subModelProperties(),
            *this
        ).ptr()
    );
    isotropyModel_.reset
    (
        IsotropyModel<FrontTrackingCloud<CloudType>>::New
        (
            this->subModelProperties(),
            *this
        ).ptr()
    );
    */
	//EA add
	surfaceTensionModel_.reset 
    (
        SurfaceTensionModel<FrontTrackingCloud<CloudType>>::New
        (
            this->subModelProperties(),
            *this
        ).ptr()
    );

    frontToFieldModel_.reset 
    (
        FrontToFieldModel<FrontTrackingCloud<CloudType>>::New
        (
            this->subModelProperties(),
            *this
        ).ptr()
    );

	volumeCorrectionModel_.reset 
    (
        VolumeCorrectionModel<FrontTrackingCloud<CloudType>>::New
        (
            this->subModelProperties(),
            *this
        ).ptr()
    );

	undulationRemovalModel_.reset 
    (
        UndulationRemovalModel<FrontTrackingCloud<CloudType>>::New
        (
            this->subModelProperties(),
            *this
        ).ptr()
    );

/*
    if( periodicalOption_ == "rotationalPeriodicalInCurvedDuct")
    {
        periodicalOptionCoeffsDict.lookup("thetaP") >> thetaP_; 
	//thetaP_ = readScalar(myDictL2.lookup("thetaP"));
        //myDictL2.readIfPresent("thetaP", thetaP_);
        periodicalOptionCoeffsDict.lookup("rCentre") >> rCentre_;
        periodicalOptionCoeffsDict.lookup("gradPAtflowDirection") >> gradPAtflowDirection_; 
        periodicalOptionCoeffsDict.lookup("Dh") >> Dh_; 

        Info << "\n   *The periodical option is selected as " << periodicalOption_ << endl;
        Info << "          The periodic angle is selected as " << thetaP_ << endl;             
        Info << "          The origin of the periodical axis is located at " << rCentre_ << endl;
        Info << "          The pressur gradient value induced in flow direction is " << gradPAtflowDirection_ << endl;
        Info << "          The length scale value used in ReDh is " << Dh_ << endl;

        //converting the thetaP to radian:
        thetaP_ = thetaP_ * mathematical::pi/180.0;
    }   
    else if( periodicalOption_ == "translationalPeriodicalInCurvedDuct")
    {
        periodicalOptionCoeffsDict.lookup("tPLength") >> tPLength_; 
        periodicalOptionCoeffsDict.lookup("rCentre") >> rCentre_;
        periodicalOptionCoeffsDict.lookup("gradPAtflowDirection") >> gradPAtflowDirection_; 
        periodicalOptionCoeffsDict.lookup("Dh") >> Dh_; 
 
        Info << "\n   *The periodical option is selected as " << periodicalOption_ << endl;
        Info << "          The periodic length is selected as " << tPLength_ << endl;
        Info << "          The origin of the periodical axis is located at " << rCentre_ << endl;
        Info << "          The pressur gradient value induced in flow direction is " << gradPAtflowDirection_ << endl;
        Info << "          The length scale value used in ReDh is " << Dh_ << endl;
    }
    else if( periodicalOption_ == "combinedRotationalPeriodicalInCurvedDuct" || periodicalOption_ == "combinedTranslationalPeriodicalInCurvedDuct" )
    {
        periodicalOptionCoeffsDict.lookup("thetaP") >> thetaP_; 
        periodicalOptionCoeffsDict.lookup("tPLength") >> tPLength_;
        periodicalOptionCoeffsDict.lookup("rCentre") >> rCentre_;
        periodicalOptionCoeffsDict.lookup("gradPAtflowDirection") >> gradPAtflowDirection_; 
        periodicalOptionCoeffsDict.lookup("Dh") >> Dh_; 

        Info << "\n   *The periodical option is selected as " << periodicalOption_ << endl;
        Info << "          The periodic angle is selected as " << thetaP_ << endl;             
        Info << "          The periodic length is selected as " << tPLength_ << endl;  
        Info << "          The origin of the periodical axis is located at " << rCentre_ << endl;
        Info << "          The pressure gradient value induced in flow direction is " << gradPAtflowDirection_ << endl;
        Info << "          The length scale value used in ReDh is " << Dh_ << endl;

        //converting the thetaP to radian:
        thetaP_ = thetaP_ * mathematical::pi/180.0;             
    }
    else if( periodicalOption_ == "straightDuctPeriodical")
    {
        periodicalOptionCoeffsDict.lookup("sDPLength") >> sDPLength_;
        periodicalOptionCoeffsDict.lookup("rCentre") >> rCentre_;
        periodicalOptionCoeffsDict.lookup("gradPAtflowDirection") >> gradPAtflowDirection_; 
        periodicalOptionCoeffsDict.lookup("Dh") >> Dh_; 

        Info << "\n   *The periodical option is selected as " << periodicalOption_ << endl;           
        Info << "          The periodic length  of straigh duct is selected as " << sDPLength_ << endl;  
        Info << "          The flow direction is parallel to the z axis" << endl;
        Info << "          The pressur gradient value induced in flow direction is " << gradPAtflowDirection_ << endl;
        Info << "          The length scale value used in ReDh is " << Dh_ << endl;             
    }    
    else if( periodicalOption_ == "none")
    {
        periodicalOptionCoeffsDict.lookup("gradPAtflowDirection") >> gradPAtflowDirection_; 
        periodicalOptionCoeffsDict.lookup("Dh") >> Dh_;     
        Info << "\n   *The periodical option is selected as " << periodicalOption_ << endl;
        Info << "          This means, there is no periodic B.C. in this problem." << endl;
        Info << "          but, we use a length scale, Dh for the ReDh calculation as: " << Dh_ << endl;
        Info << "          and the pressur gradient value induced in flow direction is: " << gradPAtflowDirection_ << endl;

    }
*/

//	flowDirectionField();
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::FrontTrackingCloud<CloudType>::FrontTrackingCloud
(
    const word& cloudName,
    const volScalarField& rho,
    const volVectorField& U,
    const volScalarField& mu,
    const dimensionedVector& g,
    const bool readFields
)
:
    CloudType(cloudName, rho, U, mu, g, false),
	//EA add
    sTensionForceFromFT_
    (
        new volVectorField::Internal
        (
            IOobject
            (
                this->name() + ":sTensionForceFromFT",
                this->db().time().timeName(),
                this->db(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            this->mesh(),
            dimensionedVector(dimForce/dimVolume,  vector::zero)
        )
    ),
    pressureJumpAtTheInterfaceFromFT_
    (
        new volVectorField::Internal
        (
            IOobject
            (
                this->name() + ":pressureJumpAtTheInterfaceFromFT",
                this->db().time().timeName(),
                this->db(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            this->mesh(),
            dimensionedVector(dimForce/dimVolume,  vector::zero)
        )
    ),
    //EA rem
    //packingModel_(nullptr),
    //dampingModel_(nullptr),
    //isotropyModel_(nullptr),
	//EA add
	surfaceTensionModel_(nullptr),
    frontToFieldModel_(nullptr),
	volumeCorrectionModel_(nullptr),
	undulationRemovalModel_(nullptr)
{
//EA - add
/*
    if (this->solution().steadyState())
    {
        FatalErrorInFunction
            << "FrontTracking modelling not available for steady state calculations"
            << exit(FatalError);
    }
*/

    setModels();

    if (readFields)
    {
        parcelType::readFields(*this);
        this->deleteLostParticles();
    }

	readInitialBubbles();
}


template<class CloudType>
Foam::FrontTrackingCloud<CloudType>::FrontTrackingCloud
(
    const word& cloudName,
    const volScalarField& rho,
    const volVectorField& U,
    const dimensionedVector& g,
    const fluidThermo& carrierThermo,
    const bool readFields
)
:
    FrontTrackingCloud(cloudName, rho, U, carrierThermo.mu(), g, readFields)
{}


template<class CloudType>
Foam::FrontTrackingCloud<CloudType>::FrontTrackingCloud
(
    FrontTrackingCloud<CloudType>& c,
    const word& name
)
:
    CloudType(c, name),
//EA add
    sTensionForceFromFT_
    (
        new volVectorField::Internal
        (
            IOobject
            (
                this->name() + ":sTensionForceFromFT",
                this->db().time().timeName(),
                this->db(),
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            c.sTensionForceFromFT_()
        )
    ),
    pressureJumpAtTheInterfaceFromFT_
    (
        new volVectorField::Internal
        (
            IOobject
            (
                this->name() + ":pressureJumpAtTheInterfaceFromFT",
                this->db().time().timeName(),
                this->db(),
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            c.pressureJumpAtTheInterfaceFromFT_()
        )
    ),
    //EA rem
    //packingModel_(c.packingModel_->clone()),
    //dampingModel_(c.dampingModel_->clone()),
    //isotropyModel_(c.isotropyModel_->clone()),
	//EA add
	surfaceTensionModel_(c.surfaceTensionModel_->clone()),
    frontToFieldModel_(c.frontToFieldModel_->clone()),
	volumeCorrectionModel_(c.volumeCorrectionModel_->clone()),
	undulationRemovalModel_(c.undulationRemovalModel_->clone())
{}


template<class CloudType>
Foam::FrontTrackingCloud<CloudType>::FrontTrackingCloud
(
    const fvMesh& mesh,
    const word& name,
    const FrontTrackingCloud<CloudType>& c
)
:
    CloudType(mesh, name, c),
//EA add
    sTensionForceFromFT_(nullptr),
    pressureJumpAtTheInterfaceFromFT_(nullptr),
    //EA rem
    //packingModel_(nullptr),
    //dampingModel_(nullptr),
    //isotropyModel_(nullptr),
	//EA add
	surfaceTensionModel_(nullptr),
    frontToFieldModel_(nullptr),
	volumeCorrectionModel_(nullptr),
	undulationRemovalModel_(nullptr)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class CloudType>
Foam::FrontTrackingCloud<CloudType>::~FrontTrackingCloud()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType>
void Foam::FrontTrackingCloud<CloudType>::storeState()
{
    cloudCopyPtr_.reset
    (
        static_cast<FrontTrackingCloud<CloudType>*>
        (
            clone(this->name() + "Copy").ptr()
        )
    );
}


template<class CloudType>
void Foam::FrontTrackingCloud<CloudType>::restoreState()
{
    this->cloudReset(cloudCopyPtr_());
    cloudCopyPtr_.clear();
}


template<class CloudType>
void Foam::FrontTrackingCloud<CloudType>::evolve()
{
    if (this->solution().canEvolve())
    {
        typename parcelType::trackingData td(*this);

        this->solve(*this, td);
    }
}

template<class CloudType>
template<class TrackCloudType>
void Foam::FrontTrackingCloud<CloudType>::motion
(
    TrackCloudType& cloud,
    typename parcelType::trackingData& td
)
{
    //const scalar time = this->db().time().value();
    //const scalar deltaTime = this->db().time().deltaTValue();

	if (this->mesh().changing())
    {
		lengthScaleOfTheMesh_ = calcLengthScaleOfTheMesh(); //EA
    	Info << "\n----------> The current lengh scale of the mesh is " << lengthScaleOfTheMesh_ << endl; //EA
	}

	//FTM subalgorithms
	//thresholdSet(); //EA9 uncomment for threshold setting
	processingThefrontMeshes();    
	surfaceTensionModel_->calculate(); //EA replace surfaceTensionForceDistribution();
	communicationsOfTheFrontAndEulerianGrid(); 
	outputTimePostProcessing(); 
	//cloudCleaning(); //EA11
	updateParcelsFromFront(cloud, td); //mapFrontToParcel(); //EA11 
	frontCloudAdvection(cloud, td);

/*
    // force calculation and tracking
    td.part() = parcelType::trackingData::tpPredictTrack;
    CloudType::move(cloud, td, this->db().time().deltaTValue());


    // Preliminary
    // ~~~~~~~~~~~

    // switch forces off so they are not applied in corrector steps
    this->forces().setCalcNonCoupled(false);
    this->forces().setCalcCoupled(false);


    // Damping
    // ~~~~~~~

    if (!isType<DampingModels::NoDamping<CloudType>>(dampingModel_()))
    {
        if (this->mesh().moving())
        {
            FatalErrorInFunction
                << "FrontTracking damping modelling does not support moving meshes."
                << exit(FatalError);
        }

        // update averages
        td.updateAverages(cloud);

        // memory allocation and eulerian calculations
        dampingModel_->cacheFields(true);

        // calculate the damping velocity corrections without moving the parcels
        td.part() = parcelType::trackingData::tpDampingNoTrack;
        CloudType::move(cloud, td, this->db().time().deltaTValue());

        // correct the parcel positions and velocities
        td.part() = parcelType::trackingData::tpCorrectTrack;
        CloudType::move(cloud, td, this->db().time().deltaTValue());

        // finalise and free memory
        dampingModel_->cacheFields(false);
    }


    // Packing
    // ~~~~~~~

    if (!isType<PackingModels::NoPacking<CloudType>>(packingModel_()))
    {
        if (this->mesh().moving())
        {
            FatalErrorInFunction
                << "FrontTracking packing modelling does not support moving meshes."
                << exit(FatalError);
        }

        // same procedure as for damping
        td.updateAverages(cloud);
        packingModel_->cacheFields(true);
        td.part() = parcelType::trackingData::tpPackingNoTrack;
        CloudType::move(cloud, td, this->db().time().deltaTValue());
        td.part() = parcelType::trackingData::tpCorrectTrack;
        CloudType::move(cloud, td, this->db().time().deltaTValue());
        packingModel_->cacheFields(false);
    }


    // Isotropy
    // ~~~~~~~~

    if (!isType<IsotropyModels::NoIsotropy<CloudType>>(isotropyModel_()))
    {
        // update averages
        td.updateAverages(cloud);

        // apply isotropy model
        isotropyModel_->calculate();
    }


    // Final
    // ~~~~~

    // update cell occupancy
    this->updateCellOccupancy();

    // switch forces back on
    this->forces().setCalcNonCoupled(true);
    this->forces().setCalcCoupled(this->solution().coupled());
*/
}


//note: func to process the front meshes.
template<class CloudType>
void Foam::FrontTrackingCloud<CloudType>::processingThefrontMeshes()
{
	volumeCorrectionModel_->calculate(); //EA 

    Info << "\n---------> The front meshes is refining/coarsening  ... " << endl;

    forAll(bDataL_,bDI)
    {//1
    	 DynamicList<pointData>& pL = bDataL_[bDI].pL;
         DynamicList<elementInfo>& eL = bDataL_[bDI].eL;

         Info << "			for bubble " << bDI << endl;

         Info << "          Before refining and coarsening, the pL,eL sizes were "
              << pL.size() << ',' << eL.size() << endl;

		 //necessary updates
		 if (!bDataL_[bDI].elPointsFlag) bDataL_[bDI].pointToElementMapping();

		 //EA-12 
		 scalar hl = fCellL(bDI); 
		 scalar& minEdge = minEdgeB_[bDI];
		 scalar& maxEdge = maxEdgeB_[bDI];
    	 minEdge = fcsRelaxation_ * factorForMinEdge_ * hl +(1-fcsRelaxation_)*minEdge;
         maxEdge = fcsRelaxation_ * factorForMaxEdge_ * hl +(1-fcsRelaxation_)*maxEdge;

         frontRefining( bDI, pL, eL, minEdge, maxEdge, maxAspectRatio_); 
         frontCoarsening( bDI, pL, eL, minEdge, maxEdge, maxAspectRatio_);

         labelList locInPointDataList = bDataL_[bDI].adjustingTheMeshDataAfterChanges(); 
		 adjustingTheCloudDataAfterChanges(locInPointDataList,bDI);  //EA11

         Info << "          After  refining and coarsening, the pL,eL sizes are "
              << pL.size() << ',' << eL.size() << endl;
    }//1
	setCRFlags();	

  	undulationRemovalModel_->calculate(); //EA 

	searchAllCellsInEeachMasket(); //EA12
}


template<class CloudType>
void Foam::FrontTrackingCloud<CloudType>::adjustingTheCloudDataAfterChanges
( 
	labelList& locInPointDataList,
	label bDI
)
{
	// remove parcels and update parcels currentIndex  
    forAllIter(typename FrontTrackingCloud<CloudType>, *this, pIter) 
    {		
		parcelType& p1 = pIter();

        label bubbleDI = p1.bubbleIndex();
		if (bubbleDI == bDI)
		{
		    label pointDI = p1.currentIndex(); //oldIndex

			if (locInPointDataList[pointDI] == -1) //(!ptD.keepPoint) // point has been removed
			{
				//removing parcel
				this->deleteParticle(p1);
			}
			else
			{
				//DynamicList<pointData>& pL = bDataL_[bubbleDI].pL;
				//pointData& ptD = pL[pointDI];
				p1.currentIndex() = locInPointDataList[pointDI]; //update parcel current Index								
			}
		}
    }
}

//--------------------------------------------------------------------------
//-----------------------------frontRefining--------------------------------
//--------------------------------------------------------------------------
	//Inputs: eL.currentIndex, eL.elementIndex, el.points, el.pointIndex
	//Outputs (updated): pL.currentIndex, pL.currentPoint, pL.PosInDomain, eL.currentIndex, eL.elementIndex, el.pointIndex
	//					 pL.elOwner, el.points 
	//Changed (not updated): pL.connectedPoints, pL.connectedPointsIndex, el.centerPosIndomain, el.elementSurfaceArea, bData.volume

//--------------------------------------------------------------------------
//-----------------------------frontCoarsening------------------------------
//--------------------------------------------------------------------------
	//Inputs: eL.currentIndex, eL.elementIndex, el.points, el.pointIndex
	//Outputs (updated): pL.currentIndex, pL.currentPoint, pL.PosInDomain, eL.currentIndex, eL.elementIndex, el.pointIndex
	//					 el.points 
	//Changed (not updated): pL.elOwner (updated), pL.connectedPoints, pL.connectedPointsIndex, el.centerPosIndomain, el.elementSurfaceArea, bData.volume


//-note:  func to refine the surface grid
template<class CloudType>
void Foam::FrontTrackingCloud<CloudType>::frontRefining
( 
	label bDI,	
	DynamicList<pointData>& pL, 
	DynamicList<elementInfo>& eL,
    scalar minEdge, 
	scalar maxEdge, 
	scalar maxAspectRatio
)
{
    label pI1, pI2, pI3, pI4;
    label c4I;
    //label c4II;
    label c4III;
    //label c5I;
    //label c5II;
    //label c5III;
    //label c6I;
    label c6II;
    //label c6III;
    //label c7I;
    label c7II;
    //label c7III;
    //label c8I;
    //label c8II;
    //label c8III;
    //label c9I, c9II, c9III;
    //label c10I, c10II, c10III;
    //label c11I, c11II, c11III;
    label I, II, III;
    label f, cNI, cNII, cNIII;

	//label cellI = -1;
	//label tetFaceI = -1;
	//label tetPtI = -1;
	//label posInList = -1;

    label numberOfELAdded = 0;

    label initialElListSize = eL.size();

    //forAll(eL, elementI)
    for(int elementI=0; elementI < initialElListSize; elementI++)
    {//2
        elementInfo& currentElement = eL[elementI];

        // order the lengths, s1 shortest, s3 longest
        scalar s1 = mag( currentElement.points[1] - currentElement.points[0]);
        scalar s2 = mag( currentElement.points[2] - currentElement.points[1]);
        scalar s3 = mag( currentElement.points[0] - currentElement.points[2]);
        label n1 = 0;
        label n2 = 1;
        label n3 = 2;

        bool refinedElement = isElementRefined(currentElement, minEdge, maxEdge, maxAspectRatio, s1, s2, s3, n1, n2, n3);

        if ( refinedElement && currentElement.keepElement )
        {//3
            label largeEdge = n3;

            // I, II, III are the point indexes of current element.
            // n1 is the index of shortest edge and p1.
            currentElement.findPointsOrder( largeEdge, I, II, III);
            point p1 = currentElement.points[I];
            point p2 = currentElement.points[II];
            point p3 = currentElement.points[III];
            pI1 = currentElement.pointIndex[I];
            pI2 = currentElement.pointIndex[II];
            pI3 = currentElement.pointIndex[III];

            // find point 5, third point of neighbElementII
            elementInfo& neighbElementII = eL[currentElement.elementIndex[II]];
            for( f=0; f<=3; f++)
            {
                if( neighbElementII.pointIndex[f] == pI2)
                {
                    //Info<< "f5 " << ' ' << f<< endl;
                    break;
                }
                if( f == 3)
                {
					FatalErrorIn
					(
						" ^^^^^^^^^^^^^ we can not find such a point."
					) << abort(FatalError);
                }
            }
            // c5I, c5II, c5III are the point indexes of neighbElementII.
            currentElement.findPointsOrder( f, cNI, cNII, cNIII);
            point p5 = neighbElementII.points[cNII];
            //c5I = cNII;
			//c5II = cNIII;
			//c5III = cNI;

            // find point 6, third point of neighbElementIII
            elementInfo& neighbElementIII = eL[currentElement.elementIndex[III]];
            for( f=0; f<=3; f++)
            {
                if( neighbElementIII.pointIndex[f] == pI3)
                {
                    //Info<< "f6 " << ' ' << f<< endl;
                    break;
                }
                if( f == 3)
                {
					FatalErrorIn
					(
						" ^^^^^^^^^^^^^ we can not find such a point."
					) << abort(FatalError);
                }
            }
            //Info<< "f6 " << ' ' << f<< endl;
            currentElement.findPointsOrder( f, cNI, cNII, cNIII);
            point p6 = neighbElementIII.points[cNII];
			//c6I = cNII;
			c6II = cNIII;
			//c6III = cNI;

            // find point 4, third point of neighbElementIII
            elementInfo& neighbElementI = eL[currentElement.elementIndex[I]];
            for( f=0; f<=3; f++)
            {
                if( neighbElementI.pointIndex[f] == pI1)
                {
                    //Info<< "f4 " << ' ' << f<< endl;
                    break;
                }
                if( f == 3)
                {
					FatalErrorIn
					(
						" ^^^^^^^^^^^^^ we can not find such a point."
					) << abort(FatalError);
                }
            }
            currentElement.findPointsOrder( f, cNI, cNII, cNIII);
            point p4 = neighbElementI.points[cNII];
            pI4 = neighbElementI.pointIndex[cNII];
			c4I = cNII;
			//c4II = cNIII;
			c4III = cNI;


            // find point 7, third point of neighbElementILeft
            elementInfo& neighbElementILeft = eL[neighbElementI.elementIndex[c4III]];
            for( f=0; f<=3; f++)
            {
                if( neighbElementILeft.pointIndex[f] == pI1)
                {
                    //Info<< "f7 " << ' ' << f<< endl;
                    break;
                }
                if( f == 3)
                {
					FatalErrorIn
					(
						" ^^^^^^^^^^^^^ we can not find such a point."
					) << abort(FatalError);
                }
            }
            currentElement.findPointsOrder( f, cNI, cNII, cNIII);
            point p7 = neighbElementILeft.points[cNII];
			//c7I = cNII;
			c7II = cNIII;
			//c7III = cNI;

            // find point 8, third point of neighbElementIRight
            elementInfo& neighbElementIRight = eL[neighbElementI.elementIndex[c4I]];

            for( f=0; f<=3; f++)
            {
                if( neighbElementIRight.pointIndex[f] == pI4)
                {
                    //Info<< "f8 " << ' ' << f<< endl;
                    break;
                }
                if( f == 3)
                {
					FatalErrorIn
					(
						" ^^^^^^^^^^^^^ we can not find such a point."
					) << abort(FatalError);
                }
            }
            currentElement.findPointsOrder( f, cNI, cNII, cNIII);
            point p8 = neighbElementIRight.points[cNII];
			//c8I = cNII;
			//c8II = cNIII;
			//c8III = cNI;

            numberOfELAdded = numberOfELAdded + 2;

		    // allocation of the new point:
		    pointData newPoint;
		    newPoint.currentIndex = pL.size();
		    newPoint.currentPoint = 0.5 * ( p1 + p2) + 0.125 * ( p3 + p4) - 0.0625 * ( p5 + p6 + p7 + p8);
			newPoint.posInDomain = newPoint.currentPoint; //EA add temp when they are the same
		    newPoint.keepPoint = true;
		    newPoint.elOwner = currentElement.currentIndex;

		    // allocation of the new element at the upper part:
		    elementInfo newElementU;
		    newElementU.currentIndex = eL.size();
		    newElementU.points[0] = newPoint.currentPoint;
		    newElementU.points[1] = p3;
		    newElementU.points[2] = p1;
		    newElementU.pointIndex[0] = newPoint.currentIndex;
		    newElementU.pointIndex[1] = pI3;
		    newElementU.pointIndex[2] = pI1;
		    newElementU.keepElement = true;

		    // allocation of the new element at the bottom part:
		    elementInfo newElementB;
		    newElementB.currentIndex = eL.size()+1;
		    newElementB.points[0] = newPoint.currentPoint;
		    newElementB.points[1] = p1;
		    newElementB.points[2] = p4;
		    newElementB.pointIndex[0] = newPoint.currentIndex;
		    newElementB.pointIndex[1] = pI1;
		    newElementB.pointIndex[2] = pI4;
		    newElementB.keepElement = true;

		    // transform currentElement[I] and neighbElementI[c4III] to point to the newPoint
		    currentElement.points[I] = newPoint.currentPoint;
		    currentElement.pointIndex[I] = newPoint.currentIndex;
		    neighbElementI.points[c4III] = newPoint.currentPoint;
		    neighbElementI.pointIndex[c4III] = newPoint.currentIndex;

		    // update the neighboures of newElementU
		    newElementU.elementIndex[0] = currentElement.currentIndex;
		    newElementU.elementIndex[1] = neighbElementIII.currentIndex;
		    newElementU.elementIndex[2] = newElementB.currentIndex;

		    // update the neighboures of newElementB
		    newElementB.elementIndex[0] = newElementU.currentIndex;
		    newElementB.elementIndex[1] = neighbElementILeft.currentIndex;
		    newElementB.elementIndex[2] = neighbElementI.currentIndex;

		    // update the neighboures of currentElement and neighbElementI
		    currentElement.elementIndex[III] = newElementU.currentIndex;
		    neighbElementI.elementIndex[c4III] = newElementB.currentIndex;

		    // update the neighboures of neighbElementIII and neighbElementILeft
		    neighbElementIII.elementIndex[c6II] = newElementU.currentIndex;
		    neighbElementILeft.elementIndex[c7II] = newElementB.currentIndex;

		    //updating the elOwner for p1
		    pL[pI1].elOwner = neighbElementILeft.currentIndex;

		    pL.append(newPoint);
		    eL.append(newElementU);
		    eL.append(newElementB);
			//adding the new parcel  //EA11
			mapPointToParcel(bDI,newPoint); //,cellI,tetFaceI,tetPtI,posInList); 
			//mapPointToParcel(bDI,newPoint,pL[pI1]); 
       }//3
   }//2
   Info << "                 After refining " << numberOfELAdded << " element(s) were added." << endl;
}


//- note: func to coarse the surface grid
template<class CloudType>
void Foam::FrontTrackingCloud<CloudType>::frontCoarsening
(
	label bDI,	
	DynamicList<pointData>& pL, 
	DynamicList<elementInfo>& eL,
	scalar minEdge, 
	scalar maxEdge, 
	scalar maxAspectRatio
)
{
    label pI1, pI2, pI3, pI4;
    label c4I;
    //label c4II;
    label c4III;
    //label c5I;
    label c5II;
    label c5III;
    label c6I;
    label c6II;
    //label c6III;
    //label c7I;
    label c7II;
    //label c7III;
    label c8I;
    label c8II;
    //label c8III;
    //label c9I, c9II, c9III;
    //label c10I, c10II, c10III;
    //label c11I, c11II, c11III;
    label I, II, III;
    label f, cNI, cNII, cNIII;

	//label cellI = -1;
	//label tetFaceI = -1;
	//label tetPtI = -1;
	//label posInList = -1;

    label numberOfELDeleted = 0;

    //label initialElListSize = eL.size();

    forAll(eL, elementI)
    {//2
        elementInfo& currentElement = eL[elementI];

        // order the lengths, s1 shortest, s3 longest
        scalar s1 = mag( currentElement.points[1] - currentElement.points[0]);
        scalar s2 = mag( currentElement.points[2] - currentElement.points[1]);
        scalar s3 = mag( currentElement.points[0] - currentElement.points[2]);
        label n1 = 0;
        label n2 = 1;
        label n3 = 2;
        bool smallElement = isElementSmall(currentElement, minEdge, maxEdge, maxAspectRatio, s1, s2, s3, n1, n2, n3);

        if ( smallElement && currentElement.keepElement)
        {//3
            numberOfELDeleted = numberOfELDeleted + 2;

            label shortEdge = n1;

            // I, II, III are the point indexes of current element.
            // n1 is the index of shortest edge and p1.
            currentElement.findPointsOrder( shortEdge, I, II, III);
            point p1 = currentElement.points[I];
            point p2 = currentElement.points[II];
            point p3 = currentElement.points[III];
            pI1 = currentElement.pointIndex[I];
            pI2 = currentElement.pointIndex[II];
            pI3 = currentElement.pointIndex[III];

            // find point 5, third point of neighbElementII
            elementInfo& neighbElementII = eL[currentElement.elementIndex[II]];
            for( f=0; f<=3; f++)
            {
                if( neighbElementII.pointIndex[f] == pI2)
                {
                    //Info<< "f5 " << ' ' << f<< endl;
                    break;
                }
                if( f == 3)
                {
					FatalErrorIn
					(
						" ^^^^^^^^^^^^^ we can not find such a point."
					) << abort(FatalError);
                }
            }
            // c5I, c5II, c5III are the point indexes of neighbElementII.
            currentElement.findPointsOrder( f, cNI, cNII, cNIII);
            point p5 = neighbElementII.points[cNII];
            //c5I = cNII;
			c5II = cNIII;
			c5III = cNI;

            // find point 6, third point of neighbElementIII
            elementInfo& neighbElementIII = eL[currentElement.elementIndex[III]];
            for( f=0; f<=3; f++)
            {
                if( neighbElementIII.pointIndex[f] == pI3)
                {
                    //Info<< "f6 " << ' ' << f<< endl;
                    break;
                }
                if( f == 3)
                {
					FatalErrorIn
					(
						" ^^^^^^^^^^^^^ we can not find such a point."
					) << abort(FatalError);
                }
            }
            //Info<< "f6 " << ' ' << f<< endl;
            currentElement.findPointsOrder( f, cNI, cNII, cNIII);
            point p6 = neighbElementIII.points[cNII];
			c6I = cNII;
			c6II = cNIII;
			//c6III = cNI;

            // find point 4, third point of neighbElementIII
            elementInfo& neighbElementI = eL[currentElement.elementIndex[I]];
            for( f=0; f<=3; f++)
            {
                if( neighbElementI.pointIndex[f] == pI1)
                {
                    //Info<< "f4 " << ' ' << f<< endl;
                    break;
                }
                if( f == 3)
                {
					FatalErrorIn
					(
						" ^^^^^^^^^^^^^ we can not find such a point."
					) << abort(FatalError);
                }
            }
            currentElement.findPointsOrder( f, cNI, cNII, cNIII);
            point p4 = neighbElementI.points[cNII];
            pI4 = neighbElementI.pointIndex[cNII];
			c4I = cNII;
			//c4II = cNIII;
			c4III = cNI;


            // find point 7, third point of neighbElementILeft
            elementInfo& neighbElementILeft = eL[neighbElementI.elementIndex[c4III]];
            for( f=0; f<=3; f++)
            {
                if( neighbElementILeft.pointIndex[f] == pI1)
                {
                    //Info<< "f7 " << ' ' << f<< endl;
                    break;
                }
                if( f == 3)
                {
					FatalErrorIn
					(
						" ^^^^^^^^^^^^^ we can not find such a point."
					) << abort(FatalError);
                }
            }
            currentElement.findPointsOrder( f, cNI, cNII, cNIII);
            point p7 = neighbElementILeft.points[cNII];
			//c7I = cNII;
			c7II = cNIII;
			//c7III = cNI;

            // find point 8, third point of neighbElementIRight
            elementInfo& neighbElementIRight = eL[neighbElementI.elementIndex[c4I]];
            for( f=0; f<=3; f++)
            {
                if( neighbElementIRight.pointIndex[f] == pI4)
                {
                    //Info<< "f8 " << ' ' << f<< endl;
                    break;
                }
                if( f == 3)
                {
					FatalErrorIn
					(
						" ^^^^^^^^^^^^^ we can not find such a point."
					) << abort(FatalError);
                }
            }
            currentElement.findPointsOrder( f, cNI, cNII, cNIII);
            point p8 = neighbElementIRight.points[cNII];
			c8I = cNII;
			c8II = cNIII;
			//c8III = cNI;

            if( neighbElementIRight.currentIndex == neighbElementII.currentIndex )
            {

                Info << " UUUUUUUUUUUUUUUUUUUUUUUUU---the right triple occured." << currentElement.currentIndex << endl;

                // allocation of the new point:
                point newPos = 0.5 * ( p1 + p2) + 0.125 * ( p3 + p4) - 0.0625 * ( p5 + p6 + p7 + p8);

                // position of p1 now is updated.
				pointData& ptD1 = pL[pI1];
				ptD1.currentPoint = newPos;
				ptD1.posInDomain = ptD1.currentPoint; //EA add, temp when they are the same 

                // transform only all positions of points of the element sharing the lable of p1 to newPos
                //cNI = c6I;
                //cNII = c6II;
                //cNIII = c6III;
                cNIII = c6I;
                cNI = c6II;
                elementInfo elTemp1 = neighbElementIII;
                neighbElementIII.points[cNI] = newPos;
                do
                {
                    elementInfo& neighbElTemp1 = eL[elTemp1.elementIndex[cNIII]];
                    for( f=0; f<=3; f++)
                    {
                       if( neighbElTemp1.pointIndex[f] == pI1 )
                       {
                           //Info<< "F I can break " << ' ' << f << endl;
                           break;
                       }
                       if( f == 3)
                       {
						   FatalErrorIn
						   (
						   " ^^^^^^^^^^^^^ we can not find such a point, pI1."
						   ) << abort(FatalError);
                       }
                    }
                    currentElement.findPointsOrder( f, cNI, cNII, cNIII);
                    neighbElTemp1.points[cNI] = newPos;
                    elTemp1 = neighbElTemp1;
                } while( elTemp1.currentIndex != neighbElementILeft.currentIndex );

                // transform all positions and indexes of points of the element sharing the p2 to newPos
                // then the p2 must be eliminated.
                neighbElementII.points[c5III] = newPos;
                neighbElementII.pointIndex[c5III] = pI1;

                // update the neighboures of neighbElementIII and neighbElementII
                neighbElementIII.elementIndex[c6II] = neighbElementII.currentIndex;
                neighbElementII.elementIndex[c5II] = neighbElementIII.currentIndex;

                // update the neighboures of neighbElementILeft and neighbElementIRight
                neighbElementILeft.elementIndex[c7II] = neighbElementIRight.currentIndex;
                neighbElementIRight.elementIndex[c8II] = neighbElementILeft.currentIndex;

                // p2 now is eliminated.
				pointData& ptD2 = pL[pI2];
                //posInList = quickSelect(pL, pI2, 0, pL.size()-1);
                //pL[posInList].keepPoint = false;
                ptD2.keepPoint = false;

                // eleminating the elements: currentElement and neighbElementI
                currentElement.keepElement = false;
                neighbElementI.keepElement = false;

            }
            else if( neighbElementILeft.currentIndex == neighbElementIII.currentIndex  )
            {

                Info << " UUUUUUUUUUUUUUUUUUUUUUUUU---the left triple occured." << currentElement.currentIndex << endl;

                // allocation of the new point:
                point newPos = 0.5 * ( p1 + p2) + 0.125 * ( p3 + p4) - 0.0625 * ( p5 + p6 + p7 + p8);

                // position of p1 now is updated.
				pointData& ptD1 = pL[pI1];
				ptD1.currentPoint = newPos;
				ptD1.posInDomain = ptD1.currentPoint; //EA add temp when they are the same 

                // transform only all positions of points of the element sharing the lable of p1 to newPos
                neighbElementIII.points[c6II] = newPos;

                // transform all positions and indexes of points of the element sharing the p2 to newPos
                // then the p2 must be eliminated.
                //cNI = c8I;
                //cNII = c8II;
                //cNIII = c8III;
                cNIII = c8I;
                cNI = c8II;
                elementInfo elTemp2 = neighbElementIRight;
                neighbElementIRight.points[cNI] = newPos;
                neighbElementIRight.pointIndex[cNI] = pI1;
                do
                {
                    elementInfo& neighbElTemp2 = eL[elTemp2.elementIndex[cNIII]];
                    for( f=0; f<=3; f++)
                    {
                        if( neighbElTemp2.pointIndex[f] == pI2 )
                        {
                            //Info<< "S I can break " << ' ' << f << endl;
                            break;
                        }
                        if( f == 3)
                        {
						    FatalErrorIn
						    (
							" ^^^^^^^^^^^^^ we can not find such a point, pI2."
						    ) << abort(FatalError);
                        }
                    }
                    currentElement.findPointsOrder( f, cNI, cNII, cNIII);
                    neighbElTemp2.points[cNI] = newPos;
                    neighbElTemp2.pointIndex[cNI] = pI1;
                    elTemp2 = neighbElTemp2;
                } while( elTemp2.currentIndex != neighbElementII.currentIndex );

                // update the neighboures of neighbElementIII and neighbElementII
                neighbElementIII.elementIndex[c6II] = neighbElementII.currentIndex;
                neighbElementII.elementIndex[c5II] = neighbElementIII.currentIndex;

                // update the neighboures of neighbElementILeft and neighbElementIRight
                neighbElementILeft.elementIndex[c7II] = neighbElementIRight.currentIndex;
                neighbElementIRight.elementIndex[c8II] = neighbElementILeft.currentIndex;

                // p2 now is eliminated.
				pointData& ptD2 = pL[pI2];
                ptD2.keepPoint = false;

                // eleminating the elements: currentElement and neighbElementI
                currentElement.keepElement = false;
                neighbElementI.keepElement = false;
            }

            else
            {
                // allocation of the new point:
                point newPos = 0.5 * ( p1 + p2) + 0.125 * ( p3 + p4) - 0.0625 * ( p5 + p6 + p7 + p8);

                // position of p1 now is updated.
				pointData& ptD1 = pL[pI1];
				ptD1.currentPoint = newPos;
				ptD1.posInDomain = ptD1.currentPoint; //EA add temp when they are the same 

                // transform only all positions of points of the element sharing the lable of p1 to newPos
                //cNI = c6I;
                //cNII = c6II;
                //cNIII = c6III;
                cNIII = c6I;
                cNI = c6II;
                elementInfo elTemp1 = neighbElementIII;
                neighbElementIII.points[cNI] = newPos;
                do
                {
                    elementInfo& neighbElTemp1 = eL[elTemp1.elementIndex[cNIII]];
                    for( f=0; f<=3; f++)
                    {
                       if( neighbElTemp1.pointIndex[f] == pI1 )
                       {
                           //Info<< "F I can break " << ' ' << f << endl;
                           break;
                       }
                       if( f == 3)
                       {
						   FatalErrorIn
						   (
						   " ^^^^^^^^^^^^^ we can not find such a point, pI1."
						   ) << abort(FatalError);
                       }
                    }
                    currentElement.findPointsOrder( f, cNI, cNII, cNIII);
                    neighbElTemp1.points[cNI] = newPos;
                    elTemp1 = neighbElTemp1;

                } while( elTemp1.currentIndex != neighbElementILeft.currentIndex );

                // transform all positions and indexes of points of the element sharing the p2 to newPos
                // then the p2 must be eliminated.
                //cNI = c8I;
                //cNII = c8II;
                //cNIII = c8III;
                cNIII = c8I;
                cNI = c8II;
                elementInfo elTemp2 = neighbElementIRight;
                neighbElementIRight.points[cNI] = newPos;
                neighbElementIRight.pointIndex[cNI] = pI1;

                do
                {
                    elementInfo& neighbElTemp2 = eL[elTemp2.elementIndex[cNIII]];
                    for( f=0; f<=3; f++)
                    {
                        if( neighbElTemp2.pointIndex[f] == pI2 )
                        {
                            //Info<< "S I can break " << ' ' << f << endl;
                            break;
                        }
                        if( f == 3)
                        {
						    FatalErrorIn
						    (
							" ^^^^^^^^^^^^^ we can not find such a point, pI2."
						    ) << abort(FatalError);
                        }
                    }
                    currentElement.findPointsOrder( f, cNI, cNII, cNIII);
                    neighbElTemp2.points[cNI] = newPos;
                    neighbElTemp2.pointIndex[cNI] = pI1;
                    elTemp2 = neighbElTemp2;

                } while( elTemp2.currentIndex != neighbElementII.currentIndex );

                // update the neighboures of neighbElementIII and neighbElementII
                neighbElementIII.elementIndex[c6II] = neighbElementII.currentIndex;
                neighbElementII.elementIndex[c5II] = neighbElementIII.currentIndex;

                // update the neighboures of neighbElementILeft and neighbElementIRight
                neighbElementILeft.elementIndex[c7II] = neighbElementIRight.currentIndex;
                neighbElementIRight.elementIndex[c8II] = neighbElementILeft.currentIndex;

                // p2 now is eliminated.
				pointData& ptD2 = pL[pI2];
                ptD2.keepPoint = false;

                // eleminating the elements: currentElement and neighbElementI
                currentElement.keepElement = false;
                neighbElementI.keepElement = false;

                // updatong the all elOwner
                //pL[pI3].elOwner = neighbElementIII.currentIndex;
                //pL[pI1].elOwner = neighbElementIII.currentIndex;
                //pL[pI4].elOwner = neighbElementIRight.currentIndex;
            }
        }//3

    }//2
    Info << "                 After coarsening " << numberOfELDeleted << " element(s) were deleted." << endl;
}


//--------------------------------------------------------------------------
//-----------------------------front REFINING AND COARSENING----------------
//----------------------------------AUX funtions----------------------------


//- note: func to detect small element. order the lengths, s1 shortest, s3 longest
template<class CloudType>
inline bool Foam::FrontTrackingCloud<CloudType>::isElementSmall
(	
	elementInfo el, 
	scalar minEdge, 
	scalar maxEdge,
    scalar maxAspectRatio, 
	scalar& s1, 
	scalar& s2,
    scalar& s3, 
	label& n1, 
	label& n2, 
	label& n3
)
{
    scalar ss = 0;
    label nn = 0;

    if( s1 >= s2)
    {
        ss = s1;
        s1 = s2;
        s2 = ss;
        nn = n1;
        n1 = n2;
        n2 = nn;
    }
    if( s2 >= s3)
    {
        ss = s2;
        s2 = s3;
        s3 = ss;
        nn = n2;
        n2 = n3;
        n3 = nn;
    }
    if( s1 >= s2)
    {
        ss = s1;
        s1 = s2;
        s2 = ss;
        nn = n1;
        n1 = n2;
        n2 = nn;
    }

    if( s1 < minEdge )
    {
        return true;
    }

    return false;
}


//--------------------------------------------------------------------------
//------------------------bubble POST PROCCESSING---------------------------
//--------------------------------------------------------------------------

//note: func to do all postprocessing jobs.
template<class CloudType>
void Foam::FrontTrackingCloud<CloudType>::outputTimePostProcessing()
{
    //Added if for checking if it shoud write output data
    if( (this->db().time().writeTime() || this->db().time().timeIndex() == 0) && Pstream::master() ) //this->db().time().outputTime()
    {
        Info << "\n----------> Performing front data post-processing and storing ... \n" << endl;

        bubblePostProcessing();
        writeBubblePositions();
        writeBubbleVelocityAndAccelaration();
        writeBubbleSphericity();
        //writeBubbleDragCoefs();
        //writeBubbleLiftCoefs();
        //writeBubbleDiameterComponents();
        printFrontsForAllTimes();
        writeDataOfTheFrontsAtTheOutputTimeIntimeName();
    }
}


//--------------------------------------------------------------------------
//-----------------Read/WRITE FRONT DATA TO FILE ---------------------------
//--------------------------------------------------------------------------

//note: func to print the front data for all times.
template<class CloudType>
void Foam::FrontTrackingCloud<CloudType>::printFrontsForAllTimes()
{
    fileName outfileName = "printFrontsForAllTimes.plt";

    // This is the directory keeping the bubble data.
    fileName FTResult = "FTResult";

    fileName outfileNamePath=this->db().time().rootPath()/this->db().time().globalCaseName()/FTResult/outfileName;

    ofstream outFile(outfileNamePath.c_str(), ios::app);
    ostream& os = outFile;

    label bubbleNum = bDataL_.size();

    const scalar time = this->db().time().value();

    // preparing the zone name according to the time value.
    OStringStream oss;
    oss << "BubbleAtTime=" << time;
    fileName zoneName = oss.str();

    label allElCounter = 0;
    label allPtCounter = 0;
    for(int bI = 0; bI < bubbleNum ; bI++)
    {//0
        bubbleData bData;
        DynamicList<pointData>& pL = bDataL_[bI].pL;
        DynamicList<elementInfo>& eL = bDataL_[bI].eL;
        allElCounter = allElCounter + eL.size();
        allPtCounter = allPtCounter + pL.size();
    }//0

    // Write header
    os  << " VARIABLES = \"X\" \"Y\" \"Z\",\"NUM\" "
	<< "\n ZONE T=\"" << zoneName << "\""
	<< "\n N=        " << allPtCounter
	<< " E=        " << allElCounter
	<< " ZONETYPE=FETriangle"
	<< "\n DATAPACKING=POINT \n";

    for(int bI = 0; bI < bubbleNum ; bI++)
    {//0
        DynamicList<pointData> pL = bDataL_[bI].pL;
        // Write vertex coords
        forAll(pL, pI)
        {
            vector curPoint = pL[pI].currentPoint;
	    os  << curPoint.x() << token::SPACE
	        << curPoint.y() << token::SPACE
	        << curPoint.z() << token::SPACE << bI << " \n";
        }

    }//0

    for(int bI = 0; bI < bubbleNum ; bI++)
    {//0
        DynamicList<elementInfo> eL = bDataL_[bI].eL;
        forAll(eL, eI)
        {
	    os  << eL[eI].pointIndex[0]+1 << token::SPACE
	        << eL[eI].pointIndex[1]+1 << token::SPACE
	        << eL[eI].pointIndex[2]+1 << token::SPACE << " \n";
        }
    }//0

    os  << " \n";

    outFile.close();
}


//note: func to write data of all bubbles
template<class CloudType>
void Foam::FrontTrackingCloud<CloudType>::writeDataOfTheFrontsAtTheOutputTimeIntimeName()
{
    fileName outfileName = "dataOfTheFrontsAtTheCurrentTime";

    mkDir
    ( 
        this->db().time().rootPath()/this->db().time().globalCaseName()/"bubbleData"/this->db().time().timeName()
    ); 

    // This is the directory keeping the bubble data. 
    fileName outfileNamePath = this->db().time().rootPath()/this->db().time().globalCaseName()/"bubbleData"/this->db().time().timeName()/outfileName;

    ofstream outFile(outfileNamePath.c_str(), ios::out);
    ostream& os = outFile;

    label bubbleNum = bDataL_.size();

    const scalar time = this->db().time().value();

    os  << "# Time = " << time << " \n";
    os  << "#" << " \n";
    os  << "#" << " \n";
    os  << "# The number of bubbles is as follows: " << " \n";
    os  << "  " << bubbleNum << " \n";

    for(int bI = 0; bI < bubbleNum ; bI++)
    {//0
        DynamicList<pointData> pL = bDataL_[bI].pL;
        DynamicList<elementInfo> eL =  bDataL_[bI].eL;

        os  << "# Data for bubble number: " << bI << " \n";
        os  << "# The number of points and elements are as follows: " << " \n";
        os  << "  " << pL.size() << ' ' << eL.size() << " \n";
        os  << "# The outerFluidDensity/density/outerFluidViscosity/viscosity/surfaceTensionCoeff/Eotvos/Morton are as follows: " << " \n";
        os  << "  " << bDataL_[bI].outerFluidDensity << ' ' << bDataL_[bI].density << ' '
                    << bDataL_[bI].outerFluidViscosity << ' ' << bDataL_[bI].viscosity << ' '
                    << bDataL_[bI].surfaceTensionCoeff << ' ' << bDataL_[bI].EotvosNumber  << ' ' << bDataL_[bI].MortonNumber << " \n";
        os  << "# The initial bubble diameter is as follows: " << " \n";
        os  << "  " << bDataL_[bI].sphereBubbleDiameter << " \n";
        os  << "# X, Y, and Z components of point positions are as follows: " << " \n";
        forAll(pL, pI)
        {
            vector curPoint = pL[pI].currentPoint;
	    os  << curPoint.x() << token::SPACE
	        << curPoint.y() << token::SPACE
	        << curPoint.z() << token::SPACE << " \n";
        }
        os  << "# all point indexes of each element are as follows: " << "\n";
        forAll(eL, eI)
        {
	    os  << eL[eI].pointIndex[0] << token::SPACE               
	        << eL[eI].pointIndex[1] << token::SPACE
	        << eL[eI].pointIndex[2] << token::SPACE << " \n";
        }

    }//0

    outFile.close();
}


//note: func for bubble post-processing.
template<class CloudType>
void Foam::FrontTrackingCloud<CloudType>::updatePlPosInDomain
(
	bubbleData& bData
)
{
	DynamicList<pointData>& pL = bData.pL;
	forAll(pL,pI) 
	{//2
    	point pos = pL[pI].currentPoint;
        pL[pI].posInDomain = pos; //calcThePositionInDomain(pos); can be updtaed for periodic fronts
    }//2
	bData.pLPosInDomainFlag = true; 
}


template<class CloudType>
void Foam::FrontTrackingCloud<CloudType>::updateElCentrePosInDomain
(
	bubbleData& bData
)
{
	DynamicList<pointData>& pL = bData.pL;
	DynamicList<elementInfo>& eL = bData.eL;
	forAll(eL,eI)
	{//2
		label p1I = eL[eI].pointIndex[0];
		label p2I = eL[eI].pointIndex[1];
		label p3I = eL[eI].pointIndex[2];
		point pos = (pL[p1I].currentPoint + pL[p2I].currentPoint + pL[p3I].currentPoint)/3.0;
	    eL[eI].centrePosInDomain = pos; //EA calcThePositionInDomain(pos); can be updtaed for periodic fronts
	}//2
	bData.elCPosInDomainFlag = true;
}


//note: func for bubble post-processing.
template<class CloudType>
void Foam::FrontTrackingCloud<CloudType>::bubblePostProcessing()
{

    //Info << "\n----------> post-processing for all bubbles is doing  ....... " << endl;

    scalar time = this->db().time().value();
    //scalar deltaTime = this->db().time().deltaTValue();

    forAll(bDataL_,bDI)
    {//1
		//necessary updates
		if (!bDataL_[bDI].elSurfaceFlag) bDataL_[bDI].updateBubbleSurface();
		if (!bDataL_[bDI].volumeFlag) bDataL_[bDI].updateBubbleVolume();

		DynamicList<pointData>& pL = bDataL_[bDI].pL;
		DynamicList<elementInfo>& eL = bDataL_[bDI].eL;

		bDataL_[bDI].sphereBubbleDiameter = pow(6.0*bDataL_[bDI].volume/mathematical::pi, 1./3.);
	   	bDataL_[bDI].sphericity = mathematical::pi * sqr( bDataL_[bDI].sphereBubbleDiameter)/bDataL_[bDI].surfaceArea;
	   	bDataL_[bDI].mass = bDataL_[bDI].density * bDataL_[bDI].volume;
	   	bDataL_[bDI].gravityForce = this->g().value() * bDataL_[bDI].mass;
	   	bDataL_[bDI].buoyantForce = this->g().value() * bDataL_[bDI].outerFluidDensity * bDataL_[bDI].volume;

	   	bDataL_[bDI].pressureForce = vector::zero;
	   	bDataL_[bDI].shearStressForce = vector::zero;
	   	bDataL_[bDI].centre = vector::zero;
	    //scalar sumSurfaceAreaX = 0.0;
	    //scalar sumSurfaceAreaY = 0.0;
	    //scalar sumSurfaceAreaZ = 0.0;

	   	forAll(eL, elementI)
		{//4
			elementInfo& elInfo = eL[elementI];

	        //vector elementSurfaceArea = elInfo.elementSurfaceArea;
	        //sumSurfaceAreaX = sumSurfaceAreaX + mag(elementSurfaceArea.x());
	        //sumSurfaceAreaY = sumSurfaceAreaY + mag(elementSurfaceArea.y());
	        //sumSurfaceAreaZ = sumSurfaceAreaZ + mag(elementSurfaceArea.z());

	        bDataL_[bDI].pressureForce = bDataL_[bDI].pressureForce - elInfo.pressure * elInfo.elementSurfaceArea;

	        bDataL_[bDI].shearStressForce = bDataL_[bDI].shearStressForce + (elInfo.shearStress & elInfo.elementSurfaceArea);

	        vector elementCentre = ( pL[elInfo.pointIndex[0]].currentPoint
	                                         + pL[elInfo.pointIndex[1]].currentPoint
	                                                    + pL[elInfo.pointIndex[2]].currentPoint )/3.0;

	        //bDataL_[bDI].centre.x() += elementCentre.x() * mag(elementSurfaceArea.x());
	        //bDataL_[bDI].centre.y() += elementCentre.y() * mag(elementSurfaceArea.y());
	        //bDataL_[bDI].centre.z() += elementCentre.z() * mag(elementSurfaceArea.z());
	        bDataL_[bDI].centre += (elementCentre & elInfo.elementSurfaceArea) * elementCentre; //EA8
	    }//4

	    //bDataL_[bDI].centre.x() = bDataL_[bDI].centre.x()/sumSurfaceAreaX;
	    //bDataL_[bDI].centre.y() = bDataL_[bDI].centre.y()/sumSurfaceAreaY;
	    //bDataL_[bDI].centre.z() = bDataL_[bDI].centre.z()/sumSurfaceAreaZ;
	    bDataL_[bDI].centre /= (4 * bDataL_[bDI].volume); //EA8

/*
		// calculating the flow direction:
		if( periodicalOption_ == "rotationalPeriodicalInCurvedDuct" || periodicalOption_ == "combinedRotationalPeriodicalInCurvedDuct" ||
		        periodicalOption_ == "combinedTranslationalPeriodicalInCurvedDuct" || periodicalOption_ == "translationalPeriodicalInCurvedDuct" )
		{
			rCentre_.z() =  bDataL_[bDI].centre.z();
			vector rStar = bDataL_[bDI].centre - rCentre_;

		    normalToflowDI_ = rStar/mag(rStar);

			flowD_.x() = -rStar.y();
			flowD_.y() =  rStar.x();
			flowD_.z() =  0;
			flowD_ = flowD_/mag(flowD_);

			normalToflowDII_.x() = 0;
			normalToflowDII_.y() = 0;
			normalToflowDII_.z() = 1;
		}
		else if( periodicalOption_ == "straightDuctPeriodical" || periodicalOption_ == "none")
		{
	        // for none option, flowD is parallel to z axis
	        //                , normalIToflowD is parallel to x axis
	        //                , normalIIToflowD is parallel to y axis
			flowD_.x() = 0;
			flowD_.y() = 0;
			flowD_.z() = 1;

			normalToflowDI_.x() = 1;
			normalToflowDI_.y() = 0;
			normalToflowDI_.z() = 0;

			normalToflowDII_.x() = 0;
			normalToflowDII_.y() = 1;
			normalToflowDII_.z() = 0;
	    }

	   	bDataL_[bDI].pressureDragForce = bDataL_[bDI].pressureForce & flowD_;
	   	bDataL_[bDI].frictionDragForce = bDataL_[bDI].shearStressForce & flowD_;
	    bDataL_[bDI].totalDragForce = bDataL_[bDI].pressureDragForce + bDataL_[bDI].frictionDragForce;

	   	bDataL_[bDI].pressureLiftForceI = bDataL_[bDI].pressureForce & normalToflowDI_;
	   	bDataL_[bDI].pressureLiftForceII = bDataL_[bDI].pressureForce & normalToflowDII_;
	   	bDataL_[bDI].frictionLiftForceI = bDataL_[bDI].shearStressForce & normalToflowDI_;
	   	bDataL_[bDI].frictionLiftForceII = bDataL_[bDI].shearStressForce & normalToflowDII_;

	    bDataL_[bDI].totalLiftForceI = bDataL_[bDI].pressureLiftForceI + bDataL_[bDI].frictionLiftForceI;
	    bDataL_[bDI].totalLiftForceII = bDataL_[bDI].pressureLiftForceII + bDataL_[bDI].frictionLiftForceII;
*/
	   	bDataL_[bDI].velocity = (bDataL_[bDI].centre - bDataL_[bDI].centreOld)/(time-timeOld_+ROOTVSMALL);
	   	bDataL_[bDI].accelaration = (bDataL_[bDI].velocity - bDataL_[bDI].velocityOld)/(time-timeOld_+ROOTVSMALL);
	  	bDataL_[bDI].velocityOld = bDataL_[bDI].velocity;
	   	bDataL_[bDI].centreOld = bDataL_[bDI].centre;
		timeOld_ = time;

/*
	    scalar sphereBubbleSurfaceArea = mathematical::pi * sqr( bDataL_[bDI].sphereBubbleDiameter);
	    scalar magVelocityDiff = mag( magFlowVelocity_ * flowD_ - bDataL_[bDI].velocity);

	   	bDataL_[bDI].dragCoeff = bDataL_[bDI].totalDragForce
		                 /(0.5 * bDataL_[bDI].outerFluidDensity * (pow(magVelocityDiff,2.0)+ROOTVSMALL) * sphereBubbleSurfaceArea);

	   	bDataL_[bDI].LiftCoeffI = bDataL_[bDI].totalLiftForceI
		                 /(0.5 * bDataL_[bDI].outerFluidDensity * (pow(magVelocityDiff,2.0)+ROOTVSMALL) * sphereBubbleSurfaceArea);

	   	bDataL_[bDI].LiftCoeffII = bDataL_[bDI].totalLiftForceII
		                 /(0.5 * bDataL_[bDI].outerFluidDensity * (pow(magVelocityDiff,2.0)+ROOTVSMALL) * sphereBubbleSurfaceArea);

	   // calculation of the droplet’s maximum extension in each coordinate direction, i.e
	   // flowD_, normalToflowDI_, and normalToflowDII_.
	   scalar maxDiameterComponent1 = 0.0; // it's in the direction of flowD_
	   scalar maxDiameterComponent2 = 0.0; // it's in the direction of normalToflowDI_
	   scalar maxDiameterComponent3 = 0.0; // it's in the direction of normalToflowDII_
	   for(int i = 0; i < pL.size() ; i++)
	   {
	       for(int j = i+1; j < pL.size() ; j++)
	       {
	           vector frontPointToPoint = pL[i].currentPoint - pL[j].currentPoint;
	           scalar diameterComponent1 = frontPointToPoint & flowD_;
	           scalar diameterComponent2 = frontPointToPoint & normalToflowDI_;
	           scalar diameterComponent3 = frontPointToPoint & normalToflowDII_;

	           if (  diameterComponent1 > maxDiameterComponent1)
	           {
	               maxDiameterComponent1 = diameterComponent1;
	           }
	           if (  diameterComponent2 > maxDiameterComponent2)
	           {
	               maxDiameterComponent2 = diameterComponent2;
	           }
	           if (  diameterComponent3 > maxDiameterComponent3)
	           {
	               maxDiameterComponent3 = diameterComponent3;
	           }
	       }
	   }
	   bDataL_[bDI].diameter[0] = maxDiameterComponent1;
	   bDataL_[bDI].diameter[1] = maxDiameterComponent2;
	   bDataL_[bDI].diameter[2] = maxDiameterComponent3;
*/
    }//1
}


//note: func to print the Bubble Positions for all times.
template<class CloudType>
void Foam::FrontTrackingCloud<CloudType>::writeBubblePositions()
{
    //Info << "\n"  << "----------> Now writeBubblePositions is done ........ " << endl;

    fileName outfileName = "BubblePositions.plt";

    // This is the directory keeping the bubble data.
    fileName FTResult = "FTResult";

    fileName outfileNamePath=this->db().time().rootPath()/this->db().time().globalCaseName()/FTResult/outfileName;

    ofstream outFile(outfileNamePath.c_str(), ios::app);
    ostream& os = outFile;

    label bubbleNum = bDataL_.size();

    const scalar time = this->db().time().value();

    if( this->db().time().timeIndex() == 0 )
    {
        os  << " VARIABLES = \"time\" ";
        for(int bI = 0; bI < bubbleNum ; bI++)
        {//0
           OStringStream ossX;
           ossX << "bubble_" << bI << "_CentreX";
           fileName bCentreX = ossX.str();

           OStringStream ossY;
           ossY << "bubble_" << bI << "_CentreY";
           fileName bCentreY = ossY.str();

           OStringStream ossZ;
           ossZ << "bubble_" << bI << "_CentreZ";
           fileName bCentreZ = ossZ.str();

           os  << ",\"" << bCentreX << "\", \"" << bCentreY << "\", \"" <<bCentreZ << "\"";
        }//0
        os  << " \n";
    }

    os  << time << setw(15);
    for(int bI = 0; bI < bubbleNum ; bI++)
    {//0
	os  << bDataL_[bI].centre.x() << setw(15)
	    << bDataL_[bI].centre.y() << setw(15)
	    << bDataL_[bI].centre.z() << token::SPACE;
    }//0
    os << " \n";

    outFile.close();
}

//note: func to print the Bubble velocity and accelaration for all times.
template<class CloudType>
void Foam::FrontTrackingCloud<CloudType>::writeBubbleVelocityAndAccelaration()
{
    //Info << "\n"  << "----------> Now writeBubbleVelocityAndAccelaration is done ........ " << endl;

    fileName outfileName = "BubbleVelocityAndAccelaration.plt";

    // This is the directory keeping the bubble data.
    fileName FTResult = "FTResult";

    fileName outfileNamePath=this->db().time().rootPath()/this->db().time().globalCaseName()/FTResult/outfileName;

    ofstream outFile(outfileNamePath.c_str(), ios::app);
    ostream& os = outFile;

    label bubbleNum = bDataL_.size();

    const scalar time = this->db().time().value();

    if( this->db().time().timeIndex() == 0 )
    {
        os  << " VARIABLES = \"time\" ";
        for(int bI = 0; bI < bubbleNum ; bI++)
        {//0

           OStringStream ossBV;
           ossBV << "bulkVelocity";
           fileName bBulkVelocity = ossBV.str();

           OStringStream ossUX;
           ossUX << "bubble_" << bI << "_CentreUX";
           fileName bCentreUX = ossUX.str();

           OStringStream ossUY;
           ossUY << "bubble_" << bI << "_CentreUY";
           fileName bCentreUY = ossUY.str();

           OStringStream ossUZ;
           ossUZ << "bubble_" << bI << "_CentreUZ";
           fileName bCentreUZ = ossUZ.str();

           OStringStream ossAX;
           ossAX << "bubble_" << bI << "_CentreAX";
           fileName bCentreAX = ossAX.str();

           OStringStream ossAY;
           ossAY << "bubble_" << bI << "_CentreAY";
           fileName bCentreAY = ossAY.str();

           OStringStream ossAZ;
           ossAZ << "bubble_" << bI << "_CentreAZ";
           fileName bCentreAZ = ossAZ.str();

           os  << ",\"" << bBulkVelocity << "\", \"" <<  bCentreUX << "\", \"" << bCentreUY << "\", \"" <<bCentreUZ << "\", "
               << "\"" << bCentreAX << "\", \"" << bCentreAY << "\", \"" <<bCentreAZ << "\"";
        }//0
        os  << " \n";
    }

    os  << time << setw(15);
    for(int bI = 0; bI < bubbleNum ; bI++)
    {//0
	os  << magFlowVelocity_ << setw(15)
	    << bDataL_[bI].velocity.x() << setw(15)
	    << bDataL_[bI].velocity.y() << setw(15)
	    << bDataL_[bI].velocity.z() << setw(15)
	    << bDataL_[bI].accelaration.x() << setw(15)
	    << bDataL_[bI].accelaration.y() << setw(15)
	    << bDataL_[bI].accelaration.z() << token::SPACE;
    }//0
    os << " \n";

    outFile.close();
}

//note: func to print the Bubble Positions for all times.
template<class CloudType>
void Foam::FrontTrackingCloud<CloudType>::writeBubbleSphericity()
{
    //Info << "\n"  << "----------> Now writeBubbleSphericity is done ........ " << endl;

    fileName outfileName = "BubbleSphericity.plt";

    // This is the directory keeping the bubble data.
    fileName FTResult = "FTResult";

    fileName outfileNamePath=this->db().time().rootPath()/this->db().time().globalCaseName()/FTResult/outfileName;

    ofstream outFile(outfileNamePath.c_str(), ios::app);
    ostream& os = outFile;

    label bubbleNum = bDataL_.size();

    const scalar time = this->db().time().value();

    if( this->db().time().timeIndex() == 0 )
    {
        os  << " VARIABLES = \"time\" ";
        for(int bI = 0; bI < bubbleNum ; bI++)
        {//0
           OStringStream ossSphericity;
           ossSphericity << "bubble_" << bI << "_Sphericity";
           fileName bSphericity = ossSphericity.str();

           OStringStream ossSurfaceArea;
           ossSurfaceArea << "bubble_" << bI << "_SurfaceArea";
           fileName bSurfaceArea = ossSurfaceArea.str();

           OStringStream ossVolume;
           ossVolume << "bubble_" << bI << "_Volume";
           fileName bVolume = ossVolume.str();

           os  << ",\"" << bSphericity << "\", \"" << bSurfaceArea << "\", \"" << bVolume << "\"";
        }//0
        os  << " \n";
    }

    os  << time << setw(15);
    for(int bI = 0; bI < bubbleNum ; bI++)
    {//0
	os  << bDataL_[bI].sphericity << setw(15)
	    << bDataL_[bI].surfaceArea << setw(15)
	    << bDataL_[bI].volume << token::SPACE;
    }//0
    os << " \n";

    outFile.close();
}

/*
//note: func to print the Bubble Positions for all times.
template<class CloudType>
void Foam::FrontTrackingCloud<CloudType>::writeBubbleDragCoefs()
{
    //Info << "\n"  << "----------> Now writeBubbleDragCoefs is done ........ " << endl;

    fileName outfileName = "BubbleDragCoefs.plt";

    // This is the directory keeping the bubble data.
    fileName FTResult = "FTResult";

    fileName outfileNamePath=this->owner().db().time().rootPath()/this->owner().db().time().globalCaseName()/FTResult/outfileName;

    ofstream outFile(outfileNamePath.c_str(), ios::app);
    ostream& os = outFile;

    label bubbleNum = bDataL_.size();

    const scalar time = this->owner().db().time().value();
    const scalar deltaTime = this->owner().db().time().deltaTValue();

    if( this->db().time().timeIndex() == 0 )
    {
        os  << " VARIABLES = \"time\" ";
        for(int bI = 0; bI < bubbleNum ; bI++)
        {//0
           OStringStream ossPDC;
           ossPDC << "bubble_" << bI << "_PressureDragCoef";
           fileName bPressureDragCoef = ossPDC.str();

           OStringStream ossFDC;
           ossFDC << "bubble_" << bI << "_FrictionDragCoef";
           fileName bFrictionDragCoef = ossFDC.str();

           OStringStream ossTDC;
           ossTDC << "bubble_" << bI << "_TotalDragCoef";
           fileName bTotalDragCoef = ossTDC.str();

           os  << ",\"" << bPressureDragCoef << "\", \"" << bFrictionDragCoef << "\", \"" << bTotalDragCoef << "\"";
        }//0
        os  << " \n";
    }

    os  << time << setw(15);
    for(int bI = 0; bI < bubbleNum ; bI++)
    {//0
        scalar denaminator = bDataL_[bI].totalDragForce/bDataL_[bI].dragCoeff + ROOTVSMALL;
	os  << bDataL_[bI].pressureDragForce/denaminator << setw(15)
	    << bDataL_[bI].frictionDragForce/denaminator << setw(15)
	    << bDataL_[bI].dragCoeff << token::SPACE;
    }//0
    os << " \n";

    outFile.close();
}

//note: func to print the Bubble Positions for all times.
template<class CloudType>
void Foam::FrontTrackingCloud<CloudType>::writeBubbleLiftCoefs()
{
    //Info << "\n"  << "----------> Now writeBubbleLiftCoefs is done ........ " << endl;

    fileName outfileName = "BubbleLiftCoefs.plt";

    // This is the directory keeping the bubble data.
    fileName FTResult = "FTResult";

    fileName outfileNamePath=this->owner().db().time().rootPath()/this->owner().db().time().globalCaseName()/FTResult/outfileName;

    ofstream outFile(outfileNamePath.c_str(), ios::app);
    ostream& os = outFile;

    label bubbleNum = bDataL_.size();

    const scalar time = this->owner().db().time().value();
    const scalar deltaTime = this->owner().db().time().deltaTValue();

    if( this->db().time().timeIndex() == 0 )
    {
        os  << " VARIABLES = \"time\" ";
        for(int bI = 0; bI < bubbleNum ; bI++)
        {//0
           OStringStream ossPLCI;
           ossPLCI << "bubble_" << bI << "_PressureLiftCoefI";
           fileName bPressureLiftCoefI = ossPLCI.str();

           OStringStream ossPLCII;
           ossPLCII << "bubble_" << bI << "_PressureLiftCoefII";
           fileName bPressureLiftCoefII = ossPLCII.str();

           OStringStream ossFLCI;
           ossFLCI << "bubble_" << bI << "_FrictionLiftCoefI";
           fileName bFrictionLiftCoefI = ossFLCI.str();

           OStringStream ossFLCII;
           ossFLCII << "bubble_" << bI << "_FrictionLiftCoefII";
           fileName bFrictionLiftCoefII = ossFLCII.str();

           OStringStream ossTLCI;
           ossTLCI << "bubble_" << bI << "_TotalLiftCoefI";
           fileName bTotalLiftCoefI = ossTLCI.str();

           OStringStream ossTLCII;
           ossTLCII << "bubble_" << bI << "_TotalLiftCoefII";
           fileName bTotalLiftCoefII = ossTLCII.str();

           os  << ",\"" << bPressureLiftCoefI << "\", \"" << bPressureLiftCoefII << "\", \""
                        << bFrictionLiftCoefI << "\", \"" << bFrictionLiftCoefII << "\", \""
                        << bTotalLiftCoefI << "\", \"" << bTotalLiftCoefII << "\"";
        }//0
        os  << " \n";
    }

    os  << time << setw(15);
    for(int bI = 0; bI < bubbleNum ; bI++)
    {//0
        scalar denaminator = bDataL_[bI].totalDragForce/bDataL_[bI].dragCoeff + ROOTVSMALL;
	os  << bDataL_[bI].pressureLiftForceI/denaminator << setw(15)
	    << bDataL_[bI].pressureLiftForceII/denaminator << setw(15)
	    << bDataL_[bI].frictionLiftForceI/denaminator << setw(15)
	    << bDataL_[bI].frictionLiftForceII/denaminator << setw(15)
	    << bDataL_[bI].LiftCoeffI << setw(15)
	    << bDataL_[bI].LiftCoeffII << token::SPACE;
    }//0
    os << " \n";

    outFile.close();
}

//note: func to print the Bubble Diameter Components.
template<class CloudType>
void Foam::FrontTrackingCloud<CloudType>::writeBubbleDiameterComponents()
{
    //Info << "\n"  << "----------> Now writeBubbleDiameterComponents is done ........ " << endl;

    fileName outfileName = "BubbleDiameterComponents.plt";

    // This is the directory keeping the bubble data.
    fileName FTResult = "FTResult";

    fileName outfileNamePath=this->owner().db().time().rootPath()/this->owner().db().time().globalCaseName()/FTResult/outfileName;

    ofstream outFile(outfileNamePath.c_str(), ios::app);
    ostream& os = outFile;

    label bubbleNum = bDataL_.size();

    const scalar time = this->owner().db().time().value();
    const scalar deltaTime = this->owner().db().time().deltaTValue();

    if( this->db().time().timeIndex() == 0 )
    {
        os  << " VARIABLES = \"time\" ";
        for(int bI = 0; bI < bubbleNum ; bI++)
        {//0
           OStringStream ossDI;
           ossDI << "bubble_" << bI << "_DI";
           fileName bDI = ossDI.str();

           OStringStream ossDII;
           ossDII << "bubble_" << bI << "_DII";
           fileName bDII = ossDII.str();

           OStringStream ossDIII;
           ossDIII << "bubble_" << bI << "_DIII";
           fileName bDIII = ossDIII.str();

           os  << ",\"" << bDI << "\", \"" << bDII << "\", \"" << bDIII << "\"";
        }//0
        os  << " \n";
    }

    os  << time << setw(15);
    for(int bI = 0; bI < bubbleNum ; bI++)
    {//0
	os  << bDataL_[bI].diameter[0] << setw(15)
	    << bDataL_[bI].diameter[1] << setw(15)
	    << bDataL_[bI].diameter[2] << token::SPACE;
    }//0
    os << " \n";

    outFile.close();
}
*/


//EA rem
/*
template<class CloudType>
void Foam::FrontTrackingCloud<CloudType>::info()
{
    CloudType::info();

    tmp<volScalarField> alpha = this->theta();

    const scalar alphaMin = gMin(alpha().primitiveField());
    const scalar alphaMax = gMax(alpha().primitiveField());

    Info<< "    Min cell volume fraction        = " << alphaMin << endl;
    Info<< "    Max cell volume fraction        = " << alphaMax << endl;

    if (alphaMax < small)
    {
        return;
    }

    scalar nMin = great;

    forAll(this->mesh().cells(), celli)
    {
        const label n = this->cellOccupancy()[celli].size();

        if (n > 0)
        {
            const scalar nPack = n*alphaMax/alpha()[celli];

            if (nPack < nMin)
            {
                nMin = nPack;
            }
        }
    }

    reduce(nMin, minOp<scalar>());

    Info<< "    Min dense number of parcels     = " << nMin << endl;
}
*/

//EA add
//--------------------------------------------------------------------------
template<class CloudType>
void Foam::FrontTrackingCloud<CloudType>::readInitialBubbles()
{
	//const scalar time = this->db().time().value();

	cloudCleaning();
	constructingInitialFront(); 
    surfaceTensionModel_->calculate(); //EA 
    communicationsOfTheFrontAndEulerianGrid();
    mapFrontToParcel();  
	outputTimePostProcessing();   
    this->writeFields();

	timeOld_ = this->db().time().value(); //before main loop in solver.C, this is the start time ( equal to Old time).
}


template<class CloudType>
void Foam::FrontTrackingCloud<CloudType>::cloudCleaning()
{
	Info << "\n---------> Removing all cloud parcels ... " << endl;

	forAllIter(typename FrontTrackingCloud<CloudType>, *this, iter1)
	{
	   parcelType& p1 = iter1();
	   this->deleteParticle(p1);
	}
}


//funvtion to read all bubbles front mesh.
template<class CloudType>
void Foam::FrontTrackingCloud<CloudType>::constructingInitialFront()
{ 
	if (initialMakingFileOfTheBubbles_ == "externalMesh" )
	{

		Info << "\n---------> At the time: " << this->db().time().value()
			 <<" , Bubble construction proccess from an external mesh  ... " << endl;

		fileName infileName = "initialFrontMesh";

		// This is the directory keeping the initial bubble data.
		fileName FTResult = "FTResult";
		fileName InfileNamePath=this->db().time().rootPath()/this->db().time().globalCaseName()/"FTResult"/infileName;
		Info << InfileNamePath;
		IFstream inFile(InfileNamePath);

		label bubbleNum;
		pointData ptI; //tempoaray pointData var
		elementInfo elI; //tempoaray elementInfo var
		label ptNum, elNum; //tempoaray vars
		scalar x, y, z; //tempoaray vars
		label poI, poII, poIII; //tempoaray vars

		string line = getLineNoComment(inFile);
		IStringStream lineStream(line);
		lineStream >> bubbleNum;
		Info << "			The total number of bubbles is " << bubbleNum << endl;

		List<fileName> printAllbubbles(bubbleNum);
		forAll(printAllbubbles,pAB)
		{
			// note, this way is special for OpenFOAM using the OStringStream object
			OStringStream oss;
			oss << "bubble" << pAB << "_AtZeroTime" << ".plt";
			fileName fName = oss.str();
			printAllbubbles[pAB] =  this->db().time().rootPath()/this->db().time().globalCaseName()/FTResult/fName;
		}
		for(int bI = 0; bI < bubbleNum ; bI++)
		{//00
			bubbleData bData;
			DynamicList<pointData>& pL = bData.pL;
			DynamicList<elementInfo>& eL = bData.eL;

			pL.reserve( 100000);
			eL.reserve( 200000);

			string line = getLineNoComment(inFile);
			IStringStream lineStream(line);
			lineStream >> ptNum >> elNum;

			Info <<"			pt Numbers:" << ptNum << "  elementNo:" << elNum << endl;

			string line2 = getLineNoComment(inFile);
			IStringStream lineStream2(line2);
			lineStream2 >> bData.outerFluidDensity >> bData.density
				        >> bData.outerFluidViscosity >> bData.viscosity
				        >> bData.surfaceTensionCoeff >> bData.EotvosNumber >> bData.MortonNumber;
			
			Info << "			density: " << bData.density <<"   Viscosity: "<< bData.viscosity << "    morton: " << bData.MortonNumber << endl ; 
			
			string line3 = getLineNoComment(inFile);
			IStringStream lineStream3(line3);
			lineStream3 >> bData.sphereBubbleDiameter;
			
			Info <<"			Spehere bubble Diameter: " << bData.sphereBubbleDiameter << endl;

			for(int I = 0; I < ptNum ; I++)
			{//11
				string line = getLineNoComment(inFile);
				{
				    IStringStream lineStream(line);
				    lineStream >> x >> y >> z;
				}
				ptI.currentPoint = point(x, y, z);
				ptI.currentIndex = I;
				ptI.keepPoint = true;
				pL.append(ptI);
			}//11

			for(int II = 0; II < elNum ; II++)
			{//12
				string line = getLineNoComment(inFile);
				{
				    IStringStream lineStream(line);
				    lineStream >> poI >> poII >> poIII;
				}
				elI.keepElement = true;
				elI.currentIndex = II;
				elI.pointIndex[0] = poI;
				elI.pointIndex[1] = poII;
				elI.pointIndex[2] = poIII;
				eL.append(elI);
			}//13

			bData.startPtIndexFromZero();
			bData.findNighboures();  
			bData.pointToElementMapping();
			bData.ptConnectedPoints();
			bData.printInitialFronts(printAllbubbles[bI]);

			bData.updateBubbleSurface();
			bData.updateBubbleVolume();
			bData.volume0 = bData.volume;
			updatePlPosInDomain(bData);
			updateElCentrePosInDomain(bData);

			bDataL_.append(bData);
		}//00
	}
	else if (initialMakingFileOfTheBubbles_ == "initialMakingFileOfTheBubbles" )
	{
		 fileName infileName = initialMakingFileOfTheBubbles_;
		 
		 Info << "\n---------> At the time: " << this->db().time().value()
			 <<" , bubble construction proccess from sphere bubble parameters  ... " << endl;

		// This is the directory keeping the bubble data.
		fileName FTResult = "FTResult";

		fileName InfileNamePath =  this->db().time().rootPath()/this->db().time().globalCaseName()/FTResult/infileName;

		IFstream inFile(InfileNamePath);

		label bubbleNum;
		scalar radin;
		label nps;
		scalar xc, yc, zc, e;// bubble centre
		scalar outerFluidDensity, density, outerFluidViscosity, viscosity;// bubble properties
		scalar sigma;// surface tension coeff
		scalar Eotvos;// Eotvos number for the current bubble
		scalar Morton;// Morton number for the current bubble

		pointData ptI; //tempoaray pointData var
		elementInfo elI; //tempoaray elementInfo var

		string line = getLineNoComment(inFile);
		IStringStream lineStream(line);
		lineStream >> bubbleNum;

		Info << "			The total number of bubbles is " << bubbleNum << endl;

		List<fileName> printAllbubbles(bubbleNum);
		forAll(printAllbubbles,pAB)
		{
			// note, this way is special for OpenFOAM using the OStringStream object
			OStringStream oss;
			oss << "bubble" << pAB << "_AtZeroTime" << ".plt";
			fileName fName = oss.str();
			printAllbubbles[pAB] =  this->db().time().rootPath()/this->db().time().globalCaseName()/FTResult/fName;
		}

		for(int bI = 0; bI < bubbleNum ; bI++)
		{//0
			bubbleData bData;
			DynamicList<pointData>& pL = bData.pL;
			DynamicList<elementInfo>& eL = bData.eL;

			pL.reserve( 100000);
			eL.reserve( 200000);

			string line = getLineNoComment(inFile);
			{
			   IStringStream lineStream(line);
			   lineStream >> radin >> nps >> xc >> yc >> zc >> e
					      >> outerFluidDensity >> density
					      >> outerFluidViscosity >> viscosity
					      >> sigma >> Eotvos >> Morton;
			}

			bData.outerFluidDensity = outerFluidDensity;
			bData.density = density;
			bData.outerFluidViscosity = outerFluidViscosity;
			bData.viscosity = viscosity;
			bData.sphereBubbleDiameter = 2.0 * radin;
			bData.surfaceTensionCoeff = sigma;
			bData.EotvosNumber = Eotvos;
			bData.MortonNumber = Morton;
// sphere nps calculator to be moved to mesh coarseining/refining
			if( frontMeshInputOption_ == "roughlyCalculated")
			{
			   scalar frontMeshLengthScale = 0;
			   nps = 0;
			   do
			   {
				   nps = nps + 1;
				   scalar bubbleSurfArea = 4.0 * mathematical::pi * pow(radin, 2.0);
				   frontMeshLengthScale = pow( 2 * bubbleSurfArea/ ( 2*nps+2*nps*(nps-1)+2*nps*nps*(4-1)), 0.5);
			   } while( 0.5*(maxEdge_+minEdge_) > frontMeshLengthScale || frontMeshLengthScale > 0.95*maxEdge_ );
			   Info << "			The front mesh LengthScale is calculated as " << frontMeshLengthScale << endl;
			}
			else if( frontMeshInputOption_ == "finerCalculated" || frontMeshInputOption_ == "manually" )
			{
			   scalar frontMeshLengthScale = 0;
			   nps = 0;
			   do
			   {
				   nps = nps + 1;
				   scalar bubbleSurfArea = 4.0 * mathematical::pi * pow(radin, 2.0);
				   frontMeshLengthScale = pow( 2 * bubbleSurfArea/ ( 2*nps+2*nps*(nps-1)+2*nps*nps*(4-1)), 0.5);
			   } while( minEdge_ > frontMeshLengthScale || frontMeshLengthScale > (1.0+howFinner_)*minEdge_ );
			   Info << "			The front mesh LengthScale is calculated as " << frontMeshLengthScale << endl;
			}
// sphere nps calculator (end)
			scalar dph = mathematical::piByTwo/double(nps);
			label nptot = 2+(4-1)*nps*nps+(nps-1)*nps+nps;
			label netot = 2*nps+2*nps*(nps-1)+2*nps*nps*(4-1);

			Info << "\n      	  Bubble" << bI << " is constructed from : " << ' ' << infileName << ' '
				 << " with radin/nps/bubblepos/e  as : " << ' ' << radin << ' '
				 << nps << ' ' << xc << ' ' << yc << ' '
				 << zc  << token::SPACE << e  << token::SPACE << endl;

		    Info << "            the outer fluid and bubble densities are : "
				  << bData.outerFluidDensity << ' ' << bData.density << endl;

		    Info << "            the outer fluid and bubble viscosities are : "
				  << bData.outerFluidViscosity << ' ' << bData.viscosity << endl;

		    Info << "            the number of points,elemenst for this bubble is : " << ' ' << nptot << ' '
				  << netot << ',' << endl;

		    //cccccccccccccccccccccccccccccccccccccccccccccccccccccc
		    //c ee is nonzero to create deformed bubbles
		    scalar ee = e;

		    //cccccccccccccccccccccccccccccccccccccccccccccccccccccc
		    //c set north and south pole----CHECK
		    scalar rad = radin-ee;

		    scalar x = xc;
		    scalar y = yc;
		    scalar z = zc - rad;
		    ptI.currentPoint = point(x, y, z);
		    ptI.currentIndex = 0;
		    ptI.keepPoint = true;
		    pL.append(ptI);

		    x = xc;
		    y = yc;
		    z = zc + rad;
		    ptI.currentPoint = point(x, y, z);
		    ptI.currentIndex = 1;
		    ptI.keepPoint = true;
		    pL.append(ptI);

		   for(int iq = 1; iq <= 4 ; iq++)
		   {//1
			 for(int i2 = 1; i2 <= nps ; i2++)
			 {//2
			   for(int i1 = 1; i1 <= nps ; i1++)
			   {//3

			   label iip = 2+(iq-1)*nps*nps+(i2-1)*nps+i1;
			   scalar phi = dph*double(i1-i2);
			   label ist = i2-1;
			   if( (i1-i2) < 0)
				   {
					 ist=i1-1;
				   }
			   scalar theta = mathematical::piByTwo*( double(iq-1) + double(ist)/( double(nps-abs(i1-i2)) + ROOTVSMALL ) );
			   rad=radin+ee*cos(2.0*phi);

			   x = xc + rad*cos(phi)*cos(theta);
			   y = yc + rad*cos(phi)*sin(theta);
			   z = zc - rad*sin(phi);
			   ptI.currentPoint = point(x, y, z);
			   ptI.currentIndex = iip - 1;
			   ptI.keepPoint = true;
			   pL.append(ptI);

			   label iie=2*i1+2*nps*(i2-1)+2*nps*nps*(iq-1);
			   label ia=iip;
			   label ib=iip+nps;
			   label ic=ib+1;
			   label id=ia+1;

			   label iqq=0;

			   if(i1 == nps)
				   {
				   iqq=iq;
				   if(iqq == 4)
					   {
					      iqq=0;
					   }
			   ic=2+iqq*nps*nps+nps+1-i2;
			   id=ic+1;
			   }

			   if(i2 == nps)
				   {
				   iqq=iq;
				   if(iqq == 4)
					   {
					      iqq=0;
					   }
			   ib=2+iqq*nps*nps+(nps+1-i1)*nps+1;
			   ic=ib-nps;
			   }

			   if((i1 == nps) && (i2 == 1))
				   {
					  id=1;
				   }

			   if((i2 == nps) && (i1 == 1))
				   {
					  ib=2;
				   }

			   elI.keepElement = true;
			   elI.currentIndex = iie-2;
			   elI.pointIndex[0] = ia-1;
			   elI.pointIndex[1] = ic-1;
			   elI.pointIndex[2] = ib-1;
			   eL.append(elI);

			   elI.keepElement = true;
			   elI.currentIndex = iie-1;
			   elI.pointIndex[0] = ia-1;
			   elI.pointIndex[1] = id-1;
			   elI.pointIndex[2] = ic-1;
			   eL.append(elI);

				}//3
			  }//2
			}//1

			bData.startPtIndexFromZero();
			bData.findNighboures();
			bData.pointToElementMapping();
			bData.ptConnectedPoints();
			bData.printInitialFronts(printAllbubbles[bI]);

			bData.updateBubbleSurface();
			bData.updateBubbleVolume();
			bData.volume0 = bData.volume;
			updatePlPosInDomain(bData);
			updateElCentrePosInDomain(bData);

			bDataL_.append(bData);
		}//0  
	}
	else if (initialMakingFileOfTheBubbles_ == "continueThePreviousRun" )
	{
		readDataOfTheFrontsAtTheCurrentTimeIntimeName();
	}
	else
	{
		FatalErrorIn ("Foam::FrontTrackingCloud<CloudType>::constructingInitialFront")
		<< "initialMakingFileOfTheBubbles_ should either be initialMakingFileOfTheBubbles, continueThePreviousRun, or externaMesh" << nl
		<< abort(FatalError);
	}

	//EA10
	Info << "\n---------> At the time: " << this->db().time().value() << endl;
	forAll(bDataL_,bDI)
	{
		Info << "			Bubble number: " << bDI << "  Initial volume = " << bDataL_[bDI].volume0 << endl;
	}

	allCellsInEeachMasket_.setSize(bDataL_.size());
	setLengthScalesNearTheFront();
	searchAllCellsInEeachMasket(); //EA12
}


//note: func to read data of all bubbles. (FROM TIME FOLDERS)
template<class CloudType>
void Foam::FrontTrackingCloud<CloudType>::readDataOfTheFrontsAtTheCurrentTimeIntimeName()
{
    Info << "\n---------> At the time: " << this->db().time().value()
         <<" , all data of the fronts is reading to continue solution  ....... " << endl;

    fileName infileName = "dataOfTheFrontsAtTheCurrentTime";

    // This is the directory keeping the bubble data.
    fileName FTResult = "FTResult";
    fileName InfileNamePath=this->db().time().rootPath()/this->db().time().globalCaseName()/"bubbleData"/this->db().time().timeName()/infileName;
    Info << InfileNamePath ;
    IFstream inFile(InfileNamePath);

    label bubbleNum;
    pointData ptI; //tempoaray pointData var
    elementInfo elI; //tempoaray elementInfo var
    label ptNum, elNum; //tempoaray vars
    scalar x, y, z; //tempoaray vars
    label poI, poII, poIII; //tempoaray vars

    string line = getLineNoComment(inFile);
    IStringStream lineStream(line);
    lineStream >> bubbleNum;

    List<fileName> printAllbubbles(bubbleNum);
    forAll(printAllbubbles,pAB)
    {
        // note, this way is special for OpenFOAM using the OStringStream object
        OStringStream oss;
		oss << "bubble" << pAB << "_AtCurrentTime" << ".plt";
        fileName fName = oss.str();
		printAllbubbles[pAB] =  this->db().time().rootPath()/this->db().time().globalCaseName()/FTResult/fName;
    }

    for(int bI = 0; bI < bubbleNum ; bI++)
    {//00
        bubbleData bData;
        DynamicList<pointData>& pL = bData.pL;
        DynamicList<elementInfo>& eL = bData.eL;

        pL.reserve( 100000);
        eL.reserve( 200000);

        string line = getLineNoComment(inFile);
        IStringStream lineStream(line);
        lineStream >> ptNum >> elNum;

        string line2 = getLineNoComment(inFile);
        IStringStream lineStream2(line2);
        lineStream2 >> bData.outerFluidDensity >> bData.density
                    >> bData.outerFluidViscosity >> bData.viscosity
                    >> bData.surfaceTensionCoeff >> bData.EotvosNumber >> bData.MortonNumber;

        string line3 = getLineNoComment(inFile);
        IStringStream lineStream3(line3);
        lineStream3 >> bData.sphereBubbleDiameter;

        for(int I = 0; I < ptNum ; I++)
        {//11
            string line = getLineNoComment(inFile);
            {
                IStringStream lineStream(line);
                lineStream >> x >> y >> z;
            }
            ptI.currentPoint = point(x, y, z);
            ptI.currentIndex = I;
            ptI.keepPoint = true;
            pL.append(ptI);
        }//11

        for(int II = 0; II < elNum ; II++)
        {//12
            string line = getLineNoComment(inFile);
            {
                IStringStream lineStream(line);
                lineStream >> poI >> poII >> poIII;
            }
            elI.keepElement = true;
            elI.currentIndex = II;
            elI.pointIndex[0] = poI;
            elI.pointIndex[1] = poII;
            elI.pointIndex[2] = poIII;
            eL.append(elI);
        }//13


		bData.startPtIndexFromZero();
		bData.findNighboures();
		bData.pointToElementMapping();
		bData.ptConnectedPoints();
		bData.printInitialFronts(printAllbubbles[bI]);

		bData.updateBubbleSurface();
		bData.updateBubbleVolume();
		bData.volume0 = bData.volume;
		updatePlPosInDomain(bData);
		updateElCentrePosInDomain(bData);

		bDataL_.append(bData);
    }//00
}


//--------------------------------------------------------------------------
//--------------Communication of the front and Eulerian grid---------------
//--------------------------------------------------------------------------
//note: func for communications of the front and Eulerian grid.
template<class CloudType>
void Foam::FrontTrackingCloud<CloudType>::communicationsOfTheFrontAndEulerianGrid()
{
    Info << "\n---------> Performing the fronts and Eulerian grid communications ... " << endl;

    //searchAllCellsInEeachMasket();
  
    frontToFieldModel_->calculate();     
}


//note: func for search all cells existed in eeach masket.
template<class CloudType>
void Foam::FrontTrackingCloud<CloudType>::searchAllCellsInEeachMasket()
{
    const fvMesh& mesh = this->mesh();
    const pointField& ctrs = mesh.cellCentres();
    const scalarField& cv = mesh.V();

    scalar hFactor = frontToFieldModel_->hFactor();

    forAll(bDataL_,bDI)
    {//10
		const DynamicList<pointData>& pL = bDataL_[bDI].pL;
        DynamicList<label>& L = allCellsInEeachMasket_[bDI];
        L.clearStorage();

        //scalar hl = lengthScaleNearTheFront_[bDI];
		//scalar hl = fCellL(bDI);
		scalar hl = frontToFieldModel_->bubbleHl(bDI);

        scalar xMin = +100000000.;
        scalar xMax = -100000000.;
        scalar yMin = +100000000.;
        scalar yMax = -100000000.;
        scalar zMin = +100000000.;
        scalar zMax = -100000000.;

        for( int i=0; i < pL.size(); i++)
        {
            // search for xMax and xMin
            if( pL[i].currentPoint.x() <= xMin )
            {
                xMin = pL[i].currentPoint.x();
            }
            if( pL[i].currentPoint.x() >= xMax )
            {
                xMax = pL[i].currentPoint.x();
            }

            // search for yMax and yMin
            if( pL[i].currentPoint.y() <= yMin )
            {
                yMin = pL[i].currentPoint.y();
            }
            if( pL[i].currentPoint.y() >= yMax )
            {
                yMax = pL[i].currentPoint.y();
            }

            // search for zMax and zMin
            if( pL[i].currentPoint.z() <= zMin )
            {
                zMin = pL[i].currentPoint.z();
            }
            if( pL[i].currentPoint.z() >= zMax )
            {
                zMax = pL[i].currentPoint.z();
            }
        }

		Info << "			For bubble:" << bDI << endl;
        Info << " 			xMin/xMax " << xMin << ' ' << xMax << endl;
        Info << " 			yMin/yMax " << yMin << ' ' << yMax << endl;
        Info << " 			zMin/zMax " << zMin << ' ' << zMax << endl;

        xMin = xMin - 1.5*hFactor * hl; //EA2 //EA4 
        xMax = xMax + 1.5*hFactor * hl; //EA2 //EA4 
        yMin = yMin - 1.5*hFactor * hl; //EA2 //EA4 
        yMax = yMax + 1.5*hFactor * hl; //EA2 //EA4 
        zMin = zMin - 1.5*hFactor * hl; //EA2 //EA4 
        zMax = zMax + 1.5*hFactor * hl; //EA2 //EA4 

        // calc the length scale near the front.
        scalar lengthSum = 0.0;
        scalar counter = 0;

        forAll(ctrs, cI)
        {//2
            if( (xMin <= ctrs[cI].x() && ctrs[cI].x() <= xMax) &&
                (yMin <= ctrs[cI].y() && ctrs[cI].y() <= yMax) &&
                (zMin <= ctrs[cI].z() && ctrs[cI].z() <= zMax)  )
            {
                L.append(cI);

                // calc the length scale near the front.
				if (fCellScaleOption_ == frontAverage)
				{
					if (frontToFieldModel_->cellsToRefine()[cI] == -1)
					{
						scalar volCell = cv[cI];
				        lengthSum += pow( volCell, 1.0/3.0);
				        counter += 1.0;
					}
				}
				//else if (fCellScaleOption_ == "gridLocal")
				//{
				//}
				else
				{
					scalar volCell = cv[cI];
		            lengthSum += pow( volCell, 1.0/3.0);
		            counter += 1.0;
				}
            }
        }//2

        // calc the length scale near the front.
        reduce(lengthSum, sumOp<scalar>());
        reduce(counter, sumOp<scalar>());
        if (counter > 0) lengthScaleNearTheFront_[bDI] = lengthSum/counter;
		Info << " 			lengthScaleNearTheFront " << lengthScaleNearTheFront_[bDI] << "\n" << endl;
    }//10
}

/*
template<class CloudType>
void Foam::FrontTrackingCloud<CloudType>::mapPointToParcel
(
	label bDI,
	pointData& ptN,
	pointData& ptC
)
{
	//parcelType* pPtr = new parcelType(ptC.parcel()); //with autoPtr<parcelType> parcel;
	parcelType* pPtr = new parcelType(*(ptC.parcel));
	
	pPtr->currentIndex() = ptN.currentIndex;
    pPtr->bubbleIndex() = bDI;

	this->addParticle(pPtr);

	// Non-owning reference from ptN to the same parcel
    //ptN.parcel.clear();       // Remove any old ownership //with autoPtr<parcelType> parcel;
    //ptN.parcel.set(pPtr);     // Just a view, not owner //with autoPtr<parcelType> parcel;
	ptN.parcel = pPtr; //dos not transfer the ownership for a regular pointer but autoPtr
}
*/

template<class CloudType>
void Foam::FrontTrackingCloud<CloudType>::mapPointToParcel
(
	label bDI,
	pointData& ptD
	//label& cellI,
	//label& tetFaceI,
	//label& tetPtI,
	//label& posInList 
)
{
	const fvMesh& mesh = this->mesh();

	const volVectorField& cellCentres = mesh.C();

	// Determine the injection position and owner cell,
    // tetFace and tetPt
    label celli = -1;
    label tetFacei = -1;
    label tetPti = -1;

	vector position = ptD.posInDomain;

 	mesh.findCellFacePt
    (
        position,
        celli,
        tetFacei,
        tetPti
    );

	label proci = -1;

    if (celli >= 0)
    {
        proci = Pstream::myProcNo();
    }

 	reduce(proci, maxOp<label>());

    // Ensure that only one processor attempts to insert this Parcel

    if (proci != Pstream::myProcNo())
    {
        celli = -1;
        tetFacei = -1;
        tetPti = -1;
    }

    // Last chance - find nearest cell and try that one - the point is
    // probably on an edge
    if (proci == -1)
    {
        celli = mesh.findNearestCell(position);

        if (celli >= 0)
        {
            position += small*(cellCentres[celli] - position);

            mesh.findCellFacePt
            (
                position,
                celli,
                tetFacei,
                tetPti
            );

            if (celli > 0)
            {
                proci = Pstream::myProcNo();
            }
        }

        reduce(proci, maxOp<label>());

        if (proci != Pstream::myProcNo())
        {
            celli = -1;
            tetFacei = -1;
            tetPti = -1;
        }
    }

	if (proci == -1)
    {
        FatalErrorInFunction
            << "Cannot find parcel injection processor. "
            << "Parcel position = " << ptD.posInDomain << nl
            << exit(FatalError);
    }

    if (celli > -1)
    {
        // Apply corrections to position for 2-D cases
        meshTools::constrainToMeshCentre(mesh, position);

        // Create a new parcel
        parcelType* pPtr = new parcelType(mesh, position, celli);

        // Check/set new parcel thermo properties
        this->setParcelThermoProperties(*pPtr, 0.0);

        // Check/set new parcel injection properties
        this->checkParcelProperties(*pPtr, 0.0, false);

        // Apply correction to velocity for 2-D cases
        meshTools::constrainDirection
        (
            mesh,
            mesh.solutionD(),
            pPtr->U()
        );

   		pPtr->currentIndex() = ptD.currentIndex;
        pPtr->bubbleIndex() = bDI;
//t		pPtr->shadowPos() = ptD.currentPoint; 

		this->addParticle(pPtr);
    }
}

/*
template<class CloudType>
void Foam::FrontTrackingCloud<CloudType>::mapPointToParcel
(
	label bDI,
	pointData& ptD
	//label& cellI,
	//label& tetFaceI,
	//label& tetPtI,
	//label& posInList 
)
{
	label cellI, tetFaceI, tetPtI, posInList = -1;
	const fvMesh& mesh = this->mesh();

	vector pos = ptD.posInDomain;
	mesh.findCellFacePt
	(
		pos,
		cellI,
		tetFaceI,
		tetPtI
	);

    if (Pstream::parRun()) 
    {//2
		// gather cellI from all proc
		labelList LcellI(Pstream::nProcs());
		LcellI[Pstream::myProcNo()] = cellI;
		Pstream::gatherList(LcellI);
		Pstream::scatterList(LcellI);

	    // gather tetFaceI from all proc
		labelList LtetFaceI(Pstream::nProcs());
		LtetFaceI[Pstream::myProcNo()] = tetFaceI;
		Pstream::gatherList(LtetFaceI);
		Pstream::scatterList(LtetFaceI);

	    // gather tetPtI from all proc
		labelList LtetPtI(Pstream::nProcs());
		LtetPtI[Pstream::myProcNo()] = tetPtI;
		Pstream::gatherList(LtetPtI);
		Pstream::scatterList(LtetPtI);

		forAll(LcellI,LcI)
		{
	 	    if( LcellI[LcI] > -1 && LtetFaceI[LcI] > -1 )
	        {
	            posInList = LcI;
	            break;
	        }
    	}

		if (  Pstream::myProcNo() == posInList )
		{
	         parcelType* pPtr = new parcelType(mesh, pos, LcellI[posInList]); //, LtetFaceI[posInList], LtetPtI[posInList]);
	         //Check/set new parcel thermo properties
	         this->setParcelThermoProperties(*pPtr, 0.0);
	         //Apply correction to velocity for 2-D cases
	         meshTools::constrainDirection(mesh, mesh.solutionD(), pPtr->U());
	         //Check/set new parcel injection properties
	         this->checkParcelProperties(*pPtr, 0.0, false);
	         //Apply correction to velocity for 2-D cases
	         meshTools::constrainDirection
			 (
				 mesh,
				 mesh.solutionD(),
				 pPtr->U()
			 );
       		 pPtr->currentIndex() = ptD.currentIndex;
             pPtr->bubbleIndex() = bDI;
//t		     pPtr->shadowPos() = ptD.currentPoint;  
     		 this->addParticle(pPtr);
			 // Non-owning reference from ptN to the same parcel
			 //ptD.parcel.clear();       // Remove any old ownership //with autoPtr<parcelType> parcel;
			 //ptD.parcel.set(pPtr);     // Just a view, not owner //with autoPtr<parcelType> parcel;
//			 ptD.parcel = pPtr; //dos not transfer the ownership for a regular pointer but autoPtr
 		}
	}//2
	
	if (!Pstream::parRun())
    {//5
		if ( cellI > -1)
		{
			parcelType* pPtr = new parcelType(mesh, pos, cellI); //, tetFaceI, tetPtI);
			//Check/set new parcel thermo properties
			this->setParcelThermoProperties(*pPtr, 0.0);
			//Apply correction to velocity for 2-D cases
			meshTools::constrainDirection(mesh, mesh.solutionD(), pPtr->U());
			// Check/set new parcel injection properties
			this->checkParcelProperties(*pPtr, 0.0, false);
			// Apply correction to velocity for 2-D cases
			meshTools::constrainDirection
			(
				mesh,
				mesh.solutionD(),
				pPtr->U()
			);
			pPtr->currentIndex() = ptD.currentIndex;
   	        pPtr->bubbleIndex() = bDI;
//			pPtr->shadowPos() = ptD.currentPoint;
			this->addParticle(pPtr);
			// Non-owning reference from ptN to the same parcel
			//ptD.parcel.clear();       // Remove any old ownership //with autoPtr<parcelType> parcel;
			//ptD.parcel.set(pPtr);     // Just a view, not owner //with autoPtr<parcelType> parcel;
//			ptD.parcel = pPtr; //dos not transfer the ownership for a regular pointer but autoPtr
		}
	}//5
}
*/


template<class CloudType>
void Foam::FrontTrackingCloud<CloudType>::mapFrontToParcel()
{
    Info << "\n---------> Mapping all fronts points to parcels ... " << endl;

    // insert correctly the particles................
    forAll(bDataL_,bDI)
    {//10
		//necessary updates
		if (!bDataL_[bDI].pLPosInDomainFlag) updatePlPosInDomain(bDataL_[bDI]);

        DynamicList<pointData>& PL = bDataL_[bDI].pL; 

		//label cellI = -1;
		//label tetFaceI = -1;
		//label tetPtI = -1;
		//label posInList = -1;

		forAll(PL, pointI) //EA check pass by ref
		{//3
			pointData& ptD = PL[pointI]; 
			mapPointToParcel(bDI,ptD); //,cellI,tetFaceI,tetPtI,posInList); 
    	}//3
    }//10
}

template<class CloudType>
template<class TrackCloudType>
void Foam::FrontTrackingCloud<CloudType>::updateParcelsFromFront
(
    TrackCloudType& cloud,
    typename parcelType::trackingData& td
)
{
	Info << "\n---------> Updating parcels from fronts ... " << endl;

	td.part() = parcelType::trackingData::tpTrackDist;
	//CloudType::move(cloud, td, this->db().time().deltaTValue());
	//this->updateCellOccupancy();	
	CloudType::motion(cloud, td);
}

template<class CloudType>
template<class TrackCloudType>
void Foam::FrontTrackingCloud<CloudType>::frontCloudAdvection
(
    TrackCloudType& cloud,
    typename parcelType::trackingData& td
)
{
	// advection
	Info << "\n----------> Performing front advection ... " << endl;
	td.part() = parcelType::trackingData::tpTrackTime;
	CloudType::motion(cloud, td);
	//CloudType::move(cloud, td, this->db().time().deltaTValue());
	//this->updateCellOccupancy();

	mapParcelToFront();
}

template<class CloudType>
void Foam::FrontTrackingCloud<CloudType>::mapParcelToFront()
{
    Info << "\n---------> Updating fronts from parcels ... " << endl;

    if (Pstream::parRun())
    {//2
		DynamicList<label>  myTranferLabelList;
		DynamicList<label>  myTranferBLabelList;
		DynamicList<vector> myTranferPointList;

		forAllIter(typename FrontTrackingCloud<CloudType>, *this, iter1)
		{
			parcelType& p1 = iter1();
			myTranferLabelList.append( p1.currentIndex());
			myTranferBLabelList.append( p1.bubbleIndex());
			myTranferPointList.append( p1.position());
		}

		DynamicList<label>  totalTranferLabelList;
		DynamicList<label>  totalTranferBLabelList;
		DynamicList<vector> totalTranferVectorList;

		// gather/scatter for myTranferLabelList
		List<DynamicList<label> > allLL(Pstream::nProcs());
		allLL[Pstream::myProcNo()] = myTranferLabelList;
		Pstream::gatherList(allLL);
		Pstream::scatterList(allLL);
		// gather/scatter for myTranferBLabelList
		List<DynamicList<label> > allBLL(Pstream::nProcs());
		allBLL[Pstream::myProcNo()] = myTranferBLabelList;
		Pstream::gatherList(allBLL);
		Pstream::scatterList(allBLL);
		// gather/scatter for myTranferPointList
		List<DynamicList<vector> > allVL(Pstream::nProcs());
		allVL[Pstream::myProcNo()] = myTranferPointList;
		Pstream::gatherList(allVL);
		Pstream::scatterList(allVL);

		// preparing the totalTranferLabelList & totalTranferVectorList
		forAll(allBLL,BLI)
		{
			DynamicList<label> tempLL = allLL[BLI];
			DynamicList<label> tempBLL = allBLL[BLI];
			DynamicList<vector> tempVL = allVL[BLI];

			forAll(tempBLL,tBLI)
			{
				totalTranferLabelList.append(tempLL[tBLI]);
				totalTranferBLabelList.append(tempBLL[tBLI]);
				totalTranferVectorList.append(tempVL[tBLI]);
			}
		}

		// updating the pL
		forAll(totalTranferBLabelList,tTBLI)
		{
	        label bubbleDI = totalTranferBLabelList[tTBLI];
	        DynamicList<pointData>& pL = bDataL_[bubbleDI].pL;

			pointData& ptD = pL[totalTranferLabelList[tTBLI]];
			ptD.currentPoint = totalTranferVectorList[tTBLI];
		}

		totalTranferLabelList.clearStorage();
		totalTranferVectorList.clearStorage();
		myTranferLabelList.clearStorage();
		myTranferPointList.clearStorage();
    }//2

    if (!Pstream::parRun())
    {//3
		forAllIter(typename FrontTrackingCloud<CloudType>, *this, iter1)
		{
			parcelType& p1 = iter1();

	        label bubbleDI = p1.bubbleIndex();
	        label pointDI = p1.currentIndex();

	        DynamicList<pointData>& pL = bDataL_[bubbleDI].pL;
			pointData& ptD = pL[pointDI];

			ptD.currentPoint = p1.position();
		}
    }//3

	setMotionFlags();
}


//EA9
/*
template<class CloudType>
void Foam::FrontTrackingCloud<CloudType>::thresholdSet()
{
    Info << "\n---------> Imposing the wall-distance treshold ... " << endl;

    forAll(bDataL_,bDI)
    {
        DynamicList<pointData>& pL = bDataL_[bDI].pL;

		forAll(pL,pI) //EA threshold alternative
		{//2
			if( pL[pI].currentPoint.z() < thresholdTh_ )
			{
				pL[pI].currentPoint.z() = thresholdTh_;
			}
		}//2
    }
}
*/

/*
//--------------------------------------------------------------------------
//-----------------------------print runtime data---------------------------
//--------------------------------------------------------------------------
//note: func to print all important runTime data. PEDI
template<class CloudType>
void Foam::FrontTrackingCloud<CloudType>::printingAllImportantRunTimeData()
{
    const fvMesh& mesh = this->mesh();
    const pointField& ctrs = mesh.cellCentres();

    // Integrate flow variables over cell set
    scalar magUbarAve = 0.0;
    scalar rhoVol = 0.0;
    const scalarField& cv = mesh.V();
    forAll(ctrs, i)
    {
        label cellI = i;//cells_[i];
        scalar volCell = cv[cellI];
        rhoVol += this->rho().primitiveFieldRef()[cellI] * volCell;
        magUbarAve += (flowDList_[cellI] & this->U().primitiveFieldRef()[cellI])
                      *this->rho().primitiveFieldRef()[cellI]*volCell;
    }

    // Collect across all processors
    reduce(magUbarAve, sumOp<scalar>());
    reduce(rhoVol, sumOp<scalar>());

    // Volume averages
    magUbarAve /= rhoVol;

    magFlowVelocity_ = magUbarAve;

    Info << "\n----------> At the time: " << this->owner().db().time().value()
         << " due to the pressure gradient value of " << gradPAtflowDirection_
         << " the average velocity at the flow direction is " << magUbarAve << endl;

    scalar basedDensity = bDataL_[0].outerFluidDensity;
    scalar basedViscosity = bDataL_[0].outerFluidViscosity;
    ReDh_ = basedDensity * magFlowVelocity_* Dh_/basedViscosity;
    Info << "\n----------> At the time: " << this->owner().db().time().value()
         << " the reynold number, ReDh_ is calculted as " << ReDh_ << endl;
}


//note: func to calculate the position in domain.
template<class CloudType>
Foam::vector Foam::FrontTrackingCloud<CloudType>::calcThePositionInDomain( point pos)
{
    vector posInDomain;

    if( periodicalOption_ == "rotationalPeriodicalInCurvedDuct")
    {//0
		rCentre_.z() = pos.z();
        vector rp = pos - rCentre_;
		scalar ptTheta = acos( rp.x()/mag(rp));
		scalar ptThetaInDomain = 0.0;
		if( rp.y() >= 0 )
		{
			ptThetaInDomain = ptTheta - thetaP_ * floor(ptTheta/thetaP_);
		}
		else if( rp.y() < 0 )
		{
			ptTheta = mathematical::twoPi - ptTheta;
			ptThetaInDomain = ptTheta - thetaP_ * floor(ptTheta/thetaP_);
		}
		// calculating the posInDoamain
		posInDomain.x() = mag(rp) * cos(ptThetaInDomain);
		posInDomain.y() = mag(rp) * sin(ptThetaInDomain);
		posInDomain.z() = pos.z();
        return posInDomain;
    }//0
    else if( periodicalOption_ == "translationalPeriodicalInCurvedDuct")
    {
		rCentre_.z() = pos.z();
        vector rp = pos - rCentre_;
		// calculating the posInDoamain
		posInDomain = pos;
		posInDomain.z() = rp.z() - tPLength_ * floor(rp.z()/tPLength_);
        return posInDomain;
    }
    else if( periodicalOption_ == "combinedRotationalPeriodicalInCurvedDuct" || periodicalOption_ == "combinedTranslationalPeriodicalInCurvedDuct" )
    {
		rCentre_.z() = pos.z();
        vector rp = pos - rCentre_;
		scalar ptTheta = acos( rp.x()/mag(rp));
		scalar ptThetaInDomain = 0.0;
		if( rp.y() >= 0 )
		{
			ptThetaInDomain = ptTheta - thetaP_ * floor(ptTheta/thetaP_);
		}
		else if( rp.y() < 0 )
		{
			ptTheta = mathematical::twoPi - ptTheta;
			ptThetaInDomain = ptTheta - thetaP_ * floor(ptTheta/thetaP_);
		}
		// calculating the posInDoamain
		posInDomain.x() = mag(rp) * cos(ptThetaInDomain);
		posInDomain.y() = mag(rp) * sin(ptThetaInDomain);
        // this is done for the translationalPeriodical part
		posInDomain.z() = rp.z() - tPLength_ * floor(rp.z()/tPLength_);
        return posInDomain;
    }
    else if( periodicalOption_ == "straightDuctPeriodical")
    {
		// calculating the posInDoamain
		posInDomain = pos;
		posInDomain.z() = pos.z() - sDPLength_ * floor(pos.z()/sDPLength_);
        return posInDomain;
    }
    else if( periodicalOption_ == "none")
    {
		// calculating the posInDoamain
		posInDomain = pos;
        return posInDomain;
    }
    else
    {
        FatalErrorIn ("Foam::vector Foam::FrontTrackingCloud<CloudType>::calcThePositionInDomain( point pos)")
            << "you should choose periodicalOption correctly" << nl
            << abort(FatalError);
    }
}

//--------------------------------------------------------------------------
//-----------------------------flow direction field-------------------------
//--------------------------------------------------------------------------
//note: func to calculate flow direction field.
template<class CloudType>
void Foam::FrontTrackingCloud<CloudType>::flowDirectionField()
{
    const fvMesh& mesh = this->owner().mesh();
    const pointField& ctrs = mesh.cellCentres();

    // calculating the flow direction:
    if( periodicalOption_ == "rotationalPeriodicalInCurvedDuct" || periodicalOption_ == "combinedRotationalPeriodicalInCurvedDuct" )
    {
        forAll(ctrs, cellI)
        {
	    rCentre_.z() =  ctrs[cellI].z();
	    vector rStar = ctrs[cellI] - rCentre_;

	    flowD_.x() = -rStar.y();
	    flowD_.y() =  rStar.x();
	    flowD_.z() =  0;
	    flowD_ = flowD_/mag(flowD_);
            this->owner().dPressureFromFT()[cellI] = flowD_ * gradPAtflowDirection_;
            flowDList_.append(flowD_);
        }
    }
    if( periodicalOption_ == "translationalPeriodicalInCurvedDuct" || periodicalOption_ == "combinedTranslationalPeriodicalInCurvedDuct")
    {
        forAll(ctrs, cellI)
        {
	    flowD_.x() = 0;
	    flowD_.y() = 0;
	    flowD_.z() = 1;
            this->owner().dPressureFromFT()[cellI] = flowD_ * gradPAtflowDirection_;
            flowDList_.append(flowD_);
        }
    }
    else if( periodicalOption_ == "straightDuctPeriodical" || periodicalOption_ == "none" )
    {
        // for none option, flowD is parallel to z axis
        //                , normalIToflowD is parallel to x axis
        //                , normalIIToflowD is parallel to y axis
        forAll(ctrs, cellI)
        {
	    flowD_.x() = 0;
	    flowD_.y() = 0;
	    flowD_.z() = 1;
            this->owner().dPressureFromFT()[cellI] = flowD_ * gradPAtflowDirection_;
            flowDList_.append(flowD_);
        }
    }
}
*/

// ************************************************************************* //
