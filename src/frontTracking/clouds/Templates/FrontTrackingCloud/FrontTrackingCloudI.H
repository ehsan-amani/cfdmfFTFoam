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

#include "NoFrontToField.H" //EA d2

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType>
inline const Foam::FrontTrackingCloud<CloudType>&
Foam::FrontTrackingCloud<CloudType>::cloudCopy() const
{
    return cloudCopyPtr_();
}


template<class CloudType> //EA
inline void
Foam::FrontTrackingCloud<CloudType>::updateDensityFromFT
(
	volScalarField& rho
)
{	
	//if (!isType<FrontToFieldModels::NoFrontToField<CloudType>>(frontToFieldModel_()))
	//if (frontToFieldModel_().type() != "none")
	if (!isType<FrontToFieldModels::NoFrontToField<FrontTrackingCloudType>>(frontToFieldModel_()))
	{
		//Info << "\n---------> At the time: " << this->db().time().value()
        // <<" , I am here2   .......\n " << endl;
		//Info << "Dynamic type: " << typeid(frontToFieldModel_()).name() << endl;
		//Info << "Expected type: " << typeid(FrontToFieldModels::NoFrontToField<CloudType>).name() << endl;

		if (this->frontToFieldModel().twoFluidFlow())
		{
			volScalarField& Id = this->frontToFieldModel().IndicatorRef();
		    dimensionedScalar oneDensity 
		    ( 
		        "oneDensity", 
		        dimensionSet(1,-3,0,0,0,0,0), 
		        scalar(1.0) 
		    );
		    
		    rho = oneDensity * ((1.0 - Id) * bDataL_[0].density + Id * bDataL_[0].outerFluidDensity);     
		}
		else
		{
			rho = this->frontToFieldModel().densityIndicator(); 
		}
	}
}

//EA added
template<class CloudType> //EA
inline void
Foam::FrontTrackingCloud<CloudType>::updateViscosityFromFT
(
	volScalarField& mu
)
{
	//if (!isType<FrontToFieldModels::NoFrontToField<CloudType>>(frontToFieldModel_()))
	//if (frontToFieldModel_().type() != "none")
	if (!isType<FrontToFieldModels::NoFrontToField<FrontTrackingCloudType>>(frontToFieldModel_()))
	{
		if (this->frontToFieldModel().twoFluidFlow())
		{
			volScalarField& Id = this->frontToFieldModel().IndicatorRef(); 
		    dimensionedScalar oneViscosity 
		    ( 
		        "oneViscosity", 
		        dimensionSet(1,-1,-1,0,0,0,0), 
		        scalar(1.0) 
		    );
		    
			mu = oneViscosity*((1.0 - Id) * bDataL_[0].viscosity + Id * bDataL_[0].outerFluidViscosity);      
		}
		else
		{		    
			mu = this->frontToFieldModel().viscosityIndicator(); 
		}
	}
}

//EA added
template<class CloudType>
inline Foam::tmp<Foam::DimensionedField<Foam::vector, Foam::volMesh>>
Foam::FrontTrackingCloud<CloudType>::sTensionForceFromFT() const
{
    return sTensionForceFromFT_();
}

//EA added
template<class CloudType>
inline Foam::DimensionedField<Foam::vector, Foam::volMesh>&
Foam::FrontTrackingCloud<CloudType>::sTensionForceFromFTRef()
{
    return sTensionForceFromFT_();
}


//EA added
template<class CloudType>
inline Foam::tmp<Foam::DimensionedField<Foam::vector, Foam::volMesh>>
Foam::FrontTrackingCloud<CloudType>::pressureJumpAtTheInterfaceFromFT() const
{
    return pressureJumpAtTheInterfaceFromFT_();
}

//EA added
template<class CloudType>
inline Foam::DimensionedField<Foam::vector, Foam::volMesh>&
Foam::FrontTrackingCloud<CloudType>::pressureJumpAtTheInterfaceFromFTRef()
{
    return pressureJumpAtTheInterfaceFromFT_();
}

template<class CloudType>
inline Foam::scalar Foam::FrontTrackingCloud<CloudType>::lengthScaleOfTheMesh() const
{
    return lengthScaleOfTheMesh_;
}

template<class CloudType>
inline Foam::List<DynamicList<label>> Foam::FrontTrackingCloud<CloudType>::allCellsInEeachMasket() const
{
    return allCellsInEeachMasket_;
}

template<class CloudType>
inline Foam::List<scalar> Foam::FrontTrackingCloud<CloudType>::lengthScaleNearTheFront() const
{
    return lengthScaleNearTheFront_;
}


//EA rem
/*
template<class CloudType>
inline const Foam::PackingModel<Foam::FrontTrackingCloud<CloudType>>&
Foam::FrontTrackingCloud<CloudType>::packingModel() const
{
    return packingModel_();
}


template<class CloudType>
inline Foam::PackingModel<Foam::FrontTrackingCloud<CloudType>>&
Foam::FrontTrackingCloud<CloudType>::packingModel()
{
    return packingModel_();
}


template<class CloudType>
inline const Foam::DampingModel<Foam::FrontTrackingCloud<CloudType>>&
Foam::FrontTrackingCloud<CloudType>::dampingModel() const
{
    return dampingModel_();
}


template<class CloudType>
inline Foam::DampingModel<Foam::FrontTrackingCloud<CloudType>>&
Foam::FrontTrackingCloud<CloudType>::dampingModel()
{
    return dampingModel_();
}


template<class CloudType>
inline const Foam::IsotropyModel<Foam::FrontTrackingCloud<CloudType>>&
Foam::FrontTrackingCloud<CloudType>::isotropyModel() const
{
    return isotropyModel_();
}


template<class CloudType>
inline Foam::IsotropyModel<Foam::FrontTrackingCloud<CloudType>>&
Foam::FrontTrackingCloud<CloudType>::isotropyModel()
{
    return isotropyModel_();
}
*/

//EA add
template<class CloudType>
inline const Foam::SurfaceTensionModel<Foam::FrontTrackingCloud<CloudType>>&
Foam::FrontTrackingCloud<CloudType>::surfaceTensionModel() const
{
    return surfaceTensionModel_();
}


template<class CloudType>
inline Foam::SurfaceTensionModel<Foam::FrontTrackingCloud<CloudType>>&
Foam::FrontTrackingCloud<CloudType>::surfaceTensionModel()
{
    return surfaceTensionModel_();
}


template<class CloudType>
inline const Foam::FrontToFieldModel<Foam::FrontTrackingCloud<CloudType>>&
Foam::FrontTrackingCloud<CloudType>::frontToFieldModel() const
{
    return frontToFieldModel_();
}


template<class CloudType>
inline Foam::FrontToFieldModel<Foam::FrontTrackingCloud<CloudType>>&
Foam::FrontTrackingCloud<CloudType>::frontToFieldModel()
{
    return frontToFieldModel_();
}


template<class CloudType>
inline const Foam::VolumeCorrectionModel<Foam::FrontTrackingCloud<CloudType>>&
Foam::FrontTrackingCloud<CloudType>::volumeCorrectionModel() const
{
    return volumeCorrectionModel_();
}


template<class CloudType>
inline Foam::VolumeCorrectionModel<Foam::FrontTrackingCloud<CloudType>>&
Foam::FrontTrackingCloud<CloudType>::volumeCorrectionModel()
{
    return volumeCorrectionModel_();
}


template<class CloudType>
inline const Foam::DynamicList<bubbleData>&
Foam::FrontTrackingCloud<CloudType>::bDataL() const
{
    return bDataL_;
}


template<class CloudType>
inline Foam::DynamicList<bubbleData>&
Foam::FrontTrackingCloud<CloudType>::bDataL()
{
    return bDataL_;
}

template<class CloudType> //EA d2
inline Foam::scalar Foam::FrontTrackingCloud<CloudType>::fCellL(label bDI)
{
	if (fCellScaleOption_ == fixed )
	{
		return h_;
	}
	else
	{
		return lengthScaleNearTheFront()[bDI];
	}
}


//EA add
//note: func to calc the length scale of mesh.
template<class CloudType>
inline Foam::scalar Foam::FrontTrackingCloud<CloudType>::calcLengthScaleOfTheMesh()
{
    const fvMesh& mesh = this->mesh();
    const pointField& ctrs = mesh.cellCentres();
    const scalarField& cv = mesh.V();

    // Integrate variables over cell set
    scalar lengthSum = 0.0;
    scalar counter = 0;
    forAll(ctrs, i)
    {
        label cellI = i;
        scalar volCell = cv[cellI];
        lengthSum += pow( volCell, 1.0/3.0);
        counter += 1.0;
    }

    // Collect across all processors
    reduce(lengthSum, sumOp<scalar>());
    reduce(counter, sumOp<scalar>());

    return lengthSum/counter;
}


template<class CloudType>
inline void Foam::FrontTrackingCloud<CloudType>::setMotionFlags()
{
	forAll(bDataL_,bDI)
    {//1
		bDataL_[bDI].elPointsFlag = false;
		bDataL_[bDI].plConnectionFlag = false;
		bDataL_[bDI].elSurfaceFlag = false;
		bDataL_[bDI].volumeFlag = false;
		bDataL_[bDI].pLPosInDomainFlag = false;
		bDataL_[bDI].elCPosInDomainFlag = false;
	}//1
}


template<class CloudType>
inline void Foam::FrontTrackingCloud<CloudType>::setCRFlags()
{
	forAll(bDataL_,bDI)
    {//1
		bDataL_[bDI].plConnectionFlag = false;
		bDataL_[bDI].elSurfaceFlag = false;
		bDataL_[bDI].volumeFlag = false;
		bDataL_[bDI].elCPosInDomainFlag = false;
	}//1
}


//-note: func to detect an element required to be refined. order the lengths, s1 shortest, s3 longest
template<class CloudType>
inline bool Foam::FrontTrackingCloud<CloudType>::isElementRefined
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
    vector e01 = ( el.points[1] - el.points[0]);
    vector e02 = ( el.points[2] - el.points[0]);

    scalar ss;
    label nn;

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

    scalar surfaceArea = mag (0.5 * (e01 ^ e02));
    //scalar elAspectRatio = 0.25 * sqrt (3.0) * sqr( ( s1 + s2 + s3) / 3.) * surfaceArea;
    scalar elAspectRatio = sqr(s1 + s2 + s3)/( sqrt (3.0) * 12.0 * surfaceArea);


    if( s3 > maxEdge )
    {
        return true;
    }

    if( s1 > minEdge && elAspectRatio > maxAspectRatio )
    {
        return true;
    }

    return false;
}


//note: func to calculate non-dimensional numbers.
template<class CloudType>
inline void Foam::FrontTrackingCloud<CloudType>::calcBubbleNonDimensionalNumbers()
{
    scalar gravity = mag (this->g().value() );
    scalar basedDensity = bDataL_[0].outerFluidDensity;
    scalar basedViscosity = bDataL_[0].outerFluidViscosity;

    forAll(bDataL_,bDI)
    {//1
        scalar db = bDataL_[bDI].sphereBubbleDiameter;
        bDataL_[bDI].surfaceTensionCoeff = basedDensity*pow(db, 2.0)*gravity/bDataL_[bDI].EotvosNumber;
        scalar sigma =bDataL_[bDI].surfaceTensionCoeff;
        bDataL_[bDI].MortonNumber = gravity*pow(basedViscosity, 4)/basedDensity/pow(sigma, 3);
    }
}


template<class CloudType>
inline void Foam::FrontTrackingCloud<CloudType>::setLengthScalesNearTheFront()
{
    lengthScaleNearTheFront_.setSize(bDataL_.size());
	minEdgeB_.setSize(bDataL_.size());
	maxEdgeB_.setSize(bDataL_.size());

    forAll(bDataL_,bDI)
    {//1
        lengthScaleNearTheFront_[bDI] = lengthScaleOfTheMesh_; //h_;
		scalar hl = fCellL(bDI); 
    	minEdgeB_[bDI] = factorForMinEdge_ * hl;
        maxEdgeB_[bDI] = factorForMaxEdge_ * hl;
    }//1
}




// ************************************************************************* //
