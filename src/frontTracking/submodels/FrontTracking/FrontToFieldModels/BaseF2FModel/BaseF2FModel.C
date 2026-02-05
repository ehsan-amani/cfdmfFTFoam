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

#include "BaseF2FModel.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::FrontToFieldModels::BaseF2FModel<CloudType>::BaseF2FModel
(
    const dictionary& dict,
    CloudType& owner, //EA d2
	const word& modelName //EA d2
)
:
    FrontToFieldModel<CloudType>(dict, owner, modelName), //EA d2
	spreadingOption_(bubbleLocal), 
	front2field_(this->coeffDict().lookup("front2field")),
	p_(16), 
    q_(16),
    r_(16),
    beta_(4),
	numberOfIndicatorFiltering_(this->coeffDict().template lookupOrDefault<label>("numberOfIndicatorFiltering", 0)), //EA d2
	filteringWeight_(this->coeffDict().template lookup<scalar>("filteringWeight")) //EA d2
{
	const word spreadingOption = this->coeffDict().lookup("spreadingOption");

    if (spreadingOption == "bubbleLocal")
    {
        spreadingOption_ = bubbleLocal;
		this->coeffDict().lookup("hFactor") >> this->hFactor_; 
    }
    else if (spreadingOption == "fixed")
    {
        spreadingOption_ = fixed;
		this->coeffDict().lookup("h") >> this->h_;
		this->coeffDict().lookup("hFactor") >> this->hFactor_; 
    }
    //else if (spreadingOption == "gridLocal")
    //{
    //    spreadingOption_ = gridLocal;
    //}
    else
    {
        FatalErrorInFunction
            << "spreadingOption must be either 'fixed', 'bubbleLocal' or 'gridLocal (future)'"
            << nl << exit(FatalError);
    }

	Info << "	Selecting spreading option " << spreadingOption << endl;
    //Info << "---------> The lengh scale of the mesh is " << this->owner().lengthScaleOfTheMesh() << endl; //EA d2
    Info << "---------> The current value of parameter h is " << this->h_ << endl;
	Info << "---------> The radius of influence is " << this->hFactor_ << "h"<< endl;
	Info << "---------> The number of indicator filtering is " << this->numberOfIndicatorFiltering_ << endl;
	Info << "---------> The filtering weight is " << this->filteringWeight_ << endl;

	adjustingTheImportantVarsForIndicators();
}


template<class CloudType>
Foam::FrontToFieldModels::BaseF2FModel<CloudType>::BaseF2FModel
(
    const BaseF2FModel<CloudType>& cm
)
:
    FrontToFieldModel<CloudType>(cm),
	spreadingOption_(cm.spreadingOption_),
	front2field_(cm.front2field_),
	p_(cm.p_), 
    q_(cm.q_),
    r_(cm.r_),
    beta_(cm.beta_),
	numberOfIndicatorFiltering_(cm.numberOfIndicatorFiltering_),
	filteringWeight_(cm.filteringWeight_)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class CloudType>
Foam::FrontToFieldModels::BaseF2FModel<CloudType>::~BaseF2FModel()
{}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType>
void Foam::FrontToFieldModels::BaseF2FModel<CloudType>::calculate()
{
	this->F2FCommunication();
	this->IndicatorsConstruction();
	this->IndicatorsFiltering();
}


template<class CloudType>
void Foam::FrontToFieldModels::BaseF2FModel<CloudType>::F2FCommunication()
{
	//Inputs: pL.currentPoint, eL.centerPosInDomain, el.elementSurfaceArea, eL.averageSurfaceTension
	//Outputs (updated): 

	simpleMatrix<scalar> mComponent(4);
    mComponent.source()[0] = 1.0;
    mComponent.source()[1] = 0.0;
    mComponent.source()[2] = 0.0;
    mComponent.source()[3] = 0.0;	

	DynamicList<label> cellsInSphereDomain;

	//DynamicList<bubbleData<CloudType>>& bDataL_= this->owner().bDataL();
	DynamicList<bubbleData>& bDataL_= this->owner().bDataL();

	const fvMesh& mesh = this->owner().mesh(); 
    const pointField& ctrs = mesh.cellCentres();

    volSymmTensorField shearStress = 2*this->owner().mu()
                                      *symm( fvc::grad(this->owner().U()) ); 
    this->cellsToRefine_.field() = 0.0;
	this->owner().sTensionForceFromFTRef().field() = vector::zero; 
    this->owner().pressureJumpAtTheInterfaceFromFTRef().field() = vector::zero; 
	this->initFields();

    forAll(bDataL_,bDI)
    {//10
		//necessary updates
		if (!bDataL_[bDI].elSurfaceFlag) bDataL_[bDI].updateBubbleSurface();
		if (!bDataL_[bDI].elCPosInDomainFlag) this->owner().updateElCentrePosInDomain(bDataL_[bDI]); 

        DynamicList<elementInfo>& eL = bDataL_[bDI].eL;
        DynamicList<label> ctrsLabel = this->owner().allCellsInEeachMasket()[bDI]; //EA d2 

		scalar hl = this->bubbleHl(bDI); //EA d2

		forAll(eL, elementI)
		{//4
			elementInfo& elInfo = eL[elementI];
			const point& pos = elInfo.centrePosInDomain;

            forAll(ctrsLabel, I)
            {
                const label& cI = ctrsLabel[I];
                scalar distance = mag (ctrs[cI] - pos);
                if ( distance <= this->hFactor_ * hl ) //EA2 //EA4
                {
                    cellsInSphereDomain.append(cI);
                }
            }

			if (!Pstream::parRun())
    		{//5
				if( front2field_ == "RKPM" ) //EA4
				{   
					for(int i=0; i < 4 ; i++)
					{//1
						for(int j=0; j < 4 ; j++)
						{//2
							label pqr = 4*i + j;
							mComponent[i][j] = 0.0;
							forAll(cellsInSphereDomain, cI)
							{
								label cellI = cellsInSphereDomain[cI];
								scalar dVolume = mesh.V()[cellI]; //EA
								vector del = ( ctrs[cellI] - pos)/(hl + ROOTVSMALL); //EA4
								scalar m = mPQR(del, hl, p_[pqr], q_[pqr], r_[pqr]) *dVolume; //EA //EA4
								mComponent[i][j] = mComponent[i][j] + m;
							}
						}//2
					}//1
					for (label i = 0; i < mComponent.n(); ++i) //regularization to avoid a singular matrix
					{
						mComponent[i][i] += SMALL;
					}
					beta_ = mComponent.solve(); 
				}     
				else if( front2field_ == "basic")
				{
					int i=0; int j=0;
					label pqr = 4*i + j;
					mComponent[i][j] = 0.0;
					forAll(cellsInSphereDomain, cI)
					{
						label cellI = cellsInSphereDomain[cI];
						scalar dVolume = mesh.V()[cellI]; //EA
						vector del = ( ctrs[cellI] - pos)/(hl + ROOTVSMALL); //EA4
						scalar m = mPQR(del, hl, p_[pqr], q_[pqr], r_[pqr]) *dVolume; //EA //EA4
						mComponent[i][j] = mComponent[i][j] + m;
					} 
					beta_[0]=1.0/mComponent[0][0]; beta_[1]=0.0;beta_[2]=0.0;beta_[3]=0.0; //EA3
				} 
			}
			else
			{
				List<scalar> mComponentL(16);

				if( front2field_ == "RKPM" ) //EA4
				{   
					for(int i=0; i < 4 ; i++)
					{//1
						for(int j=0; j < 4 ; j++)
						{//2
							label pqr = 4*i + j;
							mComponentL[pqr] = 0.0;
							mComponent[i][j] = 0.0;
							forAll(cellsInSphereDomain, cI)
							{
								label cellI = cellsInSphereDomain[cI];
								scalar dVolume = mesh.V()[cellI]; //EA
								vector del = ( ctrs[cellI] - pos)/(hl + ROOTVSMALL); //EA4
										 scalar m = mPQR(del, hl, p_[pqr], q_[pqr], r_[pqr]) * dVolume; //EA //EA4
								mComponentL[pqr] = mComponentL[pqr] + m;
							}
						}//2
					}//1

					// gather/scatter for mComponentL
					List<List<scalar> > mComponentLAll(Pstream::nProcs());
					mComponentLAll[Pstream::myProcNo()] = mComponentL;
					Pstream::gatherList(mComponentLAll);
					Pstream::scatterList(mComponentLAll);
					forAll(mComponentLAll,mCLAI)
					{
						List<scalar> tempMComponentL = mComponentLAll[mCLAI];
						for(int i=0; i < 4 ; i++)
						{//1
						    for(int j=0; j < 4 ; j++)
							{//2
								label pqr = 4*i + j;
								mComponent[i][j] = mComponent[i][j] + tempMComponentL[pqr];
							}//2
						}//1
					}
					for (label i = 0; i < mComponent.n(); ++i) //regularization to avoid a singular matrix
					{
						mComponent[i][i] += SMALL;
					}
					beta_ = mComponent.solve(); 
				}     
				else if( front2field_ == "basic")
				{
					int i=0; int j=0;
					label pqr = 4*i + j;
					mComponentL[pqr] = 0.0;
					mComponent[i][j] = 0.0;
					forAll(cellsInSphereDomain, cI)
					{
						label cellI = cellsInSphereDomain[cI];
						scalar dVolume = mesh.V()[cellI]; //EA
						vector del = ( ctrs[cellI] - pos)/(hl + ROOTVSMALL); //EA4
								 scalar m = mPQR(del, hl, p_[pqr], q_[pqr], r_[pqr]) * dVolume; //EA //EA4
						mComponentL[pqr] = mComponentL[pqr] + m;
					}

					// gather/scatter for mComponentL
					List<List<scalar> > mComponentLAll(Pstream::nProcs());
					mComponentLAll[Pstream::myProcNo()] = mComponentL;
					Pstream::gatherList(mComponentLAll);
					Pstream::scatterList(mComponentLAll);
					forAll(mComponentLAll,mCLAI)
					{
						List<scalar> tempMComponentL = mComponentLAll[mCLAI];
						int i=0; int j=0;
						label pqr = 4*i + j;
						mComponent[i][j] = mComponent[i][j] + tempMComponentL[pqr];
					}
					beta_[0]=1.0/mComponent[0][0]; beta_[1]=0.0;beta_[2]=0.0;beta_[3]=0.0; //EA3
				} 
			}//5
			
            elInfo.pressure = 0.0;
            elInfo.shearStress = symmTensor::zero;
			forAll(cellsInSphereDomain, cI)
			{
				label cellI = cellsInSphereDomain[cI];
				scalar dVolume = mesh.V()[cellI];
				vector del = ( ctrs[cellI] - pos)/(hl + ROOTVSMALL); //EA4
			    scalar wBar = windowFuncBar(del, hl, beta_); //EA4
	
//rem cache p	            elInfo.pressure += wBar * (this->owner().p().internalField()[cellI])*dVolume; //EA 
	            elInfo.shearStress += wBar * (shearStress.internalField()[cellI])*dVolume; //EA				            
	            //this->owner().sTensionForceFromFT()[cellI] = this->owner().sTensionForceFromFT()[cellI]
	            //                 + wBar * elInfo.averageSurfaceTension * mag(elInfo.elementSurfaceArea)/dVolume;
	            this->owner().sTensionForceFromFTRef()[cellI] += wBar * elInfo.averageSurfaceTension; // /dVolume; //EA
                this->owner().pressureJumpAtTheInterfaceFromFTRef()[cellI] += wBar * bDataL_[bDI].pressureJumpAtTheInterface * 
																		mag(elInfo.elementSurfaceArea); // /dVolume; //EA
				this->cellsToRefine_.primitiveFieldRef()[cellI] = -1.0; 
				this->collectFields(cellI, bDI, wBar, elInfo.elementSurfaceArea);
			}
			cellsInSphereDomain.clear();
    	}//4
    }//10
}


template<class CloudType>
void Foam::FrontToFieldModels::BaseF2FModel<CloudType>::IndicatorsFiltering()
{
    //DynamicList<bubbleData<CloudType>>& bDataL= this->owner().bDataL(); //EA d2

	if (this->twoFluidFlow_)
    {
		for(int i = 1; i <= this->numberOfIndicatorFiltering_ ; i++)
		{
			this->IndicatorRef() = filteringWeight_ * (fvc::average(fvc::interpolate(this->Indicator())))
				                             + (1.0 - filteringWeight_) * this->Indicator(); //EA d2
		}
		
		/* //EA d2 The main solver will use Indicator(), densityIndicator(), viscosityIndicator()
        dimensionedScalar oneDensity 
        ( 
            "oneDensity", 
            dimensionSet(1,-3,0,0,0,0,0), 
            scalar(1.0) 
        ); 
        
        dimensionedScalar oneViscosity 
        ( 
            "oneViscosity", 
            dimensionSet(1,-1,-1,0,0,0,0), 
            scalar(1.0) 
        );

		this->owner().alpha() = this->Indicator(); //EA d2
        
        this->owner().rho() = oneDensity * ((1.0 - this->Indicator()) * bDataL[0].density + this->Indicator() * bDataL[0].outerFluidDensity); //EA d2
        
		this->owner().mu()= oneViscosity*((1.0 - this->Indicator()) * bDataL[0].viscosity + this->Indicator() * bDataL[0].outerFluidViscosity); //EA d2  
		*/     
	}
	else
	{
		for(int i = 1; i <= this->numberOfIndicatorFiltering_ ; i++)
		{
			this->densityIndicatorRef() = filteringWeight_ * (fvc::average(fvc::interpolate(this->densityIndicator())))
				                             + (1.0 - filteringWeight_) * this->densityIndicator(); //EA d2
			this->viscosityIndicatorRef() = filteringWeight_ * (fvc::average(fvc::interpolate(this->viscosityIndicator())))
	                             + (1.0 - filteringWeight_) * this->viscosityIndicator(); //EA d2
		}

		//this->owner().rho() = this->densityIndicator(); //EA d2
        
		//this->owner().mu()= this->viscosityIndicator(); //EA d2
	}
}


//--------------------------------------------------------------------------
//-----------------------------Constructing the Indicators------------------
//---------------------------------utility function-------------------------

template<class CloudType> //EA d2
inline Foam::scalar Foam::FrontToFieldModels::BaseF2FModel<CloudType>::bubbleHl(label bDI)
{
	//return this->owner().fCellL(bDI);
	//return this->owner().lengthScaleNearTheFront()[bDI];
	if (this->spreadingOption_ == fixed )
	{
		return this->h_;
	}
	else
	{
		return this->owner().lengthScaleNearTheFront()[bDI];
	}
}

//- note: func to calculate the Peskin distribution for structured Mesh
template<class CloudType>
inline Foam::scalar Foam::FrontToFieldModels::BaseF2FModel<CloudType>::PeskinD(scalar r,scalar h)
{
   if( mag(r) < this->hFactor_ ) //EA4-2
   {
      return (1.0/(2.0*this->hFactor_*h))*( 1.0 + cos( mathematical::pi * r / this->hFactor_) ); //EA4-2
   }
   else
   {
      return 0.0;
   }
}

//- note: func to calculate the window function for structured Mesh
template<class CloudType>
inline Foam::scalar Foam::FrontToFieldModels::BaseF2FModel<CloudType>::windowFunc(vector del, scalar h)
{
    scalar dX = del.x();
    scalar dY = del.y();
    scalar dZ = del.z();
    return PeskinD(dX,h) * PeskinD(dY,h) * PeskinD(dZ,h);
}

//- note: func to calculate the window function for unstructured Mesh
template<class CloudType>
inline Foam::scalar Foam::FrontToFieldModels::BaseF2FModel<CloudType>::windowFuncBar(vector del, scalar h, List<scalar> beta)
{
    scalar dX = del.x();
    scalar dY = del.y();
    scalar dZ = del.z();
    return windowFunc(del, h) * ( beta[0] + beta[1] * dX + beta[2] * dY + beta[3] * dZ );
}

//- note: func to calculate the surface tension
template<class CloudType>
inline Foam::scalar Foam::FrontToFieldModels::BaseF2FModel<CloudType>::mPQR(vector del, scalar h, label p, label q, label r)
{
    scalar dX = del.x();
    scalar dY = del.y();
    scalar dZ = del.z();
    return pow(dX,p) * pow(dY,q) * pow(dZ,r) *  windowFunc(del, h);
}


//note: func for adjusting the important vars for indicators.
template<class CloudType>
inline void Foam::FrontToFieldModels::BaseF2FModel<CloudType>::adjustingTheImportantVarsForIndicators()
{
    p_[0] = 0;
    p_[1] = 1;
    p_[2] = 0;
    p_[3] = 0;
    p_[4] = 1;
    p_[5] = 2;
    p_[6] = 1;
    p_[7] = 1;
    p_[8] = 0;
    p_[9] = 1;
    p_[10] = 0;
    p_[11] = 0;
    p_[12] = 0;
    p_[13] = 1;
    p_[14] = 0;
    p_[15] = 0;

    q_[0] = 0;
    q_[1] = 0;
    q_[2] = 1;
    q_[3] = 0;
    q_[4] = 0;
    q_[5] = 0;
    q_[6] = 1;
    q_[7] = 0;
    q_[8] = 1;
    q_[9] = 1;
    q_[10] = 2;
    q_[11] = 1;
    q_[12] = 0;
    q_[13] = 0;
    q_[14] = 1;
    q_[15] = 0;

    r_[0] = 0;
    r_[1] = 0;
    r_[2] = 0;
    r_[3] = 1;
    r_[4] = 0;
    r_[5] = 0;
    r_[6] = 0;
    r_[7] = 1;
    r_[8] = 0;
    r_[9] = 0;
    r_[10] = 0;
    r_[11] = 1;
    r_[12] = 1;
    r_[13] = 1;
    r_[14] = 1;
    r_[15] = 2;
}

//note: func for calculating the exact closest point at the current element.
template<class CloudType>
inline Foam::scalar Foam::FrontToFieldModels::BaseF2FModel<CloudType>::LiuHeviside
( 
	scalar outerFluidProperty, 
	scalar property, 
	vector minDistance, 
	vector elementSurfaceArea, 
	scalar gama
)
{
    // indicators are built with the aid of a Heaviside function (Liu et al., 2000)

    if ( mag(minDistance) <= gama && (minDistance & elementSurfaceArea) >= 0.0 )
    {
        scalar phi = mag(minDistance);
        scalar inducedProperty =  0.5 * (outerFluidProperty - property)
                                      * (1.0 + phi/gama + sin(mathematical::pi * phi/gama)/mathematical::pi) + property;
        return inducedProperty;
    }
    else if ( mag(minDistance) <= gama && (minDistance & elementSurfaceArea) < 0.0 ) //EA d2
    {
        scalar phi = - mag(minDistance);
        scalar inducedProperty =  0.5 * (outerFluidProperty - property)
                                         * (1.0 + phi/gama + sin(mathematical::pi * phi/gama)/mathematical::pi) + property;
        return inducedProperty;
    }
    else if ( mag(minDistance) > gama && (minDistance & elementSurfaceArea) >= 0.0 ) //EA d2
    {
        return outerFluidProperty;
    }
    else // ( mag(minDistance) > gama && (minDistance & elementSurfaceArea) < 0.0 ) //EA d2
    {
        return property;
    }
}


// ************************************************************************* //
