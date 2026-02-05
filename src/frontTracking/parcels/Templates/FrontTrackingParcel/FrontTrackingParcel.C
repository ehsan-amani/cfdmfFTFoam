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

#include "FrontTrackingParcel.H"

// * * * * * * * * * * *  Protected Member Functions * * * * * * * * * * * * // //EA add

template<class ParcelType>
template<class TrackCloudType>
void Foam::FrontTrackingParcel<ParcelType>::setCellValues
(
    TrackCloudType& cloud,
    trackingData& td
)
{
    tetIndices tetIs = this->currentTetIndices();

    td.Uc() = td.UInterp().interpolate(this->coordinates(), tetIs);
}


template<class ParcelType>
template<class TrackCloudType>
void Foam::FrontTrackingParcel<ParcelType>::calc
(
    TrackCloudType& cloud,
    trackingData& td,
    const scalar dt
)
{
/*
    // Define local properties at beginning of time step
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    const scalar np0 = nParticle_;
    const scalar mass0 = mass();

    // Reynolds number
    const scalar Re = this->Re(td);


    // Sources
    //~~~~~~~~

    // Explicit momentum source for particle
    vector Su = Zero;

    // Linearised momentum source coefficient
    scalar Spu = 0.0;

    // Momentum transfer from the particle to the carrier phase
    vector dUTrans = Zero;
*/

    // Motion
    // ~~~~~~

    // Calculate new particle velocity
    this->U_ = td.Uc();
        //calcVelocity(cloud, td, dt, Re, td.muc(), mass0, Su, dUTrans, Spu);

/*
    // Accumulate carrier phase source terms
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    if (cloud.solution().coupled())
    {
        // Update momentum transfer
        cloud.UTransRef()[this->cell()] += np0*dUTrans;

        // Update momentum transfer coefficient
        cloud.UCoeffRef()[this->cell()] += np0*Spu;
    }
*/
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class ParcelType>
Foam::FrontTrackingParcel<ParcelType>::FrontTrackingParcel
(
    const FrontTrackingParcel<ParcelType>& p
)
:
    ParcelType(p),
    //UCorrect_(p.UCorrect_) //EA rem
    currentIndex_(p.currentIndex_), //EA add
    bubbleIndex_(p.bubbleIndex_) //EA add
{}


template<class ParcelType>
Foam::FrontTrackingParcel<ParcelType>::FrontTrackingParcel
(
    const FrontTrackingParcel<ParcelType>& p,
    const polyMesh& mesh
)
:
    ParcelType(p, mesh),
    //UCorrect_(p.UCorrect_) //EA rem
    currentIndex_(p.currentIndex_), //EA add
    bubbleIndex_(p.bubbleIndex_) //EA add
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class ParcelType>
template<class TrackCloudType>
void Foam::FrontTrackingParcel<ParcelType>::cellValueSourceCorrection
(
    TrackCloudType& cloud,
    trackingData& td,
    const scalar dt
)
{
    //td.Uc() += cloud.UTransRef()[this->cell()]/massCell(td); //EA10 no correction
}

template<class ParcelType>
template<class TrackCloudType>
bool Foam::FrontTrackingParcel<ParcelType>::move
(
    TrackCloudType& cloud,
    trackingData& td,
    const scalar trackTime
)
{
    typename TrackCloudType::parcelType& p =
        static_cast<typename TrackCloudType::parcelType&>(*this);

    switch (td.part())
    {
        case trackingData::tpTrackTime:
        {
            ParcelType::move(cloud, td, trackTime);

            break;
        }
        case trackingData::tpTrackDist:
        {
			label bubbleDI = p.bubbleIndex();
		    label pointDI = p.currentIndex();
		    //const DynamicList<pointData<TrackCloudType>>& pL = cloud.bDataL()[bubbleDI].pL;
		    const DynamicList<pointData>& pL = cloud.bDataL()[bubbleDI].pL;
			//const pointData<TrackCloudType>& ptD = pL[pointDI];
			const pointData& ptD = pL[pointDI];

			vector dist = ptD.currentPoint - p.position();

            moveDistance(cloud, td, dist);

            break;
        }
    }

    return td.keepParticle;
}


template<class ParcelType>
template<class TrackCloudType>
bool Foam::FrontTrackingParcel<ParcelType>::moveDistance
(
    TrackCloudType& cloud,
    trackingData& td,
    const vector s
)
{
    typename TrackCloudType::parcelType& p =
        static_cast<typename TrackCloudType::parcelType&>(*this);
    typename TrackCloudType::parcelType::trackingData& ttd =
        static_cast<typename TrackCloudType::parcelType::trackingData&>(td);

    ttd.switchProcessor = false;
    ttd.keepParticle = true;

    const scalarField& cellLengthScale = cloud.cellLengthScale();
    const scalar maxCo = cloud.solution().maxCo();

    while (ttd.keepParticle && !ttd.switchProcessor && p.stepFraction() < 1)
    {
        // Cache the current position, cell and step-fraction
        const point start = p.position();
        const scalar sfrac = p.stepFraction();

		//EA        
		// Total displacement over the time-step
        //const vector s = trackTime*U_;
		// Total time-step for the displacement (an approximation is sufficient)
		const scalar trackTime = (s & this->U_)/(magSqr(this->U_)+VSMALL);

        // Cell length scale
        const scalar l = cellLengthScale[p.cell()];

        // Deviation from the mesh centre for reduced-D cases
        const vector d = p.deviationFromMeshCentre();

        // Fraction of the displacement to track in this loop. This is limited
        // to ensure that the both the time and distance tracked is less than
        // maxCo times the total value.
        scalar f = 1 - p.stepFraction();
        f = min(f, maxCo);
        f = min(f, maxCo*l/max(small*l, mag(s)));
        if (p.moving())
        {
            // Track to the next face
            p.trackToFace(f*s - d, f);
        }
        else
        {
            // At present the only thing that sets moving_ to false is a stick
            // wall interaction. We want the position of the particle to remain
            // the same relative to the face that it is on. The local
            // coordinates therefore do not change. We still advance in time and
            // perform the relevant interactions with the fixed particle.
            p.stepFraction() += f;
        }

        const scalar dt = (p.stepFraction() - sfrac)*trackTime;

        // Avoid problems with extremely small timesteps
		//EA
        //if (dt > rootVSmall)
		if (mag(s) > VSMALL)
        {
            // Update cell based properties
            p.setCellValues(cloud, ttd);

            p.calcDispersion(cloud, ttd, dt);

            if (cloud.solution().cellValueSourceCorrection())
            {
                p.cellValueSourceCorrection(cloud, ttd, dt);
            }

            p.calc(cloud, ttd, dt);
        }

		//EA
        //p.age() += dt;

        if (p.moving() && p.onFace())
        {
            cloud.functions().postFace(p, ttd.keepParticle);
        }

        cloud.functions().postMove(p, dt, start, ttd.keepParticle);

        if (p.moving() && p.onFace() && ttd.keepParticle)
        {
            p.hitFace(f*s - d, f, cloud, ttd);
        }
    }

    return ttd.keepParticle;
}

// * * * * * * * * * * * * * * IOStream operators  * * * * * * * * * * * * * //

#include "FrontTrackingParcelIO.C"

// ************************************************************************* //
