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

\*---------------------------------------------------------------------------*/

//- note: func to update the surface area of each face
inline void Foam::elementInfo::updateElementSurfaceArea()
{
    vector e10 = ( points[1] - points[0]);
    vector e20 = ( points[2] - points[0]);
    elementSurfaceArea = 0.5 * (e10 ^ e20);
}


//- note: func to find cNI, cNII, and cNIII, according to point label of curP for current pointData.
inline void Foam::elementInfo::findPointsOrder( label f, label& cNI, label& cNII, label& cNIII )
{
    cNI = f;
    cNII = cNI +1;
    if( cNII > 2)
    {
        cNII = cNII -3;
    }
    cNIII = cNII +1;
    if( cNIII > 2)
    {
        cNIII = cNIII - 3;
    }
}


//- note: func to update the bubble surface or front surface.
//template<class CloudType>
//inline void Foam::bubbleData<CloudType>::updateBubbleSurface()
inline void Foam::bubbleData::updateBubbleSurface()
{
	if (!elPointsFlag) pointToElementMapping();

    surfaceArea = 0;
    forAll(eL, elementI)
    {
        elementInfo& elD = eL[elementI];
		elD.updateElementSurfaceArea();
        surfaceArea = surfaceArea + mag(elD.elementSurfaceArea);
    }
	elSurfaceFlag = true;
}


//- note: func to calculate the volume enclosed by a front.
//template<class CloudType>
//inline void Foam::bubbleData<CloudType>::updateBubbleVolume()
inline void Foam::bubbleData::updateBubbleVolume()
{
	if (!elPointsFlag) pointToElementMapping();

    volume = 0;
    scalar curTetrahedralVolume = 0;
    vector elementCenter = vector::zero;

    forAll(eL, elementI)
    {
        elementInfo& elD = eL[elementI];

        vector e01 = ( elD.points[1] - elD.points[0]);
        vector e02 = ( elD.points[2] - elD.points[0]);
        vector Sf = 0.5 * (e01 ^ e02);

        elementCenter = ( elD.points[0] + elD.points[1] + elD.points[2])/3.0;

        curTetrahedralVolume = ( elementCenter & Sf )/3.0;

        volume = volume + curTetrahedralVolume;
    }

    volume = mag(volume);

	volumeFlag = true;
}

// ************************************************************************* //
