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

#include "bubbleData.H"
//#include "Time.H"
//#include "localEulerDdtScheme.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //
/*
Foam::bubbleData::bubbleData
(
    const fvMesh& mesh,
    const dictionary& dict
)
:
    mesh_(mesh),
    dict_(dict),
    transient_(false),
    calcFrequency_(1),
    maxCo_(0.3),
    iter_(1),
    trackTime_(0),
    coupled_(false),
    cellValueSourceCorrection_(false),
    maxTrackTime_(0),
    resetSourcesOnStartup_(true),
    schemes_()
{
    read();
}


Foam::bubbleData::bubbleData
(
    const bubbleData& cs
)
:
    mesh_(cs.mesh_),
    dict_(cs.dict_),
    transient_(cs.transient_),
    calcFrequency_(cs.calcFrequency_),
    maxCo_(cs.maxCo_),
    iter_(cs.iter_),
    trackTime_(cs.trackTime_),
    coupled_(cs.coupled_),
    cellValueSourceCorrection_(cs.cellValueSourceCorrection_),
    maxTrackTime_(cs.maxTrackTime_),
    resetSourcesOnStartup_(cs.resetSourcesOnStartup_),
    schemes_(cs.schemes_)
{}


Foam::bubbleData::bubbleData
(
    const fvMesh& mesh
)
:
    mesh_(mesh),
    dict_(dictionary::null),
    transient_(false),
    calcFrequency_(0),
    maxCo_(great),
    iter_(0),
    trackTime_(0),
    coupled_(false),
    cellValueSourceCorrection_(false),
    maxTrackTime_(0),
    resetSourcesOnStartup_(false),
    schemes_()
{}

*/
// * * * * * * * * * * * * * * * * Destructors  * * * * * * * * * * * * * * * //
/*
//template<class CloudType>
//Foam::bubbleData<CloudType>::~bubbleData()
Foam::bubbleData::~bubbleData()
{}


Foam::elementInfo::~elementInfo()
{}
*/
// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


//- note: does it start point index from zero? if yes, this returns true.
//template<class CloudType>
//void Foam::bubbleData<CloudType>::startPtIndexFromZero()
void Foam::bubbleData::startPtIndexFromZero()
{
    forAll(pL, pI)
    {
        if ( pL[pI].currentIndex == pL.size() )
        {
            forAll(pL, pII)
            {
               pL[pII].currentIndex = pL[pII].currentIndex - 1;
               Info << "\n XXXXXXXXXXXXXXXXXXXXXXXX" << endl;
            }
            break;
        }
    }
    forAll(eL, eI)
    {
        label pEI = eL[eI].pointIndex[0];
        label pEII = eL[eI].pointIndex[1];
        label pEIII = eL[eI].pointIndex[2];
	if ( ( pEI == pL.size() ) || ( pEII == pL.size() ) || ( pEIII == pL.size() ) )
        {
               Info << "\n XXXXXXXXXXXXXXXXXXXXXXXX" << endl;
            forAll(eL, eII)
            {
                eL[eII].pointIndex[0] = eL[eII].pointIndex[0] - 1;
                eL[eII].pointIndex[1] = eL[eII].pointIndex[1] - 1;
                eL[eII].pointIndex[2] = eL[eII].pointIndex[2] - 1;
            }
            break;
        }
    }
}


//template<class CloudType>
//void Foam::bubbleData<CloudType>::findNighboures()
void Foam::bubbleData::findNighboures()
{
    forAll(eL , fI)
    {//1
        elementInfo& elfI = eL[fI];
        elfI.elementIndex[0] = -1;
        elfI.elementIndex[1] = -1;
        elfI.elementIndex[2] = -1;

        label pfI = elfI.pointIndex[0];
        label pfII = elfI.pointIndex[1];
        label pfIII = elfI.pointIndex[2];

        forAll(eL , sI)
        {//2
            elementInfo elsI = eL[sI];

            label psI = elsI.pointIndex[0];
            label psII = elsI.pointIndex[1];
            label psIII = elsI.pointIndex[2];

            if( sI != fI)
            {
                if(  (psII == pfI && psI == pfII) || (psIII == pfI && psII == pfII) || (psI == pfI && psIII == pfII)  )
                {
                    elfI.elementIndex[0] = elsI.currentIndex;
                }
                if(  (psII == pfII && psI == pfIII) || (psIII == pfII && psII == pfIII) || (psI == pfII && psIII == pfIII)  )
                {
                    elfI.elementIndex[1] = elsI.currentIndex;
                }
                if(  (psII == pfIII && psI == pfI) || (psIII == pfIII && psII == pfI) || (psI == pfIII && psIII == pfI)  )
                {
                    elfI.elementIndex[2] = elsI.currentIndex;
                }
            }
        }//2
    }//1
	elNeighborFlag = true;
}


// note: func to find connected points of each point in pL.
//template<class CloudType>
//void Foam::bubbleData<CloudType>::pointToElementMapping()
void Foam::bubbleData::pointToElementMapping()
{
    forAll(eL, II)
    {//4
        elementInfo& elI = eL[II];

        label p0 = elI.pointIndex[0];
        label p1 = elI.pointIndex[1];
        label p2 = elI.pointIndex[2];

        elI.points[0] = pL[p0].currentPoint;
        elI.points[1] = pL[p1].currentPoint;
        elI.points[2] = pL[p2].currentPoint;

        pL[p0].elOwner = elI.currentIndex;
        pL[p1].elOwner = elI.currentIndex;
        pL[p2].elOwner = elI.currentIndex;
    }//4
	elPointsFlag = true;
}


// note: func to find connected points of each point in pL.
//template<class CloudType>
//void Foam::bubbleData<CloudType>::ptConnectedPoints()
void Foam::bubbleData::ptConnectedPoints()
{
    //label posInList;
    label f, cNI, cNII, cNIII;
    label curPI, firstPI;

	if (!elNeighborFlag) findNighboures();
	if (!elPointsFlag) pointToElementMapping();

    forAll(pL, I)
    {//1
        //pointData<CloudType>& currentPt = pL[I];
		pointData& currentPt = pL[I];
        DynamicList<point>& cPoints = currentPt.connectedPoints;
        DynamicList<label>& cPointsIndex = currentPt.connectedPointsIndex;

        cPoints.clear();
        cPointsIndex.clear();

        curPI = currentPt.currentIndex;
        elementInfo elTemp = eL[currentPt.elOwner];

        for( f=0; f<=3; f++)
        {
            if( elTemp.pointIndex[f] == curPI )
            {
                break;
            }
            if( f == 3)
            {
				FatalErrorIn( " ^^^^^^^^^^^^^ we can not find such a point, curPI.") << abort(FatalError);
            }
        }
        elTemp.findPointsOrder( f, cNI, cNII, cNIII );

        firstPI = elTemp.pointIndex[cNII];
        cPoints.append( elTemp.points[cNII]);
        cPointsIndex.append( elTemp.pointIndex[cNII]);

        label counter = 0;
        do
        {
            elTemp = eL[elTemp.elementIndex[cNIII]];

            for( f=0; f<=3; f++)
            {
                if( elTemp.pointIndex[f] == curPI )
                {
                    break;
                }
                if( f == 3)
                {
		    		FatalErrorIn( " ^^^^^^^^^^^^^ we can not find such a point, curPI.") << abort(FatalError);
                }
            }
            elTemp.findPointsOrder(f, cNI, cNII, cNIII );

            cPoints.append( elTemp.points[cNII]);
            cPointsIndex.append( elTemp.pointIndex[cNII]);

            counter++;
            if( counter > 100)
            {
				FatalErrorIn( " >>>>>>>>>>>>   Sorry, some thing are going to be wrong with You !!!! ")
                << "\n point index is: " << curPI
                << "\n el owner is: " << currentPt.elOwner
                << "\n el owner neighboures: " << eL[currentPt.elOwner].elementIndex[0] << ' '
                                               << eL[currentPt.elOwner].elementIndex[1] << ' '
                                               << eL[currentPt.elOwner].elementIndex[2] << abort(FatalError);
            }
        }while( elTemp.pointIndex[cNIII] != firstPI );
    }//1
	plConnectionFlag = true;
}

//note: func to removal the front undulatuons.
//template<class CloudType>
//Foam::labelList Foam::bubbleData<CloudType>::adjustingTheMeshDataAfterChanges() //EA11
Foam::labelList Foam::bubbleData::adjustingTheMeshDataAfterChanges() //EA11
{
    //DynamicList<pointData<CloudType>> pLTemp = pL;
    DynamicList<pointData> pLTemp = pL;
    DynamicList<elementInfo> eLTemp = eL;

	labelList locInPointDataList(pL.size());
    labelList locInElementInfoList(eL.size());

    //resizing and renumbering the point and element lists.
    pL.clear();
    eL.clear();

    label ptCount = 0;
    forAll(pLTemp, I)
    {//2
        //pointData<CloudType>& ptI = pLTemp[I];
        pointData& ptI = pLTemp[I];
        locInPointDataList[ptI.currentIndex] = -1;
        if( ptI.keepPoint == true)
        {
            locInPointDataList[ptI.currentIndex] = ptCount;
            ptI.currentIndex = ptCount;
            pL.append(ptI);
            ptCount = ptCount + 1;
        }
    }//2

    label elCount = 0;
    forAll(eLTemp, II)
    {//3
        elementInfo& elI = eLTemp[II];
        locInElementInfoList[elI.currentIndex] = -1;
        if( elI.keepElement == true)
        {
            locInElementInfoList[elI.currentIndex] = elCount;
            elI.currentIndex = elCount;
            eL.append(elI);
            elCount = elCount + 1;
        }
    }//3

    pLTemp.clearStorage();
    eLTemp.clearStorage();

    forAll(eL, II)
    {//4
        elementInfo& elI = eL[II];

        label& p0 = elI.pointIndex[0];
        label& p1 = elI.pointIndex[1];
        label& p2 = elI.pointIndex[2];

        elI.pointIndex[0] = locInPointDataList[p0];
        elI.pointIndex[1] = locInPointDataList[p1];
        elI.pointIndex[2] = locInPointDataList[p2];

        elI.points[0] = pL[p0].currentPoint;
        elI.points[1] = pL[p1].currentPoint;
        elI.points[2] = pL[p2].currentPoint;

        label e0 = elI.elementIndex[0];
        label e1 = elI.elementIndex[1];
        label e2 = elI.elementIndex[2];

        elI.elementIndex[0] = locInElementInfoList[e0];
        elI.elementIndex[1] = locInElementInfoList[e1];
        elI.elementIndex[2] = locInElementInfoList[e2];

        pL[p0].elOwner = elI.currentIndex;
        pL[p1].elOwner = elI.currentIndex;
        pL[p2].elOwner = elI.currentIndex;
    }//4

/*
	forAll(eL, II)
    {//4
        elementInfo& elI = eL[II];

        label& p0 = elI.pointIndex[0];
        label& p1 = elI.pointIndex[1];
        label& p2 = elI.pointIndex[2];

        elI.pointIndex[0] = locInPointDataList[p0];
        elI.pointIndex[1] = locInPointDataList[p1];
        elI.pointIndex[2] = locInPointDataList[p2];

        label e0 = elI.elementIndex[0];
        label e1 = elI.elementIndex[1];
        label e2 = elI.elementIndex[2];

        elI.elementIndex[0] = locInElementInfoList[e0];
        elI.elementIndex[1] = locInElementInfoList[e1];
        elI.elementIndex[2] = locInElementInfoList[e2];
    }//4

	pointToElementMapping();
*/
    //ptConnectedPoints();

	return locInPointDataList; //EA11
}



//- note : func to print the front mesh in .plt format
//- there is an important point that the number
//- of points should be strated from 1  in plt format
//template<class CloudType>
//void Foam::bubbleData<CloudType>::printInitialFronts(const fileName& outfileName)
void Foam::bubbleData::printInitialFronts(const fileName& outfileName)
{
    ofstream outFile(outfileName.c_str(), ios::out);
    ostream& os = outFile;

    // Write header
    os  << " VARIABLES = \"X\" \"Y\" \"Z\",\"NUM\" "
	<< "\n ZONE T=\"Bubble\" "
	<< "\n N=        " << pL.size()
	<< " E=        " << eL.size()
	<< " ZONETYPE=FETriangle"
	<< "\n DATAPACKING=POINT \n";

    // Write vertex coords
    forAll(pL, pI)
    {
        vector curPoint = pL[pI].currentPoint;
	os  << curPoint.x() << token::SPACE
	    << curPoint.y() << token::SPACE
	    << curPoint.z() << token::SPACE << 1 << " \n";
    }

    forAll(eL, eI)
    {
	 os  << eL[eI].pointIndex[0]+1 << token::SPACE
	     << eL[eI].pointIndex[1]+1 << token::SPACE
	     << eL[eI].pointIndex[2]+1 << token::SPACE << " \n";
    }

    outFile.close();
}


/*
void Foam::bubbleData::read()
{
    // For transient runs the Lagrangian tracking may be transient or steady
    transient_ = dict_.lookupOrDefault("transient", false);

    // For LTS and steady-state runs the Lagrangian tracking cannot be transient
    if (transient_)
    {
        if (fv::localEulerDdt::enabled(mesh_))
        {
            IOWarningInFunction(dict_)
                << "Transient tracking is not supported for LTS"
                   " simulations, switching to steady state tracking."
                << endl;
            transient_ = false;
        }

        if (mesh_.steady())
        {
            IOWarningInFunction(dict_)
                << "Transient tracking is not supported for steady-state"
                   " simulations, switching to steady state tracking."
                << endl;
            transient_ = false;
        }
    }

    dict_.lookup("coupled") >> coupled_;
    dict_.lookup("cellValueSourceCorrection") >> cellValueSourceCorrection_;
    dict_.readIfPresent("maxCo", maxCo_);

    if (steadyState())
    {
        dict_.lookup("calcFrequency") >> calcFrequency_;
        dict_.lookup("maxTrackTime") >> maxTrackTime_;

        if (coupled_)
        {
            dict_.subDict("sourceTerms").lookup("resetOnStartup")
                >> resetSourcesOnStartup_;
        }
    }

    if (coupled_)
    {
        const dictionary&
            schemesDict(dict_.subDict("sourceTerms").subDict("schemes"));

        wordList vars(schemesDict.toc());
        schemes_.setSize(vars.size());
        forAll(vars, i)
        {
            // read solution variable name
            schemes_[i].first() = vars[i];

            // set semi-implicit (1) explicit (0) flag
            Istream& is = schemesDict.lookup(vars[i]);
            const word scheme(is);
            if (scheme == "semiImplicit")
            {
                schemes_[i].second().first() = true;
            }
            else if (scheme == "explicit")
            {
                schemes_[i].second().first() = false;
            }
            else
            {
                FatalErrorInFunction
                    << "Invalid scheme " << scheme << ". Valid schemes are "
                    << "explicit and semiImplicit" << exit(FatalError);
            }

            // read under-relaxation factor
            is  >> schemes_[i].second().second();
        }
    }
}


Foam::scalar Foam::bubbleData::relaxCoeff(const word& fieldName) const
{
    forAll(schemes_, i)
    {
        if (fieldName == schemes_[i].first())
        {
            return schemes_[i].second().second();
        }
    }

    FatalErrorInFunction
        << "Field name " << fieldName << " not found in schemes"
        << abort(FatalError);

    return 1;
}


bool Foam::bubbleData::semiImplicit(const word& fieldName) const
{
    forAll(schemes_, i)
    {
        if (fieldName == schemes_[i].first())
        {
            return schemes_[i].second().first();
        }
    }

    FatalErrorInFunction
        << "Field name " << fieldName << " not found in schemes"
        << abort(FatalError);

    return false;
}


bool Foam::bubbleData::solveThisStep() const
{
    return (mesh_.time().timeIndex() % calcFrequency_ == 0);
}


bool Foam::bubbleData::canEvolve()
{
    if (transient_)
    {
        trackTime_ = mesh_.time().deltaTValue();
    }
    else
    {
        trackTime_ = maxTrackTime_;
    }

    return solveThisStep();
}


bool Foam::bubbleData::output() const
{
    return mesh_.time().writeTime();
}
*/


// ************************************************************************* //
