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
#include "IOstreams.H"
#include "IOField.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

template<class ParcelType>
Foam::string Foam::FrontTrackingParcel<ParcelType>::propertyList_ =
    Foam::FrontTrackingParcel<ParcelType>::propertyList();

template<class ParcelType>
const std::size_t Foam::FrontTrackingParcel<ParcelType>::sizeofFields_
(
    sizeof(FrontTrackingParcel<ParcelType>) - sizeof(ParcelType)
);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class ParcelType>
Foam::FrontTrackingParcel<ParcelType>::FrontTrackingParcel
(
    const polyMesh& mesh,
    Istream& is,
    bool readFields
)
:
    ParcelType(mesh, is, readFields),
    //UCorrect_(Zero) //EA rem
	currentIndex_(0), //EA add
    bubbleIndex_(0) //EA add
{
    if (readFields)
    {
        if (is.format() == IOstream::ASCII)
        {
            //is >> UCorrect_; //EA rem
			currentIndex_ = readLabel(is); //EA add
            bubbleIndex_ = readLabel(is); //EA add
        }
        else
        {
            //is.read(reinterpret_cast<char*>(&UCorrect_), sizeofFields_); //EA rem
 			is.read(reinterpret_cast<char*>(&currentIndex_), sizeofFields_); //EA add
        }
    }

    is.check
    (
        "FrontTrackingParcel<ParcelType>::Collisions"
        "(const polyMesh&, Istream&, bool)"
    );
}


template<class ParcelType>
template<class CloudType>
void Foam::FrontTrackingParcel<ParcelType>::readFields(CloudType& c)
{
    bool valid = c.size();

    ParcelType::readFields(c);

/* //EA rem
    IOField<vector> UCorrect
    (
        c.fieldIOobject("UCorrect", IOobject::MUST_READ),
        valid
    );
    c.checkFieldIOobject(c, UCorrect);
*/
	//EA add
    IOField<label> currentIndex
    (
        c.fieldIOobject("currentIndex", IOobject::MUST_READ),
        valid
    );
    c.checkFieldIOobject(c, currentIndex);

	//EA add
	IOField<label> bubbleIndex
    (
        c.fieldIOobject("bubbleIndex", IOobject::MUST_READ),
        valid
    );
    c.checkFieldIOobject(c, bubbleIndex);

    label i = 0;

    forAllIter(typename CloudType, c, iter)
    {
        FrontTrackingParcel<ParcelType>& p = iter();

        //p.UCorrect_ = UCorrect[i]; //EA rem
		p.currentIndex_ = currentIndex[i]; //EA add
		p.bubbleIndex_ = bubbleIndex[i]; //EA add

        i++;
    }
}


template<class ParcelType>
template<class CloudType>
void Foam::FrontTrackingParcel<ParcelType>::writeFields(const CloudType& c)
{
    ParcelType::writeFields(c);

    label np = c.size();

    //EA rem
    //IOField<vector>
    //    UCorrect(c.fieldIOobject("UCorrect", IOobject::NO_READ), np);
	//EA add
	IOField<label>
        currentIndex(c.fieldIOobject("currentIndex", IOobject::NO_READ), np);
	IOField<label>
        bubbleIndex(c.fieldIOobject("bubbleIndex", IOobject::NO_READ), np);
	

    label i = 0;

    forAllConstIter(typename CloudType, c, iter)
    {
        const FrontTrackingParcel<ParcelType>& p = iter();

        //UCorrect[i] = p.UCorrect(); //EA rem
		//EA add
		currentIndex[i] = p.currentIndex();
		bubbleIndex[i] = p.bubbleIndex();

        i++;
    }

    //UCorrect.write(np > 0); //EA rem
	//EA add
	currentIndex.write(np > 0);
	bubbleIndex.write(np > 0);
}


// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

template<class ParcelType>
Foam::Ostream& Foam::operator<<
(
    Ostream& os,
    const FrontTrackingParcel<ParcelType>& p
)
{
    if (os.format() == IOstream::ASCII)
    {
        os  << static_cast<const ParcelType&>(p)
            //<< token::SPACE << p.UCorrect(); //EA rem
			<< token::SPACE << p.currentIndex() //EA add
			<< token::SPACE << p.bubbleIndex(); //EA add
    }
    else
    {
        os  << static_cast<const ParcelType&>(p);
        os.write
        (
            //reinterpret_cast<const char*>(&p.UCorrect_), //EA rem
			reinterpret_cast<const char*>(&p.currentIndex_), //EA add
            FrontTrackingParcel<ParcelType>::sizeofFields_
        );
    }

    os.check
    (
        "Ostream& operator<<(Ostream&, const FrontTrackingParcel<ParcelType>&)"
    );

    return os;
}


// ************************************************************************* //
