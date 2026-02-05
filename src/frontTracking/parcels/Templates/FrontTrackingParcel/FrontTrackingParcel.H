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

Class
    Foam::FrontTrackingParcel

Description
    Wrapper around parcel types to add FrontTracking modelling

SourceFiles
    FrontTrackingParcelI.H
    FrontTrackingParcelTrackingDataI.H
    FrontTrackingParcel.C
    FrontTrackingParcelIO.C

\*---------------------------------------------------------------------------*/

#ifndef FrontTrackingParcel_H
#define FrontTrackingParcel_H

#include "particle.H"
#include "labelFieldIOField.H"
#include "vectorFieldIOField.H"
#include "bubbleData.H" //EA add3

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes

template<class ParcelType>
class FrontTrackingParcel;

//EA rem2
//template<class Type>
//class AveragingMethod;

// Forward declaration of friend functions

template<class ParcelType>
Ostream& operator<<
(
    Ostream&,
    const FrontTrackingParcel<ParcelType>&
);

/*---------------------------------------------------------------------------*\
                       Class FrontTrackingParcelName Declaration
\*---------------------------------------------------------------------------*/

TemplateName(FrontTrackingParcel);


/*---------------------------------------------------------------------------*\
                         Class FrontTrackingParcel Declaration
\*---------------------------------------------------------------------------*/

template<class ParcelType>
class FrontTrackingParcel
:
    public ParcelType,
    public FrontTrackingParcelName
{
    // Private Data

        //- Size in bytes of the fields
        static const std::size_t sizeofFields_;


public:

    class trackingData
    :
        public ParcelType::trackingData
    {

    public:

        enum trackPart //EA
        {
            tpTrackTime,
            tpTrackDist,
        };


    private:

        // Private Data

            // FrontTracking Averages

                //- Volume average
                //autoPtr<AveragingMethod<scalar>> volumeAverage_; //EA rem

                //- Radius average [ volume^(1/3) ]
                //autoPtr<AveragingMethod<scalar>> radiusAverage_; //EA rem

                //- Density average
                //autoPtr<AveragingMethod<scalar>> rhoAverage_; //EA rem

                //- Velocity average
                //autoPtr<AveragingMethod<vector>> uAverage_; //EA rem

                //- Magnitude velocity sqyuared average
                //autoPtr<AveragingMethod<scalar>> uSqrAverage_; //EA rem

                //- Frequency average
                //autoPtr<AveragingMethod<scalar>> frequencyAverage_; //EA rem

                //- Mass average
                //autoPtr<AveragingMethod<scalar>> massAverage_; //EA rem


            //- Which part of the integration algorithm is taking place
            trackPart part_;


    public:

        //- Constructors

            //- Construct from components
            template<class TrackCloudType>
            inline trackingData(const TrackCloudType& cloud);


        //- Update the FrontTracking averages
        //template<class TrackCloudType> //EA rem
        //inline void updateAverages(const TrackCloudType& cloud);


        //- Access

            //- Const access to the tracking part
            inline trackPart part() const;

            //- Non const access to the tracking part
            inline trackPart& part();
    };


protected:

    // Protected data

        //- Velocity correction due to collisions [m/s]
        //vector UCorrect_; //EA rem

        //EA add- parcel Index used to match the corresponding pointData
        label currentIndex_;

        //EA add- bubble Index used to connect to the corresponding bubble 
        label bubbleIndex_;


public:

    // Static Data Members

        //- String representation of properties
        AddToPropertyList
        (
            ParcelType,
            //"(UCorrectx UCorrecty UCorrectz)" //EA rem
			" currentIndex" //EA add
          + " bubbleIndex" //EA add
        );


    // Constructors

        //- Construct from mesh, coordinates and topology
        //  Other properties initialised as null
        inline FrontTrackingParcel
        (
            const polyMesh& mesh,
            const barycentric& coordinates,
            const label celli,
            const label tetFacei,
            const label tetPti
        );

        //- Construct from a position and a cell, searching for the rest of the
        //  required topology. Other properties are initialised as null.
        inline FrontTrackingParcel
        (
            const polyMesh& mesh,
            const vector& position,
            const label celli
        );

        //- Construct from Istream
        FrontTrackingParcel
        (
            const polyMesh& mesh,
            Istream& is,
            bool readFields = true
        );

        //- Construct as a copy
        FrontTrackingParcel(const FrontTrackingParcel& p);

        //- Construct as a copy
        FrontTrackingParcel(const FrontTrackingParcel& p, const polyMesh& mesh);

        //- Construct and return a (basic particle) clone
        virtual autoPtr<particle> clone() const
        {
            return autoPtr<particle>(new FrontTrackingParcel(*this));
        }

        //- Construct and return a (basic particle) clone
        virtual autoPtr<particle> clone(const polyMesh& mesh) const
        {
            return autoPtr<particle>(new FrontTrackingParcel(*this, mesh));
        }

        //- Factory class to read-construct particles used for
        //  parallel transfer
        class iNew
        {
            const polyMesh& mesh_;

        public:

            iNew(const polyMesh& mesh)
            :
                mesh_(mesh)
            {}

            autoPtr<FrontTrackingParcel<ParcelType>> operator()(Istream& is) const
            {
                return autoPtr<FrontTrackingParcel<ParcelType>>
                (
                    new FrontTrackingParcel<ParcelType>(mesh_, is, true)
                );
            }
        };


    // Member Functions

        // Access

            //- Return const access to correction velocity
            //inline const vector& UCorrect() const; //EA del

            //- Return access to correction velocity
            //inline vector& UCorrect(); //EA del

			//EA add- Return const access to currentIndex
            inline const label& currentIndex() const;

            //EA add- Return const access to bubbleIndex
            inline const label& bubbleIndex() const;

			//EA add- Return access to currentIndex
            inline label& currentIndex();

            //EA add- Return access to bubbleIndex
            inline label& bubbleIndex();


		// Main calculation loop //EA add

            //- Set cell values
            template<class TrackCloudType>
            void setCellValues(TrackCloudType& cloud, trackingData& td);


            //- Correct cell values using latest transfer information
            template<class TrackCloudType>
            void cellValueSourceCorrection
            (
                TrackCloudType& cloud,
                trackingData& td,
                const scalar dt
            );

			//- Update parcel properties over the time interval
            template<class TrackCloudType>
            void calc
            (
                TrackCloudType& cloud,
                trackingData& td,
                const scalar dt
            );


        // Tracking

			//- Move the parcel
            template<class TrackCloudType>
            bool move
            (
                TrackCloudType& cloud,
                trackingData& td,
                const scalar trackTime
            );


            //- Move the parcel
            template<class TrackCloudType>
            bool moveDistance
            (
                TrackCloudType& cloud,
                trackingData& td,
                const vector trackDist
            );


    // Friend Functions

        // I-O

            //- Read
            template<class CloudType>
            static void readFields(CloudType& c);

            //- Write
            template<class CloudType>
            static void writeFields(const CloudType& c);


        // Ostream operator

            friend Ostream& operator<< <ParcelType>
            (
                Ostream&,
                const FrontTrackingParcel<ParcelType>&
            );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "FrontTrackingParcelI.H"
#include "FrontTrackingParcelTrackingDataI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "FrontTrackingParcel.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
