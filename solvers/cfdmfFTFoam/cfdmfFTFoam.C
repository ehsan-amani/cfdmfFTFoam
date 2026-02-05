/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2021 OpenFOAM Foundation
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

Application
    cfdmfFTFoam

Description
    Solver for N incompressible, isothermal immiscible fluids using a FTM
    (Front Tracking Method),
    with optional mesh motion and mesh topology changes including adaptive
    re-meshing. (baseCode: interFoam)

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "dynamicFvMesh.H"
#include "CMULES.H"
#include "EulerDdtScheme.H"
#include "localEulerDdtScheme.H"
#include "CrankNicolsonDdtScheme.H"
#include "subCycle.H"
//#include "immiscibleIncompressibleTwoPhaseMixture.H" //EA
#include "noPhaseChange.H"
#include "kinematicMomentumTransportModel.H"
#include "pimpleControl.H"
#include "pressureReference.H"
#include "fvModels.H"
#include "fvConstraints.H"
#include "CorrectPhi.H"
#include "fvcSmooth.H"
#include "parcelCloudList.H" //EA
#include "frontTrackingCloud.H" //EA
#include "mathematicalConstants.H" //EA e,pi,twoPi,piByTwo
#include "conservativeSmooth.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "postProcess.H"

    #include "setRootCaseLists.H"
    #include "createTime.H"
    #include "createDynamicFvMesh.H"
    #include "initContinuityErrs.H"
    #include "createDyMControls.H"
    #include "createFields.H"
    #include "createFieldRefs.H"
    //#include "createAlphaFluxes.H" //EA
    #include "initCorrectPhi.H"
    #include "createUfIfPresent.H"

    //turbulence->validate(); //EA temp

    if (!LTS)
    {
        #include "CourantNo.H"
        #include "setInitialDeltaT.H"
    }

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
    Info<< "\nStarting time loop\n" << endl;

    while (pimple.run(runTime))
    {
        #include "readDyMControls.H"

        if (LTS)
        {
            //#include "setRDeltaT.H" //EA temp needs alpha1
        }
        else
        {
            #include "CourantNo.H"
            //#include "alphaCourantNo.H" //EA
            #include "setDeltaTtemp.H" //EA temp #include "setDeltaT.H"
        }

        runTime++;

        Info<< "Time = " << runTime.timeName() << nl << endl;

/*
    	//EA //FTM update 
		Info<< "Evolving " << frontTrackingCloud.name() << endl; 
		frontTrackingCloud.evolve();  

		frontTrackingCloud.updateDensityFromFT(rho); ///rho = frontTrackingCloud.densityIndicator(); //ED d2
		frontTrackingCloud.updateViscosityFromFT(mu); //mu = frontTrackingCloud.viscosityIndicator(); //ED d2
		sTension.primitiveFieldRef() = frontTrackingCloud.sTensionForceFromFT() - frontTrackingCloud.pressureJumpAtTheInterfaceFromFT();
		//pressureJumpAtTheInterface.primitiveFieldRef() = frontTrackingCloud.pressureJumpAtTheInterfaceFromFT(); 
		if (nCFilterLoops > 0) sTension = conservativeSmooth(sTension, mesh, nCFilterLoops); 
*/
        // --- Pressure-velocity PIMPLE corrector loop
        while (pimple.loop())
        {
            if (pimple.firstPimpleIter() || moveMeshOuterCorrectors)
            {
                // Store divU from the previous mesh so that it can be mapped
                // and used in correctPhi to ensure the corrected phi has the
                // same divergence
                tmp<volScalarField> divU;
			    
                //EA
				/*
                if
                (
                    correctPhi
                 && !isType<twoPhaseChangeModels::noPhaseChange>(phaseChange)
                )
                {
                    // Construct and register divU for mapping
                    divU = new volScalarField
                    (
                        "divU0",
                        fvc::div(fvc::absolute(phi, U))
                    );
                }
				*/

                fvModels.preUpdateMesh();

				frontTrackingCloud.storeGlobalPositions(); //EA

                mesh.update();

                if (mesh.changing())
                {
                    // Do not apply previous time-step mesh compression flux
                    // if the mesh topology changed
					//EA
					/*
                    if (mesh.topoChanging())
                    {
                        talphaPhi1Corr0.clear();
                    }
					*/

                    gh = (g & mesh.C()) - ghRef;
                    ghf = (g & mesh.Cf()) - ghRef;

                    MRF.update();

                    if (correctPhi)
                    {
                        #include "correctPhi.H"
                    }

                    //mixture.correct(); //EA

                    if (checkMeshCourantNo)
                    {
                        #include "meshCourantNo.H"
                    }
                }

                divU.clear();
            }

	    	//EA //FTM update 
            if (pimple.firstPimpleIter())
            {
				Info<< "Evolving " << frontTrackingCloud.name() << endl; 
				frontTrackingCloud.evolve();  

				frontTrackingCloud.updateDensityFromFT(rho); ///rho = frontTrackingCloud.densityIndicator(); //ED d2
				frontTrackingCloud.updateViscosityFromFT(mu); //mu = frontTrackingCloud.viscosityIndicator(); //ED d2
				sTension.primitiveFieldRef() = frontTrackingCloud.sTensionForceFromFT() - frontTrackingCloud.pressureJumpAtTheInterfaceFromFT();
				//pressureJumpAtTheInterface.primitiveFieldRef() = frontTrackingCloud.pressureJumpAtTheInterfaceFromFT(); 
				if (nCFilterLoops > 0) sTension = conservativeSmooth(sTension, mesh, nCFilterLoops); 
			}
          
            fvModels.correct();
/* //EA
            surfaceScalarField rhoPhi
            (
                IOobject
                (
                    "rhoPhi",
                    runTime.timeName(),
                    mesh
                ),
                mesh,
                dimensionedScalar(dimMass/dimTime, 0)
            );
*/
            //#include "alphaControls.H" //EA
            //#include "alphaEqnSubCycle.H" //EA
	    	rhoPhi == fvc::interpolate(rho)*phi; //EA 

            //mixture.correct(); //EA

            #include "UEqn.H"

            // --- Pressure corrector loop
            while (pimple.correct())
            {
                #include "pEqn.H"
            }

			//EA temp
			/*
            if (pimple.turbCorr())
            {
                //turbulence->correct(); 
            }
			*/
        }

        runTime.write();

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
