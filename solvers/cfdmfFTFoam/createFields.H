#include "createRDeltaT.H"

Info<< "Reading field p_rgh\n" << endl;
volScalarField p_rgh
(
    IOobject
    (
        "p_rgh",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info<< "Reading field U\n" << endl;
volVectorField U
(
    IOobject
    (
        "U",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

#include "createPhi.H"

//EA
/*
Info<< "Reading transportProperties\n" << endl;
immiscibleIncompressibleTwoPhaseMixture mixture(U, phi);

autoPtr<twoPhaseChangeModel> phaseChangePtr
(
    twoPhaseChangeModel::New(mixture)
);

volScalarField& alpha1(mixture.alpha1());
volScalarField& alpha2(mixture.alpha2());

const dimensionedScalar& rho1 = mixture.rho1();
const dimensionedScalar& rho2 = mixture.rho2();


// Need to store rho for ddt(rho, U)
volScalarField rho
(
    IOobject
    (
        "rho",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT
    ),
    alpha1*rho1 + alpha2*rho2
);
rho.oldTime();
*/
volScalarField rho
(
	IOobject
	(
		"rho",
		runTime.timeName(),
		mesh,
		IOobject::MUST_READ,
		IOobject::AUTO_WRITE
	),
	mesh
);
rho.oldTime();

//EA temp
/*
// Construct incompressible turbulence model
autoPtr<incompressible::momentumTransportModel> turbulence
(
    incompressible::momentumTransportModel::New(U, phi, mixture)
);
*/

#include "readGravitationalAcceleration.H"
#include "readhRef.H"
#include "gh.H"


volScalarField p
(
    IOobject
    (
        "p",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    p_rgh + rho*gh
);

pressureReference pressureReference(p, p_rgh, pimple.dict());

if (p_rgh.needReference())
{
    p += dimensionedScalar
    (
        "p",
        p.dimensions(),
        pressureReference.refValue()
      - getRefCellValue(p, pressureReference.refCell())
    );
    p_rgh = p - rho*gh;
}

mesh.setFluxRequired(p_rgh.name());
//mesh.setFluxRequired(alpha1.name()); //EA

#include "createMRF.H"
#include "createFvModels.H"
#include "createFvConstraints.H"

//EA
volScalarField mu
(
    IOobject
    (
        "mu",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

label nCFilterLoops = runTime.controlDict().lookupOrDefault<label>("nCFilterLoops", 0);
Info<< "The number of filtering of surface tension field is " << nCFilterLoops << "\n" << endl;

// EA temp
//Info<< "Constructing clouds" << endl;
//parcelCloudList clouds(rho, U, mu, g);

Info<< nl << "Constructing frontTrackingCloud " << endl;
frontTrackingCloud frontTrackingCloud
//("frontTrackingCloud", mesh);
(
    "frontTrackingCloud",
    rho,
    U,
    mu,
    g
);

volVectorField sTension
(
    IOobject
    (
        "sTension",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedVector("0", dimForce/dimVolume, vector::zero)
);


frontTrackingCloud.updateDensityFromFT(rho); ///rho = frontTrackingCloud.densityIndicator(); //ED d2
frontTrackingCloud.updateViscosityFromFT(mu); //mu = frontTrackingCloud.viscosityIndicator(); //ED d2
sTension.primitiveFieldRef() = frontTrackingCloud.sTensionForceFromFT() - frontTrackingCloud.pressureJumpAtTheInterfaceFromFT();
//volVectorField gradP("gradP", frontTrackingCloud.dPressure());
//volVectorField pressureJumpAtTheInterface("pressureJumpAtTheInterface", frontTrackingCloud.pressureJumpAtTheInterfaceFromFT());

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
