    MRF.correctBoundaryVelocity(U);

    fvVectorMatrix UEqn
    (
        fvm::ddt(rho, U) + fvm::div(rhoPhi, U)
      + MRF.DDt(rho, U)
      //+ turbulence->divDevTau(rho, U) //EA temp
      - fvm::laplacian(mu, U) //EA temp
      - fvc::div(mu*dev2(T(fvc::grad(U)))) //EA temp
     ==
       //phaseChange.SU(rho, rhoPhi, U) //EA
     //+ fvModels.source(rho, U) //EA
	   fvModels.source(rho, U) //EA
    );

    UEqn.relax();

    fvConstraints.constrain(UEqn);

    if (pimple.momentumPredictor())
    {
        solve
        (
            UEqn
         ==
            fvc::reconstruct
            (
                (
                    //mixture.surfaceTensionForce() //EA rem
		    		(fvc::interpolate(sTension) & mesh.Sf()/mesh.magSf() ) //EA add
                  - ghf*fvc::snGrad(rho)
                  - fvc::snGrad(p_rgh)
                ) * mesh.magSf()
            )
        );

        fvConstraints.constrain(U);
    }
