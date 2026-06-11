{
    if (correctPhi)
    {
        rAU.ref() = 1.0/UEqn.A();
    }
    else
    {
        rAU = 1.0/UEqn.A();
    }

    surfaceScalarField rAUf("rAUf", fvc::interpolate(rAU()));
    volVectorField HbyA(constrainHbyA(rAU()*UEqn.H(), U, p_rgh));
    surfaceScalarField phiHbyA
    (
        "phiHbyA",
        fvc::flux(HbyA)
      + MRF.zeroFilter(fvc::interpolate(rho*rAU())*fvc::ddtCorr(U, phi, Uf))
    );
    MRF.makeRelative(phiHbyA);

    if (p_rgh.needReference())
    {
        fvc::makeRelative(phiHbyA, U);
        adjustPhi(phiHbyA, U, p_rgh);
        fvc::makeAbsolute(phiHbyA, U);
    }

    surfaceScalarField phig
    (
        (
            mixture.surfaceTensionForce()*fvc::interpolate(damper)
          +(fvc::interpolate(pVap)
           // Use smoothed interface for recoil force to avoid
           // one-cell pressure spikes on refined keyhole walls.
           *fvc::snGrad(alpha_smoothed)
           *fvc::interpolate(damper))
          - ghf*fvc::snGrad(rho*rhok)
        )*rAUf*mesh.magSf()
    );

    phiHbyA += phig;

    const Switch debugNumerics
    (
        runTime.controlDict().lookupOrDefault<Switch>
        (
            "debugNumerics",
            false
        )
    );

    const label debugEvery =
        runTime.controlDict().lookupOrDefault<label>
        (
            "debugNumericsInterval",
            25
        );

    const scalar debugSmallDeltaT =
        runTime.controlDict().lookupOrDefault<scalar>
        (
            "debugNumericsSmallDeltaT",
            5e-8
        );

    const label debugSmallDeltaTEvery =
        runTime.controlDict().lookupOrDefault<label>
        (
            "debugNumericsSmallDeltaTInterval",
            5
        );

    const bool debugThisStep =
        debugNumerics
     && (
            (
                runTime.deltaTValue() <= debugSmallDeltaT
             && debugSmallDeltaTEvery > 0
             && runTime.timeIndex() % debugSmallDeltaTEvery == 0
            )
         || (debugEvery > 0 && runTime.timeIndex() % debugEvery == 0)
        );

    if (debugThisStep)
    {
        scalarField magPhigI(mag(phig.primitiveField()));
        scalarField magPhiHbyAI(mag(phiHbyA.primitiveField()));
        scalarField magPhiI(mag(phi.primitiveField()));

        Info<< "[DBG pEqnFlux] time=" << runTime.timeName()
            << " proc=" << Pstream::myProcNo()
            << " maxMagPhig=" << gMax(magPhigI)
            << " meanMagPhig=" << gAverage(magPhigI)
            << " maxMagPhiHbyA=" << gMax(magPhiHbyAI)
            << " maxMagPhiOld=" << gMax(magPhiI)
            << " dt=" << runTime.deltaTValue()
            << endl;
    }

    // Update the pressure BCs to ensure flux consistency
    constrainPressure(p_rgh, U, phiHbyA, rAUf, MRF);

    while (pimple.correctNonOrthogonal())
    {
        fvScalarMatrix p_rghEqn
        (
            fvm::laplacian(rAUf, p_rgh) == fvc::div(phiHbyA)
        );

        p_rghEqn.setReference(pRefCell, getRefCellValue(p_rgh, pRefCell));

        p_rghEqn.solve(p_rgh.select(pimple.finalInnerIter()));

        if (pimple.finalNonOrthogonalIter())
        {
            phi = phiHbyA - p_rghEqn.flux();

            p_rgh.relax();

            U = HbyA + rAU()*fvc::reconstruct((phig - p_rghEqn.flux())/rAUf);
            U.correctBoundaryConditions();
            fvOptions.correct(U);

            if (debugThisStep && U.size())
            {
                scalarField magUI(mag(U.primitiveField()));
                label maxUCelli = 0;

                forAll(magUI, celli)
                {
                    if (magUI[celli] > magUI[maxUCelli])
                    {
                        maxUCelli = celli;
                    }
                }

                const scalar globalMaxU = gMax(magUI);

                if (magUI[maxUCelli] >= 0.95*globalMaxU)
                {
                    scalar level = -1;
                    if (mesh.foundObject<volScalarField>("cellLevel"))
                    {
                        level =
                            mesh.lookupObject<volScalarField>
                            (
                                "cellLevel"
                            )[maxUCelli];
                    }

                    tmp<surfaceScalarField> tpFlux = p_rghEqn.flux();
                    scalarField magPFluxI(mag(tpFlux().primitiveField()));

                    Info<< "[DBG U] time=" << runTime.timeName()
                        << " proc=" << Pstream::myProcNo()
                        << " cell=" << maxUCelli
                        << " C=" << mesh.C()[maxUCelli]
                        << " level=" << level
                        << " maxMagU=" << globalMaxU
                        << " localMagU=" << magUI[maxUCelli]
                        << " U=" << U[maxUCelli]
                        << " alpha=" << alpha1[maxUCelli]
                        << " T=" << T[maxUCelli]
                        << " pVap=" << pVap[maxUCelli]
                        << " recoil=" << recoilPressure[maxUCelli]
                        << " p_rgh=" << p_rgh[maxUCelli]
                        << " maxMagPressureFlux=" << gMax(magPFluxI)
                        << endl;
                }
            }
        }
    }

    #include "continuityErrs.H"

    // Correct Uf if the mesh is moving
    fvc::correctUf(Uf, U, phi);

    // Make the fluxes relative to the mesh motion
    fvc::makeRelative(phi, U);

    p == p_rgh + rho*gh*rhok;

    if (p_rgh.needReference())
    {
        p += dimensionedScalar
        (
            "p",
            p.dimensions(),
            pRefValue - getRefCellValue(p, pRefCell)
        );
        p_rgh = p - rho*gh*rhok;
    }

    pAbs = p + p0;
    pAbs.correctBoundaryConditions();

    if (!correctPhi)
    {
        rAU.clear();
    }

    Num_divU = fvc::div(phi);

}
