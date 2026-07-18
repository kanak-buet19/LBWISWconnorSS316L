{
    TSolidus = alpha1*Tsolidus1 + (1.0 - alpha1)*Tsolidus2;
    TLiquidus = alpha1*Tliquidus1 + (1.0 - alpha1)*Tliquidus2;

    if (mthd.valid())
    {
        mthd->correctIncompressibleProperties
        (
            alpha1,
            transportPropertiesMetal,
            transportPropertiesGas
        );
    }

    beta = alpha1*beta1 + (1.0 - alpha1)*beta2;
    rhok = (1.0 - max(epsilon1*(beta)*(T - TSolidus), 0.0));

    scalarField& alpha1I = alpha1;
    scalarField& TSolidusI = TSolidus;
    scalarField& TLiquidusI = TLiquidus;
    const scalar Tsolidus1Value = Tsolidus1.value();
    const scalar Tliquidus1Value = Tliquidus1.value();
    forAll(alpha1I, celli)
    {
        if (alpha1I[celli] > 0.05)
        {
            TSolidusI[celli] = Tsolidus1Value;
            TLiquidusI[celli] = Tliquidus1Value;
        }
    }

    TSolidus.correctBoundaryConditions();
    TLiquidus.correctBoundaryConditions();
    beta.correctBoundaryConditions();  // ADDED
    rhok.correctBoundaryConditions();  // ADDED

    #include "updateKappaCp.H"

    const vectorField& CI = mesh.C();
    metalaverage = fvc::average(alpha1);
    metalaverage.correctBoundaryConditions();  // ADDED

    alpha_filtered = volScalarField("alpha_filtered", alpha1);
    alpha_filtered.correctBoundaryConditions();
    scalarField& alpha_filteredI = alpha_filtered;
    // alpha_filteredI.correctBoundaryConditions();
    forAll(CI, celli)
    {
        if (alpha1I[celli] > 0.99)
        {
            alpha_filteredI[celli] = 1.0;
        }
        else
        {
            if (alpha1I[celli] < 0.01)
            {
                alpha_filteredI[celli] = 0.0;
            }
            else
            {
                alpha_filteredI[celli] = alpha1I[celli];
            }
        }
    }

    // Access the boundary fields
    volScalarField::Boundary& alpha_filteredBf = alpha_filtered.boundaryFieldRef();
    const volScalarField::Boundary& alpha1Bf = alpha1.boundaryField();

    // Loop through all boundary patches
    forAll(alpha_filteredBf, patchi)
    {
        // Get references to the patch fields
        scalarField& alpha_filteredPatch = alpha_filteredBf[patchi];
        const scalarField& alpha1Patch = alpha1Bf[patchi];

        // Apply the same logic to each face in the patch
        forAll(alpha_filteredPatch, facei)
        {
            if (alpha1Patch[facei] > 0.99)
            {
                alpha_filteredPatch[facei] = 1.0;
            }
            else if (alpha1Patch[facei] < 0.01)
            {
                alpha_filteredPatch[facei] = 0.0;
            }
            else
            {
                alpha_filteredPatch[facei] = alpha1Patch[facei];
            }
        }
    }

    alpha_filtered.correctBoundaryConditions();

    alpha_smoothed = fvc::average(alpha1);

    #include "updateInterfaceTemperature.H"

    dSigmaDtf = Marangoni_Constant;

    if (useFeSMarangoni)
    {
        scalarField& sigmaI = sigmaField.primitiveFieldRef();
        scalarField& dSigmaDtI = dSigmaDtf.primitiveFieldRef();
        const scalarField& TI = T;

        forAll(sigmaI, celli)
        {
            const scalar Tc = min(max(TI[celli], Tliquidus1.value()), Tvap.value());
            const scalar B =
                FeSAdsorptionEntropy*sulfurActivity
               *Foam::exp(FeSAdsorptionHeat/(FeSGasConstant*Tc));

            if (useFeSSigma)
            {
                sigmaI[celli] =
                    FeSPureSigma
                  - FeSPureSlope*(Tc - FeSTm)
                  - FeSGasConstant*Tc*FeSGammaS*Foam::log(1.0 + B);
            }

            dSigmaDtI[celli] =
              - FeSPureSlope
              - FeSGasConstant*FeSGammaS*Foam::log(1.0 + B)
              + FeSGammaS*FeSAdsorptionHeat/Tc*(B/(1.0 + B));
        }

        forAll(sigmaField.boundaryFieldRef(), patchi)
        {
            fvPatchScalarField& sigmaP = sigmaField.boundaryFieldRef()[patchi];
            fvPatchScalarField& dSigmaDtP = dSigmaDtf.boundaryFieldRef()[patchi];
            const scalarField& TP = T.boundaryField()[patchi];

            forAll(sigmaP, facei)
            {
                const scalar Tc = min(max(TP[facei], Tliquidus1.value()), Tvap.value());
                const scalar B =
                    FeSAdsorptionEntropy*sulfurActivity
                   *Foam::exp(FeSAdsorptionHeat/(FeSGasConstant*Tc));

                if (useFeSSigma)
                {
                    sigmaP[facei] =
                        FeSPureSigma
                      - FeSPureSlope*(Tc - FeSTm)
                      - FeSGasConstant*Tc*FeSGammaS*Foam::log(1.0 + B);
                }

                dSigmaDtP[facei] =
                  - FeSPureSlope
                  - FeSGasConstant*FeSGammaS*Foam::log(1.0 + B)
                  + FeSGammaS*FeSAdsorptionHeat/Tc*(B/(1.0 + B));
            }
        }

        if (useFeSSigma)
        {
            sigmaField.correctBoundaryConditions();
        }
        dSigmaDtf.correctBoundaryConditions();
    }
    else if (useLinearSigma)
    {
        // Linear σ(T) clamped: constant sigma_ref for T < Tliquidus,
        // linear for Tliquidus <= T <= Tvap, constant sigma_at_Tvap for T > Tvap
        sigmaField = max
        (
            min
            (
                sigma_ref + Marangoni_Constant*(T - sigma_T_ref),
                sigma_ref
            ),
            sigma_at_Tvap
        );
        sigmaField.correctBoundaryConditions();
    }

    // Read dynamicMeshDict to select the active refinement criterion
    IOdictionary dynamicMeshDict
    (
        IOobject
        (
            "dynamicMeshDict",
            runTime.constant(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );

    const word amrCriterion
    (
        dynamicMeshDict.lookupOrDefault<word>("amrCriterion", "interfaceT300")
    );

    if (amrCriterion == "interfaceGradT")
    {
        const dimensionedScalar gradTrefine
        (
            "gradTrefine",
            dimensionSet(0, -1, 0, 1, 0, 0, 0),
            dynamicMeshDict.lookupOrDefault<scalar>("gradTrefine", 5e6)
        );

        volScalarField magGradT(mag(fvc::grad(T)));
        meshRefinementIndicator =
            pos(alpha_smoothed - dimensionedScalar("alphaMin", dimless, 0.001))
           *pos(dimensionedScalar("alphaMax", dimless, 0.999) - alpha_smoothed)
           *pos(magGradT - gradTrefine);
    }
    else if (amrCriterion == "interfaceOnly")
    {
        // Refine entire alpha interface (0 < alpha < 1) plus 1 adjacent cell layer on each side
        volScalarField baseIndicator
        (
            pos(alpha_smoothed - dimensionedScalar("alphaMin", dimless, 0.001))
           *pos(dimensionedScalar("alphaMax", dimless, 0.999) - alpha_smoothed)
        );

        // Dilate by 1 layer: face-average spreads indicator to immediate neighbors
        volScalarField neighborIndicator
        (
            pos
            (
                fvc::average(fvc::interpolate(baseIndicator))
              - dimensionedScalar("dilateThresh", dimless, SMALL)
            )
        );

        meshRefinementIndicator = max(baseIndicator, neighborIndicator);
    }
    else if (amrCriterion == "interfacePersist")
    {
        // Refine current interface plus the hot metal HAZ, and retain all previously
        // refined cells. Once a cell is marked, it stays marked permanently.
        // Needed for keyhole porosity and post-laser heat-affected regions.
        volScalarField baseIndicator
        (
            pos(alpha_smoothed - dimensionedScalar("alphaMin", dimless, 0.001))
           *pos(dimensionedScalar("alphaMax", dimless, 0.999) - alpha_smoothed)
        );

        volScalarField neighborIndicator
        (
            pos
            (
                fvc::average(fvc::interpolate(baseIndicator))
              - dimensionedScalar("dilateThresh", dimless, SMALL)
            )
        );

        volScalarField hazIndicator
        (
            pos(alpha_smoothed - dimensionedScalar("alphaMin", dimless, 0.001))
           *pos(T - (Tsolidus1 - dimensionedScalar("dTrefine", dimTemperature, 300)))
        );

        // max() with existing field: accumulated markers never reset to 0
        meshRefinementIndicator =
            max(meshRefinementIndicator, max(max(baseIndicator, neighborIndicator), hazIndicator));
    }
    else if (amrCriterion == "interfaceT100")
    {
        meshRefinementIndicator =
            pos(alpha_smoothed - dimensionedScalar("alphaMin", dimless, 0.001))
           *pos(dimensionedScalar("alphaMax", dimless, 0.999) - alpha_smoothed)
           *pos(T - (Tsolidus1 - dimensionedScalar("dTrefine", dimTemperature, 100)));
    }
    else
    {
        meshRefinementIndicator =
            pos(alpha_smoothed - dimensionedScalar("alphaMin", dimless, 0.001))
           *pos(dimensionedScalar("alphaMax", dimless, 0.999) - alpha_smoothed)
           *pos(T - (Tsolidus1 - dimensionedScalar("dTrefine", dimTemperature, 300)));
    }
    meshRefinementIndicator.correctBoundaryConditions();

    // Keep amrIndicator matching the selected active criteria or default behavior
    amrIndicator = meshRefinementIndicator;
    amrIndicator.correctBoundaryConditions();

    const Switch debugNumerics
    (
        runTime.controlDict().lookupOrDefault<Switch>
        (
            "debugNumerics",
            false
        )
    );

    if (debugNumerics)
    {
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
            (
                runTime.deltaTValue() <= debugSmallDeltaT
             && debugSmallDeltaTEvery > 0
             && runTime.timeIndex() % debugSmallDeltaTEvery == 0
            )
         || (debugEvery > 0 && runTime.timeIndex() % debugEvery == 0);

        if (debugThisStep)
        {
            const scalar activeCells = gSum(amrIndicator.primitiveField());

            Info<< "[DBG AMRField] time=" << runTime.timeName()
                << " activeCells=" << activeCells
                << " totalCells=" << returnReduce(mesh.nCells(), sumOp<label>())
                << " amrMin=" << gMin(amrIndicator.primitiveField())
                << " amrMax=" << gMax(amrIndicator.primitiveField())
                << " amrMean=" << gAverage(amrIndicator.primitiveField())
                << " alphaSmoothMin=" << gMin(alpha_smoothed.primitiveField())
                << " alphaSmoothMax=" << gMax(alpha_smoothed.primitiveField())
                << " TMin=" << gMin(T.primitiveField())
                << " TMax=" << gMax(T.primitiveField())
                << endl;
        }
    }

    LatentHeat =
        alpha_filtered*LatentHeat1 + (1.0 - alpha_filtered)*LatentHeat2;

    // gradepsilon1 = fvc::grad(epsilon1);
    // gradepsilon1.correctBoundaryConditions();

    {
        const volVectorField gradAlphaFiltered(fvc::grad(alpha_filtered));
        n_filtered = gradAlphaFiltered/(mag(gradAlphaFiltered) + deltaN);
    }

    e1temp = fvc::average(epsilon1);

    const scalarField& e1tempI = e1temp;
    const scalarField& epsilon1I = epsilon1;
    scalarField& epsilon1maskI = epsilon1mask;
    forAll(CI, celli)
    {
        if (e1tempI[celli] <= 0.95)
        {
            epsilon1maskI[celli] = 0.0;
        }
        else
        {
            epsilon1maskI[celli] = epsilon1I[celli];
        }
    }

    epsilon1mask.correctBoundaryConditions();

    // e2temp = fvc::average(epsilon1mask);
    // epsilon1mask2 *= 0.0;

    // const scalarField& e2tempI = e2temp;
    // scalarField& epsilon1mask2I = epsilon1mask2;
    // forAll(CI, celli)
    // {
    //     if (e2tempI[celli] <= 0.95)
    //     {
    //         epsilon1mask2I[celli] = 0.0;
    //     }
    //     else
    //     {
    //         epsilon1mask2I[celli] = epsilon1maskI[celli];
    //     }

    // }

    // epsilon1mask2.correctBoundaryConditions();

    if (laser.powderSim())
    {
        gh = epsilon1mask*(g & mesh.C());
        ghf = fvc::interpolate(epsilon1mask)*(g & mesh.Cf());
    }
    else{
        gh = (g & mesh.C());
        ghf = (g & mesh.Cf());
    }

    // Dynamic temperature-dependent electrical resistivity lookup
    if (hasResistivityTable && resistivityProp_m.valid())
    {
        scalarField& electrical_resistivityI = electrical_resistivity;
        const scalarField& TI = T;
        const scalarField& alpha1I = alpha1;
        const scalar gasRes = 5e-7;

        forAll(electrical_resistivityI, celli)
        {
            const scalar metalRes = resistivityProp_m->value(TI[celli]);
            electrical_resistivityI[celli] = 
                max(min(alpha1I[celli], 1.0), 0.0)*metalRes
              + (1.0 - max(min(alpha1I[celli], 1.0), 0.0))*gasRes;
        }
        electrical_resistivity.correctBoundaryConditions();
    }
}
