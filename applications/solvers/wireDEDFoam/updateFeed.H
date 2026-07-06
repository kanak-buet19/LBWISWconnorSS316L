// Wire-DED (Strategy A) — update the volumetric feed source Sfeed [1/s].
//
// Sfeed is the VOF mass source (gas -> metal conversion) for the alpha
// equation. The feed region is a sphere of radius feedRadius centred on the
// laser scan point plus feedCentreOffset. The region volume Vf is measured
// each step and weighted by the available gas fraction so that, exactly,
//
//     integral( Sfeed*(1 - alpha1) dV ) = mdot / rho1
//                                                 [m^3/s of metal added]
//
// which makes the injected mass rate equal to the prescribed wire mass rate
// regardless of mesh resolution, AMR state, or overlap with existing metal
// (before any MULES clipping).
{
    Sfeed == dimensionedScalar("Sfeed", Sfeed.dimensions(), 0.0);

    const scalar tNow = runTime.value();

    // Floor for the gas-volume normalization, as a fraction of the feed-region
    // volume Vf. Caps the source coefficient so a nearly-full feed region
    // (Vgas -> 0) cannot spike: below this floor the injected rate degrades
    // below target instead of diverging.
    const scalar gasFloorFraction =
        DEDProperties.lookupOrDefault<scalar>("feedGasFloorFraction", 0.1);

    if (tNow >= feedStartTime && tNow <= feedEndTime)
    {
        scalarField& SfeedI = Sfeed.primitiveFieldRef();
        const volVectorField& C = mesh.C();

        vector feedCentre = feedPosTable(tNow) + feedCentreOffset;
        const scalar r2 = sqr(feedRadius);

        if (DEDProperties.lookupOrDefault<Switch>("feedTracksSurface", false))
        {
            const vector nominalFeedCentre = feedCentre;
            const scalar nominalSurfaceY =
                nominalFeedCentre.y() - feedCentreOffset.y();
            const scalar maxSurfaceRise =
                DEDProperties.lookupOrDefault<scalar>
                (
                    "feedMaxSurfaceRise",
                    GREAT
                );
            const scalar minTrackedSurfaceY = nominalSurfaceY - maxSurfaceRise;
            const scalar surfaceSearchRadius =
                DEDProperties.lookupOrDefault<scalar>
                (
                    "feedSurfaceSearchRadius",
                    feedRadius
                );
            const scalar rSearch2 = sqr(surfaceSearchRadius);

            scalar yTopMetal = GREAT;
            forAll(SfeedI, celli)
            {
                const vector d = C[celli] - nominalFeedCentre;
                if
                (
                    alpha1[celli] > 0.5
                 && C[celli].y() >= minTrackedSurfaceY
                 && (sqr(d.x()) + sqr(d.z())) < rSearch2
                )
                {
                    yTopMetal = min(yTopMetal, C[celli].y());
                }
            }
            reduce(yTopMetal, minOp<scalar>());

            if (yTopMetal < GREAT/2)
            {
                feedCentre.y() = yTopMetal + feedCentreOffset.y();
            }

            // Do not let a tracked surface place the source sphere outside the
            // gas domain.
            feedCentre.y() = max(feedCentre.y(), mesh.bounds().min().y() + feedRadius);
        }

        scalar Vf = 0.0;
        scalar Vgas = 0.0;
        forAll(SfeedI, celli)
        {
            if (magSqr(C[celli] - feedCentre) < r2)
            {
                SfeedI[celli] = 1.0;
                Vf += mesh.V()[celli];
                Vgas += max(scalar(1) - alpha1[celli], scalar(0))*mesh.V()[celli];
            }
        }
        reduce(Vf, sumOp<scalar>());
        reduce(Vgas, sumOp<scalar>());

        // Normalize the source to hit the target metal rate mdot/rho1 while
        // gas is plentiful, but floor the normalization volume so coeff stays
        // bounded as the feed region fills with metal (Vgas -> 0). Without the
        // floor, coeff = mdot/(rho1*Vgas) diverges and crams the whole feed
        // rate into a vanishing gas pocket, spiking the alpha/momentum/energy
        // sources and collapsing deltaT.
        const scalar Vnorm = max(Vgas, gasFloorFraction*Vf);
        if (Vnorm > VSMALL)
        {
            const scalar coeff = feedMassRate.value()/(rho1.value()*Vnorm);
            SfeedI *= coeff;
        }
        else
        {
            Info<< "Wire-DED WARNING: feed region empty (Vf = 0),"
                << " no mass injected this step" << endl;
        }

        Sfeed.correctBoundaryConditions();

        // Report both the raw source strength and the bounded alpha source.
        // With gas-volume normalization, Sfeed*(1 - alpha1) is the quantity
        // that should match the target metal volume rate.
        const scalar SfeedInt = gSum(Sfeed.primitiveField()*mesh.V());
        scalar SfeedEffectiveInt = 0.0;
        forAll(SfeedI, celli)
        {
            SfeedEffectiveInt +=
                SfeedI[celli]
               *max(scalar(1) - alpha1[celli], scalar(0))
               *mesh.V()[celli];
        }
        reduce(SfeedEffectiveInt, sumOp<scalar>());

        Info<< "Wire-DED: feed centre = " << feedCentre
            << ", feed volume = " << Vf << " m^3"
            << ", available gas volume = " << Vgas << " m^3"
            << ", integral(Sfeed dV) = " << SfeedInt
            << " m^3/s (target mdot/rho1 = "
            << feedMassRate.value()/rho1.value() << ")"
            << ", effective integral(Sfeed*(1-alpha) dV) = "
            << SfeedEffectiveInt
            << " m^3/s, effective/target = "
            << SfeedEffectiveInt/(feedMassRate.value()/rho1.value())
            << endl;
    }
}
