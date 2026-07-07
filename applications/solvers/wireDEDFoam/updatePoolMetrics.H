// Multi-probe, quasi-steady-aware melt-pool / bead geometry measurement.
//
// Reproduces the paper's measured quantities natively at runtime (no VTK):
//   melt pool  -> half-width 0.5W, penetration depth D
//   clad/bead  -> half-width 0.5W, height h, penetration depth D
//
// y increases downward. substrateSurfaceY is the original (undeformed) top
// surface. metal = alpha1>0.5; liquid = metal & T>=Tliquidus; fused =
// the monotonic "was-ever-molten" flag (paper's indication parameter).
//
// N probes (thin y-z slabs) sit along the scan path. Pool metrics PEAK as
// the beam crosses each probe (captured by running max). Clad metrics come
// from the fused zone, which only grows, so their running max is the final
// value. Quasi-steady state = the interior probes agree (low coefficient of
// variation) = translational invariance in the scan direction. That
// agreement, not any single probe, is the paper-comparable measurement.
{
    static label sampleCounter = 0;

    const label poolSampleEvery =
        DEDProperties.lookupOrDefault<label>("poolSampleEvery", 50);

    const bool doSample =
        enablePoolMetrics
     && (sampleCounter % poolSampleEvery == 0 || runTime.outputTime());
    sampleCounter++;

    if (doSample)
    {
        const scalar Tliq = Tliquidus1.value();
        const scalar Tsol = Tsolidus1.value();
        const volVectorField& C = mesh.C();

        // 1. Monotonic fused flag: once a metal cell reaches solidus, stays 1.
        forAll(fusedFlag, celli)
        {
            if (alpha1[celli] > 0.5 && T[celli] >= Tsol)
            {
                fusedFlag[celli] = 1.0;
            }
        }
        fusedFlag.correctBoundaryConditions();

        // 2. Per-probe local extremes (sentinels: -GREAT for maxima,
        //    +GREAT for minima). Single cell loop, each cell falls in at
        //    most one probe band.
        scalarList maxYliq(nProbes, -GREAT);
        scalarList minZliq(nProbes,  GREAT);
        scalarList maxZliq(nProbes, -GREAT);
        scalarList maxYfused(nProbes, -GREAT);
        scalarList minYfused(nProbes,  GREAT);
        scalarList minZfused(nProbes,  GREAT);
        scalarList maxZfused(nProbes, -GREAT);
        scalarList minYmetal(nProbes,  GREAT);

        forAll(C, celli)
        {
            const scalar x = C[celli].x();

            // find owning probe (bands are non-overlapping)
            forAll(probeX, i)
            {
                if (mag(x - probeX[i]) <= probeBandHalf)
                {
                    const scalar y = C[celli].y();
                    const scalar z = C[celli].z();
                    const bool isMetal = (alpha1[celli] > 0.5);
                    const bool isLiquid = (isMetal && T[celli] >= Tliq);
                    const bool isFused = (fusedFlag[celli] > 0.5);

                    if (isLiquid)
                    {
                        maxYliq[i] = max(maxYliq[i], y);
                        minZliq[i] = min(minZliq[i], z);
                        maxZliq[i] = max(maxZliq[i], z);
                    }
                    if (isFused)
                    {
                        maxYfused[i] = max(maxYfused[i], y);
                        minYfused[i] = min(minYfused[i], y);
                        minZfused[i] = min(minZfused[i], z);
                        maxZfused[i] = max(maxZfused[i], z);
                    }
                    if (isMetal)
                    {
                        minYmetal[i] = min(minYmetal[i], y);
                    }
                    break; // cell binned; no overlap
                }
            }
        }

        // 3. Parallel reduction of each per-probe extreme.
        forAll(probeX, i)
        {
            reduce(maxYliq[i],   maxOp<scalar>());
            reduce(minZliq[i],   minOp<scalar>());
            reduce(maxZliq[i],   maxOp<scalar>());
            reduce(maxYfused[i], maxOp<scalar>());
            reduce(minYfused[i], minOp<scalar>());
            reduce(minZfused[i], minOp<scalar>());
            reduce(maxZfused[i], maxOp<scalar>());
            reduce(minYmetal[i], minOp<scalar>());
        }

        // 4. Convert to geometry and update persistent running maxima.
        forAll(probeX, i)
        {
            if (maxYliq[i] > -GREAT/2)
            {
                peakPoolDepth[i] =
                    max(peakPoolDepth[i], maxYliq[i] - substrateSurfaceY);
                if (maxZliq[i] > -GREAT/2 && minZliq[i] < GREAT/2)
                {
                    peakPoolHalfWidth[i] =
                        max(peakPoolHalfWidth[i],
                            0.5*(maxZliq[i] - minZliq[i]));
                }
            }
            if (maxYfused[i] > -GREAT/2)
            {
                cladDepth[i] =
                    max(cladDepth[i], maxYfused[i] - substrateSurfaceY);
            }
            if (minYfused[i] < GREAT/2)
            {
                cladHeight[i] =
                    max(cladHeight[i], substrateSurfaceY - minYfused[i]);
            }
            if (maxZfused[i] > -GREAT/2 && minZfused[i] < GREAT/2)
            {
                cladHalfWidth[i] =
                    max(cladHalfWidth[i], 0.5*(maxZfused[i] - minZfused[i]));
            }
            if (minYmetal[i] < GREAT/2)
            {
                beadMetalHeight[i] =
                    max(beadMetalHeight[i], substrateSurfaceY - minYmetal[i]);
            }
        }
    }

    // 5. CSV output at write times (master only). All lengths in mm.
    if (enablePoolMetrics && runTime.outputTime() && Pstream::master())
    {
        // runTime.globalPath() gets intercepted by the parallel fileHandler
        // and still resolves inside processorN/. Navigate up explicitly
        // instead so the CSVs land at the case root regardless of
        // decomposition.
        const fileName caseRoot =
            Pstream::parRun() ? (runTime.path()/"..") : runTime.path();

        mkDir(caseRoot/"post-processing-data");

        // 5a. Per-probe running values (full audit trail).
        {
            fileName f =
                caseRoot/"post-processing-data"/"pool_metrics_perProbe.csv";
            const bool header = !isFile(f);
            std::ofstream os(f.c_str(), std::ios_base::app);
            if (os.is_open())
            {
                if (header)
                {
                    os << "Time_s,probe,probeX_mm,"
                       << "poolDepth_mm,poolHalfWidth_mm,"
                       << "cladDepth_mm,cladHeight_mm,cladHalfWidth_mm,"
                       << "beadMetalHeight_mm\n";
                }
                forAll(probeX, i)
                {
                    os << runTime.value() << "," << i
                       << "," << probeX[i]*1e3
                       << "," << peakPoolDepth[i]*1e3
                       << "," << peakPoolHalfWidth[i]*1e3
                       << "," << cladDepth[i]*1e3
                       << "," << cladHeight[i]*1e3
                       << "," << cladHalfWidth[i]*1e3
                       << "," << beadMetalHeight[i]*1e3 << "\n";
                }
                os.close();
            }
        }

        // 5b. Steady summary: mean/std/CoV across interior probes.
        const label i0 = interiorProbeMargin;
        const label i1 = nProbes - 1 - interiorProbeMargin;
        const label nInt = i1 - i0 + 1;

        if (nInt >= 1)
        {
            // mean/std helper over interior probes for one metric list
            auto meanStd = [&](const scalarList& v, scalar& m, scalar& s)
            {
                m = 0.0;
                for (label i = i0; i <= i1; ++i) m += v[i];
                m /= scalar(nInt);
                s = 0.0;
                for (label i = i0; i <= i1; ++i) s += sqr(v[i] - m);
                s = (nInt > 1) ? Foam::sqrt(s/scalar(nInt)) : 0.0;
            };

            scalar mPD, sPD, mPW, sPW, mCD, sCD, mCH, sCH, mCW, sCW;
            meanStd(peakPoolDepth,     mPD, sPD);
            meanStd(peakPoolHalfWidth, mPW, sPW);
            meanStd(cladDepth,         mCD, sCD);
            meanStd(cladHeight,        mCH, sCH);
            meanStd(cladHalfWidth,     mCW, sCW);

            const scalar covCD = mCD > SMALL ? sCD/mCD : GREAT;
            const scalar covCH = mCH > SMALL ? sCH/mCH : GREAT;
            const scalar covCW = mCW > SMALL ? sCW/mCW : GREAT;
            const scalar covCladMax = max(covCD, max(covCH, covCW));

            // Steady once the clad geometry is populated at every interior
            // probe and translationally invariant (low CoV).
            bool allCladPopulated = true;
            for (label i = i0; i <= i1; ++i)
            {
                if (cladHeight[i] <= SMALL) { allCladPopulated = false; break; }
            }
            const label steady =
                (allCladPopulated && covCladMax < steadyCoVThreshold) ? 1 : 0;

            fileName f =
                caseRoot/"post-processing-data"/"pool_metrics_steady.csv";
            const bool header = !isFile(f);
            std::ofstream os(f.c_str(), std::ios_base::app);
            if (os.is_open())
            {
                if (header)
                {
                    os << "Time_s,nInterior,"
                       << "poolDepth_mm,poolDepth_std,"
                       << "poolHalfWidth_mm,poolHalfWidth_std,"
                       << "cladDepth_mm,cladDepth_std,"
                       << "cladHeight_mm,cladHeight_std,"
                       << "cladHalfWidth_mm,cladHalfWidth_std,"
                       << "cladCoVmax,STEADY\n";
                }
                os << runTime.value() << "," << nInt
                   << "," << mPD*1e3 << "," << sPD*1e3
                   << "," << mPW*1e3 << "," << sPW*1e3
                   << "," << mCD*1e3 << "," << sCD*1e3
                   << "," << mCH*1e3 << "," << sCH*1e3
                   << "," << mCW*1e3 << "," << sCW*1e3
                   << "," << covCladMax << "," << steady << "\n";
                os.close();
            }

            Info<< "[PoolMetrics] clad(mm) 0.5W=" << mCW*1e3
                << " h=" << mCH*1e3 << " D=" << mCD*1e3
                << " | pool(mm) 0.5W=" << mPW*1e3 << " D=" << mPD*1e3
                << " | cladCoV=" << covCladMax
                << (steady ? " STEADY" : "") << endl;
        }
    }
}
