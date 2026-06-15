/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2020 OpenCFD Ltd.
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
    laserbeamFoam

Group
    grpMultiphaseSolvers

Description
    Ray-Tracing heat source implementation with two phase incompressible VoF
    description of the metallic substrate and shielding gas phase,
    with optional mesh motion and mesh topology changes including adaptive
    re-meshing.
Authors

    Tom Flint, UoM.
    Philip Cardiff, UCD.
    Gowthaman Parivendhan, UCD.
    Joe Robson, UoM.
    Petar Cosic, UCD
    Simon Rodriguez, UCD

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include <fstream>
#include "dynamicFvMesh.H"
#include "isoAdvection.H"
#include "CMULES.H"
#include "EulerDdtScheme.H"
#include "localEulerDdtScheme.H"
#include "CrankNicolsonDdtScheme.H"
#include "subCycle.H"
#include "immiscibleIncompressibleTwoPhaseMixture.H"
#include "incompressibleInterPhaseTransportModel.H"
#include "turbulentTransportModel.H"
#include "pimpleControl.H"
#include "fvOptions.H"
#include "CorrectPhi.H"
#include "fvcSmooth.H"
#include "dynamicRefineFvMesh.H"

#include "Polynomial.H"
#include "thermoScalarProperty.H"
#include "laserHeatSource.H"
#include "mthdModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::addNote
    (
        "Solver for two incompressible, isothermal immiscible fluids"
        " using VOF phase-fraction based interface capturing.\n"
        "With optional mesh motion and mesh topology changes including"
        " adaptive re-meshing."
    );

    #include "postProcess.H"

    #include "addCheckCaseOptions.H"
    #include "setRootCaseLists.H"
    #include "createTime.H"
    #include "createDynamicFvMesh.H"
    #include "initContinuityErrs.H"
    #include "createDyMControls.H"
    #include "createFields.H"
    #include "MULES/createAlphaFluxes.H"
    #include "initCorrectPhi.H"
    #include "createUfIfPresent.H"

    if (interfaceTrackingScheme == "MULES")
    {
        if (!LTS)
        {
            #include "MULES/CourantNo.H"
            #include "setInitialDeltaT.H"
        }
    }
    else if (interfaceTrackingScheme == "isoAdvector")
    {
        #include "isoAdvector/porousCourantNo.H"
        #include "setInitialDeltaT.H"
    }

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
    Info<< "\nStarting time loop\n" << endl;

    while (runTime.run())
    {
        #include "readControls.H"
        #include "readDyMControls.H"

        if (interfaceTrackingScheme == "MULES")
        {
            if (LTS)
            {
                #include "MULES/setRDeltaT.H"
            }
            else
            {
                #include "MULES/CourantNo.H"
                #include "MULES/alphaCourantNo.H"
                #include "MULES/setDeltaT.H"
            }
        }
        else if (interfaceTrackingScheme == "isoAdvector")
        {
            #include "isoAdvector/porousCourantNo.H"
            #include "isoAdvector/porousAlphaCourantNo.H"
            #include "isoAdvector/setDeltaT.H"
        }

        ++runTime;

        Info<< "Time = " << runTime.timeName() << nl << endl;

        // Automated track ID calculation
        label currentTrack = floor(runTime.value()/trackDuration) + 1;

        // --- Pressure-velocity PIMPLE corrector loop
        while (pimple.loop())
        {

            if (interfaceTrackingScheme == "MULES")
            {
                #include "MULES/firstIter.H"
                #include "MULES/alphaControls.H"
                #include "MULES/alphaEqnSubCycle.H"
            }
            else if (interfaceTrackingScheme == "isoAdvector")
            {
                #include "isoAdvector/firstIter.H"
                #include "isoAdvector/alphaControls.H"
                #include "isoAdvector/alphaEqnSubCycle.H"
            }

            #include "updateProps.H"

            // Update the laser deposition field
            laser.updateDeposition
            (
                alpha_filtered, n_filtered, electrical_resistivity
            );

            mixture.correct();

            if (pimple.frozenFlow())
            {
                continue;
            }

            #include "UEqn.H"
            if (mthd.valid())
            {
                mthd->solve(phi, U);
            }
            #include "TEqn.H"

            // --- Pressure corrector loop
            while (pimple.correct())
            {
                #include "pEqn.H"
            }

            if (pimple.turbCorr())
            {
                turbulence->correct();
            }
        }

        // Update the melt history
        const volScalarField& alphaMetal =
            mesh.lookupObject<volScalarField>("alpha.metal");
        condition = pos(alphaMetal - 0.5) * pos(epsilon1 - 0.5);
        meltHistory += condition;

        forAll(meltTrackID, celli)
        {
            if (condition[celli] > 0.5)
            {
                meltTrackID[celli] = currentTrack;
            }
            else if (meltTrackID[celli] > 0.5)
            {
                meltTrackID[celli] = floor(meltTrackID[celli] + 0.5);
            }
        }
        meltTrackID.correctBoundaryConditions();

        // --- Calculate Solidification Metrics ---
        Info << "[DEBUG-Solidification] Updating fields G, R, CR_spatial, CR_temporal..." << endl;
        
        // Ensure directory exists for post-processing
        if (Pstream::master())
        {
            mkDir(runTime.path()/"post-processing-data");
        }

        const volScalarField metalMask(pos(alphaMetal - 0.5));

        // 1. Thermal Gradient
        G = mag(fvc::grad(T)) * metalMask;
        G.correctBoundaryConditions();

        // 2. Solidification Rate R_solid (with automatically calculated scanDirection)
        interpolationTable<vector> laserPosTable(laser.subDict("timeVsLaserPosition"));
        scalar tCurrent = runTime.value();
        vector posCurrent = laserPosTable(tCurrent);
        vector posNext = laserPosTable(tCurrent + 1e-6);
        vector dir = posNext - posCurrent;
        scalar magDir = mag(dir);
        vector normScanDir = vector::zero;
        if (magDir > 1e-12)
        {
            normScanDir = dir / magDir;
        }

        dimensionedScalar smallGradT("smallGradT", dimensionSet(0, -1, 0, 1, 0), 1e-6);
        volScalarField cosTheta((fvc::grad(T) & normScanDir) / (mag(fvc::grad(T)) + smallGradT));
        R_solid = V_scan * max(cosTheta, 0.0) * metalMask;
        R_solid.correctBoundaryConditions();

        // 3. Spatial Cooling Rate (non-negative)
        CR_spatial = G * R_solid;
        CR_spatial.correctBoundaryConditions();

        // 4. Temporal Cooling Rate (absolute value magnitude)
        CR_temporal = mag(fvc::ddt(T)) * metalMask;
        CR_temporal.correctBoundaryConditions();

        Info << "[DEBUG-Solidification] Solidification fields updated." << endl;

        // 5. Native CSV Logging for specific points from melt pool bottom
        Info << "[DEBUG-Solidification] Locating melt pool bottom..." << endl;
        label bottomCelli = -1;
        scalar maxMeltPoolY = -1.0;

        forAll(T, celli)
        {
            // check liquidus/solidus region in metal
            if (alphaMetal[celli] > 0.5 && T[celli] >= 1658.0)
            {
                scalar cellY = mesh.C()[celli].y();
                if (cellY > maxMeltPoolY)
                {
                    maxMeltPoolY = cellY;
                    bottomCelli = celli;
                }
            }
        }

        scalar globalMaxMeltPoolY = maxMeltPoolY;
        reduce(globalMaxMeltPoolY, maxOp<scalar>());
        point bottomPoint(-1e9, -1e9, -1e9);

        if (bottomCelli != -1 && mag(maxMeltPoolY - globalMaxMeltPoolY) < 1e-12)
        {
            bottomPoint = mesh.C()[bottomCelli];
        }

        // Parallel reduction to share the bottom point coordinates across all processors
        reduce(bottomPoint.x(), maxOp<scalar>());
        reduce(bottomPoint.y(), maxOp<scalar>());
        reduce(bottomPoint.z(), maxOp<scalar>());

        if (globalMaxMeltPoolY >= 0.0)
        {
            Info << "[DEBUG-Solidification] Melt pool bottom found at: " << bottomPoint << endl;
            
            scalarList distances(4);
            distances[0] = 20.0;
            distances[1] = 80.0;
            distances[2] = 140.0;
            distances[3] = 200.0;

            scalarList globalT(4, 0.0);
            scalarList globalG(4, 0.0);
            scalarList globalR(4, 0.0);
            scalarList globalCRS(4, 0.0);
            scalarList globalCRT(4, 0.0);

            Info << "[DEBUG-Solidification] Finding closest cells for distances..." << endl;
            forAll(distances, i)
            {
                scalar d = distances[i];
                // Y increases downwards, so subtracting distance goes upwards towards surface
                point targetP(bottomPoint.x(), bottomPoint.y() - d * 1e-6, bottomPoint.z());
                
                label closestCelli = -1;
                scalar minDistanceSq = 1e10;
                
                forAll(mesh.C(), celli)
                {
                    scalar distSq = magSqr(mesh.C()[celli] - targetP);
                    if (distSq < minDistanceSq)
                    {
                        minDistanceSq = distSq;
                        closestCelli = celli;
                    }
                }
                
                scalar globalMinDistSq = minDistanceSq;
                reduce(globalMinDistSq, minOp<scalar>());
                
                if (closestCelli != -1 && mag(minDistanceSq - globalMinDistSq) < 1e-20)
                {
                    globalT[i] = T[closestCelli];
                    globalG[i] = G[closestCelli];
                    globalR[i] = R_solid[closestCelli];
                    globalCRS[i] = CR_spatial[closestCelli];
                    globalCRT[i] = CR_temporal[closestCelli];
                }
                
                reduce(globalT[i], sumOp<scalar>());
                reduce(globalG[i], sumOp<scalar>());
                reduce(globalR[i], sumOp<scalar>());
                reduce(globalCRS[i], sumOp<scalar>());
                reduce(globalCRT[i], sumOp<scalar>());
            }

            if (Pstream::master() && runTime.outputTime())
            {
                Info << "[DEBUG-Solidification] Writing to cooling_rate_comparison.csv..." << endl;
                fileName csvPath = runTime.path()/"post-processing-data"/"cooling_rate_comparison.csv";
                bool writeHeader = !isFile(csvPath);
                
                std::ofstream csvFile;
                csvFile.open(csvPath.c_str(), std::ios_base::app);
                
                if (csvFile.is_open())
                {
                    if (writeHeader)
                    {
                        csvFile << "Time_s,MeltPoolBottom_Y_m,MeltPoolBottom_Z_m,"
                                << "T_20um,G_20um,R_20um,CRS_20um,CRT_20um,"
                                << "T_80um,G_80um,R_80um,CRS_80um,CRT_80um,"
                                << "T_140um,G_140um,R_140um,CRS_140um,CRT_140um,"
                                << "T_200um,G_200um,R_200um,CRS_200um,CRT_200um\n";
                    }
                    
                    csvFile << runTime.value() << "," << bottomPoint.y() << "," << bottomPoint.z();
                    forAll(distances, i)
                    {
                        csvFile << "," << globalT[i]
                                << "," << globalG[i]
                                << "," << globalR[i]
                                << "," << globalCRS[i]
                                << "," << globalCRT[i];
                    }
                    csvFile << "\n";
                    csvFile.close();
                }
                else
                {
                    Info << "[WARNING-Solidification] Failed to open CSV file for writing!" << endl;
                }
            }
        }
        else
        {
            Info << "[DEBUG-Solidification] No active melt pool found. Skipping CSV log." << endl;
        }

        Info<< "TMax = " << gMax(T.primitiveField())
            << ", pVapMax = " << gMax(pVap.primitiveField()) << endl;

        runTime.write();

        // Write ray paths to VTK files
        if (runTime.outputTime())
        {
            laser.writeRayPathsToVTK();
            laser.writeRayPathVTKSeriesFile();
        }

        runTime.printExecutionTime(Info);
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
