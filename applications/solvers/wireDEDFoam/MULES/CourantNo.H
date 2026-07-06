/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
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

Global
    CourantNo

Description
    Calculates and outputs the mean and maximum Courant Numbers.

\*---------------------------------------------------------------------------*/

scalar CoNum = 0.0;
scalar meanCoNum = 0.0;

{
    scalarField sumPhi
    (
        fvc::surfaceSum(mag(phi))().primitiveField()
    );

    CoNum = 0.5*gMax(sumPhi/mesh.V().field())*runTime.deltaTValue();

    meanCoNum =
        0.5*(gSum(sumPhi)/gSum(mesh.V().field()))*runTime.deltaTValue();

    const Switch debugNumerics
    (
        runTime.controlDict().lookupOrDefault<Switch>
        (
            "debugNumerics",
            false
        )
    );

    if (debugNumerics && sumPhi.size())
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

        const scalar debugCoFraction =
            runTime.controlDict().lookupOrDefault<scalar>
            (
                "debugNumericsCoFraction",
                0.95
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
            scalarField coByDt(sumPhi/mesh.V().field());
            label maxCelli = 0;

            forAll(coByDt, celli)
            {
                if (coByDt[celli] > coByDt[maxCelli])
                {
                    maxCelli = celli;
                }
            }

            const scalar localCo =
                0.5*coByDt[maxCelli]*runTime.deltaTValue();

            if (localCo >= debugCoFraction*CoNum || localCo > maxCo)
            {
                scalar level = -1;
                if (mesh.foundObject<volScalarField>("cellLevel"))
                {
                    level =
                        mesh.lookupObject<volScalarField>
                        (
                            "cellLevel"
                        )[maxCelli];
                }

                Info<< "[DBG Courant] time=" << runTime.timeName()
                    << " proc=" << Pstream::myProcNo()
                    << " cell=" << maxCelli
                    << " C=" << mesh.C()[maxCelli]
                    << " V=" << mesh.V()[maxCelli]
                    << " level=" << level
                    << " dt=" << runTime.deltaTValue()
                    << " localCo=" << localCo
                    << " globalCo=" << CoNum
                    << " sumPhi=" << sumPhi[maxCelli]
                    << " magU=" << mag(U[maxCelli])
                    << " U=" << U[maxCelli]
                    << " alpha=" << alpha1[maxCelli]
                    << " T=" << T[maxCelli]
                    << " pVap=" << pVap[maxCelli]
                    << " recoil=" << recoilPressure[maxCelli]
                    << " p_rgh=" << p_rgh[maxCelli]
                    << endl;
            }
        }
    }
}

Info<< "Courant Number mean: " << meanCoNum
    << " max: " << CoNum << endl;

// ************************************************************************* //
