/*--------------------------------*- C++ -*----------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Version:  10
     \\/     M anipulation  |
\*---------------------------------------------------------------------------*/
FoamFile
{
    format      ascii;
    class       dictionary;
    location    "constant";
    object      physicalProperties.water;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

viscosityModel  constant;

nu               6.7e-7;

rho              7600; // original 7578

elec_resistivity  1e-6; //original 1.25e-6




// // Thermal conductivity table [W/(m·K)]
// table_kappa_ v1
// (
//    (300    9.0)
//    (400    11.0)
//   (500    13.0)
//    (600    18.0)
//   (800    21.0)
//    (1000   27.0)
//   (1100   27.0)
//    (1300   28.0)
//    (1533   32.0)
//    (1609   32.0)
//    (2500   32.0)
//    (5000   32.0) // changed from 31 
// );

table_kappa
(
    (300    8.55)
    (400    10.45)
    (500    12.35)
    (600    17.10)
    (800    19.95)
    (1000   25.65)
    (1100   25.65)
    (1300   29.79)
    (1531   30.86)
    (1612   31.81)
    (2500   32.98)
    (5000   32.98)
);

// table_kappa
// (
//     (300  9.0)
//     (400  11.0)
//     (500  13.0)
//     (600  18.0)
//     (800  22.0)
//     (1000 22.0)
//     (1100 22.0)
//     (1300 23.0)    // Reduced gradient - was 28.0
//     (1531 23.0)    // Reduced gradient - was 29.0  
//     (1612 32.0)    // Reduced gradient - was 30.0
//     (2500 32.0)    // Reduced gradient - was 31.0
//     (5000 32.0)    // Flatter at high temps
// );
// Specific heat capacity table [J/(kg·K)]
table_cp
(
    (298   425)
    (400   447)
    (500   468)
    (600   489)
    (700   510)
    (800   531)
    (1100  620)
    (1200  620)
    (1300  620)
    (1400  620)
    (2500  630)   // Hold constant beyond 1400 K
    (5000  630)   // increase cp ---wider width
);
   
  Tsolidus 1533;
  Tliquidus 1609;
    LatentHeat 250e3;
    beta    1.3e-5;


// ************************************************************************* //