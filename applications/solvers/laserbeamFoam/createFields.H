#include "createRDeltaT.H"

const dimensionedScalar deltaN = 1e-8/pow(average(mesh.V()), 1.0/3.0);

// Read track properties from a dedicated file
IOdictionary trackProperties
(
    IOobject
    (
        "trackProperties",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    )
);

const scalar trackDuration(readScalar(trackProperties.lookup("trackDuration")));
const label meshRefinementBufferLayers
(
    trackProperties.lookupOrDefault<label>("meshRefinementBufferLayers", 2)
);

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

Info<< "Reading field Temperature\n" << endl;
volScalarField T
(
    IOobject
    (
        "T",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

#include "createPhi.H"


Info<< "Reading transportProperties\n" << endl;
immiscibleIncompressibleTwoPhaseMixture mixture(U, phi);

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


// Mass flux
surfaceScalarField rhoPhi
(
    IOobject
    (
        "rhoPhi",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    fvc::interpolate(rho)*phi
);


IOdictionary transportProperties
(
    IOobject
    (
        "transportProperties",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ_IF_MODIFIED,
        IOobject::NO_WRITE
    )
);

const word interfaceTrackingScheme
(
    transportProperties.lookup("interfaceTrackingScheme")
);

if
(
    interfaceTrackingScheme != "MULES"
 && interfaceTrackingScheme != "isoAdvector"
)
{
    FatalErrorIn("laserbeamFoam.C")
        << "Unknown interface tracking scheme "
        << interfaceTrackingScheme
        << ". Available options are 'MULES' and 'isoAdvector'"
        << exit(FatalError);
}

dictionary& transportPropertiesMetal = transportProperties.subDict("metal");
dictionary& transportPropertiesGas = transportProperties.subDict("gas");


const dimensionedScalar Tsolidus1
(
    "Tsolidus",
    dimensionSet(0, 0, 0, 1, 0),
    transportPropertiesMetal//.lookup("Tsolidus")
);
// Info<<Tsolidus1<<endl;
const dimensionedScalar Tsolidus2
(
    "Tsolidus",
    dimensionSet(0, 0, 0, 1, 0),
    transportPropertiesGas//.lookup("Tsolidus")
);
// Info<<Tsolidus2<<endl;

const dimensionedScalar Tliquidus1
(
    "Tliquidus",
    dimensionSet(0, 0, 0, 1, 0),
    transportPropertiesMetal//.lookup("Tliquidus")
);

const dimensionedScalar Tliquidus2
(
    "Tliquidus",
    dimensionSet(0, 0, 0, 1, 0),
    transportPropertiesGas//.lookup("Tliquidus")
);

const dimensionedScalar LatentHeat1
(
    "LatentHeat",
    dimensionSet(0, 2, -2, 0, 0),
    transportPropertiesMetal//.lookup("LatentHeat")
);

const dimensionedScalar LatentHeat2
(
    "LatentHeat",
    dimensionSet(0, 2, -2, 0, 0),
    transportPropertiesGas//.lookup("LatentHeat")
);

const dimensionedScalar beta1
(
    "beta",
    dimensionSet(0, 0, 0, -1, 0),
    transportPropertiesMetal//.lookup("beta")
);

const dimensionedScalar beta2
(
    "beta",
    dimensionSet(0, 0, 0, -1, 0),
    transportPropertiesGas//.lookup("beta")
);

const dimensionedScalar Marangoni_Constant
(
    "Marangoni_Constant",
    dimensionSet(1, 0, -2, -1, 0),
    transportProperties//.lookup("dsigmadT")
);

// Temperature-dependent surface tension σ(T) = sigma + Marangoni_Constant*(T - Tliquidus1)
// Backwards compatible: if Marangoni_Constant is zero, uses constant sigma (original behaviour).
const dimensionedScalar sigma_ref
(
    "sigma_ref",
    dimensionSet(1, 0, -2, 0, 0, 0, 0),
    transportProperties.lookupOrDefault<scalar>("sigma", 1.0)
);

const dimensionedScalar sigma_T_ref("sigma_T_ref", dimTemperature, Tliquidus1.value());

const Switch useTdepSigma
(
    transportProperties.lookupOrDefault<Switch>("variableSurfaceTension", false)
);

volScalarField sigmaField
(
    IOobject
    (
        "sigmaField",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    sigma_ref
);

const dimensionedScalar p0
(
    "p0",
    dimensionSet(1, -1, -2, 0, 0, 0, 0),
    transportProperties//.lookup("p0")
); // atmospheric pressure

const dimensionedScalar p_amb
(
    transportProperties.found("p_amb")
  ? dimensionedScalar
    (
        "p_amb",
        p0.dimensions(),
        transportProperties.lookup("p_amb")
    )
  : dimensionedScalar("p_amb", p0.dimensions(), p0.value())
); // ambient pressure for recoil model

const dimensionedScalar recoil_coeff
(
    "recoil_coeff",
    dimless,
    transportProperties.lookupOrDefault<scalar>("recoil_coeff", 0.54)
);

const dimensionedScalar Tvap
(
    "Tvap",
    dimensionSet(0, 0, 0, 1, 0),
    transportProperties//.lookup("Tvap")
);

const dimensionedScalar Tvap_atm
(
    transportProperties.found("Tvap_atm")
  ? dimensionedScalar
    (
        "Tvap_atm",
        Tvap.dimensions(),
        transportProperties.lookup("Tvap_atm")
    )
  : dimensionedScalar("Tvap_atm", Tvap.dimensions(), Tvap.value())
);

// sigma floor: σ held constant at this value for T > Tvap
const dimensionedScalar sigma_at_Tvap
(
    "sigma_at_Tvap",
    sigma_ref.dimensions(),
    sigma_ref.value() + Marangoni_Constant.value()*(Tvap.value() - sigma_T_ref.value())
);

const dimensionedScalar recoilTransitionWidth
(
    "recoilTransitionWidth",
    dimTemperature,
    transportProperties.lookupOrDefault<scalar>("recoilTransitionWidth", 100.0)
);

const word recoilPressureModel
(
    transportProperties.lookupOrDefault<word>
    (
        "recoilPressureModel",
        "modifiedReducedPressure"
    )
);

// Cap on temperature fed to the recoil pressure CC formula.
// T can rise freely in the energy equation; only pVap evaluation is capped.
// Default GREAT = no cap (original behaviour).
const dimensionedScalar T_recoil_cap
(
    "T_recoil_cap",
    dimTemperature,
    transportProperties.lookupOrDefault<scalar>("T_recoil_cap", GREAT)
);

const dimensionedScalar Mm
(
    "Mm",
    dimensionSet(1, 0, 0, 0, -1),
    transportProperties//.lookup("Mm")
);

const dimensionedScalar LatentHeatVap
(
    "LatentHeatVap",
    dimensionSet(0, 2, -2, 0, 0),
    transportProperties//.lookup("LatentHeatVap")
);

// Lee model mass-transfer intensity factor [1/s].
// Volumetric evaporative cooling: S_E = -r * rho1 * alpha1 * Lv * (T - Tvap)/Tvap
// Derived from Hertz-Knudsen kinetics (Chen et al., 2020, IJHMT).
//   0       = disabled (only interface Qv via Hertz-Knudsen)
//   1e5     = moderate clamping (~300 K overshoot at 200 W / 50 µm spot)
//   1e7     = tight clamping (~30 K overshoot, near-instant equilibrium)
const dimensionedScalar LeeCoeff
(
    "LeeCoeff",
    dimensionSet(0, 0, -1, 0, 0),
    transportProperties.lookupOrDefault<scalar>("LeeCoeff", 0.0)
);

const dimensionedScalar R
(
    "R",
    dimensionSet(1, 2, -2, -1, -1, 0, 0),
        8.314
);

volScalarField em
(
    IOobject
    (
        "em",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("zero", dimless, 0.0)
);

const dimensionedScalar SB
(
    "SB",
    dimensionSet(1, 0, -3, -4, 0, 0, 0),
    scalar(5.67e-8)
); // Stefan-Boltzmann Constant

[[maybe_unused]]const scalar emL(readScalar(transportProperties.lookup("emL")));
[[maybe_unused]]const scalar emS(readScalar(transportProperties.lookup("emS")));
const dimensionedScalar TRef
(
    "TRef",
    dimTemperature,
    transportProperties.lookup("TRef")
);

const dimensionedScalar T_ambient
(
    transportProperties.found("T_ambient")
  ? dimensionedScalar
    (
        "T_ambient",
        TRef.dimensions(),
        transportProperties.lookup("T_ambient")
    )
  : dimensionedScalar("T_ambient", TRef.dimensions(), TRef.value())
);

const dimensionedScalar h_convection
(
    "h_convection",
    dimensionSet(1, 0, -3, -1, 0, 0, 0),
    transportProperties.lookupOrDefault<scalar>("h_convection", 0.0)
);

const thermoScalarProperty kappaProp_m
(
    transportPropertiesMetal,
    "table_kappa",
    "poly_kappa",
    "metal kappa"
);

const Polynomial<8> polykappa_g(transportPropertiesGas.lookup("poly_kappa"));

const thermoScalarProperty cpProp_m
(
    transportPropertiesMetal,
    "table_cp",
    "poly_cp",
    "metal cp"
);

const Polynomial<8> polycp_g(transportPropertiesGas.lookup("poly_cp"));

volScalarField cp
(
    IOobject
    (
        "cp",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("cp", dimensionSet(0, 2, -2, -1, 0), 0.0)
);

volScalarField kappa
(
    IOobject
    (
        "kappa",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("kappa",dimensionSet(1, 1, -3, -1, 0),0.0)
);

#include "updateKappaCp.H"


volScalarField TSolidus
(
    IOobject
    (
        "TSolidus",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    alpha1*Tsolidus1 + (1.0 - alpha1)*Tsolidus2
);


volScalarField TLiquidus
(
    IOobject
    (
        "TLiquidus",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    alpha1*Tliquidus1 + (1.0 - alpha1)*Tliquidus2
);

volScalarField DC
(
    IOobject
    (
        "DC",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("DC", dimensionSet(1,-3,-1,0,0), 1.0e14)
);

volScalarField metalaverage
(
    IOobject
    (
        "metalaverage",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::NO_WRITE
    ),
    fvc::average(alpha1)
);

volScalarField epsilon1
(
    IOobject
    (
        "epsilon1",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    max(min((T - TSolidus)/(TLiquidus - TSolidus), 1.0), 0.0)
);

volVectorField n_filtered
(
    IOobject
    (
        "n_filtered",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedVector("n_filtered", dimless, vector::zero)//,
    // zeroGradientFvPatchVectorField::typeName
);

volScalarField epsilon1mask
(
    IOobject
    (
        "epsilon1mask",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("epsilon1mask",dimless,0.0)
);

volScalarField epsilon1mask2
(
    IOobject
    (
        "epsilon1mask2",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("epsilon1mask2", dimless, 0.0)
);

volVectorField gradT
(
    IOobject
    (
        "gradT",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedVector("gradT", dimensionSet(0, -1, 0, 1, 0), vector::zero)
);

volScalarField LatentHeat
(
    IOobject
    (
        "LatentHeat",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    alpha1*LatentHeat1 + (1.0 - alpha1)*LatentHeat2
);

volScalarField TRHS
(
    IOobject
    (
        "TRHS",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("TRHS", dimensionSet(1, -1, -3, 0, 0), 0.0)
);

volScalarField Tcorr
(
    IOobject
    (
        "Tcorr",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("Tcorr", dimensionSet(0, 0, 0, 1, 0), 0.0)
);

volScalarField beta
(
    IOobject
    (
        "beta",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    alpha1*beta1 + (1.0-alpha1)*beta2
);

volScalarField rhok
(
    IOobject
    (
        "rhok",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("rhok", dimless, 0.0)
);

volScalarField Num_divU
(
    IOobject
    (
        "Num_divU",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("Num_divU", dimensionSet(0,0,-1,0,0), 0.0)
);

volVectorField Marangoni
(
    IOobject
    (
        "Marangoni",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedVector("Marangoni", dimensionSet(1, -2, -2, 0, 0), vector::zero)
);

const volVectorField zeroJcrossB
(
    IOobject
    (
        "zeroJcrossB",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedVector("zeroJcrossB", dimensionSet(1, -2, -2, 0, 0), vector::zero)
);

const volScalarField zeroJoulesource
(
    IOobject
    (
        "zeroJoulesource",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("zeroJoulesource", dimensionSet(1, -1, -3, 0, 0), 0.0)
);

// Recoil Pressure
volScalarField pVap
(
    IOobject
    (
        "pVap",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("0", dimensionSet(1,-1,-2,0,0,0,0), 0.0)
);

volScalarField recoilPressure
(
    IOobject
    (
        "recoilPressure",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    pVap
);

// Evaporation Cooling
volScalarField Qv
(
    IOobject
    (
        "Qv",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("0", dimensionSet(1,0,-3,0,0,0,0), 0.0)
);


volScalarField damper
(
    IOobject
    (
        "damper",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("damper", dimless, 1.0),
    zeroGradientFvPatchScalarField::typeName
);

volScalarField thermalDamper
(
    IOobject
    (
        "thermalDamper",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("thermalDamper", dimless, 1.0),
    zeroGradientFvPatchScalarField::typeName
);

volScalarField alpha_filtered
(
    IOobject
    (
        "alpha_filtered",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("alpha_filtered", dimensionSet(0, 0, 0, -0, 0), 0.0)
);

volScalarField alpha_smoothed
(
    IOobject
    (
        "alpha_smoothed",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("alpha_smoothed", dimensionSet(0, 0, 0, -0, 0), 0.0)
);

volScalarField amrIndicator
(
    IOobject
    (
        "amrIndicator",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("amrIndicator", dimless, 0.0)
);

volScalarField e1temp
(
    IOobject
    (
        "e1temp",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("e1temp", dimensionSet(0, 0, 0, -0, 0), 0.0)
);
e1temp = fvc::average(epsilon1);

volScalarField rhoCp
(
    IOobject
    (
        "rhoCp",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("rhoCp", dimensionSet(1, -1, -2, -1, 0), 0.0)
);
rhoCp = (rho*cp);
// volScalarField rhoCp(rho*cp);
rhoCp.oldTime();


// surfaceScalarField rhophicp(fvc::interpolate(cp)*rhoPhi);

surfaceScalarField rhophicp
(
    IOobject
    (
        "rhophicp",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("rhophicp", dimensionSet(1, 2, -3, -1, 0), 0.0)
);
rhophicp = fvc::interpolate(cp)*rhoPhi;











const dimensionedScalar DarcyConstantlarge
(
    "DarcyConstantlarge", dimensionSet(1,-3,-1,0,0), 1.0e6
);

const dimensionedScalar DarcyConstantsmall
(
    "DarcyConstantsmall", dimless, 1.0e-12
);



// volVectorField gradAlpha(fvc::grad(alpha1));

volVectorField gradAlpha
(
    IOobject
    (
        "gradAlpha",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedVector("gradAlpha", dimensionSet(0, -1, 0, -0, 0), vector::zero)
);
gradAlpha = fvc::grad(alpha1);

// volVectorField nHatM(gradAlpha/(mag(gradAlpha) + deltaN));


volVectorField nHatM
(
    IOobject
    (
        "nHatM",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedVector("nHatM", dimensionSet(0, 0, 0, -0, 0), vector::zero)
);
nHatM = gradAlpha/(mag(gradAlpha) + deltaN);


// This variable is >=0 if the cell has been liquid
volScalarField meltHistory
(
    IOobject
    (
        "meltHistory",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("meltHistory", dimensionSet(0, 0, 0, -0, 0), 0.0)
);

// Peak temperature reached by each cell over the simulation.
volScalarField TmaxHistory
(
    IOobject
    (
        "TmaxHistory",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    T
);

// This variable stores the latest laser track ID that melted each cell
volScalarField meltTrackID
(
    IOobject
    (
        "meltTrackID",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("meltTrackID", dimless, 0.0)
);

// Mesh-only refinement marker. Keep amrIndicator unchanged for diagnostics.
volScalarField meshRefinementIndicator
(
    IOobject
    (
        "meshRefinementIndicator",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("meshRefinementIndicator", dimless, 0.0)
);

// This variable is 1 if (alpha.metal >= 0.5 and epsilon1 >=0.5), 0 otherwise
volScalarField condition
(
    IOobject
    (
        "condition",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("condition", dimensionSet(0, 0, 0, -0, 0), 0.0)
);




typedef incompressibleInterPhaseTransportModel
    <
        immiscibleIncompressibleTwoPhaseMixture
    > transportModelType;

autoPtr<transportModelType> turbulence
(
    new transportModelType
    (
        rho, U, phi, rhoPhi, mixture
    )
);

#include "readGravitationalAcceleration.H"
#include "readhRef.H"
#include "gh.H"





autoPtr<mthdModel> mthd(nullptr);

if (mthdModel::electromagneticsEnabled(mesh, "transportProperties"))
{
    mthd.reset(new mthdModel(mesh, "transportProperties"));
}

// Check if temperature-dependent resistivity tables are provided
bool hasResistivityTable = 
    transportPropertiesMetal.found("table_resistivity")
 || transportPropertiesMetal.found("poly_resistivity");

volScalarField electrical_resistivity
(
    IOobject
    (
        "electrical_resistivity",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    hasResistivityTable
  ? dimensionedScalar("res", dimless, 5e-7)
  : 1.0 / mthdModel::readElectricalConductivity(transportPropertiesMetal)
);

autoPtr<thermoScalarProperty> resistivityProp_m(nullptr);
if (hasResistivityTable)
{
    resistivityProp_m.reset
    (
        new thermoScalarProperty
        (
            transportPropertiesMetal,
            "table_resistivity",
            "poly_resistivity",
            "metal electrical resistivity"
        )
    );
    Info<< "Using temperature-dependent electrical resistivity" << endl;
}
else
{
    Info<< "Using constant electrical resistivity from dictionary." << endl;
}




// Create laser heat source
laserHeatSource laser(mesh);

#include "updateYDim.H"
#include "updateProps.H"

// volScalarField depositionTemp
// (
//     IOobject
//     (
//         "depositionTemp",
//         runTime.timeName(),
//         mesh,
//         IOobject::NO_READ,
//         IOobject::NO_WRITE
//     ),
//     mesh,
//     dimensionedScalar("depositionTemp", dimensionSet(1, -1, -2, -1, 0), 0.0),
//     zeroGradientFvPatchScalarField::typeName
// );


// depositionTemp=laser.deposition();
// depositionTemp.correctBoundaryConditions();

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
    p_rgh + rho*gh*rhok
);

label pRefCell = 0;
scalar pRefValue = 0.0;
setRefCell
(
    p,
    p_rgh,
    pimple.dict(),
    pRefCell,
    pRefValue
);

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

volScalarField pAbs
(
    IOobject
    (
        "pAbs",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    p + p0
);

mesh.setFluxRequired(p_rgh.name());
mesh.setFluxRequired(alpha1.name());

#include "createMRF.H"
#include "createFvOptions.H"


#include "isoAdvector/createPorosity.H"
isoAdvection advector(alpha1, phi, U);

// --- Solidification parameters and fields ---
Info<< "Reading transportProperties for V_scan and scanDirection..." << endl;
const dimensionedScalar V_scan
(
    "V_scan",
    dimVelocity,
    transportProperties.lookupOrDefault<scalar>("V_scan", 0.0)
);


Info<< "Declaring solidification fields..." << endl;
volScalarField G
(
    IOobject
    (
        "G",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("G", dimensionSet(0, -1, 0, 1, 0), 0.0)
);

volScalarField R_solid
(
    IOobject
    (
        "R_solid",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("R_solid", dimVelocity, 0.0)
);

volScalarField CR_spatial
(
    IOobject
    (
        "CR_spatial",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("CR_spatial", dimensionSet(0, 0, -1, 1, 0), 0.0)
);

volScalarField CR_temporal
(
    IOobject
    (
        "CR_temporal",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("CR_temporal", dimensionSet(0, 0, -1, 1, 0), 0.0)
);
