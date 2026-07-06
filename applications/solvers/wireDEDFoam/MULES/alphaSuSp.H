// Wire-DED (Strategy A): real Su/Sp/divU internal fields so the alpha
// equation carries the volumetric feed mass source.
//   Su  [1/s]  explicit source: gas -> metal conversion in the feed region
//   Sp  [1/s]  implicit source (unused here, kept zero)
//   divU [1/s] kept identically zero (matches the original zeroField
//              behaviour; the one-fluid model stays incompressible/divU-free)
volScalarField::Internal Su
(
    IOobject("Su", runTime.timeName(), mesh),
    mesh,
    dimensionedScalar(dimless/dimTime, Zero)
);

volScalarField::Internal Sp
(
    IOobject("Sp", runTime.timeName(), mesh),
    mesh,
    dimensionedScalar(dimless/dimTime, Zero)
);

volScalarField::Internal divU
(
    IOobject("divU", runTime.timeName(), mesh),
    mesh,
    dimensionedScalar(dimless/dimTime, Zero)
);

// Inject the feed mass source (converts gas to metal where Sfeed > 0).
// Split as Su + Sp*alpha1 = Sfeed*(1 - alpha1): the implicit Sp = -Sfeed
// makes the source self-limiting (stops as a cell fills, alpha1 -> 1), so
// alpha1 stays bounded in [0,1] instead of overshooting. In the gas feed
// zone (alpha1 ~ 0) the injected rate is still ~Sfeed = mdot/rho1.
Su += Sfeed();
Sp -= Sfeed();
