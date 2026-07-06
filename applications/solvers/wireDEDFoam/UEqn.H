MRF.correctBoundaryVelocity(U);


if (laser.powderSim())
{
    DC =
        DarcyConstantlarge*Foam::pow((1.0-epsilon1mask),2)
       /(Foam::pow(epsilon1mask, 3) + DarcyConstantsmall);
}
else
{
    DC =
        DarcyConstantlarge*Foam::pow((1.0 - epsilon1), 2)
       /(Foam::pow(epsilon1, 3) + DarcyConstantsmall);
}

DC.correctBoundaryConditions();  // ADDED

gradAlpha = fvc::grad(alpha1);
gradAlpha.correctBoundaryConditions();  // ADDED - CRITICAL!

nHatM = gradAlpha/(mag(gradAlpha) + deltaN);
nHatM.correctBoundaryConditions();  // ADDED

Marangoni =
    epsilon1*dSigmaDtf*(gradT - nHatM*(nHatM & gradT))*mag(gradAlpha);

const scalar Tliquidus1Value = Tliquidus1.value();
const scalarField& TI = T;
vectorField& MarangoniI = Marangoni;
forAll(MarangoniI, celli)
{
    if(TI[celli] < Tliquidus1Value)
    {
        MarangoniI[celli] = vector::zero;
    }
}

Marangoni.correctBoundaryConditions();

if (damperSwitch)
{
    // To account for high density difference
    damper = 2.0*rho/(rho1 + rho2);
    damper.correctBoundaryConditions();  // ADDED
}

{
    const volScalarField T_rc(min(T, T_recoil_cap));
    pVap = recoil_coeff*p0*Foam::exp(LatentHeatVap*Mm*((T_rc - Tvap)/(R*T_rc*Tvap)));
}
pVap.correctBoundaryConditions();  // ADDED


if (mthd.valid())
{
    mthd->updateLorentzForce();
}

const volVectorField& JcrossB =
    mthd.valid()
  ? mthd->JcrossB()
  : zeroJcrossB;

fvVectorMatrix UEqn
(
    fvm::ddt(rho, U) + fvm::div(rhoPhi, U)
  + MRF.DDt(rho, U)
  + turbulence->divDevRhoReff(rho, U)
  + fvm::Sp(DC, U)
  - Marangoni*damper
  - JcrossB*damper
 ==
    fvOptions(rho, U)
    // Wire-DED (Strategy A): momentum of the injected feedstock.
    // Scaled by (1-alpha1) to match the bounded mass source, so momentum is
    // added only in proportion to the metal actually injected into gas.
    // S_U = mdot_vol * u_feed = rho1 * Sfeed*(1-alpha1) * feedVelocity  [N/m^3]
  + rho1*Sfeed*max(scalar(1) - alpha1, scalar(0))*feedVelocity
);

UEqn.relax();

fvOptions.constrain(UEqn);

if (pimple.momentumPredictor())
{
    solve
    (
        UEqn
     ==
        fvc::reconstruct
        (
            (
                (
                    useTdepSigma
                  ? fvc::interpolate(sigmaField/sigma_ref)*mixture.surfaceTensionForce()
                  : mixture.surfaceTensionForce()
                )*fvc::interpolate(damper)
                + (fvc::interpolate(pVap)
                    *fvc::snGrad(alpha1)
                    *fvc::interpolate(damper))
              - ghf*fvc::snGrad(rho*rhok)
              - fvc::snGrad(p_rgh)
            ) * mesh.magSf()
        )
    );

    fvOptions.correct(U);
}
