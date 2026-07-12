// Temperature at the interface midplane, used to evaluate the temperature-
// dependent evaporation fluxes (recoil pressure in UEqn.H, Hertz-Knudsen
// evaporative cooling in TEqn.H).
//
// Both fluxes scale with exp(Lv*Mm*(T - Tvap)/(R*T*Tvap)). Evaluating that
// with the local cell temperature is the classical continuum-surface-flux
// approach, and it breaks down across a diffuse VOF interface: gas-side cells
// carry rho*cp some four orders of magnitude below the metal, so they overheat
// under the deposited laser power and then feed a runaway exponential into the
// recoil pressure. Weighting the temperature by the metal fraction and
// smoothing it across the interface band makes those cells adopt the
// liquid-side temperature instead.
//
// Meier et al., "Improved accuracy of continuum surface flux models for metal
// additive manufacturing melt pool simulations", arXiv:2401.12114 (2024).

if (interfaceFluxTemperature)
{
    const dimensionedScalar wFloor("wFloor", dimless, 1e-3);

    volScalarField w(max(alpha1, wFloor));
    volScalarField wT(w*T);

    for (label i = 0; i < nInterfaceTSmooth; i++)
    {
        wT = fvc::average(fvc::interpolate(wT));
        w = fvc::average(fvc::interpolate(w));
    }

    TInterface = wT/max(w, wFloor);
}
else
{
    TInterface = T;
}

TInterface.correctBoundaryConditions();
