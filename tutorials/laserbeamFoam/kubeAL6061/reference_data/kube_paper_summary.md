# Kube Al6061 Paper Summary: Baseline Stationary Laser Case

Source files:

- `kube_al6061_paper.pdf`
- `supplementary_files/43246_2024_584_MOESM2_ESM.pdf`
- `supplementary_files/43246_2024_584_MOESM3_ESM.pdf`

Paper: Lovejoy Mutswatiwa et al., Communications Materials, 2024.

This summary only records the baseline stationary laser melting case.

## Experimental Case

- Material: bare Al6061 plate.
- Sample size: 20 mm length, 12 mm height, 1.5 mm thickness.
- Process type: stationary single-spot laser melting.
- Laser: continuous-wave ytterbium fiber laser, IPG YLR-500-AC.
- Laser wavelength: 1070 nm.
- Laser power: 350 W.
- Laser spot diameter: 80 um.
- Laser-on time: 3.34 ms.
- Scan speed: none reported; the laser was stationary.
- X-ray imaging location: Sector 32-ID-B, Advanced Photon Source, Argonne National Laboratory.
- X-ray imaging rate: 50 kHz.
- X-ray spatial resolution: 2 um/pixel.
- X-ray field of view: 1.8 mm x 1 mm.

## Simulation Case

- CFD platform: Flow-3D.
- Domain in paper body: 1 mm2.
- Domain in supplementary material: 0.5 mm x 0.5 mm x 0.8 mm.
- Mesh size: 4 um.
- Simulation finish time: 1.3 ms.
- Laser-on time in simulation: 0.8 ms.
- Laser power in simulation: 350 W.
- Laser spot diameter in simulation: 80 um.
- Execution time reported in supplementary material: about 15 hours.
- Fluid model: Newtonian, laminar flow, with temperature-dependent density.
- Bubble and phase-change model: constant-pressure bubbles with vaporization enabled.
- Evaporation pressure coefficients: A = 31000 g/(cm s2), B = 8.
- Rising pressure magnification: 0.2.
- Surface tension: 0.45 N/m.
- Surface tension temperature coefficient in Flow-3D: 0.001 g/s2/K.
- Specific heat ratio for aluminum vapor-pressure model: 2.1.
- Gravity: -981 cm/s2 in the Z direction.

The simulation laser duration was shorter than the experiment: 0.8 ms in Flow-3D versus 3.34 ms in the X-ray experiment.

## Experimental Data Available For Validation

The main validation data for the baseline case comes from high-speed synchrotron X-ray imaging.

- Time-resolved melt pool morphology.
- Time-resolved keyhole/vapor depression morphology.
- Melt pool depth.
- Keyhole depth.
- Melt pool width.
- Melt pool and keyhole aspect ratios.
- Bubble motion and keyhole-induced pore behavior.
- Final solidification structure observed in X-ray imaging.
- Optical image of the stationary laser-generated melt pool.
- EBSD grain map of the baseline melt pool.

The paper compares simulated and experimental melt pool geometry, keyhole morphology, and melt pool aspect ratio. It reports that the baseline case shows a deep and narrow keyhole in both the X-ray data and the Flow-3D simulation.

## Supplementary Movie Data

- Supplementary Movie 1: high-speed synchrotron X-ray sequence of a stationary Al6061 melt pool showing conduction and keyhole modes.
- Supplementary Movie 4: Flow-3D keyhole dynamics for the baseline case.
- Supplementary Movie 6: Flow-3D melt pool dynamics for the baseline case.
- Supplementary Movie 8: Flow-3D melt pool velocity vectors for the baseline case.
- Supplementary Movie 10: Flow-3D pressure distribution for the baseline case.

## Notes For Reproducing In LaserbeamFoam

- Use a stationary laser source; no scan speed is given for this baseline case.
- The strongest direct validation quantities are melt pool width, melt pool depth, keyhole shape/depth, aspect ratio, and time-resolved X-ray morphology.
- EBSD and optical images are useful for qualitative comparison of the final solidified melt pool, not as primary thermal-fluid calibration targets.
