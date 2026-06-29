NIST Data Publication:
Asynchronous AM Bench 2022 Challenge Data: Real-time, simultaneous absorptance
and high-speed x-ray imaging
Version 1.0.0
DOI: https://doi.org/10.18434/mds2-2525

Authors:
  Brian J Simonds
    National Institute of Standards and Technology
  Jack  Tanner
    National Institute of Standards and Technology
  Alexandra  Artusio-Glimpse
    National Institute of Standards and Technology
  Paul A. Williams
    National Institute of Standards and Technology
  Niranjan  Parab
    National Institute of Standards and Technology
  Cang  Zhao
    National Institute of Standards and Technology
  Tao  Sun
    National Institute of Standards and Technology

Contact:
  Brian Simonds
    brian.simonds@nist.gov

Description:

The absolute laser absorption was measured simultaneously with X-ray imaging during laser melting of Ti-6Al-4V solid metal. The data included here are the time-resolved absolute absorbed power and the X-ray images acquired at the same time, along with timing data for synchronization. Also included is information about the experimental configuration including applied laser power, laser beam spatial profile, and the experimental setup. A text document is included that describes all files. There are datasets from two experimental configurations: 1) A stationary laser of pulse duration 2.0 ms, and 2) a scanned laser with velocity 700 mm/s.

This collection is a supplement to:
  Simonds, B. J., Tanner, J., Artusio-Glimpse, A., Williams, P. A., Parab, N., Zhao, C., & Sun, T. (2021). The causal relationship between melt pool geometry and energy absorption measured in real time during laser-based manufacturing. Applied Materials Today, 23, 101049. doi:10.1016/j.apmt.2021.101049
  Khairallah, S. A., Sun, T., & Simonds, B. J. (2021). Onset of periodic oscillations as a precursor of a transition to pore-generating turbulence in  
laser melting. Additive Manufacturing Letters, 1, 100002. doi:10.1016/j.addlet.2021.100002
  Simonds, B. J., Tanner, J., Artusio-Glimpse, A., Williams, P. A., Parab, N., Zhao, C., & Sun, T. (2020). Simultaneous high-speed x-ray transmission imaging and absolute dynamic absorptance measurements during high-power laser-metal processing. Procedia CIRP, 94, 775â779. doi:10.1016/j.procir.2020.09.135
  Derimow, N., Schwalbach, E.J., Benzing, J.T., Kilgore, J.P., Artusio-Glimpse, A.B., Hrabe, N., & Simonds, B.J. (2021). In situ absorption synchrotron measurements, predictive modeling, microstructural analysis, and scanning probe measurements of laser melted Ti-6Al-4V single tracks for additive manufacturing applications. Journal of Alloys and Compounds, 900, 163494, doi:10.1016/j.jallcom.2021.163494


-------------------
General Information
-------------------

This data was collected: 12 June 2019 - 13 June 2019

--------------
Data Use Notes
--------------

This data is publicly available according to the NIST statements of
copyright, fair use and licensing; see https://www.nist.gov/director/copyright-fair-use-and-licensing-statements-srd-data-and-software

You may cite the use of this data as follows:
Simonds, Brian J, Tanner, Jack, Artusio-Glimpse, Alexandra, Williams, Paul A., Parab, Niranjan, Zhao, Cang, Sun, Tao (2022), Asynchronous AM Bench 2022 Challenge Data: Real-time, simultaneous absorptance and high-speed x-ray imaging, Version 1.0.0, National Institute of Standards and Technology,
https://doi.org/10.18434/mds2-2525 (Accessed: [give download date])

----------
References
----------

[1] Simonds, B. J., Tanner, J., Artusio-Glimpse, A., Williams, P. A., Parab, N., Zhao, C., & Sun, T. (2021). The causal relationship between melt pool geometry and energy absorption measured in real time during laser-based manufacturing. Applied Materials Today, 23, 101049. doi:10.1016/j.apmt.2021.101049

[2] Khairallah, S. A., Sun, T., & Simonds, B. J. (2021). Onset of periodic oscillations as a precursor of a transition to pore-generating turbulence in laser melting. Additive Manufacturing Letters, 1, 100002. doi:10.1016/j.addlet.2021.100002

[3] Simonds, B. J., Tanner, J., Artusio-Glimpse, A., Williams, P. A., Parab, N., Zhao, C., & Sun, T. (2020). Simultaneous high-speed x-ray transmission imaging and absolute dynamic absorptance measurements during high-power laser-metal processing. Procedia CIRP, 94, 775-779. doi:10.1016/j.procir.2020.09.135

[4] Simonds, B.J., Sowards, J., Hadler, J., Pfeif, E., Wilthan, B., Tanner, J., Harris, C., Wiliams, P., Lehman, J. (2018). Time-resolved Absorptance and Melt Pool Dynamics during Intense Laser Irradiation of a Metal. Physical Review Applied, 10, 044061. doi:10.1103/PhysRevApplied.10.044061 

[5] Derimow, N., Schwalbach, E.J., Benzing, J.T., Kilgore, J.P., Artusio-Glimpse, A.B., Hrabe, N., Simonds, B.J. (2021). In situ absorption synchrotron measurements, predictive modeling, microstructural analysis, and scanning probe measurements of laser melted Ti-6Al-4V single tracks for additive manufacturing applications. Journal of Alloys and Compounds, 900, 163494, doi:10.1016/j.jallcom.2021.163494

[6] Parab, N.D., Escano, L.I., Fezzaa, K., Everhart, W., Rollett, A.D., Chen, L., Sun, T. (2018). Ultrafast X-ray imaging of laser –metal additive manufacturing processes, Journal of Synchrotron Radiation, 25, 1467–1477, doi: 10.1107/S160 05775180 09554 


-------------
Data Overview
-------------

Files included in this publication:


  2525_README.txt
	This document.  

  Sample Diagram.pdf

    	This file contains information about the experimental setup. Page 1 is a schematic of the experimental setup including beam diameter information. See references [1] and [3] for additional information. Page 2 contains a detailed schematic of the sample and sample mount, including dimensions. 

  Beam Profile_5p5micron pixels_Normalized Integral.csv

    	This is a matrix of beam profile data that was measured at the beam focus. The spatial dimensions of the matrix is 5.5 micrometers. The data are normalized such that the integral of the entire dataset is equal to 1. Therefore, multiplying these data by the laser power gives the correct laser irradiance. An example Python script is given as “Beam Profile Imaging.ipynb” This data was obtained at focal plane of laser. However, please note that all experimental data was obtained with the metal sample surface 2.8 mm below the focal plan. At this location, the calculated laser beam spot size (1/e^2) was 122.5 mircometers. 

  Beam Profile Imaging.ipynb

	This is a script for visualizing the beam profile data given in the "Beam Profile_5p5micron pixels_Normalized Integral.csv". The user should adjust the “LaserPower” parameter to be equal to the value specified for the absorption data. The user also needs to input an appropriate value for the data file location.
  
  NIST SRM 654b Ti64 data sheet.pdf

	Composition datasheet for NIST standard reference material (SRM) 654b, a Ti-6Al-4V alloy. The data presented in this record were obtained from material taken from this SRM feedstock. Additional information found at https://www-s.nist.gov/srmors/view_detail.cfm?srm=654B
 
  Data Synchronization Example.ipynb

	This is a sample script for comparing raw and processed Xray images along with the absorption up to the time that the images were obtained. The time is selected by the user by a choice of image frame to display.

  Absorption_Uncertainty_Analysis.pdf
	
	This document gives detailed information on the uncertainty analysis of the absorption experiment. It includes derivations of the expanded uncertainties of the absolute reflected laser power, absolute absorbed laser power, and the laser coupling efficiency. Also included are sample calculations of these quantities for the data provided in "Spot on Bare Metal_Calibrated Absorption Data.csv."  

  Spot on Bare Metal_Calibrated Absorption Data.csv

	File containing time-dependent absorption and camera timing data. There are nine data columns.
	   1: Numerical index. 
	   2: ‘Time’ – in seconds. 
	   3: ‘InputLaser’ – Temporal profile of processing laser pulse. The sharp peak for first few hundred nanoseconds is an artifact and can be ignored. 
	   4: ‘AbsoluteAbsorption’ – The calculated absolute absorbed power in units of Watts. See references [1] and [4] for additional information on how this is calculated. 
	   5: ‘AbsAbsorptionUncertainty’ – Absolute expanded uncertainty in the absorbed power in units of Watts. See "Absorption_Uncertainty_Analysis.pdf" for how this is calculated. 
	   6: ‘RelativeAbsorption’ – Percent absorption. 
	   7: ‘CameraTrigger’ – Binary data that indicates when the camera was triggered to start capturing images. A value of 1 means the camera has been triggered.  
	   8: ‘FrameTrigger’ – Binary data that indicates when the camera is capturing a frame. A value of 1 means an images is being captured with the leading edge marking the start of frame exposure. Both 'CameraTrigger' and 'FrameTrigger' must be equal to 1 for a frame to be captured. 
	   9: ‘FrameNumber’ – The frame number is given at the start of each frame. These numbers correspond to those of the raw image files and all Xray processed data.

  Spot on Bare Metal_XrayImages_Raw.zip

	Unprocessed (raw) Xray images. The last three digits of the file name indicate the frame number, which is consistent with the column 'FrameNumber' in "Spot on Bare Metal_Calibrated Absorption Data.csv."

  Spot on Bare Metal_XrayImages_processed_captioned.zip

	Xray images that have been processed according the reference [3], which have also been timestamped.   
    
  Spot on Bare Metal_Xray with Absorption.avi

	An example movie made from the frames in “Spot on Bare Metal_XrayImages_processed_wAbsorption”.

  Spot on Bare Metal_XrayImages_processed_wAbsorption.zip

	Image frames from "Spot on Bare Metal_Xray with Absorption.avi" containing processed Xray images along with a time series plot of the measured laser absorption up to the time that the Xray image as captured.

  Scan on Bare Metal_Calibrated Absorption Data.csv

    	File containing time-dependent absorption and camera timing data. There are nine data columns.
	   1: Numerical index. 
	   2: ‘Time’ – in seconds. 
	   3: ‘InputLaser’ – Temporal profile of processing laser pulse. The sharp peak for first few hundred nanoseconds is an artifact and can be ignored. 
	   4: ‘AbsoluteAbsorption’ – The calculated absolute absorbed power in units of Watts. See references [1] and [4] for additional information on how this is calculated. 
	   5: ‘AbsAbsorptionUncertainty’ – Absolute expanded uncertainty in the absorbed power in units of Watts. See "Absorption_Uncertainty_Analysis.pdf" for how this is calculated. 
	   6: ‘RelativeAbsorption’ – Percent absorption. 
	   7: ‘CameraTrigger’ – Binary data that indicates when the camera was triggered to start capturing images. A value of 1 means the camera has been triggered.  
	   8: ‘FrameTrigger’ – Binary data that indicates when the camera is capturing a frame. A value of 1 means an images is being captured with the leading edge marking the start of frame exposure. Both 'CameraTrigger' and 'FrameTrigger' must be equal to 1 for a frame to be captured. 
	   9: ‘FrameNumber’ – The frame number is given at the start of each frame. These numbers correspond to those of the raw image files and all Xray processed data.

  Scan on Bare Metal_XrayImages_Raw.zip

	Unprocessed (raw) Xray images. The last three digits of the file name indicate the frame number, which is consistent with the column 'FrameNumber' in "Spot on Bare Metal_Calibrated Absorption Data.csv."

  Scan on Bare Metal_XrayImages_processed_captioned.zip

	Xray images that have been processed according the reference [3], which have also been time stamped.

  Scan on Bare Metal_Xray with Absorption.avi

	An example movie made from the frames in “Scan on Bare Metal_XrayImages_processed_wAbsorption”.

  Scan on Bare Metal_XrayImages_processed_wAbsorption.zip

	Image frames for "Scan on Bare Metal_Xray with Absorption.avi" containing processed Xray images along with a time series plot of the measured laser absorption up to the time that the Xray image as captured.


---------------
Version History
---------------

1.0.0 
  initial release
1.1.0 (this version)
  Laser beam size correction

--------------------------
Methodological Information
--------------------------

Specific experimental details for the data given in this record can be found in reference [1], with a brief description given here. For general information about the development of integrating sphere radiometry for laser processing, see reference [4]. For general information about the Xray imaging technique used here, please see reference [6].

All data presented here were carried out at the 32-ID-B beamline of Advanced Photon Source at Argonne National Laboratory. Page 1 in "Sample Diagram.pdf" illustrates the experimental configuration. To measure the absorptivity, a custom-made, calibrated integrating sphere was placed over the metal sample which collected the backscattered laser light during laser processing.  The amount of absorbed laser power is determined from an energy balance calculation between the measured input and scattered laser light (zero light transmission). The sphere used in these measurements had an elongated inlet aperture at the top to allow for a scanned process laser. A small channel along the sphere base accommodated clear passage of the 2 mm wide x-ray beam. The process laser had an angle of incidence of 7° relative to the sample surface normal. This ensured that the strong backscattered light from the initially solid surface would be captured within the sphere. The backscattered light was measured by a photodiode that was fiber-coupled to the sphere surface; it was bandpass filtered at 1070 nm. The photodiode voltage was measured by a high-speed oscilloscope with a 40 ns time resolution. A calibration procedure using a well-characterized scattering surface in place of the experimental target produces a coefficient that converts the photodiode signal to an absolute amount of scattered light power. This procedure was repeated between laser injections so that changes to the internal surface of the sphere due to weld spatter could be accounted for. Over the course of these experiments, the loss within the sphere changed only slightly, increasing on average by 0.6 % per laser exposure. The dynamic absorptance is determined by calculating the difference between the backscattered and input laser powers and then dividing by the input laser power. Normalizing this to 100 gives the percent ab- sorption. A complete description of the absorption measurement equations and the uncertainty analysis is give in "Absorption_Uncertainty_Analysis.pdf."

The high-speed Xray imaging experiments were performed similar to those described in reference [6]. All videos in the datasets released here were captured at 50,000 frames per second with a 2.5 μs exposure time. Image analysis was performed according to the procedure in reference [3].

The process laser was a 1070 nm wavelength Yb-doped fiber. The beam profile at the sample surface for all simultaneous X-ray and absorption measurements was calculated to have a diameter (1/e^2 value) of 122.5 μm ± 3.0 μm. This is based on beam profile measurements made at multiple heights, with measured minimum spot size (beam waist) of 49.5 μm ± 5 μm. The beam profile data measured at beam focus are provided (see "Beam Profile_5p5micron pixels_Normalized Integral.csv"), which is not the exact location of the sample surface. The X-ray and absorption measurements were performed with the sample 2.8 mm below the beam waist, resulting in the calculated 122.5 μm spot size. The beam profile is approximately Gaussian over the range that was measured with only a slight ellipticity of greater than 92 %.  A galvanometer scanner that maintains a 7° angle of incidence of the process laser beam to the sample surface is used for controlling the laser duration and scanning. The scanned challenge data ("Scan on Bare Metal_*") had a scan speed of 700 mm/s where a "Sky-write" scan strategy was used. This ensured that the scan speed was constant the entire time the laser was on. The temporal laser pulse information is given with the tabulated absorptivity data and was approximately 2 ms long. The sample and integrating sphere apparatus were placed in an inert gas chamber with two Kapton windows that allowed for transmission of the X-ray beam. This chamber was vacuum pumped and backfilled with argon to room pressure.

The materials used were Ti-6Al-4V samples machined from NIST Standard Reference Material 654b using wire electric discharge machining. These were thinned to approximately 300 μm by polishing of the sides to which the Xrays were incident. This thickness was necessary for adequate Xray transmission at the high frame rate (50,000 frames per second) Xray imaging used here. The polished sides also created a highly specular surface which is necessary for good contrast of keyhole and melt pool in the Xray images. The laser incident surface was also polished to a specular finish. Please see "Sample Diagram.pdf" for a sample diagram and more information.

