Backup made before reducing connor316AMR to the small one-pulse test case.

Contains the previous case setup:
- constant/: laser parameters, time-vs-position/power, dynamic mesh, transport
- system/: original blockMeshDict, setFieldsDict, controlDict, fvSchemes/fvSolution
- initial/: initial fields copied by Allrun

Restore example from connor316AMR:
cp -a backups/pre_small_one_pulse/constant .
cp -a backups/pre_small_one_pulse/system .
cp -a backups/pre_small_one_pulse/initial .
