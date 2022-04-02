# Source Code Updates

## Modifications to SWAT 
We make updates to the reservoir module of the Soil & Water Assessment Tool by modifying the following files:
* `res.f`: add a 6th option for reservoir operations. The new option calculates daily optimal releases at all system reservoirs based on system state information and optimized radial basis function parameters.
* `modparm.f`: add necessary new reservoir variables 

And adding the following file:
* `interp_lin.f`: code for linear interpolation. Used to calculate head based on reservoir volume. 

## Compiling Source Code for Linux
This is a fork of [an existing repo](https://github.com/WatershedModels/SWAT).

Run `compile.sh` to compile the executable with the updated source code. 
