The current approach is specific to the Omo River Basin. However, this approach could be generalized to any context. 

The following information is hard-coded in the res.f file, but could be transferred to the .res I/O input files:
* The number of RBFs
* The number of inputs 
* The number of outputs 


More steps:
* Define flexible framework for pulling in information for all reservoirs at once (e.g., volume, evol, pvol).
* Create a standard input file for 
  * RBF parameters (in this example, called ra.txt)
  * volume-storage-level curves for each reservoir 
