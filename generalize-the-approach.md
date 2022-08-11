The current approach is specific to the Omo River Basin. However, this approach could be generalized to any context. 

The following information is hard-coded in the res.f file, but could be transferred to the .res input file so users would not have to edit the source code to change the information:
* The number of RBFs
* The number of inputs 
* The number of outputs 


Additional required steps for generalizability:
* The process requires start of day reservoir storages for all reservoirs in the system. The original source code deals with each reservoir individually. In this study, we manually pull this information for all the reservoirs when dealing with the first reservoir, so we can later use that information when calculating reservoir releases for the remaining reservoirs. This could be done more flexibly by defining an array of the reservoir numbers and then pulling the reqired info (e.g., volume, evol, pvol). 
* Create a standard input file for 
  * RBF parameters (in this example, called ra.txt)
  * volume-storage-level curves for each reservoir 
