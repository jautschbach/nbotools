convert-47-to-namelists.py is a python script that converts an NBO 'file47'
input to Fortran namelist input, to be read by the provided parse-namelist
program. The latter code will check basis set conventions in file47 and 
attempt to recreate the overlap matrix. If the parse-namelist output prints 
warnings about the overlap matrix not agreeing with what's in file47, it 
is likely not safe to use NBOTools or nbo2cube 
(https://github.com/jautschbach/nbo2cube).

Basis set conventions in file47 data created by Gaussian, Orca, 
Molden2AIM (https://github.com/zorkzou/Molden2AIM) and
and molcasto47 (https://github.com/jautschbach/molcasto47) have been
confirmed to work. 
