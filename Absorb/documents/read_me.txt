
ABSORB code  

INFORMATION

This folder \documents contains the most recent copies of manuals, and published papers.

The folder \src contains the fortran code for the public version of the standalone Absorb program.
It should not be publically released, but can be copied for use in the cfml

There are also private, commercial versions with confidential features developed for specific clients. 
The private features are contained in fortran code files that are called, for example, company.for. 
They contain the same subroutines as those in the file normal.for in the public version.

The private commercial versions are not in this folder and will not be released.


CODE

The entire processing of a crystal model and a dataset is controlled by a series of subroutine 
calls from absorb70_main.for . The methods and execution are also described in the programmers' manual.

24/04/2019 Ross Angel