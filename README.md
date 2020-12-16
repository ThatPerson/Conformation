# Conformation

Program to generate different conformations of a molecule. Can input a file containing a molecule (currently only XYZ files), assign bonds (given uniform bond length), and then perform rotations about bonds.

Compile as `gcc main.c -lm -o conf -O3`, and run as either `./conf script.r` (to run a script) or `./conf` to open interactive mode.

## Scripting

The following commands are available

- open FILENAME
  Opens xyz file at FILENAME
- bond BL
  Generates bonds between atoms where distance < BL (angstroms)
- rotate A B THETA
  Rotates atoms connected to B about the A-B axis an amount theta (radians)
- output FILENAME
  Writes output as xyz file to FILENAME.
  
## Interactive

- open FILENAME
   opens xyz file FILENAME
- bond BL
   sets up bonds between all atoms with distance < BL angstroms
- graph N
   prints out graph starting at atom N
- print
   prints out entire system
- rotate A B THETA
   applies a rotation to atoms connected to B about the A-B axis
   an amount theta radians.
- output FILENAME
   outputs as XYZ file into FILENAME
- run SCRIPT
   runs script in SCRIPT.
- exit
   exits program
   
## OpenBabel

`obabel -:"<smiles>" -oxyz -O<filename.xyz> --gen3d`

Can be used to generate xyz files to test this with.


