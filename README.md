# CmbRecTool

The Library has been written by Fortran90 to reconstruct lensing potential, cosmic bi-refrimgence, and patchy reionization from cosmic microwave background anisotropies (CMB) in full and flat sky. The Library also includes subrotuines for delensing, bi-spectrum calculation, and so on. Installing this library at NERSC is straightforward, but the library also works at other environment. 

Python modules are also available at https://github.com/toshiyan/cmblensplus

# Installation

  0) The code assumes "ifort" (intel fortran compiler) as a fortran compiler. 

  1) Go to "pub" directory, and install each pubclic package (FFTW, Healpix, LAPACK and cfitsio). The clibrary utilizes static links to those packages.

  2) Go to the top directory, and type "./MALEALL.sh install". 
  
  3) You will find the Library at "lib" and "mod" directories. 

  4) To use some example codes, you also need "pylib/makefile.py" to create a Makefile if no Makefle exists. 

# Test

You can test the library by using the example codes at "examples" directory. At each subdirectory, you need to compile each code 
by just typying "make" (you do need to edit each Makefile). Then you can run each code by "./exe" or "./exe params.ini" where exe is 
an executable file. 


# References

In flat sky analysis, the verification of the source code is given in 

  - https://arxiv.org/abs/1209.0091 : the temperature lensing reconstruction, 
  - https://arxiv.org/abs/1310.2372 : the polarization lensing reconstruction, 
  - https://arxiv.org/abs/1612.07855 : the cosmic bi-refringence reconstruction, 
  - https://arxiv.org/abs/1703.00169 : the CMB delensing
  - https://arxiv.org/abs/1706.05133 : the bi-spectrum

For full sky analysis, the library also supports the subroutines for the lensing reconstruction, patchy reionization reconstruction and delensing as shown in 

  - https://arxiv.org/abs/1405.6568 : the temperature/polarization lensing reconstruction, and delensing
  - https://arxiv.org/abs/1711.00058 : the patchy reionization reconstruction
  - https://arxiv.org/abs/1812.10635 : the bi-spectrum

The algorithm written in nldd/xxx.f90 is described in doc/note.pdf


# Acknowledgment

The library software uses the following public codes: FFTW, HEALPix, LAPACK, CFITSIO and also part of CAMB subroutines. 

# Contact

  namikawa at slac.stanford.edu

