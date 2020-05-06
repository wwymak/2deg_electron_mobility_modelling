***********README.TXT************************************************************
Contact Kanti for details/bugs/suggestions!
<kdasgupta@phy.iitb.ac.in>
<kd241@cam.ac.uk> 
*********************************************************************************

This file is a simple documentation of the program that calculates 2DEGmobility at 
given densities using the wafer paramters provided in the input file.

The source files are all written in "C" and do not use any external libraries like 
NAG, IMSL etc. The numerical integration subroutines are all part of the program,
they are slightly modified versions of "Numerical recipes in C". 

I have not yet devoted much time to error checking....for example there are very
minimal (or none) checks against inconsistent units/mistyped/absurd input values!!!!


--------------------------------------------------------------------------------
1. IF YOU ARE USING A COMPILED VERSION (e.g. a windows/dos ".exe")
--------------------------------------------------------------------------------

The program reads input instructions from the file "input_for_mobility.txt". 
This file is organised in a section wise manner, with latex-like sectioning.
Anything not included in \begin{SECTIONNAME} & \end{SECTIONNAME} will be 
ignored. They are effectively comment lines where you can write your notes
and comments. For example I use the same file to store the best fit parameters
that I have got for a particular wafer.

If you need to include comments within a section then start the 
line with "#". 

e.g. 

\begin{OutputSpecifications}  # This section must exist
#outputunits=SI
outputunits=cm2byVs
#outputfile=OUTPUT_110nm_HET[V656].dat
outputfile=OUTPUT_compare_HET_QW.dat
\end{OutputSpecifications}


to retain the # for some reason use '\#', e.g. "filename_with_a\#.dat"
The GaAs-AlGaAs heterointerface is at z=0 (for heterostructure)
The QW is between (z=0,z=L); the impurity profiles should be set accordingly 
The surface will have negative co-ordinates - this is a convention carried over
from what Ando-Fowler-Stern used. 

To run the program you only need to type the name of the executable "2DEGmobility.exe",
it takes no arguments, all the input will come from the file "input_for_mobility.txt".


The names of the sections are fairly self explanatory and the units are mentioned in comment 
lines. In some cases you have a choice of units. For example the mobility itself may be desired 
in cm^2/Vs or SI - depending on how you want/like to plot them or how the experimental data 
points are organised.

If you omit a section then the scattering rate corresponding to that process will not be included. 
The outputfile will consist of columns which give the mobility limit imposed by that process. 
The final column is the total scattering rate calculated using "Mathiessen's rule".

Notice that the program can handle umultiple layers of delta doping or modulation doping or
impurity profile (background doping). For two delta doped layers you can use:

\begin{DeltaDoping}            # This section is optional 
# units : distances in nanometers & charge in per mt^2
-50.00   7e14
-70.00   1e15
\end{DeltaDoping}
 
This means that at 50nm below the surface there is a delta-doped layer with charge density of 
7e14/m^2 and at 70 nm depth there is another layer. There is no limit to the number of layers
you can put in. 

Similarly for two slabs of modulation doping you use

\begin{ModulationDoping}       # This section is optional 
#units : distances in nanometers & background in per mt^3
-80.0  -40    1.0e23
60   100      1.0e23
\end{ModulationDoping}

there can be as many slabs as one needs. 

The background/impurity concentration is treated exactly as modulation 
doping, the only difference being that the typical concentrations would 
be 2-3 orders of magnitude lower.

e.g.
\begin{ImpurityProfile}
#units : distances in nanometers & background in per mt^3
-310.0  -300.0   4.6e19
-300.0     0.1   4.6e19
  0.1    25.0   3.0e19
  0.1   400.0   4.6e19
\end{ImpurityProfile}
 
These could be appropriate for a quantum well 25nm in width, 300nm below the surface.


The interface roughness is described by an amplitude and correlation length, the 
distribution is assumed to be gaussian. 
e.g.

\begin{InterfaceRoughness}
roughness_amplitude_in_nanometer = 0.1;
roughness_correl_length_in_nanometer=15;
N_depletion_charge_density=1e13; 
\end{InterfaceRoughness}
 
The N_depletion_charge_density is almost irrelevant for GaAs/AlGaAs, it is carried 
over from the original formulation for MOSFETs, where it can be quite significant.  




-------------------------------------------------------------------------------------
2. IF YOU ARE COMPILING THE FILES
------------------------------------------------------------------------------------

You need 4 source files:

2degmobility.c
inputfile_read_subroutines.c
2degchi_adaptive_subroutine.c
formfactor2deg_subroutine.c

For windows/DOS compilers...
----------------------------
1. Open a project, name it whatever you like.
2. Add the files to the project.
3. Depending on the compiler you are using, "compile" and then do "link" 
   OR use "build all" to create the executable directly.

For Unix/Linux
---------------
gcc 2degmobility.c inputfile_read_subroutines.c 2DEGchi_adaptive_subroutine.c formfactor2deg_subroutine.c -lm -o OUTPUT_FILENAME

 
You can also do it using a Makefile, but then it is a small enough project...
Remember to use "-lm"  option to include the math libraries.

For Mac
---------

I have never used one.....somebody can fill in this section

If you use the prog...you can add your comments/updates/bug-reports......

Kanti
Jan 2011
