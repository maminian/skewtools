**README last updated: 15 August 2017**

## Introduction.

![preview](web_images/pipevid0180.png)

This is a collection of codes and scripts for 
analyzing a passive tracer in 
laminar pipe flow. The code was produced 
as a part of my PhD dissertation work, 
but it is in the process of being generalized 
enough to handle a variety of different situations. 

A passive tracer is a solute, 
pollutant, field of particles, or otherwise, 
which is carried by a fluid flow, 
but doesn't influence the actual flow behavior. 
While in principle this never occurs, it is 
a very reasonable assumption for many situations, 
as long as the bulk density or temperature 
of the fluid doesn't change noticeably (which 
would affect the flow) and the particles are 
sufficiently small relative to the length scales 
of interest.

This project is primarily interested in the 
behavior of the tracer itself; 
how its statistics (mean, median, skewness) 
change in time, and how they depend on the 
type of pipe flow they are carried by. In this 
case, the laminar flow field is found 
beforehand, then simulations are run by 
releasing a large number of particles, and 
allowing them to diffuse and be carried by the fluid.

The simulations are done using a Monte Carlo method, 
which essentially drops many (on the order of millions to 
billions) of particles, allowing them to diffuse and 
be carried by the flow, and studies their 
collective statistics throughout.

This code was written mainly for my own use, so 
work needs to be done to make it more modular,
readable, etc. As progress gets made in this 
area, fewer "core" requirements (as described below) 
will be necessary to compile a minimal 
working example and visualize it.

Most (if not all) of these prerequisites should 
be installed with your operating system's package 
manager (for example, apt-get in Ubuntu). 
For the python packages, I recommend using 
the "pip" package manager. An example setup 
is described below.

I have not tried installing/running in 
OSX or Windows, so you will need to figure 
that out on your own. If you do, I would 
appreciate any feedback or notes for 
your installation process.


## Prerequisites for compiling and running the executables.


The simulation code is written entirely in Fortran 90, 
and i/o is handled with HDF5, a data container which 
allows mixed output (doubles, integers, strings, etc) 
in a single file with a directory structure. This is 
both compact (i.e., one file per simulation) and 
being self-documenting (data in the file
comes with plain-text descriptions).

For all this, you need:

 - A modern fortran compiler (e.g., gfortran or ifort)
 - An installation of HDF5. This must include the binary "h5fc", 
   which is a wrapper compiler. (hdf5-tools package in ubuntu)
 - An installation of LAPACK that you can link. (liblapack-dev in ubuntu)
   In particular, dgemv.f is used, if you want to just copy the 
   source file (can be found online) and link it properly in the makefile.
 - Python2.7 is recommended to automate the process of running 
   the executable multiple times and/or on a cluster.


## Prerequisites for analyzing the output of a simulation.


In principle, here, you can use any tool that can 
import HDF files, as this is the output of the simulation.
This is done in python2.7 in this project. Again, 
I recommend installing python-pip (apt-get install python-pip) 
then installing these via "pip install <package>".

 - h5py
 - numpy
 - matplotlib
 - scipy, only necessary for some older scripts
 - ipython (recommended)


## Compiling and running a test example.


Currently, the simulation executables are broken up by 
the class of cross-sectional geometry, and each has its own 
separate make command:

 - channel (make channel_mc)
 - rectangular duct (make duct_mc)
 - elliptical pipe (make ellipse_mc) 
 - equilateral triangle (make triangle_mc)
 - "racetrack" pipe (make racetrack_mc)

After compiling, modify the parameter file parameters_mc.txt 
with the problem parameters you wish. It is then possible 
to start a simulation via a command like

     ./channel_mc parameters_mc.txt output.h5

where the third argument is the name of the file where 
the simulation data will be written (note, you need 
to change the RNG seed in parameters_mc.txt to an 
integer if you go this route).

**However**, it is highly recommended you instead 
use "batch_submit.py" to call one of the executables. 
This allows for a few things:

1. Automatically running multiple simulations; i.e., 
   running the same code multiple times with different 
   seeds for the RNG
2. Running the code on one of UNC's clusters:
  * "bsub" assumes an LSF system (on the kure/killdevil clusters)
  * "longleaf" assumes a SLURM system (on the longleaf cluster)
  * "local" assumes you're running it on your own computer

  The cluster approach submits N separate jobs, with different 
  initializations of the seed

3. Automatically creating folders and naming simulation files.

Only the first few lines of batch_submit.py which contain 
parameter values (filenames, switches, etc) should be changed.

Once this is set up, start the simulation in the terminal with 

```
python batch_submit.py
```

If you run the job locally, you should see something like this:

>
>================================================================================
> 
>Geometry: channel
> 
>                             Uniform initial data.
> 
>                                          Peclet:  1.000E+01
> 
>                             Number of particles: 1000000
>                                   Time interval: ( 0.000E+00 ,  1.000E+00 )
>                   Number of requested timesteps: 1001
>                    Number of internal timesteps: 1002
>                       Largest internal timestep:  1.000E-03
> 
>                 Mersenne Twister seed: 7270442
> 
>================================================================================
> 
>Simulating...   6%   Time remaining:   2.8min
>


This provides some information about the simulation being 
run, and is helpful to verify that the parameters you specified 
were processed correctly. A crude progress meter has been implemented
and will update dynamically until the simulation is completed.


## Examining the output of a test example.

Let's assume the simulation finished successfully, and the output 
is written to output.h5. The scripts/ folder contains a large number of 
scripts for visualizing the moments of the statistics. 
These are mostly "one-purpose-only" scripts, unfortunately, and require 
editing depending on your needs. 

The main exception to this is the skewtools.py file, which has a large number of 
functions which are used commonly for our purposes. These 
include things such as collecting a common subset of data from a 
large number of data files, averaging data, interpolation and comparison 
of different numerical results, along with some miscellaneous 
analytical formulas necessary for validating the numerics. 

Here is a simple python example to look at the skewness of the 
cross-sectionally averaged distribution over time:

```python
import scripts.skewtools as st
from matplotlib import pyplot

t,sk = st.importDatasets('output.h5','Time','Avgd_Skewness')

pyplot.plot(t,sk)
pyplot.show()
```

You should see a plot window of the skewness as a function of time.

Depending on the values specified in the input parameter file, 
different statistics are calculated, and the output file output.h5 
will contain different arrays. One way to look at all the arrays 
contained in an hdf file while in python is to do:

```python
import h5py

of = h5py.File('output.h5','r')
of.keys()
```

It is also possible to use the command-line tool h5dump that should come with the hdf5-tools package. For example, running this at the terminal:

```
h5dump -H output.h5 
```

gives a (relatively cryptic) output of all the arrays stored in the file.


## Contact.


If you have questions, don't hesitate to contact me via 
email at manuchehr.aminian@gmail.com. Further tutorials to 
demonstrate the capabilities of the code are planned for the future.

## Relevant research publications
- Aminian, M., Bernardi, F., Camassa, R., Harris, D., McLaughlin, R. [“How boundaries shape chemical delivery in microfluidics”](http://science.sciencemag.org/content/354/6317/1252), Science, Vol. 354, Iss. 6317, 2016, pp. 1252-1256.
- Aminian, M., Bernardi, F., Camassa, R., McLaughlin, R. [“Squaring the circle: Geometric skewness and symmetry breaking for passive scalar transport in ducts and pipes”](https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.115.154503), Phys. Rev. Lett., Vol. 115, Iss. 15, 2015.

