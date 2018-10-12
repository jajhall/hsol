# hsol
Original high performance simplex solver of Qi Huangfu

Welcome to HSOL 1.0
===================

HSOL is a high performance serial and parallel solver for large scale
sparse linear programming (LP) problems. It is written in C++ with
OpenMP directives. It is based on the dual revised simplex method,
exploiting parallelism using either "parallel minor iterations" (PAMI)
or "single iteration parallelism" (SIP). A full technical reference to
PAMI and SIP is

Parallelizing the dual revised simplex method Q. Huangfu and J. A. J. Hall
Mathematical Programming Computation, 10 (1), 119-142, 2018.
DOI: 10.1007/s12532-017-0130-5

http://www.maths.ed.ac.uk/hall/HuHa13/

HSOL has been developed and tested on various linux installations
using both the GNU (g++) and Intel (icc) C++ compilers.

Download
--------

Unzip the file hsol_1.0 into your folder of choice

Compilation
-----------

A serial executable is obtained by using your favourite C++ compiler
to compile and link the downloaded .cpp files. For example, using g++

g++ *.cpp

The solver has also been tested successfully with the Intel "icc"
compiler and with various levels of optimization (for g++ and icc),
and -O2 is recommended. Since they exist as comments, the OpenMP
directives are ignored by default.

A parallel executable is obtained invoking the OpenMP directives. This
is achieved using GNU g++ thus

g++ -fopenmp *.cpp

and with icc  thus

icc -openmp *.cpp

Run-time options
----------------

In the following discussion, the name of the executable file generated
is assumed to be "hsol"

HSOL can only read plain text MPS files, and the following command
solves the model in ml.mps

hsol ml.mps

Run-time options -pami and -sip direct hsol to use PAMI or SIP. 

When compiled with the OpenMP directives invoked, the number of
threads used at run time is the value of the environment variable
OMP_NUM_THREADS. For example, to use HSOL with PAMI and eight
threads to solve ml.mps execute

export OMP_NUM_THREADS=8
hsol -pami ml.mps

If OMP_NUM_THREADS is not set, either because it has not been set or
due to executing the command

unset OMP_NUM_THREADS

then all available threads will be used.

Observations
------------

When compiled without the OpenMP directives, or if run with
OMP_NUM_THREADS=1, hsol is serial. The -sip run-time option will not
affect performance. The -pami run-time option will cause hsol to use
serial minor iterations and, although this could lead to better
performance on some problems, performance will typically be
diminished.

When compiled with the OpenMP directives and OMP_NUM_THREADS>1 or
unset, hsol will use multiple threads if the -pami or -sip run-time
option is specified. If OMP_NUM_THREADS is unset, hsol will try to use
all available threads so performance may be very slow. Although the
best value will be problem and architecture dependent,
OMP_NUM_THREADS=8 is typically a good choice. Although hsol is slower
when run in parallel than in serial for some problems, it is typically
faster, with the -pami option usually faster than the -sip option.
