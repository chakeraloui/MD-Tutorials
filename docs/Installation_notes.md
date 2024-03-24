# install GROMACS on Linux
1. Update Your Packages
```bash
sudo apt update 
sudo apt upgrade
```
 2. Install Prerequisites
Let’s start by installing cmake.
```bash
sudo apt install cmake
```
Now we can proceed to install build-essential:
```bash
sudo apt install build-essential
```
3. Download and Unpack your Desired GROMACS Version
the first thing you should do is download the tarball of this specific version of GROMACS
```bash
wget https://ftp.gromacs.org/gromacs/gromacs-2022.3.tar.gz
tar xvf gromacs-2022.3.tar.gz
```
4. Now create a build directory and move into it.

```bash
cd gromacs-2022.3
mkdir build
cd build
```
5. Now, it’s time to compile the source code and install GROMACS:
```bash
cmake .. -DGMX_BUILD_OWN_FFTW=ON -DREGRESSION_DOWNLOAD=ON
make 
make check 
sudo make install
```
Note

 This is the stage where one can customize its GROMACS installation as needed. For instance, if you want to get the mpi version of GROMACS you can simply add this flag -DGMX_MPI=on.

By default GROMACS will be installed in /usr/local/gromacs. If you want to install it in a specific directory you can add this to the cmake command -DCMAKE_INSTALL_PREFIX=/path/to/dir/.

```bash
cmake .. -DGMX_BUILD_OWN_FFTW=ON -DREGRESSION_DOWNLOAD=ON -DGMX_MPI=on -DCMAKE_INSTALL_PREFIX=/path/to/dir
make 
make check 
sudo make install
```
For more details on how to customize your installation you can refer to the [official guide](https://manual.gromacs.org/current/install-guide/index.html).

 Activate GROMACS

Activate GROMACS by using the following command every time you open a new terminal session:
```bash
source /usr/local/gromacs/bin/GMXRC
```
Please note that if you opted to install it in a different directory, make sure to source from that specific location.




## Enabling Optimizations for Your Simulation

The optimizations in this section are not yet enabled by default. The new code paths have been verified by the standard GROMACS regression tests, but still lack substantial “real-world” testing. 

GROMACS should be built using its internal threadMPI library instead of any external MPI library. 

for example: to build  Gromacs 2024 on linux ubuntu eauipped whith two rtx 3090 GPU we use 

```bash
cmake .. -DGMX_BUILD_OWN_FFTW=ON -DREGRESSIONTEST_DOWNLOAD=ON -DGMX_GPU=CUDA -DGMX_USE_OPENCL=off -DGMX_CUDA_TARGET_SM=75
```
 At runtime, the optimizations can be fully enabled by setting the following three environment variables to any non-NULL value in your shell (shown for bash shell here).

For halo exchange communications between PP tasks, use the following command:

 ```bash
export GMX_GPU_DD_COMMS=true

```
For communications between PME and PP tasks, use the following command:

```bash
export GMX_GPU_PME_PP_COMMS=true

```
To enable the update and constraints part of the timestep for multi-GPU:

```
export GMX_FORCE_UPDATE_DEFAULT_GPU=true

```
The combination of these settings triggers all optimizations, including dependencies such as GPU-acceleration of buffer operations.

When running on a single GPU, only GMX_FORCE_UPDATE_DEFAULT_GPU is required (where for single-GPU only, this can alternatively be enabled by adding the -update gpu option to the mdrun command).

In both single and multi GPU cases, it is also necessary to assign the three classes of force calculations to the GPU through the following options to mdrun:

```
-nb gpu -bonded gpu -pme gpu

```
On multi-GPU, the -npme 1 option is also required to limit PME to a single GPU. 

The new GPU update and constraints code path is only supported in combination with domain decomposition of the PP tasks across multiple GPUs when update groups are used. This means constraining all bonds is not supported, except for small molecules.

To facilitate this, we recommend that you convert bonds with hydrogen bonds to constraints, rather than all bonds. To do this, in the .mdp input file to the GROMACS preprocessor grompp, change the following line:
```bash
constraints = all-bonds
```
The new line should read as follows:
```bash
constraints = h-bonds
```
For more information, see the mdp options section of the User Guide.

Here are the full set of mdrun options that we used when running the 4-GPU performance comparisons:
```bash
gmx mdrun -v -nsteps 100000 -resetstep 90000 -noconfout \
-ntmpi 4 -ntomp 10 \
-nb gpu -bonded gpu -pme gpu -npme 1 \
-nstlist 400
```
The first line instructs GROMACS to run 100,000 steps for this relatively short benchmark test. The timing counters are reset at step 90,000 to avoid the initialization costs being included in the timing, as these are not typically important for long production runs. These specific values have been chosen to allow GROMACS time enough to perform PME tuning automatically, where the size of the PME grid is tuned to give an optimal load balance between PP and PME tasks. Similarly, the -noconfout option instructs GROMACS not to write output files at the end of this short benchmarking run, to avoid artificially high I/O overheads.

The second line specifies that four thread-MPI tasks should be used (one per GPU), with 10 OpenMP threads per thread-MPI task. A total of 40 OpenMP threads are in use to match the number of physical CPU cores in the server.
The third line offloads all force calculations to the GPU, as described earlier.

The final -nstlist 400 option instructs GROMACS to update the neighbor list with a frequency of 400 steps. This option can be adapted without any loss of accuracy when using the Verlet scheme. We found this value to give the best performance through experimentation.

for more information visit this [tutorial](https://developer.nvidia.com/blog/creating-faster-molecular-dynamics-simulations-with-gromacs-2020/)
 
# Download and install Chimera 
1. Here's what you need to do to get chimera working on your Windows system: 

1.1.If you want Chimera to be used by multiple users:

  Install Chimera as an administrator.

1.2. Download [Chimera release](https://www.cgl.ucsf.edu/chimera/download.html);

We recommend that you use the latest production release.
 A small minority of browsers may download the file as chimera-get.py. If yours does, rename the file to chimera-installer.exe and then run the downloaded executable file. It should install everything you need. 
1.3. Run the chimera installer:

If you choose to save the chimera installer, rather than run it immediately, run it now. 
    
1.4. Double check that shortcuts were created:
If requested, after the installation there should be an additional shortcut for chimera on the desktop. There should also be a UCSF Chimera menu in the Start\Programs menu with various shortcuts.


2. Quick Installation instructuctions for Ubuntu. 
 
 
 Download [Chimera release](https://www.cgl.ucsf.edu/chimera/download.html);
 
  We need to make this file executable. In the directory with downloaded file type this:

```bash
 chmod +x chimera-1.13-linux_x86_64.bin 
./chimera-1.17.3-linux_x86_64.bin 
```

When asked for installation location delete what is written and type this: 

```
~/chimera/
 yes
 1
```

That's all. 



# Download and install FoldX

Download statically linked binary
go to [http://foldx.crg.es](http://foldx.crg.es)
login
select a file from the download section

Copy the binary to a system-wide location, for example:
```bash
sudo cp foldx_2.5.2.linux /usr/local/bin/foldx
```
Make sure it is executable
```bash
sudo chmod 755 /usr/local/bin/foldx
```
Test

First test that you can start Fold-X itself (sometimes there are license-issues):
 ```bash
foldx
```

# Download and install Gnuplot

from rpm or debian

install the gnuplot package with your favorite package manager

from source:

 download from http://www.gnuplot.info/download.html

unpack the tarball and then build, test, and install it:
```bash
cd gnuplot-4.2.0 ; ./configure ; make
make check
make install
```

# Download and install xmgrace
from source: 
from rpm or debian

install the xmgrace package with your favorite package manager

# Download and install PROPKA3
The easiest way to install a release of PROPKA3 is from the PyPI archive with the command
```bash
pip install --upgrade propka
```
   