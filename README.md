 #                                  MDSimulation Tutorial




# Hardware/Software requirements for the tutorial

- REQUIRED: Access to a Linux-like machine, ideally Ubuntu OS;
- REQUIRED: Gromacs v. 5.1.4, [intallations instructions](http://www.gromacs.org/Downloads/Installation_Instructions);
- REQUIRED: VMD v. 1.9.4, [intallations instructions](https://www.biostars.org/p/196147/);
- OPTIONAL: Chimera v. 1.13, [intallations instructions](http://www.cgl.ucsf.edu/chimera/download.html);
- OPTIONAL: Anaconda Python, [intallations instructions](https://docs.anaconda.com/anaconda/install/);
- OPTIONAL: FoldX, [download here](http://foldxsuite.crg.eu/academic-license-info)
- OPTIONAL: Propka, [installation instructions](https://github.com/jensengroup/propka-3.1)
...



#  How to build and run GROMACS

GROMAX installation instructions [here](http://www.gromacs.org/Documentation/Installation_Instructions_5.0). At first, download latest version of GROMAX.
Installation requires administration rights. This is the way of local installation:

```bash
tar xfz gromacs-5.1.4.tar.gz
cd gromacs-5.1.4
mkdir build
cd build
cmake .. -DGMX_BUILD_OWN_FFTW=ON -DCMAKE_INSTALL_PREFIX=/home/usr/dir_where_gromacx_will_be
make
make check
make install
source /usr/local/gromacs/bin/GMXRX 
``` 

(In my case - source /home/aloui/Soft/bin/GMXRC)

Everytime you will use GROMAX, you need to put the last string (with your address)

## Enabling Optimizations for Your Simulation

The optimizations in this post are not yet enabled by default. The new code paths have been verified by the standard GROMACS regression tests, but still lack substantial “real-world” testing.

GROMACS should be built using its internal threadMPI library instead of any external MPI library. For more information, see the Multi-GPU Optimization section earlier.

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

```bash
export GMX_FORCE_UPDATE_DEFAULT_GPU=true

```
The combination of these settings triggers all optimizations, including dependencies such as GPU-acceleration of buffer operations.

When running on a single GPU, only GMX_FORCE_UPDATE_DEFAULT_GPU is required (where for single-GPU only, this can alternatively be enabled by adding the -update gpu option to the mdrun command).

In both single and multi GPU cases, it is also necessary to assign the three classes of force calculations to the GPU through the following options to mdrun:

```bash
-nb gpu -bonded gpu -pme gpu

```
On multi-GPU, the -npme 1 option is also required to limit PME to a single GPU. 

Preparing a system for GROMACS

 this process could apply to a general system where you start from a .pdb file and end up with a set of input files, ready to run a GROMACS simulation.



# 1.System preparation

The initial phase of an MD (Molecular Dynamics) simulation involves the preparation of the system for analysis. This entails establishing an initial configuration of the system before commencing the experiment. This decision must be made meticulously as it significantly impacts the outcome of the simulation. Beginning with an inaccurate structure could result in substantial errors in subsequent stages.

The structure can be derived from experimental or computational techniques, or a combination of both. Ideally, the goal is to obtain a starting structure that closely resembles equilibrium configurations. This facilitates the reduction of the time required to stabilize the system during the equilibration phase.

Once the initial configuration is obtained, several additional steps are necessary to refine the molecular system. This editing process encompasses various actions:

    Defining a simulation box: Initially, a virtual box containing the molecule is established. Periodic boundary conditions (PBC) are then implemented to simulate bulk properties. This involves generating replicas of the box in all directions, enabling the simulation of a large system from a smaller box.


## 1.1.download PDB File

In this section we will describe how to choose and prepare the input files for a system before to running a Gromacs simulation.

As an example we will be looking at 5pep. This is an enzyme used in digestion. The PDB entry can be found here and the PDB file can be downloaded from the protein data bank using:


```bash
wget http://files.rcsb.org/download/5pep.pdb
```
NOTE: The file can contain additional molecules that should be removed before proceeding further.



## 1.2.Validation of the system

Before applying the Preparation step in GROMACS, we first need to validate and preprocess the system.


The following step involves assessing whether adjustments are necessary for water and ion conditions.

Initially, we'll retain the bound water molecules. Subsequently, you'll need to inspect the ions using either the VMD program or by examining the .pdb file. Once you've identified the ions, you must decide whether to keep or delete them. In the selected structure, only Mn and Cl ions are present. Since Mn is a divalent ion, it should be deleted due to the reasons explained earlier. On the other hand, Cl, being monovalent and situated in regions of high positive potential, will be retained.

Following this, it's essential to examine the conformation of side chains.

Basic information regarding side chain conformation is provided here.

The determination of side chain conformation isn't always accurate in the original .pdb structure. For instance, discerning the location of nitrogen or oxygen atoms at the end of asparagine and glutamine is typically challenging due to their similar electron densities.

The residues that require examination include asparagine, glutamine, and histidine. Asn and Gln residues hold significant structural importance in proteins due to their side-chain amide groups, which can act as both hydrogen bond acceptors and donors, thereby stabilizing the protein structure. Moreover, these residues, along with His, are frequently found on protein surfaces or within enzyme active sites, where the hydrogen bonds they form play crucial roles in stabilizing protein-protein or protein-substrate interactions.

To assess side chain conformation, we'll employ the FoldX program. FoldX utilizes an empirical force field for protein design, enabling the determination of the energetic impact of point mutations and the interaction energy of protein complexes (including Protein-DNA).

We'll utilize the RepairPDB command, which identifies residues with unfavorable torsion angles, VanderWaals' clashes, or high total energy, and rectifies them. Initially, it identifies all Asn, Gln, and His residues and rotates them by 180 degrees to prevent incorrect rotamer assignment in the structure. This rotation is necessary because the electron density of Asn and Gln carboxamide groups is nearly symmetrical, and determining the correct placement requires calculating interactions with surrounding atoms. The same principle applies to His residues.
```bash
FoldX --command=RepairPDB --pdb=RP.pdb
```
 FoldX uses output-file as a tag to label different outputs from different commands in batch runs. After running RepairPDB you'll get two files to look at. Given pdb=PDB.pdb the output files are:

    PDB_Repair.fxout -> energies of the repaired residues
    PDB_Repair.pdb -> repaired Pdb

## 1.3.Assessing the protonation states of the protein and preparing the protein 

There are several methods to assess the protonation state of the residues in a protein (for further information, here). For example, the H++ server and PROPKA software can help you in the future. We are only going to use PROPKA .
Install propka 
```bash
pip install propka
```

After installation get back to your working directory and type this in terminal:

```bash
propka3 protein.pdb 
```

PROPKA returns the calculated pKa values for tritable residues in the protein. These are aminoacids, whose sidechain contains a chemical moiety that can change its protonation state depending on the pH of the environment they are found in. You have acidic residues (typically ASP and GLU, but in extreme environments TYR, SER and CYS can act as acids) and basic residues (LYS, ARG and HIS).

Protein forcefield encode these protonation states as different residue names (i.e. GLU and GLH for the two states of glutamic acid). You can find more information about these residues, and how its pKa depends on the chemical environment here .

When preparing a protein, we have to compare each one of the predicted pKa values with the chosen pH: if the pKa is below the chosen pH the residue should be deprotonated (that is acidic residues will be negatively charged and basic residues will be neutral), if above it will be protonated (acidic residues will be neutral and basic residues will be positively charged).

    TIP:

    This step is especially relevant when you are simulating catalytic cycles, where the reaction mechanism often involves proton transfer events. Please check the literature on your target protein before you start to run your simulations.


## 1.4.Removing alternate locations if existed

If the system contains atoms with alternate locations, then they need to be removed. The  system used in this tutorial does not have any alternate locations so you do not need to apply this step.



These are "alternative locations", meaning that in high resolution structure, you may observe different conformations in the electron density. Depending on the programs you are using, you either need to set some parameter to tell the program to ignore all but the most highly occupied conformation, or you need to pre-process the PDB files to remove the secondary locations. Many structure visualisation programs contain options to select, and by extension, selectively delete alternative locations - check the manual for the program you are using
* In VMD, you can use the PDB plugin to do so. https://www.ks.uiuc.edu/Research/vmd/plugins/molfile/pdbplugin.html
* You can also use PDBtools in Phenix https://www.phenix-online.org/documentation/reference/pdbtools.html
* In Rosetta, a simple python script can be used to clean up a pdb file:
tools/protein_tools/scripts/clean_pdb.py   - Prepare PDBs for Rosetta by cleaning and renumbering residues.
https://www.rosettacommons.org/docs/latest/application_documentation/tools/Tools
* In PyMOL, you can use the removal.py script: https://pymolwiki.org/index.php/Removealt or by simply specifying
remove not (alt ' '+'A')
alter all, alt=' '




## 1.5.Selecting the Appropriate Force Field

Next, we must select the force field to simulate the system.

It's crucial to understand that the choice of force field can profoundly impact the simulation's outcome. We recommend consulting publications on various force fields to determine the most suitable one for your system.

For this tutorial, we will utilize the CHARMM  force field. To do this, we'll visit the [MacKerell lab website](http://mackerell.umaryland.edu/charmm_ff.shtml#gromacs), a reliable source for up-to-date CHARMM force field files that are compatible with GROMACS.


Download this file:
> charmm36-jul2022.ff.tgz

or just type:
```bash
wget http://mackerell.umaryland.edu/download.php?filename=CHARMM_ff_params_files/charmm36_ljpme-jul2022.ff.tgz
```
Once the Force Field software is downloaded, extract its contents into the current working directory. This completes the installation process.  



## 1.6.Removing Water Molecules

Our next step involves removing the existing water molecules from the system.

It's important to recognize that this action may not be suitable for all cases. Refer to the "How to delete existing crystal waters outside of the active site" section for guidance on selectively deleting water molecules while preserving those essential for system functionality, such as those tightly bound or functional in the active site.
```bash
grep -v 'HOH' 5pep.pdb > 5pep_protein.pdb
```
alternative method using PYMOL
```pymol

fetch ..pdb


display>sequence check
select sel, resn HOH
remove sel

```
alternative method using VMD
```vmd 
set selprotein [atomselect top protein] 	 
$selprotein writepdb common/ubq_ww_eq.pdb
```

Relaxation:

The system is now prepared for relaxation. Before initiating dynamics, it's essential to verify that the system is free of steric clashes or inappropriate geometry. This is achieved through a process known as energy minimization (EM)

## 1.7.Creating a Gromacs topology (PDB2GMX)


We are now ready to generate a GROMACS topology for the system. The GROMACS "pdb2gmx" command is employed for this purpose, converting a PDB coordinate file into a GROMACS topology file while also creating a processed structure file in the GROMACS format (.gro).

Then you can run:

```bash

gmx pdb2gmx -f 5pep.pdb -o 5pep.gro -water tip3p

```


You will be prompted to select the forcefield you would like to use. GROMACS comes with a number of AMBER and GROMOS forcefields, as well as a CHARMM and an OPLS-AA option. You will also need to specify your water model (choices included are TIP models, and the SPC and SPC/E models).

As it was discussed above we will choose CHARMM36 FF. Type 
```
1
```

Then the selection of water model  will be offered. we will choose TIP3P. Type

```
1
```


If you did everything right you will get this:

    You have successfully generated a topology from: 5pep.pdb. The Charmm36-jul2017 force field and the tip3p water model are used.


Now we will put our dimer into a periodic box using the following command
```bash
gmx editconf -f conf.gro -c -d 1.0 -bt cubic -o 5pep-box.gro
```
The above command centers the protein in the box (-c), and places it at least 1.0 nm from the box edge (-d 1.0). The box type is defined as a cube (-bt cubic). 


Warning:
The size of the box needs to be selected with great care. When using periodic boundary conditions we need to make sure that we respect the minimum image convention.

If you now look at the new 5pep-box.gro file you should see the box dimensions have changed - the box is now cubic.


## 1.7.Solvate the system


 Now we need to add solvent to the system , we are going to use water as a solvent.

Type the following command:
```bash
gmx solvate -cp 5pep-box.gro -cs spc216.gro -o 5pep-solv.gro -p topol.top
```
-cs spc216.gro determines the water structure file The configuration of the protein (-cp) is contained in the output of the previous editconf step, and the configuration of the solvent (-cs) is part of the standard GROMACS installation.

## 1.8.Neutralize the system

 

The last step of our system preparation phase consists of neutralizing the system by adding ions. This is accomplished via the genion command to replace solvent molecules with ions. By adding ions we are able to  neutralise any charge in your system; and it allows you to simulate systems with similar salt concentrations to their real-world equivalents.


Adding ions is done in two parts: first, you need to use the grompp tool to generate a .tpr file to be used when adding ions,
```bash
# create and empty mdp file 
touch mdrun.mdp
gmx grompp -f mdrun.mdp -c 5pep-solv.gro -p topol.top -o ions.tpr 
```
Now that the ".tpr" file has been generated, "genion" can be employed to neutralize the system's charge. This process involves reducing the system's charge by replacing certain components with anions and cations. To execute this, run the following command:


```bash
gmx genion -s ions.tpr -o system.gro -p topol.top -pname NA -nname CL -neutral -conc 0.15
```

When it is promted choose "SOL" (17). It will replace solvent molecules with positive ions. -pname and -nname defines names of positive and negative charges -np defines the amount of the positive charges -nn defines the amount of negative charges
# 2.Simulation
The second phase of this tutorial consists of actually simulating the system we prepared.

The simulation procedure involves multiple steps.

## 2.1.Energy minimisation
This procedure aims to minimize the potential energy of the system by adjusting the atomic coordinates. It will stabilize the overall structure and avoid steric clashes.

Before running a “useful” simulation, we will want to equilibrate our system. In the previous session, we added water and ions to our system to random positions within our simulation box. This may mean that some of their positions or orientations are not physically likely. For example when randomly adding water molecules to our simulation box, we may have ended up with a number of molecules whose oxygen atoms are close to one another, or with an ion in close proximity to a part of our protein with a similar charge. By equilibrating our system before running it, we can push our simulation towards a more physically realistic state.

There are many ways of equilibrating a system. You can:

    1.Run an equilibration simulation that starts at a temperature of 0 K and increases slowly to the temperature that you’re wanting to simulate;
    2.Run an equilibration simulation with a very small timestep to allow your system to “relax” into a more realistic configuration
    3.Reduce the potential energy of a system by using an energy-minimisation techniques;
    4.etc.


 Creating portable binary run file (or .tpr file) whith those parameters.
    integrator = md
    dt = 0.002 (ps)
    nsteps = 50000

In this particular example, the default parameters will be used: md integrator algorithm, a time step of 2fs, and a total of 50,000 md steps (100ps).

fistly, download   minim.mdp File, it contains the input parameters for energy minimization.

secondly, we will group all  the information for our simulation into a simulation input file (.tpr file). We can do this by running:


```bash

gmx grompp -f minim.mdp -c 1kx5_solv_ions.gro -p topol.top -o em.tpr

```
We are now ready to invoke mdrun to carry out the EM:
```bash
gmx mdrun -v -deffnm em
```


-v makes mdrun verbose, it makes mdrun write every step it makes on the screen. -deffnm - defines file names

There are two crucial factors to consider when evaluating the success of energy minimization (EM). 
1. Firstly, the potential energy (Epot) should be negative, typically ranging from 105 to 106 for a simple protein in water, depending on the system's size and the number of water molecules. 

2. Secondly, the maximum force (Fmax) is essential, with the target set in minim.mdp as "emtol = 1000.0," indicating a desired Fmax of no more than 1000 kJ mol-1 nm-1.



Now let's do some analysis. The em.edr file contains all of the energy terms that GROMACS collects during EM. We are going to analyze them via GROMACS energy module.
```bash
gmx energy -f em.edr -o potential.xvg
```
Once you launched the command you will be asked to choose which property you want to obtain.
From here, you can select the property you want to extract by simply typing the corresponding number followed by a 0.

To obtain an xvg file for the potential you can type:
```
14 0
```
Let’s assume that you generated a potential.xvg file, that is an XMGrace file that has the timestep in the first column and the total potential energy of the system in the second column. We can plot it to check that the energy has decreased and is now at or near a minimum.ofcorse we supposed that you have Grace installed.

	```bash

xmgrace volume.xvg
```
To save the image in the .png format you can enter the following:

```bash
xmgrace -nxy volume.xvg -hdevice PNG -hardcopy -printfile volume.png

```


Now that we have equilibrated our system, we can move on to simulating it in more physical conditions. We will run a simulation in the isobaric-isothermic ensemble (fixed pressure and temperature).


## 2.2.Running NVT equilibration
After successful energy minimization, confirming the system's geometry and solvent orientation, we move on to equilibrating the solvent and ions surrounding the protein before initiating real dynamics. This initial phase is performed within an NVT ensemble (constant Number of particles, Volume, and Temperature), also known as "isothermal-isochoric" or "canonical" ensemble. The duration of this procedure depends on the system's composition, but in NVT, the system's temperature should stabilize at the desired value over time.

In this stage, we will constrain the water molecules and allow only the movement of water and ions. We will utilize the “nvt.mdp” file, which contains the input parameters necessary for NVT equilibration.


```bash	
gmx grompp -f nvt.mdp -c em.gro -r em.gro -p topol.top -o nvt.tpr
```
Now you can run the simulation with this command

```bash
gmx mdrun -v -deffnm nvt
```
Once the simulation is over we can use the gmx energy to extract properties from the run.

	
```bash
 gmx energy -f nvt.edr -o nvt_temperature.xvg
```
When it is promted choose (temperature) 15 then 0  to exit the program.
```bash
15 0
```
Plot the xvg file.

```bash
xmgrace nvt_temperature.xvg
```

NB.As you can see the temperature is quite stable around the selected value (300K). In general, if the temperature has not yet stabilized, additional time will be required.

To verify that the system is properly equilibrated we can check the RMSD value via the gmx rms module.
```bash
gmx rms -f nvt.trr -s nvt.tpr -o nvt_rmsd.xvg
```
When it is promted choose (Backbone) 4 then 4 (Backbone)
```bash
4 4
```
This will output an xvg file (nvt_rmsd.xvg) with the backbone to backbone RMSD for the protein. Let’s plot it.
```bash
xmgrace nvt_rmsd.xvg
```
The RMSD quickly converges to a stable value signaling that the system has equilibrated at the desired temperature. We can move to the next phase
## 2.3.Running NPT  equilibration

Following the NVT equilibration step, which stabilized the system's temperature, the next phase involves equilibrating the system's pressure under an NPT ensemble. This ensemble maintains constant values for the Number of particles, Pressure, and Temperature. During this step, only water and ions are allowed to move.

Similar to the procedure for NVT equilibration, we will execute "grompp" and "mdrun." However, this time, we will include the "-t" flag to incorporate the checkpoint file from the NVT equilibration. This file contains all the necessary state variables required to continue our simulation.
```bash
gmx grompp -f npt.mdp -c nvt.gro -t nvt.cpt -p topol.top -o npt.tpr
```
```bash
gmx mdrun -deffnm npt
```


For analyzing:
```bash
gmx energy -f npt.edr -o pressure.xvg
```
Type "17 0" at the prompt to select the pressure of the system and exit.


You can see that the protein widely fluctuates around the pressure value we selected in the mdp file (1 bar). Don’t worry, large pressure fluctuations are completely normal for nanoscale simulations of aqueous solutions.

Once again, we can check the RMSD value to control that the system has properly equilibrated at constant pressure.
```bash
printf "4 4" | gmx rms -f npt.trr -s npt.tpr -o npt_rmsd.xvg
```
plot it whith 
```bash
xmgrace npt_rmsd.xvg
```
As you see, the RMSD values are very stable over time, indicating that the system is well-equilibrated at the desired pressure.



We will also look at the density

gmx energy -f npt.edr -o density.xvg
Type "23 0" at the prompt to select the density of the system and exit.

Then we will unfix protein. Create npt1.mdp from npt.mdp commenting the second string in it.


# 3.MD Calculation
The production run will still be carried out within the NPT ensemble using the "V-rescale" thermostat and the "Parrinello-Rahman" barostat.

For the purpose of demonstration, we are just going to perform a 1 ns simulation. Keep in mind that in a real experiment you will need to run hundreds of ns to extract useful information from your system.

```bash
gmx grompp -f md.mdp -c npt.gro -t npt.cpt -p topol.top -o md.tpr
```
```bash
gmx mdrun -v -deffnm md
```
The production run will be a little bit longer than the previous one. If you want to run it in background use:
```bash
nohup gmx mdrun -v -deffnm md &
```
This will append the output in a temporary file named nohup.out

You can check the progress of the simulation with the following command:
```bash
tail -f nohup.out
```
When the simulation is over you will still find the relevant informations in the usual files (md.edr, md.log, md.xtc)

Congratulations! you have now successfully conducted your MD simulations.