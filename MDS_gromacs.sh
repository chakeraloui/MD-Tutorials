#!/bin/bash -e

# Dr Chaker Aloui, molecular dynamics simulation, gromacs 
tput bold
echo -e "Chaker ALOUI \n"
echo -e "Running gromacs MDS ...\n"

#configPath="./configs"
#dataPath="./data" # mettre les models pdb ici
#prot2simulate="./data/model.pdb" # a changer a chaque fois
#outPutDir="simulation_date" # à modifier à chaque analyse
#mkdir $outPutDir

#pip install propka
# Vérifie si 'propka' est installé en utilisant pip list
if pip list | grep -F propka > /dev/null; then
    echo "Le package 'propka' est déjà installé."
else
    echo "Installation du package 'propka'..."
    pip install propka
    if [ $? -eq 0 ]; then
        echo "'propka' a été installé avec succès."
    else
        echo "L'installation de 'propka' a échoué."
    fi
fi

# script a partir d un tuto internet:
# Lysozyme PDB code 1AKI http://www.mdtutorials.com/gmx/lysozyme/01_pdb2gmx.html

# telecharcher  le fichier pdb (wget https://files.rcsb.org/download/1AKI.pdb) 
# Preparer les fichiers de parametres rep configs/
#wget http://www.mdtutorials.com/gmx/lysozyme/Files/ions.mdp
#wget http://www.mdtutorials.com/gmx/lysozyme/Files/minim.mdp
#wget http://www.mdtutorials.com/gmx/lysozyme/Files/nvt.mdp
#wget http://www.mdtutorials.com/gmx/lysozyme/Files/npt.mdp
#wget http://www.mdtutorials.com/gmx/lysozyme/Files/md.mdp

###GROMAX TUTORAIL#####
#REMOVE ALL THE WATER CRYSTAL


############TOPOLOGY##########################
##1##pdb2gmx to generate topolgy, restrain and processed pdb file used to generate a file that can be used by GROMACS

#RepairPDB:It is highly recommended to repair your structures before you do any modelling with FoldX. RepairPDB identify those residues which have bad torsion angles, or VanderWaals' clashes, or total energy, and repairs them. The minimal configuration file for RepairPDB is:


#It can be run from the command line:
echo "Downloading $1.pdb "
wget http://files.rcsb.org/download/$1.pdb

echo ""
FoldX --command=RepairPDB --pdb=$1.pdb


echo "Assessing the protonation states of the protein.."
propka3 $1\_clean.pdb

echo "Fetching  the CHARMM force field.."

wget http://mackerell.umaryland.edu/download.php?filename=CHARMM_ff_params_files/charmm36_ljpme-jul2022.ff.tgz

echo "Removing Water Molecules"

grep -v HOH $1 > $1\_clean_protein.pdb
echo "Creating a Gromacs topology (PDB2GMX)"
gmx pdb2gmx \  
-f $1\_clean_protein.pdb \
-o $1\_processed.gro \
-water  tip3p #or spce or tip4p



#Syntax

    #Call the program with gmx
    #Select the pdb2gmx module
    #Call the -f flag and provide the starting PDB structure protein.pdb
    #Select the -o flag and decide how you want to name the output file
    #The -water flag allows us to select the water model we want to use

#Note. there is  three water models, namely, TIP3P, SPC/E et TIP4P

 echo "topolgy file generated ..."


#defines the size and shape of the box for our simulation  using editconf.
echo "Now we will put our dimer into a cubic box"

gmx editconf \
-f $1\_processed.gro \
-o $1\_newbox.gro \ 
-c \
-d 1.0 \ 
-bt cubic
 

#Syntax

   # Call the program with gmx
    #Select the editconf command
    #Select the -f flag and provide the starting structure (protein.gro)
    #Call the -o flag and decide how you want to name the output file (protein_box.gro)
    #The -c places our system in the center of the box
    #Choose the -d flag and place our system at the selected distance from the box (at least 1.0 nm in our case)
    #Decide the shape of the simulation box with the -bt flag

#Warning: The size of the box needs to be selected with great care.   


# Adding Solvent
echo "adding solvation"
gmx solvate \
-cp $1\_newbox.gro \
-cs spc216.gro \
-o $1\_solv.gro \
-p topol.top

#Syntax

    #Call the program with gmx
    #Select the solvate command
    #Select the -cp flag and provide the starting structure for the solute (protein_box.gro)
    #Specify the solvent with the -cs flag. The default solvent is Simple Point Charge (spc) water
    #Call the -o flag and decide how you want to name the output file (protein_solv.gro)
    #Finally, we use the -p flag and provide the previously generated topology file


#In the final step of the system preparation, we use the genion command 
#to replace solvent molecules with ions. By adding ions we are able to neutralize the system.
# Adding Ions (wget http://www.mdtutorials.com/gmx/lysozyme/Files/ions.mdp)  
 
#Generate .tpr file
echo "creating ions.tpr file ."
touch ions.mdp
gmx grompp \
-f ions.mdp \
-c $1\_solv.gro \
-p topol.top \
-o ions.tpr 
-po mdout.mdp


echo "Adding ions"
gmx genion \
-s ions.tpr \
-o $1\_solv_ions.gro \
-p topol.top \
-pname NA \
-nname CL 
-conc 0.15
-neutral
#or -conc 0.1


#Syntax

    #Call the program with gmx
    #Select the genion command
    #Select the -s flag and provide the .tpr file for the ions (ions.tpr)
    #Call the -o flag and decide how you want to name the output file (system.gro)
    #Use the -p flag and provide the topology file (topol.top)
    #Choose the -pname flag to decide the name of the positive ions name
    #The -nname flag is used to decide the name of the negative charge ions
    #The -neutral flag tells GROMACS to neutralize the net charge of the system

#Warning: To use this command we need to generate a tpr file for the ions. This is done via the grompp command
echo " system preparation end  successfully"
#  Energy Minimization (wget http://www.mdtutorials.com/gmx/lysozyme/Files/minim.mdp) 
# em.tpr production: we get a .tpr file with coordinates (system.gro), topology (topol.top), and parameters for minimization (simulation.mdp). Now we can use this file to run the simulation.
echo "we have our system ready to be simulated ..."
gmx grompp \
-f minim.mdp \
-c $1\_solv_ions.gro \
-r $1\_solv_ions.gro \
-t system.cpt \
-p topol.top \
-o em.tpr
-po mdout.mdp
#Syntax

    #Call the program with gmx
    #Select the grompp command
    #Select the -f flag and provide the mdp file with the parameters for the run (simulation.mdp)
    #The -c flag is used to provide the coordinates of atoms
    #The -r flag is needed if we are using position restraints. Mostly during the equilibration steps of your simulation.
    #In some cases we may want to continue our simulation starting from a previous one. To this aim, we can use the -t flag coupled with a previously generated cpt file.
    #Use the -p flag and provide the topology file (topol.top)
    #Call the -o flag and decide how you want to name the output file (simulation.tpr)




#launch the simulation with the mdrun command
gmx mdrun -v -deffnm em -nb gpu


#Syntax
        
        
            
#Call the program with gmx
#Select the mdrun command
#The -v makes the command verbose. It is useful to make clearer what the program is actually doing while running.
#The -deffnm option is followed by the prefix of the tpr file we are using. This will generate files having the same prefix.

#  (simulation.tpr $\Rightarrow$ -deffnm simulation)

#After completion you will generally get five files:

    #simulation.log containing information about the run
    #simulation.gro with the final structure of the system
    #simulation.edr with energies that GROMACS collects during the simulation
    #simulation.xtc a “light” weight trajectory with the coordinates of our system in low precision
   # simulation.trr with the higher precision trajectory of positions, velocities, and forces during the simulation

#You can use these files to perform a variety of different analysis.
#Now let's do some analysis. The em.edr file contains all of the energy terms that GROMACS collects during EM. We are going to analyze them via GROMACS energy module.
printf "14 0" | gmx energy -f em.edr -o potential.xvg



##6##NVT equilibration
echo "minimising energy:NVT equilibration."

gmx grompp -f nvt.mdp -c em.gro -r em.gro -p topol.top -o nvt.tpr


gmx mdrun -deffnm nvt
printf "15 0" | gmx energy -f nvt.edr -o nvt_temperature.xvg
echo "checking the RMSD value To verify that the system is properly equilibrated"
printf "4 4" | gmx rms -f nvt.trr -s nvt.tpr -o nvt_rmsd.xvg
##7##NPT equilibration
echo "minimising energy:NPT equilibration."
gmx grompp -f npt.mdp -c nvt.gro -r nvt.gro -t nvt.cpt -p topol.top -o npt.tpr

gmx mdrun -deffnm npt -nb gpu

printf "17 0" | gmx energy -f npt.edr -o pressure.xvg
echo "check the RMSD value to control that the system has properly equilibrated at constant pressure"
printf "4 4" | gmx rms -f npt.trr -s npt.tpr -o npt_rmsd.xvg
echo "starting simulation"

gmx grompp -f md.mdp -c npt.gro -t npt.cpt -p topol.top -o md.tpr -po mdout.mdp
gmx mdrun -v -deffnm md -nb gpu

# fin de la simulation
echo -e "Fin de la simulation \n"
#*******************************************************
#*******************************************************
#*******************************************************
# Step Nine: Analysis


#https://tutorials.gromacs.org/md-intro-tutorial.html à tester notebook jupyter
# xmgrace potential.xvg


