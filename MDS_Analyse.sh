#!/usr/bin/bash
###GROMAX TUTORAIL#####
##################Analysis###################

gmx trjconv -s md.tpr -f md.xtc -o md_noPBC.xtc -pbc mol -center

gmx rms -s md.tpr -f md_noPBC.xtc -o rmsd.xvg -tu ns


###calculate RMSD relative to the crystal structure, we could issue the following:

gmx rms -s em.tpr -f md_noPBC.xtc -o rmsd_xtal.xvg -tu ns