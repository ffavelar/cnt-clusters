# grompp and molmaker.py must be in $PATH

for MOLECULE in hydroxyectoine # cholesteryloleate glycerinetrioleate oleicacid palmiticacid palmytoyloleolylphosphatidylcoline
do

	molmaker.py 	-ff martini_v2.2P.itp 	-i $MOLECULE.itp 			-o $MOLECULE-raw.gro

	editconf -f $MOLECULE-raw.gro -bt cubic -d 1.3 -o $MOLECULE-box.gro
	grompp -f minimization.mdp -c $MOLECULE-box.gro -p $MOLECULE.top -o $MOLECULE.tpr
	mdrun -v -deffnm $MOLECULE
	
	editconf -f $MOLECULE.gro -bt cubic -d 0.1 -o ../$MOLECULE.gro

	vmd ../$MOLECULE.gro

	rm *#
	rm step*.pdb

done
