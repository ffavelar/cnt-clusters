# grompp must be in $PATH

molmaker.py 	-ff martini_v2.0.itp 	-i cholesteryloleate.itp 			-o cholesteryloleate-raw.gro
molmaker.py 	-ff martini_v2.0.itp 	-i glycerinetrioleate.itp			-o glycerinetrioleate-raw.gro
molmaker.py 	-ff martini_v2.0.itp 	-i oleicacid.itp  				-o oleicacid-raw.gro
molmaker.py 	-ff martini_v2.0.itp 	-i palmiticacid.itp  				-o palmiticacid-raw.gro
molmaker.py 	-ff martini_v2.0.itp 	-i palmytoyloleolylphosphatidylcoline.itp	-o palmytoyloleolylphosphatidylcoline-raw.gro

rm *#
rm step*.pdb

