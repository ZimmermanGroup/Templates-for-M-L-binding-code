module unload Babel
#module load Babel

export PATH=/home/paulzim/openbabel241/bin/:$PATH
export PYTHONPATH=/home/paulzim/openbabel241/lib/
export LD_LIBRARY_PATH=/home/paulzim/openbabel241/lib:$LD_LIBRARY_PATH


~joshkamm/miniconda3/bin/python constr_conformer.py $1 $2 $3

# use josh's python 
# $mol = molecule
# output = conf_id
# followed by atoms you want to unfreeze ex: atoms 18, 37, 38, 41 can be rotated - unfreeze atoms involved in torsion you want to rotate
done
