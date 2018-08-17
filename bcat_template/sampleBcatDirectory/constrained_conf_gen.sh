module unload Babel
#module load Babel

export PATH=/home/paulzim/openbabel241/bin/:$PATH
export PYTHONPATH=/home/paulzim/openbabel241/lib/
export LD_LIBRARY_PATH=/home/paulzim/openbabel241/lib:$LD_LIBRARY_PATH


~joshkamm/miniconda3/bin/python constr_conformer.py $1 $2 $3


