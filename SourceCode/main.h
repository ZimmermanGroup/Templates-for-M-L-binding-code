#include <iostream>
#include <fstream>
#include <stdio.h> 
#include <sys/stat.h>
    
#include "utils.h"
#include "mopac.h"
#include "xtb.h"
#include "stringtools.h"
#include "align.h"   
     
#include <openbabel/mol.h>
#include <openbabel/obconversion.h>
#include <openbabel/math/align.h>
      
using namespace std;
using namespace OpenBabel;
      
void xyz_read(int natoms, string* anames, double* coords, string xyzfile);
void xyz_read_last(int natoms1, double* coords, string xyzfile);
int get_natoms(string filename);
int get_charge(string filename);
void get_all_xyz(int natoms, string* anames, vector<double*> &xyzs, string xyzfile);
int get_unique_conf(int nstruct, int* unique);
void align_and_opt(int natoms1, int natoms2, string* anames, string* anamesm, string* anamest, int* anumbers, int* anumbersm, int charget, int nstruct, int* unique, vector<double*> xyzall, double* xyzm, int align, int rotation);
void write_all_xyz(int natoms, string* anames, double* E, vector<double*> xyzs, string xyzfile_string);
void write_all_xyz(int natoms, string* anames, int nstruct, double* E, double** xyzs, string xyzfile_string);
void write_gsm(int natoms, string* anames, int charge, int nstruct, double* E, double** xyzs, int nadd, int* adds);
void do_gsm(int nstruct);


