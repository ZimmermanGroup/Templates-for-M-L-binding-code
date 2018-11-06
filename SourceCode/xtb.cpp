#include "xtb.h"
using namespace std;
#include "constants.h"

//GRADNORM was 1.5
#define GRADNORM 5.

void XTB::write_ic_input(ofstream& inpfile, int anum, ICoord ic){

  printf(" ERROR: XTB write_ic_input not implemented \n");
  exit(1);

  return;
}

// standard opt after checking for existing file
double XTB::opt_check(string filename)
{
#if SKIPXTB
  return 0.;
#endif

  struct stat sts;
  if (stat(filename.c_str(), &sts) != -1)
  {
    //printf("  XTB already done \n");
    energy = read_output(filename);
    xyz_read(filename);
    return energy;
  }
   
  return opt(filename);
}

double XTB::opt() {

  string filename = "scratch/testxtb.mop";

  energy = opt(filename);
 
  return energy;
}

void XTB::opt_write() {

  string filename = "scratch/testxtb.mop";

  opt_write(filename);
 
  return;
}

double XTB::opt(string filename, ICoord icoords) {

  printf(" IC opt for XTB not implemented! \n");
  exit(1);

#if 0
  printf(" WARNING: bypassing ICs! \n");
  return opt(filename);
#endif

   //printf(" in mopac/opt() w/IC freeze \n");

  energy0 = energy = 0;

#if SKIPXTB
  printf(" skipping mopac opt! \n");
#endif

  ofstream inpfile;
  string inpfile_string = filename;
  ofstream xyzfile;
  string xyzfile_string = "scratch/testxtb.xyz";
#if !SKIPXTB
  inpfile.open(inpfile_string.c_str());
  inpfile.setf(ios::fixed);
  inpfile.setf(ios::left);
  inpfile << setprecision(6);

  xyzfile.open(xyzfile_string.c_str());
  xyzfile.setf(ios::fixed);
  xyzfile.setf(ios::left);
  xyzfile << setprecision(6);
  xyzfile << " " << natoms << endl << endl;

  for (int i=0;i<natoms;i++)
  //if (!skip[i])
  {
    xyzfile << " " << anames[i] << " " << xyz[3*i+0] << " " << xyz[3*i+1] << " " << xyz[3*i+2] << endl;
  }

  string cmd = "/export/zimmerman/adewyer/bin/xtb/xtb "+filename+" -opt ";
  system(cmd.c_str());
#endif

  energy = read_output(filename);
 
  //printf(" initial energy: %1.4f final energy: %1.4f \n",energy0,energy); 

  // need to retrieve final geometry, write to xyz
  xyz_read(inpfile_string);
  xyz_save(inpfile_string+".xyz");

  xyzfile.close();
  inpfile.close();

  if (abs(energy)<0.00001)
  {
    printf(" energy zero, XTB failed \n");
    return 10000;
  }

  return energy;
}


void XTB::opt_write(string filename, ICoord icoords) {

  printf(" ERROR: opt_write for ic's not implemented \n");
  exit(1);

  return;
}


double XTB::grad(string filename) 
{
  printf(" in XTB/grad() \n");
 
  energy0 = energy = 0;

#if SKIPXTB
  printf(" skipping xtb grad! \n");
#endif

  ofstream inpfile;
  string inpfile_string = "scratch/"+filename;
  string outfile0 = filename+".out";
  string outfile = "scratch/"+filename+".xyz";

  //printf("  filename: %s inpfile: %s outfile: %s \n",filename.c_str(),inpfile_string.c_str(),outfile.c_str());
#if !SKIPXTB
  inpfile.open(inpfile_string.c_str());
  inpfile.setf(ios::fixed);
  inpfile.setf(ios::left);
  inpfile << setprecision(6);

  inpfile << " " << natoms << endl << endl;

  for (int i=0;i<natoms;i++)
  //if (!skip[i])
  {
    inpfile << " " << anames[i] << " " << xyz[3*i+0] << " " << xyz[3*i+1] << " " << xyz[3*i+2] << endl;
  }

  inpfile << endl << "$set" << endl;
  inpfile << "charge " << charge << endl;
  
  inpfile << "$end" << endl;


  inpfile.close();

  string cmd = "cd scratch; /export/zimmerman/adewyer/bin/xtb/xtb "+filename+" -opt > "+outfile0+"; mv xtbopt.xyz "+filename+".xyz";
  //printf(" cmd: %s \n",cmd.c_str());
  system(cmd.c_str());
#endif

  energy = read_output("scratch/"+filename);
 
  // need to retrieve final geometry, write to xyz
  xyz_read(outfile);
  xyz_save(inpfile_string+".xyz");


  if (abs(energy)<0.00001)
  {
    printf(" energy zero, xtb failed \n");
    return 10000.;
  }

  return energy;
}

// standard opt
double XTB::opt(string filename) 
{
  //printf(" in XTB/opt() \n");
#if 0
  printf(" frozen:");
  for (int i=0;i<natoms;i++)
    printf(" %i",frzlistb[i]);
  printf("\n");
#endif
 
  energy0 = energy = 0;

#if SKIPXTB
  printf(" skipping xtb opt! \n");
#endif

  ofstream inpfile;
  string inpfile_string = "scratch/"+filename;
  string outfile0 = filename+".out";
  string outfile = "scratch/"+filename+".xyz";

  //printf("  filename: %s inpfile: %s outfile: %s \n",filename.c_str(),inpfile_string.c_str(),outfile.c_str());
#if !SKIPXTB
  inpfile.open(inpfile_string.c_str());
  inpfile.setf(ios::fixed);
  inpfile.setf(ios::left);
  inpfile << setprecision(6);

  inpfile << " " << natoms << endl << endl;

  for (int i=0;i<natoms;i++)
  //if (!skip[i])
  {
    inpfile << " " << anames[i] << " " << xyz[3*i+0] << " " << xyz[3*i+1] << " " << xyz[3*i+2] << endl;
  }

  inpfile << endl << "$set" << endl;
  inpfile << "charge " << charge << endl;
  if (nfrz>0)
  {
    inpfile << "fix";
    for (int j=0;j<natoms;j++)
    if (frzlistb[j])
      inpfile << " " << j+1;
    inpfile << endl;
  }
  inpfile << "$end" << endl;


  inpfile.close();

  string cmd = "cd scratch; /export/zimmerman/adewyer/bin/xtb/xtb "+filename+" -opt > "+outfile0+"; mv xtbopt.xyz "+filename+".xyz";
  //printf(" cmd: %s \n",cmd.c_str());
  system(cmd.c_str());
#endif

  energy = read_output("scratch/"+filename);
 
  // need to retrieve final geometry, write to xyz
  xyz_read(outfile);
  xyz_save(inpfile_string+".xyz");


  if (abs(energy)<0.00001)
  {
    printf(" energy zero, xtb failed \n");
    return 10000.;
  }

  return energy;
}

void XTB::opt_write(string filename) {

   //printf(" in mopac/opt() \n");

  energy0 = energy = 0;

#if SKIPXTB
  printf(" skipping mopac opt! \n");
#endif

  ofstream inpfile;
  string inpfile_string = filename;
#if !SKIPXTB
  inpfile.open(inpfile_string.c_str());
  inpfile.setf(ios::fixed);
  inpfile.setf(ios::left);
  inpfile << setprecision(6);

  inpfile << " " << natoms << endl << endl;

  for (int i=0;i<natoms;i++)
  //if (!skip[i])
  {
    inpfile << " " << anames[i] << " " << xyz[3*i+0] << " " << xyz[3*i+1] << " " << xyz[3*i+2] << endl;
  }

  inpfile << "$set" << endl;
  inpfile << "charge " << charge << endl;
  inpfile << "$end" << endl;

  inpfile.close();
#endif

  return;
}

double XTB::read_output(string filename) {

  //double energy = -1;
  energy = 10000;

  string oname = filename+".out";
  ifstream output(oname.c_str(),ios::in);
  string line;
  vector<string> tok_line;
  int ne = 0;
  while(!output.eof()) 
  { 
    getline(output,line);
    //cout << " RR " << line << endl;
    if (line.find("total E")!=string::npos)
    {
      cout << " found E: " << line << endl;
      tok_line = StringTools::tokenize(line, " \t");
      if (ne++==0)
        energy0 = atof(tok_line[3].c_str());
      else
        energy = atof(tok_line[3].c_str());
    }
  }

  return energy;
}

void XTB::alloc(int natoms_i) {
 
//  printf(" in mopac/alloc() \n");

  charge = 0;

  natoms = natoms_i;
  anumbers = new int[natoms];
  anames = new string[natoms];

  xyz0 = new double[3*natoms];
  xyz = new double[3*natoms];

  nfrz = 0;
  frzlist = new int[natoms];
  frzlistb = new int[natoms];

  nskip = 0;
  skip = new int[natoms];
  for (int i=0;i<natoms;i++)
    skip[i] = 0;

  return;
}

void XTB::init(int natoms_i, int* anumbers_i, string* anames_i, double* xyz_i) {
 
//  printf(" in mopac/init() \n");

  charge = 0;

  natoms = natoms_i;
  anumbers = new int[natoms];
  anames = new string[natoms];

  xyz0 = new double[3*natoms];
  xyz = new double[3*natoms];

  for (int i=0;i<natoms;i++)
    anumbers[i] = anumbers_i[i];
  for (int i=0;i<natoms;i++)
    anames[i] = anames_i[i];
  for (int i=0;i<3*natoms;i++)
    xyz0[i] = xyz[i] = xyz_i[i];  

  nfrz = 0;
  frzlist = new int[natoms_i]();
  frzlistb = new int[natoms_i]();

  nskip = 0;
  skip = new int[natoms];
  for (int i=0;i<natoms;i++)
    skip[i] = 0;
  for (int i=0;i<natoms;i++)
  if (anames[i]=="X")
  {
    nskip++;
    skip[i] = 1;
  }

  return;
}

void XTB::freemem() {

  delete [] xyz0;
  delete [] xyz;
  delete [] anumbers;
  delete [] anames;

  delete [] frzlist;
  delete [] frzlistb;

  delete [] skip;

  return;
}

void XTB::freeze(int* frzlist_new, int nfrz_new, int nfrz0_new) {

  nfrz0 = nfrz0_new;
  nfrz = nfrz_new;
  for (int i=0;i<natoms;i++)
    frzlistb[i] = 0;
  for (int i=0;i<nfrz;i++)
    frzlistb[frzlist_new[i]] = 1;
  for (int i=0;i<nfrz;i++)
    frzlist[i] = frzlist_new[i];

#if 0
  printf(" freeze list: ");
  for (int i=0;i<nfrz;i++)
    printf("%i ",frzlist[i]);
  printf("\n");
#endif
  
  return;
}

void XTB::freeze_d(int* frzlist_new) {

  for (int i=0;i<natoms;i++)
    frzlist[i] = frzlist_new[i];
  for (int i=0;i<natoms;i++)
    frzlistb[i] = frzlist_new[i];

  nfrz = 0;
  for (int i=0;i<natoms;i++)
  if (frzlistb[i])
    nfrz++;

#if 1
  printf(" freeze list: ");
  for (int i=0;i<natoms;i++)
    printf("%i ",frzlistb[i]);
  printf("\n");
#endif
  
  return;
}

void XTB::reset(int natoms_i, int* anumbers_i, string* anames_i, double* xyz_i) {
 
//  printf(" in mopac/reset() \n");

#if 0
  if (natoms!=natoms_i)
  {
    printf(" mopac reset failed due to different # of atoms \n");
    return;
  }
#endif
  natoms = natoms_i;

  for (int i=0;i<natoms;i++)
    anumbers[i] = anumbers_i[i];
  for (int i=0;i<natoms;i++)
    anames[i] = anames_i[i];
  for (int i=0;i<3*natoms;i++)
    xyz0[i] = xyz[i] = xyz_i[i];  

  nskip = 0;
  for (int i=0;i<natoms;i++)
  if (anames[i]=="X")
  {
    nskip++;
    skip[i] = 1;
  }
  else
    skip[i] = 0;

  nfrz = 0;

  return;
}

void XTB::set_charge(int c0)
{
  charge = c0;
  return;
}

void XTB::xyz_save(string filename){

  ofstream xyzfile;
  xyzfile.open(filename.c_str());
  xyzfile.setf(ios::fixed);
  xyzfile.setf(ios::left);
  xyzfile << setprecision(6);

   xyzfile << " " << natoms << endl;
   xyzfile << " " << endl;
   for (int i=0;i<natoms;i++)
   {
     xyzfile << "  " << anames[i];
     xyzfile << " " << xyz[3*i+0] << " " << xyz[3*i+1] << " " << xyz[3*i+2];
     xyzfile << endl;
   }

  xyzfile.close();
  return;
}

void XTB::xyz_read(string filename)
{
  printf(" in xyz_read for file: %s \n",filename.c_str());
  string oname = filename;
  ifstream output(oname.c_str(),ios::in);
  string line;
  vector<string> tok_line;
  int count = 0;
  int i = 0;

  while(!output.eof()) 
  { 
    getline(output,line);
    //cout << " RR " << line << endl;
    if (count == 2)
    {
      while (skip[i] && i<natoms)
        i++;
      if (i==natoms) break;
      if (!StringTools::cleanstring(line)) break;
      vector<string> tok_line = StringTools::tokenize(line, " \t");
      xyz[3*i+0]=atof(tok_line[1].c_str());
      xyz[3*i+1]=atof(tok_line[2].c_str());
      xyz[3*i+2]=atof(tok_line[3].c_str());
//      cout << tok_line[0] << " " << tok_line[1] << " " << tok_line[2] << " " << endl; 
      i++;
    }
    count++;
    if (i>=natoms) break;
  }

  output.close();
 
  return;
}   

