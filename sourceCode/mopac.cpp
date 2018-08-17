#include "mopac.h"
//#include "utils.h"
using namespace std;
#include "constants.h"

//NOTE: max time for MOPAC: 12M

//GRADNORM was 1.5
#define GRADNORM 5.
#define USE_AUX 1

void Mopac::write_ic_input(ofstream& inpfile, int anum, ICoord ic){

  //printf("\n\n in write_ic_input() for atom %i \n",anum);

//  printf(" nbonds: %i \n",ic.nbonds);
 
  for (int i=0;i<nfrz0;i++)
  {
    if (anum==frzlist[i] || anum < 4)
    {
      //printf(" this atom was moved or is in the first 3, doing xyz \n");
      inpfile << xyz[3*anum+0] << " 0 " << xyz[3*anum+1] << " 0 " << xyz[3*anum+2] << " 0 " << endl;
      return;
    }
  }

//  printf(" atom attached to moved atom, freezing bond \n"); 

  int b1,c1,d1;
  for (int i=0;i<ic.nbonds;i++)
  {
    b1 = -1;
    if (ic.bonds[i][0] == anum || ic.bonds[i][1] == anum)
    for (int j=0;j<nfrz0;j++)
    {
      //printf(" found %i, checking for bond partner %i \n",anum,frzlist[j]);
      if (ic.bonds[i][0]==frzlist[j]
        || ic.bonds[i][1]==frzlist[j])
      {
        b1 = frzlist[j];
        break;
      }
    }

    if (b1>-1)
    {
      //printf(" found attachment: %i %i \n",b1,anum);
      break;
    }
  }

//could change this to "xyz list" if statement
  c1 = -1; d1 = -1;
  for (int i=0;i<natoms;i++)
  {
    if (c1==-1)
    for (int j=0;j<i;j++)
    {
      if (!frzlistb[i] && !frzlistb[j])
      {
        c1 = i; d1 = j;
        //printf(" found angle/tor possibility: %i %i \n",c1,d1);
        break;
      }
    }
  }

  if (b1==-1 || c1==-1 || d1==-1 || b1==d1)
  {
    printf(" failed to find IC for atom %i, writing xyz \n",anum);
    inpfile << xyz[3*anum+0] << " 0 " << xyz[3*anum+1] << " 0 " << xyz[3*anum+2] << " 0 " << endl;
  }
  else
  {
    double rv,anglev,torv;
    rv = ic.distance(anum,b1);
    anglev = ic.angle_val(anum,b1,c1);
    torv = ic.torsion_val(anum,b1,c1,d1);
    //printf(" r: %1.3f angle: %1.3f tor: %1.3f \n",rv,anglev,torv);
    inpfile << " " << rv << " 0 " << anglev << " 1 " << torv << " 1 " << b1+1 << " " << c1+1 << " " << d1+1 << endl;
  }


  return;
}

// standard opt after checking for existing file
double Mopac::opt_check(string filename)
{
#if SKIPMOPAC
  return 0.;
#endif

  struct stat sts;
  if (stat(filename.c_str(), &sts) != -1)
  {
    //printf("  Mopac already done \n");
    energy = read_output(filename);
#if USE_AUX
    xyz_read_aux(filename);
#else
    xyz_read(filename);
#endif
    return energy;
  }
   
  return opt(filename);
}

double Mopac::opt() {

  string filename = "scratch/testmopac.mop";

  energy = opt(filename);
 
  return energy;
}

void Mopac::opt_write() {

  string filename = "scratch/testmopac.mop";

  opt_write(filename);
 
  return;
}

void Mopac::opt_header(ofstream& inpfile) {

  double gnorm = GRADNORM;
//  inpfile << " PM7 C.I.=2 " << endl; 
//  inpfile << " PM7 NOSYM T=12M " << endl;
//  inpfile << " PM7 NOSYM T=12M GNORM=0.7 " << endl;
#if !UNRESTRICTED
  inpfile << " PM6 NOSYM T=12M GNORM=" << gnorm << " CHARGE=" << charge << " AUX LOCALIZE " << endl;
#else
#if !TRIPLET
  inpfile << " PM6 NOSYM T=12M GNORM=" << gnorm << " UHF CHARGE=" << charge << " AUX LOCALIZE " << endl;
#else
  inpfile << " PM6 NOSYM T=12M GNORM=" << gnorm << " UHF TRIPLET CHARGE=" << charge << " AUX LOCALIZE " << endl;
#endif
#endif
  inpfile << "   MOPAC run " << endl;
  inpfile << "   ready to go? " << endl;

  return;
}

double Mopac::opt(string filename, ICoord icoords) {

#if 0
  printf(" WARNING: bypassing ICs! \n");
  return opt(filename);
#endif

   //printf(" in mopac/opt() w/IC freeze \n");

  energy0 = energy = 0;

#if SKIPMOPAC
  printf(" skipping mopac opt! \n");
#endif

  ofstream inpfile;
  string inpfile_string = filename;
  ofstream xyzfile;
  string xyzfile_string = "scratch/testmopac.xyz";
#if !SKIPMOPAC
  inpfile.open(inpfile_string.c_str());
  inpfile.setf(ios::fixed);
  inpfile.setf(ios::left);
  inpfile << setprecision(6);

  xyzfile.open(xyzfile_string.c_str());
  xyzfile.setf(ios::fixed);
  xyzfile.setf(ios::left);
  xyzfile << setprecision(6);
  xyzfile << " " << natoms << endl << endl;

  opt_header(inpfile);

  for (int i=0;i<natoms;i++)
  if (!skip[i])
  {
    if (!frzlistb[i])
      inpfile << " " << anames[i] << " " << xyz[3*i+0] << " 1 " << xyz[3*i+1] << " 1 " << xyz[3*i+2] << " 1 " << endl;
    else 
    {
      inpfile << " " << anames[i] << " ";
      write_ic_input(inpfile,i,icoords);
    }
    xyzfile << " " << anames[i] << " " << xyz[3*i+0] << " " << xyz[3*i+1] << " " << xyz[3*i+2] << endl;
  }
  
  string cmd = "/export/zimmerman/adewyer/mopac/MOPAC2016.exe "+filename;
//  string cmd = "/tmp/MOPAC2016.exe "+filename;
  system(cmd.c_str());
#endif

  energy = read_output(filename);
 
  //printf(" initial energy: %1.4f final energy: %1.4f \n",energy0,energy); 

  // need to retrieve final geometry, write to xyz
#if USE_AUX
  xyz_read_aux(inpfile_string);
#else
  xyz_read(inpfile_string);
#endif
  xyz_save(inpfile_string+".xyz");

  xyzfile.close();
  inpfile.close();

  if (abs(energy)<0.00001)
  {
    printf(" energy zero, mopac failed \n");
    return 10000;
  }

  return energy;
}


void Mopac::opt_write(string filename, ICoord icoords) {

#if 0
  printf(" WARNING: bypassing ICs! \n");
  return opt(filename);
#endif

  printf(" in mopac/opt_write() w/IC freeze ");
  printf(" filename: %s \n",filename.c_str());

  energy0 = energy = 0;

  ofstream inpfile;
  string inpfile_string = filename;
  ofstream xyzfile;
  string xyzfile_string = "scratch/testmopac.xyz";
  inpfile.open(inpfile_string.c_str());
  inpfile.setf(ios::fixed);
  inpfile.setf(ios::left);
  inpfile << setprecision(6);

  xyzfile.open(xyzfile_string.c_str());
  xyzfile.setf(ios::fixed);
  xyzfile.setf(ios::left);
  xyzfile << setprecision(6);
  xyzfile << " " << natoms << endl << endl;

  opt_header(inpfile);

  for (int i=0;i<natoms;i++)
  if (!skip[i])
  {
    if (!frzlistb[i])
      inpfile << " " << anames[i] << " " << xyz[3*i+0] << " 1 " << xyz[3*i+1] << " 1 " << xyz[3*i+2] << " 1 " << endl;
    else
    {
      inpfile << " " << anames[i] << " ";
      write_ic_input(inpfile,i,icoords);
    }
    xyzfile << " " << anames[i] << " " << xyz[3*i+0] << " " << xyz[3*i+1] << " " << xyz[3*i+2] << endl;
  }

  xyzfile.close();
  inpfile.close();

  return;
}


// standard opt
double Mopac::opt(string filename) 
{
  //printf(" in mopac/opt() \n");
#if 0
  printf(" frozen:");
  for (int i=0;i<natoms;i++)
    printf(" %i",frzlistb[i]);
  printf("\n");
#endif
 
  energy0 = energy = 0;

#if SKIPMOPAC
  printf(" skipping mopac opt! \n");
#endif

  ofstream inpfile;
  string inpfile_string = filename;
  ofstream xyzfile;
  string xyzfile_string = "scratch/testmopac.xyz";
#if !SKIPMOPAC
  inpfile.open(inpfile_string.c_str());
  inpfile.setf(ios::fixed);
  inpfile.setf(ios::left);
  inpfile << setprecision(6);

  xyzfile.open(xyzfile_string.c_str());
  xyzfile.setf(ios::fixed);
  xyzfile.setf(ios::left);
  xyzfile << setprecision(6);
  xyzfile << " " << natoms << endl << endl;

  opt_header(inpfile);

  for (int i=0;i<natoms;i++)
  if (!skip[i])
  {
    if (!frzlistb[i])
      inpfile << " " << anames[i] << " " << xyz[3*i+0] << " 1 " << xyz[3*i+1] << " 1 " << xyz[3*i+2] << " 1 " << endl;
    else
      inpfile << " " << anames[i] << " " << xyz[3*i+0] << " 0 " << xyz[3*i+1] << " 0 " << xyz[3*i+2] << " 0 " << endl;
    xyzfile << " " << anames[i] << " " << xyz[3*i+0] << " " << xyz[3*i+1] << " " << xyz[3*i+2] << endl;
  }
//  xyzfile << " " << natoms << endl << endl;
   string cmd = "/export/zimmerman/adewyer/mopac/MOPAC2016.exe "+filename;
//  string cmd = "/tmp/MOPAC2016.exe "+filename;
  system(cmd.c_str());
#endif

  energy = read_output(filename);
 
  //printf(" initial energy: %1.4f final energy: %1.4f \n",energy0,energy); 

  // need to retrieve final geometry, write to xyz
#if USE_AUX
  xyz_read_aux(inpfile_string);
#else
  xyz_read(inpfile_string);
#endif
  xyz_save(inpfile_string+".xyz");

  xyzfile.close();
  inpfile.close();

  if (abs(energy)<0.00001)
  {
    printf(" energy zero, mopac failed \n");
    return 10000;
  }

  return energy;
}

void Mopac::opt_write(string filename) {

   //printf(" in mopac/opt() \n");

  energy0 = energy = 0;

#if SKIPMOPAC
  printf(" skipping mopac opt! \n");
#endif

  ofstream inpfile;
  string inpfile_string = filename;
  ofstream xyzfile;
  string xyzfile_string = "scratch/testmopac.xyz";
#if !SKIPMOPAC
  inpfile.open(inpfile_string.c_str());
  inpfile.setf(ios::fixed);
  inpfile.setf(ios::left);
  inpfile << setprecision(6);

  xyzfile.open(xyzfile_string.c_str());
  xyzfile.setf(ios::fixed);
  xyzfile.setf(ios::left);
  xyzfile << setprecision(6);
  xyzfile << " " << natoms << endl << endl;

  opt_header(inpfile);

  for (int i=0;i<natoms;i++)
  if (!skip[i])
  {
    if (!frzlistb[i])
      inpfile << " " << anames[i] << " " << xyz[3*i+0] << " 1 " << xyz[3*i+1] << " 1 " << xyz[3*i+2] << " 1 " << endl;
    else
      inpfile << " " << anames[i] << " " << xyz[3*i+0] << " 0 " << xyz[3*i+1] << " 0 " << xyz[3*i+2] << " 0 " << endl;
    xyzfile << " " << anames[i] << " " << xyz[3*i+0] << " " << xyz[3*i+1] << " " << xyz[3*i+2] << endl;
  }
  xyzfile.close();
  inpfile.close();
#endif

  return;
}

double Mopac::read_output(string filename) {

  //double energy = -1;
  energy = 10000;

  string oname = filename+".out";
  ifstream output(oname.c_str(),ios::in);
  string line;
  vector<string> tok_line;
  while(!output.eof()) 
  { 
    getline(output,line);
//    cout << " RR " << line << endl;
    if (line.find("CYCLE:     1")!=string::npos)
    {
//      cout << "Initial: " << line << endl;
      tok_line = StringTools::tokenize(line, " \t");
//      energy0=atof(tok_line[10].c_str());
    }
    if (line.find("FINAL HEAT")!=string::npos)
    {
     // cout << "Final: " << line << endl;
      tok_line = StringTools::tokenize(line, " \t");
      energy=atof(tok_line[5].c_str());
    }
  }

  return energy;
}

void Mopac::alloc(int natoms_i) {
 
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

void Mopac::init(int natoms_i, int* anumbers_i, string* anames_i, double* xyz_i) {
 
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
  frzlist = new int[natoms_i];
  frzlistb = new int[natoms_i];

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

void Mopac::freemem() {

  delete [] xyz0;
  delete [] xyz;
  delete [] anumbers;
  delete [] anames;

  delete [] frzlist;
  delete [] frzlistb;

  delete [] skip;

  return;
}

void Mopac::freeze(int* frzlist_new, int nfrz_new, int nfrz0_new) {

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

void Mopac::freeze_d(int* frzlist_new) {

  for (int i=0;i<natoms;i++)
    frzlist[i] = frzlist_new[i];
  for (int i=0;i<natoms;i++)
    frzlistb[i] = frzlist_new[i];

#if 1
  printf(" freeze list: ");
  for (int i=0;i<natoms;i++)
    printf("%i ",frzlistb[i]);
  printf("\n");
#endif
  
  return;
}

// what is reset doing? //
void Mopac::reset(int natoms_i, int* anumbers_i, string* anames_i, double* xyz_i) {
 
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

void Mopac::set_charge(int c0)
{
  charge = c0;
  return;
}

void Mopac::xyz_save(string filename){

  ofstream xyzfile;
  string xyzfile_string = "xyzfile.txt";
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

void Mopac::xyz_read(string filename)
{
//  printf(" in xyz_read \n");
  string oname = filename+".out";
  ifstream output(oname.c_str(),ios::in);
  string line;
  vector<string> tok_line;
  int count = 0;
  int i = 0;

  while(!output.eof()) 
  { 
    getline(output,line);
//    cout << " RR " << line << endl;
    if (count == 2)
    {
      while (skip[i] && i<natoms)
        i++;
      if (i==natoms) break;
      if (!StringTools::cleanstring(line)) break;
      vector<string> tok_line = StringTools::tokenize(line, " \t");
      xyz[3*i+0]=atof(tok_line[2].c_str());
      xyz[3*i+1]=atof(tok_line[3].c_str());
      xyz[3*i+2]=atof(tok_line[4].c_str());
//      cout << tok_line[0] << " " << tok_line[1] << " " << tok_line[2] << " " << endl; 
      i++;
    }

    if (line.find("CARTESIAN COORDINATES")!=string::npos && count == 0)
    {
//      cout << "Initial: " << line << endl;
      count++;
    }
    else if (line.find("CARTESIAN COORDINATES")!=string::npos && count == 1)
    {
      count++;
      getline(output,line);
    }

  }

  output.close();
 
  return;
}   

void Mopac::xyz_read_aux(string filename){ 
   
  string oname = filename+".aux";
  //printf(" in xyz_read_aux: %s \n",oname.c_str());
  ifstream output(oname.c_str(),ios::in);
  string line;
  vector<string> tok_line;
  int count = 0;
  int i = 0;
  while(!output.eof()) 
  { 
    getline(output,line);
    //cout << " RR " << line << endl; fflush(stdout);
    if (line.find("ATOM_X_OPT")!=string::npos)
    {
      count++;
      getline(output,line);
    }
    if (count>0 && i<natoms)
    {
      //printf(" xyz_read1: skip[%i]: %i natoms: %i \n",i,skip[i],natoms);
      while (skip[i] && i<natoms)
        i++;
      if (i>=natoms) break;
      //printf(" xyz_read2: skip[%i]: %i natoms: %i \n",i,skip[i],natoms); fflush(stdout);

      if (!StringTools::cleanstring(line)) break;
      vector<string> tok_line = StringTools::tokenize(line, " \t");
      if (tok_line.size()<3) 
      {
        printf(" error: couldn't find XYZ in aux: %s  lines_read: %2i \n",filename.c_str(),i);
        cout << " RR: " << line << endl;
       // break;
        exit(1);
      }
      xyz[3*i+0]=atof(tok_line[0].c_str());
      xyz[3*i+1]=atof(tok_line[1].c_str());
      xyz[3*i+2]=atof(tok_line[2].c_str());
      //cout << tok_line[0] << " " << tok_line[1] << " " << tok_line[2] << " " << endl; 
      i++;
    }
    if (i>=natoms) break;

  }

  output.close();
  return;
}   


#if 0
//currently accepts data from geom1, but not actual class instantiation
void do_mp1(ICoord geom0, ICoord geom1, int i, int rem1, int mov1) {

  printf(" in do_mp1: %i %i %i \n",i,rem1,mov1);
    geom1.reset(geom0.natoms,geom0.anames,geom0.anumbers,geoms0.coords);
    geom1.ic_create();
    geom1.nbonds = geoms[id].nbonds;
//    geom1.bonds[rem1][0]=addsubmp[i][0];
    geom1.bonds[rem1][0]=i;
    geom1.bonds[rem1][1]=mov1;
    geom1.ic_create_nobonds(); // add to this, to do imptor if needed!!
    geom1.mm_init();
    geom1.update_ic();
//    geom1.print_ic();

  return;
}
#endif



