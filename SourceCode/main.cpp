#include "main.h"
void check_bonding(int nstruct, ICoord ic1, int natoms1, string* anames, int* anumbers, vector<double*> xyzalla, int* unique);
void procedure_1(string xyzfile, string targetfile, int rotation, int align, int constrain, string atoms);
void procedure_2(int nconf, string xyzfile, string targetfile, string atoms);


#define DIST_THRESH 0.2

#define USE_XTB 1

void opt_xtb(int charge, int natoms, string* anames, int* anumbers, vector<double*> xyzall, double* E, int type)
{
  int N3 = natoms*3;
  int nstruct = xyzall.size();

  int* frz = new int[natoms]();

  int* done = new int[nstruct]();

// #pragma omp parallel for
  for (int i=0;i<nstruct;i++)
  if (!done[i])
  {
    XTB xopt;
    xopt.alloc(natoms);
    xopt.set_charge(charge);
    xopt.freeze_d(frz);

    string nstr = StringTools::int2str(i,3,"0");
    string xfilename;
    if (type==1)
      xfilename = "xopt"+nstr+".xyz";
    else if (type==2)
      xfilename = "xoptb"+nstr+".xyz";
    else
      xfilename = "xoptc"+nstr+".xyz";

    xopt.reset(natoms,anumbers,anames,xyzall[i]);
    E[i] = xopt.opt_check(xfilename);
    for (int j=0;j<N3;j++)
      xyzall[i][j] = xopt.xyz[j];

    done[i] = 1;
    xopt.freemem();
  }

#if 1
  printf("  conformer energies:");
  for (int i=0;i<nstruct;i++)
    printf(" %5.2f",E[i]);
  printf("\n");
#endif

  delete [] done;
  delete [] frz;

  return;
}

void opt_mopac(int charge, int natoms, string* anames, int* anumbers, vector<double*> xyzall, double* E, int type)
{
  int N3 = natoms*3; //total coordinates
  int nstruct = xyzall.size(); //total number of structures
  printf("nstruct opt_mopac = %2i \n", nstruct);
  int* frz = new int[natoms];
  for (int i=0;i<natoms;i++)
  {
    frz[i] = 0;
  }
  Mopac mopt;
  mopt.alloc(natoms);
  mopt.set_charge(charge);
  // AD - Feb 16, 2017
  // pointer called freeez which is an array of size natoms
  // prints out a line of zeros in progress.log
  mopt.freeze_d(frz);

  int* done = new int[nstruct];
  for (int i=0;i<nstruct;i++)
  {
    done[i] = 0;
  }
  for (int i=0;i<nstruct;i++)
  if (!done[i])
  {
    // AD - Feb 16, 2017
    // what is the difference between these three types
    // 1 = mopt = ligand only
    // 2 = moptb = ligand + metal
    // 3 = ? 
    string nstr = StringTools::int2str(i,3,"0");
    string mfilename;
    if (type==1)
      mfilename = "scratch/mopt"+nstr+".xyz";
    else if (type==2)
      mfilename = "scratch/moptb"+nstr+".xyz";
    else
      mfilename = "scratch/moptc"+nstr+".xyz";
  
//    ICoord ic1,ic2;
//    ic1.init(natoms,anumbers,anames,xyzall[i]);
//    ic1.ic_create();

    mopt.reset(natoms,anumbers,anames,xyzall[i]);
    E[i] = mopt.opt_check(mfilename);
    for (int j=0;j<N3;j++)
      xyzall[i][j] = mopt.xyz[j];
      
//    ic2.init(natoms,anumbers,anames,xyzall[i]);
//    ic2.ic_create();

//    int intact = compare_ic(ic1,ic2);

    done[i] = 1;
  }

#if 1
  printf("  conformer energies:");
  for (int i=0;i<nstruct;i++)
    printf(" %5.2f",E[i]);
  printf("\n");
#endif

  delete [] done;
  delete [] frz;

  return;
}

int generate_conformers_and_opt(int nconf, string filename, double* &E, vector<double*> &xyzall, int constrain, string atoms)
{
//  string xyzfile = filename;
  
  printf("  generating %4i conformers of %s \n",nconf,filename.c_str());
  string confile = "all.xyz";

  //constrained on (AD 5-31-18)
  //run josh's code
  //need to read in frozen atom values instead of "1 2 3 4"
  
  if(constrain==1)
  {
    //constrained on
    string cmd = "./constrained_conf_gen.sh "+filename+" "+confile+" "+atoms;
    //cout << "cmd string: " << cmd << endl;
    system(cmd.c_str());
  }
  else
  {
    //constrained off
    string nconfstr = StringTools::int2str(nconf,1,"0");
    //added in --ecutoff line (AD 3-16-17)
    string cmd = "obabel -ixyz "+filename+" -oxyz -O "+confile+" --confab --conf "+nconfstr+ " --ecutoff 25 --rcutoff 0.5";
    system(cmd.c_str());
  }
   
  //appends original "test.xyz" file to the end of all.xyz so that it is sampled (AD 2/20/17)
  ifstream infile;
  fstream outfile;
  infile.open("test.xyz");
  outfile.open("all.xyz", ios::out | ios::app);
  while(infile.good())
  {
    string line;
    getline(infile, line);
    outfile << line << endl;
  }
  outfile.close();
  infile.close();  

  //CHARGE1 = test.xyz charge, ligand charge - structure that conformational search is done on.
  string cfilename = "CHARGE1";
  int charge = get_charge(cfilename);
  int natoms = get_natoms(filename);
  int N3 = natoms*3;
  printf("  there are %2i atoms \n",natoms);

  string* anames = new string[natoms];
  get_all_xyz(natoms,anames,xyzall,confile);
  
  int nstruct = xyzall.size();
  cout << "nstruct = " << nstruct << endl;
  int* anumbers = new int[natoms];
  for (int i=0;i<natoms;i++)
    anumbers[i] = PTable::atom_number(anames[i]);


  E = new double[nstruct];
  for (int i=0;i<nstruct;i++) E[i] = 0.;

//mopac opt all.xyz
//write structures to all2.xyz 
#if USE_XTB
  opt_xtb(charge,natoms,anames,anumbers,xyzall,E,1);
#else
  opt_mopac(charge,natoms,anames,anumbers,xyzall,E,1);
#endif
  write_all_xyz(natoms,anames,E,xyzall,"all2.xyz");

#if 0
  printf("   showing all structures \n");
  for (int i=0;i<nstruct;i++)
  {
    printf(" %2i \n\n",natoms);
    for (int j=0;j<natoms;j++)
      printf(" %2s %8.5f %8.5f %8.5f \n",anames[j].c_str(),xyzall[i][3*j+0],xyzall[i][3*j+1],xyzall[i][3*j+2]);
  }
#endif

  delete [] anumbers;
  delete [] anames;

  return nstruct;
}

int main(int argc, char* argv[])
{
  printf("\n\n in main() \n");
  string xyzfile = "test.xyz";
  string targetfile = "target.xyz";
  switch (argc)
  {
    case 2:
      //printf(" case 2 \n");
      xyzfile = argv[1];
      break;
    default:
      break;
  }
  
  // Read in constrained file (AD 5-31-18)
  int rotation = 0;
  int align = 0;
  int constrain = 0;
  int atom= 0; 
  string atoms; 
  string file="bcat_inp";
  ifstream infile;
  infile.open(file.c_str());
  if (!infile)
  {
    printf("Error opening %s \n", file.c_str());
    exit(-1);
  }

  //pass file to stringtools + grab lines
  string tag="BCat Info";
  bool found=StringTools::findstr(infile, tag);
  if (!found)
  {
    printf("Could not find tag for default info");
    exit(-1);
  }
  string line, templine, tagname;

  //following section parses file
  bool reading=true;
  while(reading)
  {
    reading=false;
    getline(infile, line);
    vector<string> tok_line = StringTools::tokenize(line, " ,\t");
    templine=StringTools::newCleanString(tok_line[0]);
    tagname=StringTools::trimRight(templine);
    if (tagname=="ALIGN")
    {
      align=atoi(tok_line[1].c_str());
      reading=true;
    }
    if (tagname=="ROTATION")
    {
      rotation=atoi(tok_line[1].c_str());
      reading=true;
    }
    if (tagname=="CONSTRAIN")
    {
      constrain=atoi(tok_line[1].c_str());
      reading=true;
    }
    if (tagname=="ATOMS")
    {
      reading=true;
      int count=0;
      string done="done";
      bool is_done = false;
      while (is_done == false)
      {
        string item=tok_line[count];
        atom=atoi(tok_line[count].c_str());
        if(item == done)
        {
          is_done = true;
        }
        else
        {
          if (count != 0)
          {
            atoms.append(item);
            atoms.append(" ");   
          }
        }
        count++; 
      }
    }
  }
  infile.close();

  int nconf = 250000;
  //cout << "PROCEDURE 1 INPUT" << endl;
  //cout << "*****************" << endl;
  //cout << "xyzfile: " << xyzfile << " is of type " << typeid(xyzfile).name() <<  endl;
  //cout << "targetfile: " << targetfile << " is of type " << typeid(targetfile).name() <<  endl;
  //cout << "align: " << align << " is of type " << typeid(align).name() << endl;
  //cout << "rotation: " << rotation << " is of type " << typeid(rotation).name() << endl;
  //cout << "constrain: " << constrain << " is of type " << typeid(constrain).name() << endl;
  //cout << "atoms: " << atoms << " is of type " << typeid(atoms).name() << endl;
  //cout << "*****************" << endl;
  //cout << "END OF INPUT" << endl;
  procedure_1(xyzfile,targetfile,rotation,align,constrain,atoms);
  
  //turned off procedure 2 (AD 3-9-17)
  //procedure_2(nconf,xyzfile,targetfile);

  return 0;
}

void procedure_2(int nconf, string xyzfile, string targetfile, string atoms)
{
  // AD - Feb 16, 2017
  // takes xyz file and runs conformer generation, and inital mopac opt
  ICoord ic1; 
  ic1.init(xyzfile);
  printf(" initial bonds \n");
  ic1.print_bonds();
  int constrain=0;
  vector<double*> xyzall;
  double* E;
  double* E2;
  // AD - Feb 16, 2017
  //why is nconf here so large in comparison to nconf for procedure 1?
  int nstruct = generate_conformers_and_opt(nconf,xyzfile,E,xyzall,constrain, atoms);

  return;
}

void procedure_1(string xyzfile, string targetfile, int rotation, int align, int constrain, string atoms)
{

  cout << endl;
  cout << "PROCEDURE 1 VARIABLE" << endl;
  cout << "********************" << endl;
  cout << "xyz file: " <<  xyzfile << endl;
  cout << "target file: " << targetfile << endl;
  cout << "rotation: " << rotation << endl;
  cout << "align: " << align << endl;
  cout << "constrain: " << constrain << endl;
  cout << "atoms: " << atoms << endl;
  cout << "*******************" << endl;
  cout << "END OF VARIABLES" << endl;

  // AD - Feb 16, 2017
  // Runs gsm to push ligand and metal center together
  // sums up charges to get overall charge 
  ICoord ic1; 
  ic1.init(xyzfile);
  printf(" initial bonds \n");
  ic1.print_bonds();

  vector<double*> xyzall; //structures from conf search/opt
  double* E;
  double* E2; 
  int nstruct = generate_conformers_and_opt(2000,xyzfile,E,xyzall,constrain,atoms);
  
  string cfilename = "CHARGE1";
  int charge1 = get_charge(cfilename);
  cfilename = "CHARGE2";
  int charge2 = get_charge(cfilename);
  int charget = charge1 + charge2;

  int natoms1 = get_natoms(xyzfile);
  int natoms2 = get_natoms(targetfile);
  string* anames = new string[natoms1];
  string* anamesm = new string[natoms2];
  int* anumbers = new int[natoms1];
  int* anumbersm = new int[natoms2];
  double* xyz0 = new double[3*natoms1]; //xyz of test.xyz
  double* xyzm = new double[3*natoms2]; //xyz of target.xyz
  xyz_read(natoms1,anames,xyz0,xyzfile);
  xyz_read(natoms2,anamesm,xyzm,targetfile);
  for (int i=0;i<natoms1;i++)
    anumbers[i] = PTable::atom_number(anames[i]);
  for (int i=0;i<natoms2;i++)
    anumbersm[i] = PTable::atom_number(anamesm[i]);

  int natomst = natoms1 + natoms2;
  int N3t = natomst*3;
  string* anamest = new string[natomst];
  int* anumberst = new int[natomst];
  for (int i=0;i<natoms1;i++)
    anamest[i] = anames[i];
  for (int i=0;i<natoms2;i++)
    anamest[natoms1+i] = anamesm[i];
  for (int i=0;i<natomst;i++)
    anumberst[i] = PTable::atom_number(anamest[i]);

  printf("  target structure: \n");
  for (int i=0;i<natoms2;i++)
    printf(" %2s %8.5f %8.5f %8.5f \n",anamesm[i].c_str(),xyzm[3*i+0],xyzm[3*i+1],xyzm[3*i+2]);

  
  int* unique = new int[nstruct];
  // AD - Feb 16, 2017
  // checks if structures are unique by reading in array of structures - unique is a pointer towards an array of structures
  // reads in all2.xyz to check unique structures
  get_unique_conf(nstruct,unique);

  cout << "pre align_and_opt nstruct: " << nstruct << endl;
  //reads in unique structures, and all structures generated as well as metal center to align structures & optimizes them
  align_and_opt(natoms1,natoms2,anames,anamesm,anamest,anumbers,anumbersm,charget,nstruct,unique,xyzall,xyzm,align,rotation);
  
  int n;
  if(align!=0)
  { 
    if(align==1)
    {
      nstruct=nstruct*6;
    }
  }
  if(rotation!=0)
  {
    n=(rotation+1);
    cout << "n post rotation = " << n << endl;
    nstruct=nstruct*n;
    cout << "nstruct post rotation = " << nstruct << endl;
  }
  cout << "post align_and_opt nstruct: " << nstruct << endl;
  //retrieve firstnode.xyz files from GSM
  vector<double*> xyzalla;
  for (int i=0;i<nstruct;i++)
  {
    double* xyz1 = new double[3*natomst];
    xyzalla.push_back(xyz1);
  }
  for (int i=0;i<nstruct;i++)
  {
    string nstr = StringTools::int2str(i,4,"0");
    string gsmfilename = "scratch/firstnode.xyz"+nstr;
    xyz_read_last(natomst,xyzalla[i],gsmfilename);
  }

//prints out a x value of first atom from final node of firstnode.xyz
//** FIGURE OUT HOW TO PRINT ALL POINTS**//
//  for (int i=0;i<nstruct;i++)
//  {
//    printf("Printing out xyzalla[%d] \n", i);
//    {
//      for (int j=0;j<N3t;j++)
//      {
//        if (j%3 == 0)
//        {
//          cout << endl;
//        }
//        cout << xyzalla[i][j] << " ";
//        if (j == (N3t-1))
//        {
//          cout << endl << endl;
//        }
//      }
//    }
//  }
  //cout << "E: " << E << endl;

  if((rotation!=0)||(align!=0))
  {
    E2 = new double[nstruct]; 
#if USE_XTB
    opt_xtb(charget,natomst,anamest,anumberst,xyzalla,E2,2);    
#else    
    opt_mopac(charget,natomst,anamest,anumberst,xyzalla,E2,2);
#endif    
    write_all_xyz(natomst,anamest,E2,xyzalla,"all4.xyz");
  }
  else
  {
#if USE_XTB
    opt_xtb(charget,natomst,anamest,anumberst,xyzalla,E,2);
#else
    opt_mopac(charget,natomst,anamest,anumberst,xyzalla,E,2);
#endif
    write_all_xyz(natomst,anamest,E,xyzalla,"all4.xyz");
  }
  printf(" done generating complexes \n");

  int* unique2 = new int[nstruct];
  
  check_bonding(nstruct,ic1,natoms1,anames,anumbers,xyzalla,unique2);
  printf("\n");

  // AD - Feb 16, 2017
  // Section of code chooses lowest E structures to add to all5.xyz
  // Not sure how it decides what the energy cutoff for strctures is
  // ** How does it decide?
  vector<pair<double,int> > bestc;
  pair<double,int> newc;
  pair<double,int> newc1;
  if((rotation!=0)||(align!=0))
  {
    printf("  conformer energies E2:");
    for (int i=0;i<nstruct;i++)
      printf(" %5.2f",E2[i]);
      printf("\n");
  }
  else
  {
    printf("  conformer energies E:");
    for (int i=0;i<nstruct;i++)
      printf(" %5.2f",E[i]);
      printf("\n");
  }
 
//  cout << "  n struct post check bonding: " << nstruct << endl; 
  
  int u=0;
  while(unique2[u])
  {
//    cout << "u= " << u << endl;
    if((rotation!=0)||(align!=0))
    {
      pair<double,int> newc1(E2[u],u);
      bestc.push_back(newc1);
    } 
    else
    {
      pair<double,int> newc(E[u],u);
      bestc.push_back(newc);
    }
    u++;
  }  
  sort(bestc.begin(),bestc.end());
  int nstructf = bestc.size();
  cout << "nstructf = " << nstructf << endl;
  printf("  best structures: \n");
  for (int i=0;i<nstructf;i++)
  {
    printf("   %2i: %8.5f \n",bestc[i].second+1,bestc[i].first);
  }

  double* E1 = new double[nstructf];
  double** xyzb = new double*[nstructf];
  for (int i=0;i<nstructf;i++)
  {
    xyzb[i] = new double[N3t];
  }
  for (int i=0;i<nstructf;i++)
  {
    int i1 = bestc[i].second;
    if((rotation!=0)||(align!=0))
    {
      E1[i] = E2[i1];
    }
    else
    {
      E1[i] = E[i1];
    }
    double* xyz1 = xyzalla[i1];
    for (int j=0;j<N3t;j++)
      xyzb[i][j] = xyz1[j];
  }

  write_all_xyz(natomst,anamest,nstructf,E1,xyzb,"all5.xyz");

  return;
}


void check_bonding(int nstruct, ICoord ic1, int natoms1, string* anames, int* anumbers, vector<double*> xyzalla, int* unique)
{
  printf("\n  now checking bonding \n");
  for (int i=0;i<nstruct;i++)
  {
    ICoord ic2;
    ic2.isOpt = 0;
    ic2.init(natoms1,anames,anumbers,xyzalla[i]);
    //ic2.print_bonds();

    if (ic2.nbonds!=ic1.nbonds)
    {
      printf("  struct %2i bonds changed, eliminating \n",i+1);
      unique[i] = 0;
    }
    else
    {
     //here check that bonding is the same
      int nf = 0;
      //for (int j=0;j<ic1.nbonds;j++)
    }
  }

  return;
}


int read_adds(int* adds, string addfile)
{
  int nadd = 0;

  ifstream infile;
  infile.open(addfile.c_str());
  if (!infile){
    printf(" Error: couldn't open file: %s \n",addfile.c_str());
    exit(-1);
  } 
  
  string line;
  while (!infile.eof())
  {
    getline(infile,line);
    vector<string> tok_line = StringTools::tokenize(line, " \t");
    if (tok_line.size()<2) break;

    adds[2*nadd+0] = atoi(tok_line[0].c_str())-1;
    adds[2*nadd+1] = atoi(tok_line[1].c_str())-1;
    //printf(" ADD found: %2i %2i \n",adds[2*nadd+0]+1,adds[2*nadd+1]+1);
    nadd++;
  }

  return nadd;
}


void align_and_opt(int natoms1, int natoms2, string* anames, string* anamesm, string* anamest, int* anumbers, int* anumbersm, int charget, int nstruct, int* unique, vector<double*> xyzall, double* xyzm, int align, int rotation)
{ 
  int nstruct_new=nstruct; //structures for multiple alignments and rotations
  int s=0; //counter for multiple rotations and alignments
  int natomst = natoms1 + natoms2;
  int* adds = new int[8]; //no more than 8 coord
  int nadd = read_adds(adds,"ADD");
  for (int i=0;i<nadd;i++)
  {
    adds[2*i+1] += natoms1;
    printf(" ADD found: %2i %2i \n",adds[2*i+0]+1,adds[2*i+1]+1);
  }
  cout << "nstruct align_opt: " << nstruct << endl;
  int N3t = (natoms1+natoms2)*3;
  if(align!=0)
  {
    if(align==1)
    {
      nstruct_new=nstruct_new*6; 
    }
  }
  if(rotation!=0)
  {
    int ns = rotation+1;
    nstruct_new=nstruct_new*ns;
  }
  //cout << "nstruct align_opt: " << nstruct << endl;
  //cout << "nstruct align_opt_2: " << nstruct_new << endl;
  double** xyzalign = new double*[nstruct_new];
  for (int i=0;i<nstruct_new;i++)
    xyzalign[i] = new double[N3t];
  
  printf("  now creating initial.xyz \n");
  for (int i=0;i<nstruct;i++)
  {
    printf(" nstruct = %2i", nstruct);
    printf(" working on structure %2i \n",i+1);
    Align align1;
    align1.init(natoms1,anames,anumbers,xyzall[i],natoms2,anamesm,anumbersm,xyzm);

//adewyer 5/3/18 - rotations
//adewyer 8-8-18 - alignments
    cout << "align = " << align << endl;
    if(align!=0)
    { 
      if(align==1)
      {
        for(int malign=0;malign<6;malign++)
        {  
          for(int mrot=0;mrot<=rotation;mrot++)
          {
            cout << "mrot = " << mrot << endl;
            cout << "rotation = " << rotation << endl;
            align1.add_align(nadd,adds,mrot,malign,rotation);
            for (int j=0;j<N3t;j++)
              xyzalign[s][j] = align1.xyza[j];
            printf("s in main.cpp = %2i \n", s);
            s++;
          }
        }
      }
    } 
    else
    {
      for(int mrot=0;mrot<=rotation;mrot++)
      {
        cout << "mrot = " << mrot << endl;
        cout << "rotation = " << rotation << endl;
        align1.add_align(nadd,adds,mrot,align,rotation);
        for (int j=0;j<N3t;j++)
        {
          xyzalign[s][j] = align1.xyza[j];
        }
        printf("s in main.cpp = %2i \n", s);
        s++;
      }
    }
  }  
  printf("\n");
  if(rotation!=0 || align!=0) 
  {
    write_all_xyz(natomst,anamest,nstruct_new,NULL,xyzalign,"all3.xyz");
    write_gsm(natomst,anamest,charget,nstruct_new,NULL,xyzalign,nadd,adds);
  }
  else
  {
    write_all_xyz(natomst,anamest,nstruct,NULL,xyzalign,"all3.xyz");
    write_gsm(natomst,anamest,charget,nstruct,NULL,xyzalign,nadd,adds);
  }
    
#if 0
  printf(" printing aligned structures \n");
  for (int i=0;i<nstruct;i++)
  {
    printf(" %2i \n \n",natomst);
    for (int j=0;j<natomst;j++)
      printf(" %2s %8.5f %8.5f %8.5f \n",anamest[j].c_str(),xyzalign[i][3*j+0],xyzalign[i][3*j+1],xyzalign[i][3*j+2]);
  }
#endif
  if((rotation!=0)||(align!=0))
  {
    nstruct=nstruct_new;
  }
  do_gsm(nstruct);

//moved delete lines outside of for loop so that xyzalign does not get erased
  for (int i=0;i<nstruct;i++)
    delete [] xyzalign[i];
  delete [] xyzalign;
  return;
}

void do_gsm(int nstruct)
{
//  printf("  skipping GSM step \n");
//  return;

  string cmd;
  for (int i=0;i<nstruct;i++)
  {
    string nstr = StringTools::int2str(i,4,"0");
    string istr = StringTools::int2str(i,1,"0");
    string filename = "scratch/firstnode.xyz"+nstr;
    struct stat sts;
    if (stat(filename.c_str(), &sts) != -1)
    {
      printf("  string %2i already done \n",i);
    }
    else
    {
      cmd = "./gfstringq.exe "+istr+" > scratch/paragsm"+nstr;
      system(cmd.c_str());
    }
  }
}

void get_all_xyz(int natoms, string* anames, vector<double*> &xyzs, string xyzfile)
{    
  ifstream infile;
  infile.open(xyzfile.c_str());
  if (!infile){
    printf(" Error: couldn't open XYZ file: %s \n",xyzfile.c_str());
    exit(-1);
  } 
  
  string line;
  int nf = 0;
  int nf2 = 0;
  while (!infile.eof())
  { 
    int nf2 = 0;
    getline(infile, line);
    int natoms1 = atoi(line.c_str());
    if (natoms1==natoms)
    {
      double* coords = new double[3*natoms];
      getline(infile, line);
  
      for (int i=0;i<natoms;i++)
      {
        getline(infile, line);
        int length=StringTools::cleanstring(line);
        vector<string> tok_line = StringTools::tokenize(line, " \t");
        if (nf==0)
          anames[i]=tok_line[0];
        coords[3*i+0]=atof(tok_line[1].c_str());
        coords[3*i+1]=atof(tok_line[2].c_str());
        coords[3*i+2]=atof(tok_line[3].c_str());
      }
      xyzs.push_back(coords);
      nf++;
      nf2++;
    } //if adding new geom
    else
    {
      cout << "nf2 = " << nf2 << endl;
      printf(" done reading after %2i structures \n",nf);
      break;
    }
  }
  infile.close();
 
  return;
}   

void xyz_read(int natoms, string* anames, double* coords, string xyzfile)
{  
  ifstream infile;
  infile.open(xyzfile.c_str());
  if (!infile){
    printf(" Error: couldn't open XYZ file \n");
    exit(-1);
  } 
  
  string line;
  bool success=true;
  success=getline(infile, line);
  if (success){
    int length=StringTools::cleanstring(line);
    natoms=atoi(line.c_str());
  }
  printf(" natoms: %i \n",natoms);
  
  success=getline(infile, line);
  
  //cout <<"  -Reading the atomic names...";
  for (int i=0;i<natoms;i++){
    success=getline(infile, line);
    int length=StringTools::cleanstring(line);
    vector<string> tok_line = StringTools::tokenize(line, " \t");
    anames[i]=tok_line[0];
    coords[3*i+0]=atof(tok_line[1].c_str());
    coords[3*i+1]=atof(tok_line[2].c_str());
    coords[3*i+2]=atof(tok_line[3].c_str());
  }
  
  infile.close();

#if 0
  printf(" XYZ: \n");
  for (int i=0;i<natoms;i++)
    printf(" %s %8.6f %8.6f %8.6f \n",anames[i].c_str(),coords[3*i+0],coords[3*i+1],coords[3*i+2]);
#endif

  printf(" done reading XYZ \n"); fflush(stdout);
 
  return;
}   

void xyz_read_last(int natoms, double* coords, string xyzfile)
{ 
  ifstream infile;
  infile.open(xyzfile.c_str());
  if (!infile){
    printf(" Error: couldn't open XYZ file: %s \n",xyzfile.c_str());
    exit(-1);
  } 
  
  string line;
  int nf = 0;
  while (!infile.eof())
  {
    getline(infile, line);
    int natoms1 = atoi(line.c_str());
    if (natoms1==natoms)
    {
      getline(infile, line);  
      for (int i=0;i<natoms;i++)
      {
        getline(infile, line);
        int length=StringTools::cleanstring(line);
        vector<string> tok_line = StringTools::tokenize(line, " \t");
        coords[3*i+0]=atof(tok_line[1].c_str());
        coords[3*i+1]=atof(tok_line[2].c_str());
        coords[3*i+2]=atof(tok_line[3].c_str());
      }
      nf++;
    } //if adding new geom
    else
    {
      printf(" done reading after %2i structures \n",nf);
      break;
    }
  }
  infile.close();
 
  return;
}   

int get_charge(string filename)
{
  return get_natoms(filename);
}

int get_natoms(string filename)
{
  ifstream infile;
  infile.open(filename.c_str());
  if (!infile){
    printf("  couldn't find file %s \n",filename.c_str());
    exit(-1);
  }

  string line;
  getline(infile, line);
  int length=StringTools::cleanstring(line);
  int natoms = atoi(line.c_str());

  infile.close();

  return natoms;
}

int get_unique_conf(int nstruct, int* unique)
{
  int ns2 = nstruct * nstruct;
  int* similar = new int[ns2];
  for (int i=0;i<ns2;i++) similar[i] = 0;

  OBConversion obc1,obc2;
  obc1.SetInFormat("xyz");
  obc2.SetInFormat("xyz");
  OBMol mol1;
  OBMol mol2;

  int nf1 = 0;
  bool notatend1 = obc1.ReadFile(&mol1,"all2.xyz");
  while (notatend1)
  {
    //printf("  Molecular Weight: %5.3f \n",mol1.GetMolWt());

    int nf2 = 0;
    bool notatend2 = obc2.ReadFile(&mol2,"all2.xyz");
    while (notatend2)
    {
      OBAlign align1 = OBAlign(mol1,mol2);
      align1.Align();
      double dist = align1.GetRMSD();
      //printf(" difference (%2i-%2i): %10.6f \n",nf1,nf2,dist);
      if (dist<DIST_THRESH) similar[nf1*nstruct+nf2] = 1;

      mol2.Clear();
      notatend2 = obc2.Read(&mol2);
      nf2++;
    }
    mol1.Clear();
    notatend1 = obc1.Read(&mol1);
    nf1++;
  }
  
  for (int i=0;i<nstruct;i++) unique[i] = 1;
  for (int i=0;i<nstruct;i++)
  for (int j=0;j<i;j++)
  if (similar[i*nstruct+j])
    unique[j] = 0;

  int nf = 0;
  for (int i=0;i<nstruct;i++)
  if (unique[i])
    nf++;
  
  printf("\n similarity matrix: \n");
  for (int i=0;i<nstruct;i++)
  {
    for (int j=0;j<nstruct;j++)
      printf(" %i",similar[i*nstruct+j]);
    printf("\n");
  }
  printf(" unique list:");
  for (int i=0;i<nstruct;i++)
    printf(" %i",unique[i]);
  printf("\n\n");

  return nf;
}

void write_all_xyz(int natoms, string* anames, double* E, vector<double*> xyzs, string xyzfile_string)
{
  int nstruct = xyzs.size();
  printf("nstruct = %i", nstruct);
  double** xyzs1 = &xyzs[0];
  write_all_xyz(natoms,anames,nstruct,E,xyzs1,xyzfile_string);
  return;
}

void write_all_xyz(int natoms, string* anames, int nstruct, double* E, double** xyzs, string xyzfile_string)
{
  ofstream xyzfile;
  xyzfile.open(xyzfile_string.c_str());
  xyzfile.setf(ios::fixed);
  xyzfile.setf(ios::left);
  xyzfile << setprecision(6);

  char* sbuff = new char[10000];
  for (int i=0;i<nstruct;i++)
  {
    if (E==NULL)
      xyzfile << natoms << endl << endl;
    else
      xyzfile << natoms << endl << E[i] << endl;
    for (int j=0;j<natoms;j++)
    {
      sprintf(sbuff," %2s %10.6f %10.6f %10.6f \n",anames[j].c_str(),xyzs[i][3*j+0],xyzs[i][3*j+1],xyzs[i][3*j+2]);
      xyzfile << sbuff;
    }
  }

  xyzfile.close();

  delete [] sbuff;

  printf(" done writing \n");

  return;
}

void write_gsm(int natoms, string* anames, int charge, int nstruct, double* E, double** xyzs, int nadd, int* adds)
{
  char* sbuff = new char[10000];

  for (int i=0;i<nstruct;i++)
  {
    string nstr = StringTools::int2str(i,4,"0");

    string xyzfile_string = "scratch/initial"+nstr+".xyz";
    ofstream xyzfile;
    xyzfile.open(xyzfile_string.c_str());
    xyzfile.setf(ios::fixed);
    xyzfile.setf(ios::left);
    xyzfile << setprecision(8);

    xyzfile << natoms << endl << charge << endl;
    for (int j=0;j<natoms;j++)
    {
      sprintf(sbuff," %2s %10.6f %10.6f %10.6f \n",anames[j].c_str(),xyzs[i][3*j+0],xyzs[i][3*j+1],xyzs[i][3*j+2]);
      xyzfile << sbuff;
    }
    xyzfile.close();

    string isofile_string = "scratch/ISOMERS"+nstr;
    ofstream isofile;
    isofile.open(isofile_string.c_str());
    isofile << "NEW" << endl;
    for (int j=0;j<nadd;j++)
      isofile << "BOND " << adds[2*j+0]+1 << " " << adds[2*j+1]+1 << endl;
    isofile.close();

    double dist1 = 2.0;
    string forcefile_string = "scratch/FORCE"+nstr;
    ofstream forcefile;
    forcefile.open(forcefile_string.c_str());
    for (int j=0;j<nadd;j++)
      forcefile << adds[2*j+0]+1 << " " << adds[2*j+1]+1 << " " << dist1 << " 0.15 " << endl;
    forcefile.close();
  }


  delete [] sbuff;

  printf(" done writing \n");

  return;
}


