#ifndef DESCRIBE_H_
#define DESCRIBE_H_

#include "describe.hpp"

// Describe constructor
Describe::Describe(string nc_fname, string parm_fname)
{
  ReadNetCDF(nc_fname);
//  cout << "netcdf success" << endl;
  ReadParm(parm_fname);
//  cout << "parm sucess" << endl;
  GetAtomIndices();
}

// parser for the netCDF file
void Describe::ReadNetCDF(string fname)
{
  // instantiates NcFile class
  NcFile dataFile(fname, NcFile::read);

  // save coordinate data to NcVar class
  NcVar data = dataFile.getVar("coordinates");

  // vector of NcDim classes of the data
  vector<NcDim> data_dims = data.getDims();

  // vector of dimensions of the data
  vector<int> dims;

  // populating dimensions vector
  // also calculating the length of the array required
  int data_size = 1;
  for(auto i : data_dims)
  {
      int dim_size = i.getSize();
      dims.push_back(dim_size);
      data_size *= dim_size;
  }

  // instantiating float array of size "data_size" where the data is held
  float *coords;
  coords = new float [data_size];

  // data is saved as a 1-D array
  data.getVar(coords);

  // converting the 1-D array to its proper shape as a vector (number of frames, number of atoms, xyz coordinates)
  vector<vector<vector<float>>> coord_vector;
  for(unsigned int i=0; i<dims[0]; i++)
  {
    vector<vector<float>> frame;
    for(unsigned int j=0; j<dims[1]; j++)
    {
      vector<float> xyz;
      for(unsigned int k=0; k<dims[2]; k++)
      {
        unsigned int current_ind = (i * dims[1] * dims[2]) + (j * dims[2]) + k;
        xyz.push_back(coords[current_ind]);
      }
      frame.push_back(xyz);
    }
    coord_vector.push_back(frame);
  }

  n_frames = dims[0];
  n_atoms = dims[1];

  // deallocating the long 1-D array
  delete[] coords;

  // saving the large vector
  trajectory_data.coords = coord_vector;
}

// parser for the parm7 file
void Describe::ReadParm(string fname)
{
  ifstream infile;
  infile.open(fname);

  // exits if parm file is not found
  if(!infile)
  {
    cout << "Error: Unable to open 'parm7' file." << endl;
    exit(1);
  }

  string buff; 
  int buff_int;

  // ignores the first four lines
  for(unsigned int i=0; i<4; i++)
  {
    getline(infile, buff);
  }
  
  stringstream ss;

  // saves FLAG POINTERS
  string pointers_str = ParmReadBlock(&infile, 32);
  ss << pointers_str;
  ss >> buff_int;
  if(buff_int != n_atoms)
  {
    cout << buff_int << " not equal to " << n_atoms << endl;
    exit(1); 
  }
  ss >> n_types >> nbonh >> mbona >> ntheth >> mtheta;

  for(unsigned int i=0; i<5; i++)
  {
    ss >> buff_int;
  }

  ss >> n_residues;

  for(unsigned int i=0; i<3; i++)
  {
    ss >> buff_int;
  }

  ss >> n_bonds;
  ss >> n_angles;
  ss >> n_torsions;
  ss >> n_solty;

  for(unsigned int i=0; i<13; i++)
  {
    ss >> buff_int;
  }

  // saves ATOM NAMES
  string atom_str = ParmReadBlock(&infile, n_atoms);

  // saves atom names to proper data type
  for(unsigned int i=0; i<n_atoms; i++)
  {
    buff = atom_str.substr(0,4);
    atom_str = atom_str.substr(4);
    trajectory_data.atoms.push_back(buff);
  }

  // saves CHARGES
  string charge_str = ParmReadBlock(&infile, n_atoms);

  // saves charges to proper data type
  for(unsigned int i=0; i<n_atoms; i++)
  {
    string buff_string = charge_str.substr(0,16);
    charge_str.erase(0,16);
    float buff_float = stof(buff_string);
    trajectory_data.charges.push_back(buff_float / 18.223);
  }

  // ignores ATOMIC NUMBERS
  buff = ParmReadBlock(&infile, n_atoms);

  // save ATOM MASSES
  string mass_str = ParmReadBlock(&infile, n_atoms);

  // saves masses to proper data type
  ss << mass_str;
  for(unsigned int i=0; i<n_atoms; i++)
  {
    float buff_float;
    ss >> buff_float;
    trajectory_data.masses.push_back(buff_float);
  }

  // clears stringstream for reassignment
  ss.str("");
  ss.clear();

  // ignores ATOM TYPE INDICES
  buff = ParmReadBlock(&infile, n_atoms);
  ss << buff;
  for(unsigned int i=0; i<n_atoms; i++)
  {
    ss >> buff_int;
    if(buff_int > n_types)
    {
      n_types = buff_int;
    }
  }

  // clears stringstream for reassignment
  ss.str("");
  ss.clear();

  // ignores NUMBER OF EXCLUDED ATOMS
  buff = ParmReadBlock(&infile, n_atoms);

  // ignores NONBONDED PARM INDICES
  buff = ParmReadBlock(&infile, n_types*n_types);

  // saves RESIDUE NAMES
  buff = ParmReadBlock(&infile, n_residues);
  ss << buff;
  for(unsigned int i=0; i<n_residues; i++)
  {
    string res_name = buff.substr(0,4);
    buff = buff.substr(4);
    trajectory_data.residues.push_back(res_name);
  }

  trajectory_data.residues.front() = 'N' + trajectory_data.residues.front().substr(0,3);

  // clears stringstream for reassignment
  ss.str("");
  ss.clear();

  // saves RESIDUE POINTERS
  buff = ParmReadBlock(&infile, n_residues);
  ss << buff;
  for(unsigned int i=0; i<n_residues; i++)
  {
    int index;
    ss >> index;
    residue_pointers.push_back(index);
  }

  // clears stringstream for reassignment
  ss.str("");
  ss.clear();

  // ignores BOND FORCE CONSANTS
  buff = ParmReadBlock(&infile, n_bonds);

  // ignores BOND EQUILIBRIUM VALUES
  buff = ParmReadBlock(&infile, n_bonds);

  // ignores ANGLE FORCE CONSTANTS
  buff = ParmReadBlock(&infile, n_angles);

  // ignores ANGLE EQUILIBRIUM VALUES
  buff = ParmReadBlock(&infile, n_angles);

  // ignores DIHEDRAL FORCE CONSTANTS
  buff = ParmReadBlock(&infile, n_torsions);

  // ignores DIHEDRAL PERIODICITY
  buff = ParmReadBlock(&infile, n_torsions);

  // ignores DIHEDRAL PHASE
  buff = ParmReadBlock(&infile, n_torsions);

  // ignores SCEE SCALE FACTOR
  buff = ParmReadBlock(&infile, n_torsions);

  // ignores SCNB SCALE FACTOR
  buff = ParmReadBlock(&infile, n_torsions);

  // ignores SOLTY
  buff = ParmReadBlock(&infile, n_solty);

  buff_int = (n_types * (n_types+1)) / 2;

  // ignores LJ A COEFFICIENTS
  buff = ParmReadBlock(&infile, buff_int);

  // ignores LJ B COEFFICIENTS
  buff = ParmReadBlock(&infile, buff_int);

  // saves BONDS WITH HYDROGEN
  buff = ParmReadBlock(&infile, 3*nbonh);
  ss << buff;

  for(unsigned int i=0; i<nbonh; i++)
  {
    int atom1, atom2;
    ss >> atom1 >> atom2 >> buff_int;
    atom1 = atom1 / 3;
    atom2 = atom2 / 3;
    vector<int> buff_intv;
    if(trajectory_data.atoms[atom2][0] == 'H')
    {
      buff_intv.push_back(atom1);
      buff_intv.push_back(atom2);
    }
    else
    {
      buff_intv.push_back(atom2);
      buff_intv.push_back(atom1);
    }
    bonh.push_back(buff_intv);
  }

  // THE REST OF THE PARM FILE IS IGNORED

  // closes filestream
  infile.close();
}

// reads a block of data and returns the data in a string
string Describe::ParmReadBlock(ifstream* parm_stream, int divisor)
{
  int n_columns, buff_int;
  char buff_char;
  string dum, buff, block;

  getline(*parm_stream, buff);
  *parm_stream >> dum;

  size_t pos = dum.find('(');
  dum = dum.substr(pos+1);

  stringstream ss(dum);
  ss >> n_columns >> buff_char >> buff_int;
  getline(*parm_stream, buff);

  int n_iter;
  if(divisor % n_columns > 0)
  {
    n_iter = (divisor / n_columns) + 1;
  }
  else
  {
    n_iter = (divisor / n_columns);
  }

  for(unsigned int i=0; i<n_iter; i++)
  {
    getline(*parm_stream, buff);
    block += buff;
  }
  return block;
}

// creates arrays of atom indices for each type of atom
void Describe::GetAtomIndices(void)
{
  for(unsigned int i=0; i<n_atoms; i++)
  {
    char buff_char = trajectory_data.atoms[i][0];
    if(buff_char == 'C')
    {
      Cs.push_back(i);
    }
    else if(buff_char == 'N')
    {
      Ns.push_back(i);
    }
    else if(buff_char == 'O')
    {
      Os.push_back(i);
    }
    else if(buff_char == 'H')
    {
      Hs.push_back(i);
    }
    else if(buff_char == 'S')
    {
      Ss.push_back(i);
    }
  }
}

// finds and counts hydrogen bonds in each frame
void Describe::HydrogenBonds(string num_hbond_fname, string hbond_fname)
{
//  vector<vector<int>> NHs, OHs;
  int n_hbonds = 0;
  for(unsigned int i=0; i<bonh.size(); i++)
  {
    if(IsIn(bonh[i][0], Ns))
    {
      NHs.push_back(bonh[i]);
    }
    else if(IsIn(bonh[i][0], Os))
    {
      OHs.push_back(bonh[i]);
    }
  }

  num_hbond_out.open(num_hbond_fname);
  num_hbond_out << "# Hydrogen Bonds" << endl;

  hbond_out.open(hbond_fname);
  hbond_out << "# Hydrogen Bonds" << endl;

  num_hbond_out << setw(10) << left << "#Frame" << '\t' << left << "Count" << endl;

  for(unsigned int frame=0; frame<n_frames; frame++)
  {
    hbond_out << "#Frame " << frame+1 << endl;
    int NHO_count = CountHBonds(NHs, Os, frame);
    int OHN_count = CountHBonds(OHs, Ns, frame);
    int frame_count = NHO_count + OHN_count;
    hbond_out << endl;
    num_hbond_out << setw(10) << left << frame+1 << '\t' << left << frame_count << endl;
  }
  num_hbond_out.close();
}

// counts hydrogen bonds
// takes in the frame number and atom indices of each hydrogen bond donor and each hydrogen bond acceptor
// returns the count of hydrogen bonds in the frame
int Describe::CountHBonds(vector<vector<int>> donor, vector<int> acceptor, int frame)
{
  /* hydrogen bonds are defined using a distance and angle cutoff
   * distance should be less than 3.5 angstroms
   * angle should be larger than 120 degrees
   */

  int count=0;
  for(unsigned int i=0; i<donor.size(); i++)
  {
    int donor_ind = donor[i][0];
    int H_ind = donor[i][1];
    for(unsigned int j=0; j<acceptor.size(); j++)
    {
      vector<float> first, second;
      int acceptor_ind = acceptor[j];

      first = trajectory_data.coords[frame][donor_ind];
      second = trajectory_data.coords[frame][acceptor_ind];

      float dist = Distance(first, second);
      if(dist < 3.5)
      {
        vector<float> vertex = trajectory_data.coords[frame][H_ind];
        float angle = Angle(first, vertex, second);
        if(angle > 120)
        {
          hbond_out << setw(8) << left << GetResidue(donor_ind) << '\t' << trajectory_data.atoms[donor_ind] 
            << '\t' << setw(8) << left << GetResidue(H_ind) << '\t' << trajectory_data.atoms[H_ind] 
            << '\t' << setw(8) << left << GetResidue(acceptor_ind) << '\t' << trajectory_data.atoms[acceptor_ind] 
            << fixed << setw(11) << setprecision(6) << setfill(' ') << dist << '\t' << angle
            << endl;
          count++;
        }
      }
    }
  }
  return count;
}

// takes in the index of an atom and returns the residue it is in
string Describe::GetResidue(int ind)
{
  vector<int> translated;
  string out;
  ind = ind + 1;
  translated = Translate(residue_pointers, -1*ind);
  int res_ind = FindLargestNegativeInteger(translated);

  if(res_ind > 0)
  {
    out = trajectory_data.residues[res_ind].substr(0,3) + to_string(res_ind+1);
  }
  else
  {
    out = trajectory_data.residues[res_ind].substr(1,3) + to_string(res_ind+1);
  }

  return out;
}

// finds and counts ionic interactions in each frame
void Describe::IonicInteractions(string num_ionic_fname, string ionic_fname)
{
  /* Ionic interactions are defined based on a distance cutoff on the presumed charged groups in the side chains
   * The center of mass of the atoms in a charged group is taken to be the position of the charge
   * A distance of less than 6 angstroms is counted as an ionic interaction
   * NOTE: Definition of charged groups are based on chemical intuition. Actual charges calculated in the force field are not used.
   */

  if(negatives.empty())
  {
    FindChargedResidues();
  }

  ionic_out.open(ionic_fname);
  ionic_out << "# Ionic Interactions" << endl;

  num_ionic_out.open(num_ionic_fname);
  num_ionic_out << "# Ionic Interactions" << endl;

  num_ionic_out << setw(10) << left << "#Frame" << '\t' << left << "Count" << endl;

  for(unsigned int frame=0; frame<n_frames; frame++)
  {
    ionic_out << "#Frame " << frame+1 << endl;
    vector<vector<float>> positive_groups, negative_groups;
    int count=0;

    positive_groups = CalculateChargedGroupsCOM(positives, positives_identity, frame);
    negative_groups = CalculateChargedGroupsCOM(negatives, negatives_identity, frame);

    for(unsigned int i=0; i<positive_groups.size(); i++)
    {
      for(unsigned int j=0; j<negative_groups.size(); j++)
      {
        float dist = Distance(positive_groups[i], negative_groups[j]);
        if(dist < 6.0)
        {
          ionic_out << trajectory_data.residues[positives[i]].substr(0,3) << setw(5) << left << positives[i]+1
            << '\t' << trajectory_data.residues[negatives[j]].substr(0,3) << setw(5) << left << negatives[j]+1 
            << fixed << setw(11) << setprecision(6) << setfill(' ') << dist << endl;
          count++;
        }
      }
    }

    num_ionic_out << setw(10) << left << frame+1 << '\t' << left << count << endl;
  }

  ionic_out.close();
  num_ionic_out.close();
}

// calculates the center of mass of a charged group
// takes in residue pointer of charged residues, atom indices of charged group atoms, and frame number 
// returns a list of the centers of mass
vector<vector<float>> Describe::CalculateChargedGroupsCOM(vector<int> a, vector<vector<int>> b, int frame)
{
  vector<vector<float>> groups;
  for(unsigned int j=0; j<a.size(); j++)
  {
    int first_atom = residue_pointers[a[j]] - 1;
    vector<vector<float>> charged_group;
    vector<float> charged_group_masses;

    for(unsigned int atom = first_atom + b[j][0]; atom < first_atom + b[j][1] + 1; atom++)
    {
      charged_group.push_back(trajectory_data.coords[frame][atom]);
      charged_group_masses.push_back(trajectory_data.masses[atom]);
    }

    vector<float> group_com = CenterOfMass(charged_group, charged_group_masses);
    groups.push_back(group_com);
  }

  return groups;
}

// locates charged residues in the peptide
void Describe::FindChargedResidues(void)
{
  vector<string> positive, negative;
  vector<int> buff_vector;
  vector<vector<int>> positive_side_chains, negative_side_chains;

  // atoms of ARG charged group
  positive.push_back("ARG ");
  buff_vector = {15,21};
  positive_side_chains.push_back(buff_vector);

  positive.push_back("NARG");
  buff_vector = Translate(buff_vector, 2);
  positive_side_chains.push_back(buff_vector);

  buff_vector.clear(); // clear vector

  // atoms of LYS charged group
  positive.push_back("LYS ");
  buff_vector = {16,19};
  positive_side_chains.push_back(buff_vector);

  positive.push_back("NLYS");
  buff_vector = Translate(buff_vector, 2);
  positive_side_chains.push_back(buff_vector);

  buff_vector.clear(); // clear vector

  // atoms of protonated HIS charged group
  positive.push_back("HIP ");
  buff_vector = {8,13};
  positive_side_chains.push_back(buff_vector);

  positive.push_back("NHIP");
  buff_vector = Translate(buff_vector, 2);
  positive_side_chains.push_back(buff_vector);

  buff_vector.clear(); // clear vector

  // atoms of ASP charged group
  negative.push_back("ASP ");
  buff_vector = {7,9};
  negative_side_chains.push_back(buff_vector);

  negative.push_back("NASP");
  buff_vector = Translate(buff_vector, 2);
  negative_side_chains.push_back(buff_vector);

  buff_vector.clear(); // clear vector

  // atoms of GLU charged group
  negative.push_back("GLU ");
  buff_vector = {10,12};
  negative_side_chains.push_back(buff_vector);
  buff_vector.clear();

  negative.push_back("NGLU");
  buff_vector = Translate(buff_vector, 2);
  negative_side_chains.push_back(buff_vector);

  buff_vector.clear(); // clear vector

  for(unsigned int i=0; i<trajectory_data.residues.size(); i++)
  {
    int ind;
    if(IsIn(trajectory_data.residues[i], positive))
    {
      ind = Where(trajectory_data.residues[i], positive);
      positives_identity.push_back(positive_side_chains[ind]);
      positives.push_back(i);
    }
    else if(IsIn(trajectory_data.residues[i], negative))
    {
      ind = Where(trajectory_data.residues[i], negative);
      negatives_identity.push_back(negative_side_chains[ind]);
      negatives.push_back(i);
    }
  }
}

// finds and counts the number of aromatic-aromatic (pi-pi) interactions
void Describe::AromaticAromatic(string num_aromatic_fname, string aromatic_fname)
{
  /* Pi-pi interactions are defined through a distance cutoff between centers of mass of phenyl rings
   * Centers of mass of phenyl rings of aromatic side chains should be between 4.5 and 7 angstroms apart
   */

  if(aromatics.empty())
  {
    FindAromaticResidues();
  }

  aromatic_out.open(aromatic_fname);
  aromatic_out << "# Aromatic-Aromatic (pi-pi) Interactions" << endl;

  num_aromatic_out.open(num_aromatic_fname);
  num_aromatic_out << "# Aromatic-Aromatic (pi-pi) Interactions" << endl;

  num_aromatic_out << setw(10) << left << "#Frame" << '\t' << left << "Count" << endl;

  for(unsigned int frame=0; frame<n_frames; frame++)
  {
    aromatic_out << "#Frame " << frame+1 << endl;
    vector<vector<float>> aromatic_groups;
    int count=0;

    aromatic_groups = CalculateAromaticGroupsCOM(aromatics, aromatics_identity, frame);

    for(unsigned int i=0; i<aromatic_groups.size(); i++)
    {
      for(unsigned int j=i+1; j<aromatic_groups.size(); j++)
      {
        float dist = Distance(aromatic_groups[i], aromatic_groups[j]);
        if(dist > 4.5 && dist < 7.0)
        {
          aromatic_out << trajectory_data.residues[aromatics[i]].substr(0,3) << setw(5) << left << aromatics[i]+1
            << '\t'    << trajectory_data.residues[aromatics[j]].substr(0,3) << setw(5) << left << aromatics[j]+1 
            << '\t' << fixed << setw(11) << setprecision(6) << setfill(' ') << dist << endl;
          count++;
        }
      }
    }

    num_aromatic_out << setw(10) << left << frame+1 << '\t' << left << count << endl;
  }

  aromatic_out.close();
  num_aromatic_out.close();
}

// locates the aromatic residues in the peptide
void Describe::FindAromaticResidues(void)
{
  vector<string> aromatic;
  vector<int> buff_vector;
  vector<vector<int>> aromatic_side_chains;

  // atoms of PHE phenyl ring 
  aromatic.push_back("PHE ");
  buff_vector = {7,8,10,12,14,16};
  aromatic_side_chains.push_back(buff_vector);

  aromatic.push_back("NPHE");
  buff_vector = Translate(buff_vector, 2);
  aromatic_side_chains.push_back(buff_vector);

  buff_vector.clear(); // clear vector

  // atoms of TRP phenyl ring
  aromatic.push_back("TRP ");
  buff_vector = {12,13,15,17,19,21};
  aromatic_side_chains.push_back(buff_vector);

  aromatic.push_back("NTRP");
  buff_vector = Translate(buff_vector, 2);
  aromatic_side_chains.push_back(buff_vector);

  buff_vector.clear(); // clear vector

  // atoms of TYR phenyl ring
  aromatic.push_back("TYR ");
  buff_vector = {7,8,10,12,15,17};
  aromatic_side_chains.push_back(buff_vector);

  aromatic.push_back("NTYR");
  buff_vector = Translate(buff_vector, 2);
  aromatic_side_chains.push_back(buff_vector);

  buff_vector.clear(); // clear vector

  for(unsigned int i=0; i<trajectory_data.residues.size(); i++)
  {
    int ind;
    if(IsIn(trajectory_data.residues[i], aromatic))
    {
      ind = Where(trajectory_data.residues[i], aromatic);
      aromatics_identity.push_back(aromatic_side_chains[ind]);
      aromatics.push_back(i);
    }
  }
}

// calculates centers of mass of phenyl rings
vector<vector<float>> Describe::CalculateAromaticGroupsCOM(vector<int> a, vector<vector<int>> b, int frame)
{
  vector<vector<float>> groups;
  for(unsigned int j=0; j<a.size(); j++)
  {
    int first_atom = residue_pointers[a[j]] - 1;
    vector<vector<float>> aromatic_group;
    vector<float> aromatic_group_masses;

    vector<int> c = Translate(b[j], first_atom);

//    cout << trajectory_data.residues[a[j]] << '\t';
    for(auto atom : c)
    {
      aromatic_group.push_back(trajectory_data.coords[frame][atom]);
      aromatic_group_masses.push_back(trajectory_data.masses[atom]);
      
//      cout << trajectory_data.atoms[atom] << '\t';
    }
//    cout << endl;

    vector<float> group_com = CenterOfMass(aromatic_group, aromatic_group_masses);
    groups.push_back(group_com);
  }

  return groups;
}

// finds and counts cation-pi interactions
void Describe::CationPi(string num_cationpi_fname, string cationpi_fname)
{
  /* Cationic groups are defined as described in ionic interactions
   * Pi groups are defined are described in pi-pi interactions
   */

  if(aromatics.empty())
  {
    FindAromaticResidues();
  }

  if(positives.empty())
  {
    FindChargedResidues();
  }

  cationpi_out.open(cationpi_fname);
  cationpi_out << "# Cation-Pi Interactions" << endl;

  num_cationpi_out.open(num_cationpi_fname);
  num_cationpi_out << "# Cation-Pi Interactions" << endl;

  num_cationpi_out << setw(10) << left << "#Frame" << '\t' << left << "Count" << endl;

  for(unsigned int frame=0; frame<n_frames; frame++)
  {
    cationpi_out << "#Frame " << frame+1 << endl;
    vector<vector<float>> aromatic_groups, positive_groups;
    int count=0;

    aromatic_groups = CalculateAromaticGroupsCOM(aromatics, aromatics_identity, frame);
    positive_groups = CalculateChargedGroupsCOM(positives, positives_identity, frame);

    for(unsigned int i=0; i<positive_groups.size(); i++)
    {
      for(unsigned int j=0; j<aromatic_groups.size(); j++)
      {
        float dist = Distance(positive_groups[i], aromatic_groups[j]);
        if(dist < 6.0)
        {
          cationpi_out << trajectory_data.residues[positives[i]].substr(0,3) << setw(5) << left << positives[i]+1 
            << '\t'    << trajectory_data.residues[aromatics[j]].substr(0,3) << setw(5) << left << aromatics[j]+1
            << '\t' << fixed << setw(11) << setprecision(6) << setfill(' ') << dist << endl;
          count++;
        }
      }
    }

    num_cationpi_out << setw(10) << left << frame+1 << '\t' << left << count << endl;
  }

  cationpi_out.close();
  num_cationpi_out.close();
}

// finds and counts aromatic-sulphur interactions
void Describe::AromaticSulphur(string num_aromaticsulphur_fname, string aromaticsulphur_fname)
{
  if(aromatics.empty())
  {
    FindAromaticResidues();
  }

  if(cysteines.empty())
  {
    FindCysteineResidues();
  }

  aromaticsulphur_out.open(aromaticsulphur_fname);
  aromaticsulphur_out << "# Aromatic-Sulphur Interactions" << endl;

  num_aromaticsulphur_out.open(num_aromaticsulphur_fname);
  num_aromaticsulphur_out << "# Aromatic-Sulphur Interactions" << endl;

  num_aromaticsulphur_out << setw(10) << left << "#Frame" << '\t' << left << "Count" << endl;

  for(unsigned int frame=0; frame<n_frames; frame++)
  {
    aromaticsulphur_out << "#Frame " << frame+1 << endl;
    vector<vector<float>> aromatic_groups;
    int count=0;

    aromatic_groups = CalculateAromaticGroupsCOM(aromatics, aromatics_identity, frame);

    for(unsigned int i=0; i<aromatic_groups.size(); i++)
    {
      for(unsigned int j=0; j<cysteines.size(); j++)
      {
        int cys_first_atom = residue_pointers[cysteines[j]] - 1;
        int sulphur = cys_first_atom + cysteines_identity[j];

        float dist = Distance(aromatic_groups[i], trajectory_data.coords[frame][sulphur]);
        if(dist < 5.3)
        {
          aromaticsulphur_out << trajectory_data.residues[aromatics[i]].substr(0,3) << setw(5) << left << aromatics[i]+1 
            << '\t'           << trajectory_data.residues[cysteines[j]].substr(0,3) << setw(5) << left << cysteines[j]+1
            << '\t' << fixed << setw(11) << setprecision(6) << setfill(' ') << dist << endl;
          count++;
        }
      }
    }

    num_aromaticsulphur_out << setw(10) << left << frame+1 << '\t' << left << count << endl;
  }

  aromaticsulphur_out.close();
  num_aromaticsulphur_out.close();
}

void Describe::FindCysteineResidues(void)
{
  vector<string> cysteine;
  vector<int> cysteine_side_chains;

  cysteine.push_back("CYS ");
  cysteine_side_chains.push_back(7);

  cysteine.push_back("NCYS");
  cysteine_side_chains.push_back(9);

  cysteine.push_back("CYX ");
  cysteine_side_chains.push_back(7);

  cysteine.push_back("NCYX");
  cysteine_side_chains.push_back(9);

  for(unsigned int i=0; i<trajectory_data.residues.size(); i++)
  {
    int ind;
    if(IsIn(trajectory_data.residues[i], cysteine))
    {
      ind = Where(trajectory_data.residues[i], cysteine);
      cysteines_identity.push_back(cysteine_side_chains[ind]);
      cysteines.push_back(i);
    }
  }
}

// finds and counts disulfide bonds
void Describe::Disulfide(string num_disulfide_fname, string disulfide_fname)
{
  if(cysteines.empty())
  {
    FindCysteineResidues();
  }

  disulfide_out.open(disulfide_fname);
  disulfide_out << "# Disulfide Bonds" << endl;

  num_disulfide_out.open(num_disulfide_fname);
  num_disulfide_out << "# Disulfide Bonds" << endl;

  num_disulfide_out << setw(10) << left << "#Frame" << '\t' << left << "Count" << endl;

  for(unsigned int frame=0; frame<n_frames; frame++)
  {
    disulfide_out << "#Frame " << frame+1 << endl;
    vector<vector<float>> aromatic_groups;
    int count=0;

    for(unsigned int i=0; i<cysteines.size(); i++)
    {
      for(unsigned int j=i+1; j<cysteines.size(); j++)
      {
        int cys_first_atom_1 = residue_pointers[cysteines[i]] - 1;
        int sulphur_1 = cys_first_atom_1 + cysteines_identity[i];

        int cys_first_atom_2 = residue_pointers[cysteines[j]] - 1;
        int sulphur_2 = cys_first_atom_2 + cysteines_identity[j];

        float dist = Distance(trajectory_data.coords[frame][sulphur_1], trajectory_data.coords[frame][sulphur_2]);
        if(dist < 2.2)
        {
          disulfide_out << trajectory_data.residues[cysteines[i]].substr(0,3) << setw(5) << left << cysteines[i]+1 
            << '\t'     << trajectory_data.residues[cysteines[j]].substr(0,3) << setw(5) << left << cysteines[j]+1
            << '\t' << fixed << setw(11) << setprecision(6) << setfill(' ') << dist << endl;
          count++;
        }
      }
    }

    num_disulfide_out << setw(10) << left << frame+1 << '\t' << left << count << endl;
  }

  disulfide_out.close();
  num_disulfide_out.close();
}

// finds and counts hydrophobic interactions
void Describe::Hydrophobic(string num_hydrophobic_fname, string hydrophobic_fname)
{
  /* Hydrophobic interactions are defined using a distance cutoff of C-H bonds
   * 
   */

//  vector<vector<int>> CHs;
  vector<int> buff_vector = {-100,-100};

  CHs.push_back(buff_vector);

  for(unsigned int i=0; i<bonh.size(); i++)
  {
    if(IsIn(bonh[i][0], Cs) && trajectory_data.atoms[bonh[i][0]] != "CA  " && bonh[i][0] != CHs.back()[0]) 
    {
      CHs.push_back(bonh[i]);
    }
  }

  CHs.erase(CHs.begin());

  hydrophobic_out.open(hydrophobic_fname);
  hydrophobic_out << "# Hydrophobic Interactions" << endl;

  num_hydrophobic_out.open(num_hydrophobic_fname);
  num_hydrophobic_out << "# Hydrophobic Interactions" << endl;

  num_hydrophobic_out << setw(10) << left << "#Frame" << '\t' << left << "Count" << endl;

  for(unsigned int frame=0; frame<n_frames; frame++)
  {
    hydrophobic_out << "#Frame " << frame+1 << endl;
    int count=0;

    for(unsigned int i=0; i<CHs.size(); i++)
    {
      int one_ind = CHs[i][0];

      for(unsigned int j=i+1; j<CHs.size(); j++)
      {
        int two_ind = CHs[j][0];

        float dist = Distance(trajectory_data.coords[frame][one_ind], trajectory_data.coords[frame][two_ind]);
        if(dist < 5.0 && dist > 0)
        {
          hydrophobic_out << setw(8) << left << GetResidue(one_ind) << trajectory_data.atoms[one_ind]
            << '\t'       << setw(8) << left << GetResidue(two_ind) << trajectory_data.atoms[two_ind]
            << '\t' << fixed << setw(11) << setprecision(6) << setfill(' ') << dist << endl;
          count++;
        }
      }
    }

    num_hydrophobic_out << setw(10) << left << frame+1 << '\t' << left << count << endl;
  }

  hydrophobic_out.close();
  num_hydrophobic_out.close();
}

void Describe::HydrophobicOverHydrophilic(string num_hoverh_fname, string hoverh_fname)
{
  hoverh_out.open(hoverh_fname);
  hoverh_out << "# Hydrophobic over Hydrophilic" << endl;

  num_hoverh_out.open(num_hoverh_fname);
  num_hoverh_out << "# Hydrophobic over Hydrophilic" << endl;

  num_hoverh_out << setw(10) << left << "#Frame" << '\t' << left << "Count" << endl;
  for(unsigned int i=0; i<n_frames; i++)
  {
    hoverh_out << "#Frame " << i+1 << endl;

    vector<float> frame_com = CenterOfMass(trajectory_data.coords[i], trajectory_data.masses);
    float nonpolar_count, polar_count, dist, ratio;

    for(unsigned int j=0; j<CHs.size(); j++)
    {
      int ind = CHs[j][0];
      dist = Distance(trajectory_data.coords[i][j], frame_com);
      if(dist<4.0)
      {
        nonpolar_count++;
      }
    }

    nonpolar_count = nonpolar_count / CHs.size();

    for(unsigned int j=0;j<NHs.size(); j++)
    {
      int ind = NHs[j][0];
      dist = Distance(trajectory_data.coords[i][j], frame_com);
      if(dist>4.0)
      {
        polar_count++;
      }
    }

    for(unsigned int j=0; j<OHs.size(); j++)
    {
      int ind = OHs[j][0];
      dist = Distance(trajectory_data.coords[i][j], frame_com);
      if(dist>4.0)
      {
        polar_count++;
      }
    }

    polar_count = polar_count / (NHs.size() + OHs.size());

    ratio = nonpolar_count / polar_count;

    num_hoverh_out << setw(10) << left << i+1 << '\t' << left << ratio << endl;
  }

  hoverh_out.close();
  num_hoverh_out.close();
}

#endif
