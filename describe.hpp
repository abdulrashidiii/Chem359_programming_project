#include <netcdf>
#include <iostream>
#include <sstream>
#include <fstream>
#include <string>
#include <vector>
#include <cmath>
#include <iomanip>
#include <algorithm>
#include "vecop.hpp"

using namespace std;
using namespace netCDF;
using namespace netCDF::exceptions;

struct traj{
  vector<vector<vector<float>>> coords; // vector of trajectory xyz coordinates
  vector<string> atoms; // vector of atom names
  vector<string> residues; // vector of residue names
  vector<float> charges; // vector of atom charges
  vector<float> masses; // vector of masses
};

class Describe
{
protected:
  traj trajectory_data; // trajectory data

  int n_frames; // number of frames
  int n_atoms; // number of atoms
  int n_types; // number of types
  int nbonh; // number of bonds containing hydrogen
  int mbona; // number of bonds not containing hydrogen
  int ntheth; // number of angles containing hydrogen
  int mtheta; // number of angles not containing hydrogen
  int n_residues; // number of residues
  int n_bonds; // number of unique bond types
  int n_angles; // number of unique angle types
  int n_torsions; // number of unique torsion types
  int n_solty; // number of solty types
  vector<vector<int>> bonh; // bond containing hydrogen
  vector<int> residue_pointers; // atom indices of first atoms in each residue

  vector<int> Os, Ns, Hs, Cs, Ss; // indices of each specific atom
  vector<vector<int>> CHs, NHs, OHs; // indices of atoms involved in covalent bonds with hydrogen atoms

  vector<int> positives, negatives; // indices of charged residues
  vector<vector<int>> positives_identity, negatives_identity; // indices of charged groups in charge residues

  vector<int> aromatics; // indices of aromatic residues
  vector<vector<int>> aromatics_identity; // indices of aromatic atoms in aromatic residues
  vector<int> aromatic_atoms; // indices of aromatic atoms

  vector<int> cysteines; // indices of cysteine residues
  vector<int> cysteines_identity; // indices of S atoms in cysteine residues

  // output file streams and file names
  ofstream num_hbond_out, hbond_out;
  string num_hbond_fname, hbond_fname;
  ofstream num_ionic_out, ionic_out;
  string num_ionic_fname, ionic_fname;
  ofstream num_aromatic_out, aromatic_out;
  string num_aromatic_fname, aromatic_fname;
  ofstream num_cationpi_out, cationpi_out;
  string num_cationpi_fname, cationpi_fname;
  ofstream num_aromaticsulphur_out, aromaticsulphur_out;
  string num_aromaticsulphur_fname, aromaticsulphur_fname;
  ofstream num_disulfide_out, disulfide_out;
  string num_disulfide_fname, disulfide_fname;
  ofstream num_hydrophobic_out, hydrophobic_out;
  string num_hydrophobic_fname, hydrophobic_fname;
  ofstream num_hoverh_out, hoverh_out;
  string num_hoverh_fname, hoverh_fname;

  // NetCDF parser
  void ReadNetCDF(string);

  // prmtop parser
  void ReadParm(string);

  // reads a block of information in prmtop files
  string ParmReadBlock(ifstream*, int);

  // populates atomic index arrays
  void GetAtomIndices(void);

  // returns residue where an atom belongs
  string GetResidue(int);

  // counts the number of hydrogen bonds
  int CountHBonds(vector<vector<int>>, vector<int>, int);

  // finds charged residues
  void FindChargedResidues(void);

  // finds aromatic residues
  void FindAromaticResidues(void);

  // finds cysteine residues
  void FindCysteineResidues(void);

  // calculates center of mass of aromatic groups
  vector<vector<float>> CalculateAromaticGroupsCOM(vector<int>, vector<vector<int>>, int);

  // calculates center of mass of charged groups
  vector<vector<float>> CalculateChargedGroupsCOM(vector<int>, vector<vector<int>>, int);

public:
  // class constructor
  Describe(string, string);

  // functions for each structural descriptor
  void HydrogenBonds(string, string);
  void IonicInteractions(string, string);
  void AromaticAromatic(string, string);
  void CationPi(string, string);
  void AromaticSulphur(string, string);
  void Disulfide(string, string);
  void Hydrophobic(string, string);
  void HydrophobicOverHydrophilic(string, string);
};
