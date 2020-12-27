#include "describe.hpp"
#include <iostream>
#include <string>

using namespace std;


int main(int argc, char* argv[])
{
  string traj_fname, parm_fname;

  if(argc != 3)
  {
    cerr << "Usage: " << argv[0] << " \"NetCDF file\" " << " \"prmtop/parm7 file\" " << endl;
    exit(1);
  }
  else
  {
    traj_fname = argv[1];
    parm_fname = argv[2];

    cout << "NetCDF file: " << traj_fname << endl;
    cout << "Parm7 file: " << parm_fname << endl;

    Describe test(traj_fname, parm_fname);
    test.HydrogenBonds("results/num_hbond.dat", "results/hbond.dat");
    test.IonicInteractions("results/num_ionic.dat", "results/ionic.dat");
    test.AromaticAromatic("results/num_aromatic.dat", "results/aromatic.dat");
    test.CationPi("results/num_cationpi.dat", "results/cationpi.dat");
    test.AromaticSulphur("results/num_aromaticsulphur.dat", "results/aromaticsulphur.dat");
    test.Disulfide("results/num_disulfide.dat", "results/disulfide.dat");
    test.Hydrophobic("results/num_hydrophobic.dat", "results/hydrophobic.dat");
    test.HydrophobicOverHydrophilic("results/num_hoverh.dat", "results/hoverh.dat");
  }
}
