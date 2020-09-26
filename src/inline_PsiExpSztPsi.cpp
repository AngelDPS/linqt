// C & C++ libraries
#include <iostream> /* for std::cout mostly */
#include <string>   /* for std::string class */
#include <fstream>  /* for std::ofstream and std::ifstream functions classes*/
#include <stdlib.h>
#include <chrono>


#include "kpm_noneqop.hpp" //Message functions
#include "chebyshev_moments.hpp"
#include "sparse_matrix.hpp"
#include "quantum_states.hpp"
#include "chebyshev_solver.hpp"

int main(int argc, char *argv[])
{
  const std::string
    S_NTIME= argv[1],
    S_TMAX = argv[2];

  const int numTimes = atoi( S_NTIME.c_str() );
  const double tmax = stod( S_TMAX );

  chebyshev::MomentsTD chebMoms(1, numTimes);

  SparseMatrixType Sz;
  Sz.SetID("Sz");
  Sz.setDimensions(2, 2);
  vector<int> Sz_colrowidx = {0, 1};
  vector<complex<double>> Sz_data = {-1.0*0.1, 1.0*0.1};
  Sz.ConvertFromCOO(Sz_colrowidx, Sz_colrowidx, Sz_data);
  std::array<double, 2> spectral_bounds = {-1.0*0.1, 1.0*0.1};

  //CONFIGURE THE CHEBYSHEV MOMENTS
  chebMoms.SystemLabel("Sz");
  chebMoms.BandWidth ( (spectral_bounds[1]-spectral_bounds[0])*1.0 );
  chebMoms.BandCenter( (spectral_bounds[1]+spectral_bounds[0])*0.5 );
  chebMoms.TimeDiff( tmax/(numTimes-1) );
  chebMoms.SetAndRescaleHamiltonian(Sz);
  chebMoms.Print();

  chebMoms.ResetTime();
  vector<complex<double>>
    PhiL = {1.0, 1.0},
    PhiR = {1.0, 1.0};

  std::string outputName = "PsiExp_miSztw_Psi.dat";
  std::ofstream outputfile( outputName.c_str() );
  
  while ( chebMoms.CurrentTimeStep() != chebMoms.MaxTimeStep() )
    {
      const auto n = chebMoms.CurrentTimeStep();
      chebMoms(0, n) += linalg::vdot(PhiL, PhiR);
      chebMoms.IncreaseTimeStep();
      chebMoms.Evolve(PhiL);
      outputfile << n * chebMoms.TimeDiff() << " " << chebMoms(0, n).real() << std::endl;
    }
  outputfile.close();
}
