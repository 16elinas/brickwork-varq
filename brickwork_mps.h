#ifndef BRICKMPS_H
#define BRICKMPS_H

#include <array>
#include <chrono>
#include <cmath>
#include <complex>
#include <fstream>
#include <iostream>
#include <map>
#include <random>
#include <string>
#include <vector>
#include "itensor/all.h"
using namespace itensor;

// an array of unitary q^2 by q^2 matrices
// 1-indexed in the outer map
// 0-indexed in the inner maps (matrices)
typedef std::map<int, std::map<int, std::map<int, std::complex<double> > > >
    UnitaryData;

// 0-indexed
struct entanglement_data {
  std::vector<std::vector<double>> svals;
  std::vector<double> entropies;
} ;

class BrickMPS {
  std::map<int, ITensor> qudits;
  std::string name;
  int q;

  public:
    BrickMPS(std::string nam, int len, int qdim, bool prod);

    MPS to_mps();

    int get_len();
    std::string get_name();
    ITensor& get_tensor(int pos);
    std::vector<double> get_singvals(int b);
    entanglement_data get_entanglement();
    ITensor apply_unitary(ITensor T1, ITensor T2, std::map<int, std::map<int, std::complex<double>>> u);

    BrickMPS& measure(std::map<int, Index> inds_to_meas, std::map<int, Index> new_phys_inds);
    BrickMPS& truncate(int maxdim);
    BrickMPS& truncate(float maxerror);
    BrickMPS& iterate(int iter, UnitaryData udata);
};

#endif
