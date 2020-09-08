#include "brickwork_mps.h"
#include <stdlib.h>
#include "itensor/all.h"

using namespace itensor;
using namespace std::chrono;

int main(int argc, char *argv[]) {
  int n = strtol(argv[1], NULL, 10); // length of MPS
  int m = strtol(argv[2], NULL, 10); // length of the other side of the circuit
  // m is also the number of iterations of SEBD
  int q = strtol(argv[3], NULL, 10); // local qudit dimension
  int method = strtol(argv[4], NULL, 10); // 1 for error, 2 for maxdim
  float err = atof(argv[5]); // maximum error allowed when truncating
  int maxdim = strtol(argv[6], NULL, 10);
  std::string unitary_file = argv[7]; // where the random unitary gates are held
  std::string out_file = argv[8]; // where the entropies are written to
  int circuit_type = strtol(argv[9], NULL, 10); // 1 for brickwork, 2 for cluster
  int orth_type = strtol(argv[10], NULL, 10);

  UnitaryData udata;

  BrickMPS circuit("circuit1", n, q, false);
  entanglement_data ent = circuit.get_entanglement();
  entanglement_data ent_prev;

  std::ofstream ofile;
  ofile.open(out_file);
  ofile << "n = " << n << std::endl;
  ofile << "m = " << m << std::endl;
  ofile << "q = " << q << std::endl;
  ofile << "cutoff error = " << err << std::endl;
  ofile << "gates file = " << unitary_file << std::endl;

  ofile << "iteration " << 0 << std::endl;
  ofile << "entanglement entropies" << std::endl;
  for (int i = 0; i < n - 1; i++) {
    ofile << ent.entropies[i] << ",";
  }
  ofile << std::endl;


  std::cout << "entanglement across middle: " << ent.entropies[(n / 2) - 1] << std::endl;
  auto sample_tensor = circuit.get_tensor(n / 2);

  std::cout << "initial circuit: ";
  for (int i = 1; i <= n; i++) {
    // PrintData(circuit.get_tensor(i));
  }


  std::ifstream myfile;
  myfile.open(unitary_file,std::ios::in);
  std::string line;
  std::istringstream is;
  std::complex<double> c;

  // print the number of unitary matrices per iteration
  std::cout << "Unitary number: " <<  (n + 1 + int(2 * (n / 8.0))) << std::endl;

  // define the time variables
  duration<double> it_time, trunc_time;
  steady_clock::time_point begin, end;

  int num_gates;
  if (circuit_type == 1) {
    num_gates = n + 1 + int(2 * (n / 8.0));
  } else {
    num_gates = 2 * n + 1;
  }
  // for each iteration get the gates and do the measurement
  for (int it = 1; it <= m; it++) {
    // get unitary matrix from the file
    for (int p = 1; p <= num_gates; p++) {
        //std::cout << "Unitary Number: "<< p << std::endl;
        for (int l2 = 0; l2 < std::pow(q, 2); l2++) {
            for (int l = 0; l < std::pow(q, 2); l++) {
                std::getline(myfile,line,';');
                std::istringstream is(line);
                is >> c;
                udata[p][l2][l]=c;
                //std::cout << "("<<l2<<", "<<l<<"): "<< c << std::endl;
            }
        }
    }
    println("iteration " + std::to_string(it));
    ofile << "iteration " << it << std::endl;

    // do the actual iteration of sebd
    // and print the time it takes to iterate
    begin = steady_clock::now();
    circuit.iterate(it, udata, circuit_type, orth_type);
    end = steady_clock::now();
    it_time = duration_cast<duration<double>>(end - begin);

    ent_prev = circuit.get_entanglement();
    // println("singular values before truncation: ");
    // std::cout << std::endl << "singular values in the middle: " << std::endl;
    // for (int i = 0; i < ent_prev.svals[n/2].size(); i++) {
    //   std::cout << ent_prev.svals[int(n / 2)][i] << ",";
    //   // ofile << ent.entropies[i] << ",";
    // }
    // std::cout << std::endl;

    begin = steady_clock::now();
    if (method == 1) {
      circuit.truncate_err(err);
    } else {
      circuit.truncate_maxdim(maxdim);
    }
    end = steady_clock::now();
    trunc_time = duration_cast<duration<double>>(end - begin);

    // get and output the entanglement entropy
    ent = circuit.get_entanglement();


    println("entanglement entropies: ");
    ofile << "entanglement entropies" << std::endl;
    for (int i = 0; i < n - 1; i++) {
      std::cout << ent.entropies[i] << " - ";
      ofile << ent.entropies[i] << ",";
    }
    // sort the singular values
    std::sort (ent.svals[int(n/2)].begin(), ent.svals[int(n/2)].end());
    // println("singular values: ");
    ofile << std::endl << "singular values in the middle: " << std::endl;
    for (int i = 0; i < ent.svals[n/2].size(); i++) {
      // std::cout << ent.svals[int(n / 2)][i] << ",";
      ofile << ent.svals[int(n / 2)][i] << ",";
      // ofile << ent.entropies[i] << ",";
    }
    std::cout << std::endl;
    ofile << std::endl;

    //output the time it took to do everything:
    std::cout << "it time: " << it_time.count() << std::endl;
    std::cout << "truncation time: " << trunc_time.count() << std::endl;
    ofile << "it time: " << std::endl << it_time.count() << std::endl;
    ofile << "truncation time: " << std::endl << trunc_time.count() << std::endl;

    ofile << "entanglement entropies before truncation" << std:: endl;
    for (int i = 0; i < n - 1; i++) {
      ofile << ent_prev.entropies[i] << ",";
    }
    // sort the singular values
    std::sort (ent_prev.svals[int(n/2)].begin(), ent_prev.svals[int(n/2)].end());
    // println("singular values: ");
    ofile << std::endl << "singular values in the middle: " << std::endl;
    for (int i = 0; i < ent_prev.svals[n/2].size(); i++) {
      // std::cout << ent.svals[int(n / 2)][i] << ",";
      ofile << ent_prev.svals[int(n / 2)][i] << ",";
      // ofile << ent.entropies[i] << ",";
    }
    ofile << std::endl;
  }

  ofile.close();
  myfile.close();

  return 0;
}
