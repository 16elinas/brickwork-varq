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

  UnitaryData udata;

  BrickMPS circuit("circuit1", n, q, false);
  entanglement_data ent = circuit.get_entanglement();

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

  // for each iteration get the gates and do the measurement
  for (int it = 1; it <= m; it++) {
    // get unitary matrix from the file
    for (int p = 1; p <= n + 1 + int(2 * (n / 8.0)); p++) {
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
    circuit.iterate(it, udata);
    end = steady_clock::now();
    it_time = duration_cast<duration<double>>(end - begin);

    begin = steady_clock::now();
    if (method == 1) {
      circuit.truncate(err);
    } else {
      circuit.truncate(maxdim);
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
    println("singular values: ");
    ofile << std::endl << "singular values in the middle: " << std::endl;
    for (int i = 0; i < ent.svals[n/2].size(); i++) {
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
  }

  ofile.close();
  myfile.close();

  return 0;
}
