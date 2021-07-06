#include "brickwork_mps.h"
#include <stdlib.h>
#include "itensor/all.h"

using namespace itensor;
using namespace std::chrono;

/////////////////////////////////////////////////////////
// Run the SEBD algorithm on an ITensor MPS representing
// the qubits, using parameters from
// the command line, and write the entanglement
// data of the circuit to a file at each step
////////////////////////////////////////////////////////
int main(int argc, char *argv[]) {
  // command-line arguments in order of how they should be entered
  int n = strtol(argv[1], NULL, 10); // length of MPS
  int m = strtol(argv[2], NULL, 10); // length of the other side of the circuit
  // m is also the number of iterations of SEBD
  int q = strtol(argv[3], NULL, 10); // local qudit dimension
  int q2 = std::pow(q, 2);
  int method = strtol(argv[4], NULL, 10); // 1 for error, 2 for maxdim
  float err = atof(argv[5]); // maximum error allowed when truncating
  int maxdim = strtol(argv[6], NULL, 10);
  std::string unitary_file = argv[7]; // name of file where the random unitary gates are written
  std::string out_file = argv[8]; // where the entropies are written to
  int circuit_type = strtol(argv[9], NULL, 10); // 1 for brickwork, 2 for cluster, 3 for more-entangling bw, 4 for less-entangling
  int orth_type = strtol(argv[10], NULL, 10); // type of orthogonalization just before measurement

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

  std::cout << "q = " << q << std::endl;
  std::cout << "method = " << method << std::endl;
  std::cout << "cutoff error = " << err << std::endl;
  std::cout << "maxdim = " << maxdim << std::endl;

  std::cout << "entanglement across middle: " << ent.entropies[(n / 2) - 1] << std::endl;
  auto sample_tensor = circuit.get_tensor(n / 2);


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
        for (int l2 = 0; l2 < q2; l2++) {
            for (int l = 0; l < q2; l++) {
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

    begin = steady_clock::now();
    circuit.truncate(err, maxdim, method, ofile);
    end = steady_clock::now();
    trunc_time = duration_cast<duration<double>>(end - begin);


    //output the time it took to do everything:
    std::cout << "it time: " << it_time.count() << std::endl;
    std::cout << "truncation time: " << trunc_time.count() << std::endl;
    ofile << "it time: " << std::endl << it_time.count() << std::endl;
    ofile << "truncation time: " << std::endl << trunc_time.count() << std::endl;
  }

  ofile.close();
  myfile.close();

  delete &circuit;

  return 0;
}
