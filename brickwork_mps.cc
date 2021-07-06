#include "brickwork_mps.h"


// Constructor for the BrickMPS class
// len is the length of the MPS, the number of qudits on the chain
// q is the local dimension of the qudit
// prod is whether we want a product state or not
// name should be unique
BrickMPS::BrickMPS(std::string nam, int len, int qdim, bool prod) {
  name = nam;
  q = qdim;

  if (prod == false) {
    // i is the physical index
    auto i = Index(q, "i=1,Site,it=" + name);
    auto right = Index(q, "bond=1,Link");

    qudits[1] = ITensor(i, right);
    // make the 'left'-most tensor an equal diagonal one
    for (int j = 1; j <= q; j++) {
      qudits[1].set(i(j), right(j), 1 / sqrt(q));
    }

    // make the rest of the tensors be some entangled state
    for (int ind = 2; ind <= len - 1; ind++) {
      auto left = right;
      auto new_right = Index(q, "bond=" + std::to_string(ind) + ",Link");
      i = Index(q, "i=" + std::to_string(ind) + ",Site,it=" + name);
      qudits[ind] = ITensor(i, left, new_right);

      for (int ival = 1; ival <= q; ival++) {
        for (int j = 1; j <= q; j++) {
          if (ival % 2 == 0 && j % 2 == 0) {
            qudits[ind].set(i(ival), left(j), new_right(ival), -1 / sqrt(q));
          } else {
            qudits[ind].set(i(ival), left(j), new_right(ival), 1 / sqrt(q));
          }

        }
      }
      qudits[ind] = qudits[ind] / norm(qudits[ind]) * sqrt(q);
      right = new_right;
    }

    // make the rightmost tensor
    i = Index(q, "i=" + std::to_string(len) + ",Site,it=" + name);
    auto left = right;
    qudits[len] = ITensor(i, left);
    for (int ival = 1; ival <= q; ival++) {
      for (int lval = 1; lval <= q; lval++) {
        if (ival == q && lval == q) {
          qudits[len].set(i(ival), left(lval), -1 / sqrt(q));
        } else {
          qudits[len].set(i(ival), left(lval), 1 / sqrt(q));
        }
      }
    }

  } else {
    // for a product state
    for (int ind = 1; ind <= len; ind++) {
      auto i = Index(q, "i=" + std::to_string(ind) + ",Site,it=" + name);
      qudits[ind] = ITensor(i);

      qudits[ind].set(i(1), 1);
   }
  }

}

int BrickMPS::get_len() { return (this->qudits).size(); }

ITensor& BrickMPS::get_tensor(int pos) { return (this->qudits)[pos]; }

// Turn the BrickMPS into an iTensor MPS
MPS BrickMPS::to_mps() {
  int N = this->get_len();
  auto sites = SiteSet(N, q);
  auto psi = MPS(sites);
  for (int pos = 1; pos <= N; pos++) {
    psi.set(pos, qudits[pos]);
  }
  return psi;
}


// output the entanglement data of the current state of the circuit
// in the form of entanglement_data which has a list of the
// total von Neumann entropies across the bonds and also
// the lists of all the singular values
// the lists are zero indexed
// assumes the MPS is already orthogonalized and normalized
entanglement_data BrickMPS::get_entanglement() {
  int N = this->get_len();
  entanglement_data out;
  std::vector<std::vector<double>> svals;
  std::vector<double> entropies;

  auto psi = this->to_mps();
  // psi.orthogonalize();
  // psi.normalize();
  for (int i = 1; i <= this->get_len(); i++) {
    qudits[i] = psi(i);
  }

  // variables for each iteration
  double ent;
  std::vector<double> bond_svals;


  std::map<int, ITensor> lambdas;

  ITensor conjd = conj(prime(qudits[1], "Link"));

  lambdas[1] = qudits[1] * conjd;

  for (int m = 2; m <= N; m++) {
    conjd = conj(prime(qudits[m], "Link"));
    lambdas[m] = lambdas[m - 1] * conjd * qudits[m];
  }

  for (int m = 1; m <= N - 1; m++) {
    int D = dim(lambdas[m].inds().index(1));
    ent = 0;
    bond_svals.clear();
    // println("position: " + std::to_string(m));
    for (int p = 1; p <= D; p++) {
      double val = real(lambdas[m].eltC(p, p));
      bond_svals.push_back(val);
      // println("singular value: " + std::to_string(val));
      ent += -val * log2(val);
    }
    svals.push_back(bond_svals);
    entropies.push_back(ent);
  }

  out.svals = svals;
  out.entropies = entropies;

  return out;
}


// output the entanglement data of the current state of the circuit
// in the form of entanglement_data which has a list of the
// total von Neumann entropies across the bonds and also
// the lists of all the singular values
// the lists are zero indexed
// assumes the MPS is already orthogonalized and normalized
entanglement_data BrickMPS::get_entanglement(MPS& mps) {
  int N = length(mps);
  entanglement_data out;
  std::vector<std::vector<double>> svals;
  std::vector<double> entropies;

  // variables for each iteration
  double ent;
  std::vector<double> bond_svals;

  std::map<int, ITensor> lambdas;

  ITensor conjd = conj(prime(mps(1), "Link"));

  lambdas[1] = mps(1) * conjd;

  for (int m = 2; m <= N; m++) {
    conjd = conj(prime(mps(m), "Link"));
    lambdas[m] = lambdas[m - 1] * conjd * mps(m);
  }

  for (int m = 1; m <= N - 1; m++) {
    int D = dim(lambdas[m].inds().index(1));
    ent = 0;
    bond_svals.clear();
    // println("position: " + std::to_string(m));
    for (int p = 1; p <= D; p++) {
      double val = real(lambdas[m].eltC(p, p));
      bond_svals.push_back(val);
      // println("singular value: " + std::to_string(val));
      ent += -val * log2(val);
    }
    svals.push_back(bond_svals);
    entropies.push_back(ent);
  }

  out.svals = svals;
  out.entropies = entropies;

  return out;
}


// apply a unitary gate between two ITensors and return the product
// assumes that T1 and T2 only have 1 physical (Site) index at a time
// u is a q^2 by q^2 unitary matrix
ITensor BrickMPS::apply_unitary(ITensor T1, ITensor T2, std::map<int, std::map<int, std::complex<double>>> u) {
  // find the physical indices of t1 and t2
  Index i, j;
  for (auto& index : inds(T1)) {
    if (hasTags(index, "Site")) {
      i = index;
    }
  }
  for (auto& index : inds(T2)) {
    if (hasTags(index, "Site")) {
      j = index;
    }
  }
  // these are the other 2 indices of U (the 'output' indices)
  auto ii = prime(i);
  auto jj = prime(j);

  // turn the matrix into a 4-index tensor, each index with dimension q
  auto U = ITensor(i, j, ii, jj);

  int uind1, uind2;

  for (int ival = 1; ival <= q; ival++) {
    for (int jval = 1; jval <= q; jval++) {
      for (int iival = 1; iival <= q; iival++) {
        for (int jjval = 1; jjval <= q; jjval++) {
          uind1 = (iival-1) * q + (jjval-1);
          uind2 = (ival-1) * q + (jval-1);
          U.set(i(ival), j(jval), ii(iival), jj(jjval), u[uind1][uind2]);
        }
      }
    }
  }


   auto out = T1 * T2;
   out *= U;
   out.noPrime();
   return out;
}


// measure the indices which are to be measured and orthogonalize the MPS
// inds_to_meas: a 1-indexed array of the indices to be measured
// new_phys_inds: the indices to be left alone
// orth_type: different types of orthogonalization before the measurement
//            1 - do nothing
//            2 - use ITensor's "position" method on the first position
//            3 - use ITensor's "orthogonalize" method
// what this method does for each index:
//   1. pick a measurement result, based on magnitude of projected qudit
//   2. project the qudit onto that result, for the index to be measured and normalize
//   3. contract into the next qudit to prepare
BrickMPS& BrickMPS::measure(std::map<int, Index> inds_to_meas, std::map<int, Index> new_phys_inds, int orth_type) {
  // put it in canonical form
  // orhogonalization: 1 is none, 2 is position, 3 is orthogonalize
  if (orth_type == 2) {
    auto copyMPS = this->to_mps();
    println("positioning");
    copyMPS.position(1);
    for (int i = 1; i <= this->get_len(); i++) {
      qudits[i] = copyMPS(i);
    }
  } else if (orth_type == 3) {
    auto copyMPS = this->to_mps();
    println("orthogonalizing");
    copyMPS.orthogonalize();
    for (int i = 1; i <= this->get_len(); i++) {
      qudits[i] = copyMPS(i);
    }
  }

  // projectors for each position handled separately
  std::vector<ITensor> projs;

  // record the outcome string
  string z;
  int N = this->get_len();
  println("measuring");

  // make an MPS where the physical indices are combined, into ones of dimension q^2
  // this is for ease of operations and gauging
  auto sites = SiteSet(N, q*q);
  auto combinedMPS = MPS(sites);


  // get the probability of each measurement outcome

  // treat the first position differently
  // probs is a vector of length q of the probabilities of gettnig a measurement
  // result <= to i, where i is the index in the vector
  std::vector<float> probs;
  float prob;
  float prev_prob = 0;

  // create the projector
  auto n = new_phys_inds[1];
  auto nn = prime(n);
  auto m = inds_to_meas[1];
  auto mm = prime(m);

  auto [C, c] = combiner({n, m}, {"Tags", "Site,i="+std::to_string(1)});
  auto [CC, cc] = combiner({nn, mm}, {"Tags", "prout,Site,i="+std::to_string(1)});

  auto combinedSite = qudits[1] * C;
  combinedMPS.set(1, combinedSite);

  // leave the new physical index alone
  auto identity = ITensor(n, nn);
  for (int j = 1; j <= q; j++) {
    identity.set(n=j, nn=j, 1);
  }

  // make the projector onto each measurement result
  for (int i = 0; i < q; i++) {
    auto proj =  ITensor(m, mm);
    // set the right index to one to project onto the result
    // the +1s are because ITensors are 1-indexed
    proj.set(m(i + 1), mm(i + 1), 1);
    auto total_proj = proj * identity;
    total_proj = C * total_proj * CC;
    projs.push_back(total_proj);
  }


  std::cout << "pos = 1 ";
  for (int i = 0; i < q; i++) {
    prob = std::pow(norm(combinedMPS(1) * projs[i]), 2);
    probs.push_back(prob + prev_prob);
    prev_prob += prob;
  }

  // pick a measurement result using a random_device
  std::random_device rand_dev;
  std::mt19937 generator(rand_dev());
  std::uniform_real_distribution<double> sample(0.0, 1.0);

  double result = sample(generator);
  // store the contracted LHS tensors in contracted
  ITensor contracted;
  // temporary storage for projected tensor
  ITensor projd;
  for (int i = 0; i < q; i++) {
    // pick i=1
    if (result <= probs[i]) {
      z = z.append(std::to_string(i));

      projd = combinedMPS(1) * projs[i];
      // PrintData(projs[1][i]);

      // create the projector for the "official" mps
      auto projop = ITensor(inds_to_meas[1]);
      projop.set(inds_to_meas[1](i + 1), 1);
      qudits[1] = (qudits[1] * projop);

      // normalize by the square root of the measurement probability
      // the norm doesn't have to be 1 and probably shouldn't be
      if (i > 0) {
        // println("prob: " + std::to_string(probs[i] - probs[i - 1]));
        projd = projd / sqrt((probs[i] - probs[i - 1]));
        qudits[1] = qudits[1] / sqrt((probs[i] - probs[i - 1]));
      } else {
        // println("prob: " + std::to_string(probs[i]));
        projd = projd / sqrt(probs[i]);
        qudits[1] = qudits[1] / sqrt(probs[i]);
      }
      combinedMPS.set(1, projd);
      // println("right after projecting onto result: ");
      // PrintData(combinedMPS(1));
      // println("norm of qudits[1]: " + std::to_string(norm(qudits[1])));
      // PrintData(qudits[1]);
      break;
    }
    if (i == (q - 1)) {
      projd = combinedMPS(1) * projs[i];
      // PrintData(projs[1][i]);

      // create the projector for the "official" mps
      auto projop = ITensor(inds_to_meas[1]);
      projop.set(inds_to_meas[1](i + 1), 1);
      qudits[1] = (qudits[1] * projop);

      // normalize by the square root of the measurement probability
      // the norm doesn't have to be 1 and probably shouldn't be
        // println("prob: " + std::to_string(probs[i]));
      projd = projd / sqrt(probs[i]);
      qudits[1] = qudits[1] / sqrt(probs[i]);
      combinedMPS.set(1, projd);
    }
  }

  // contract with itself and store
  // to store the conjugated tensor
  ITensor conjd = conj(prime(projd, "Link"));

  contracted = projd * conjd;
  ITensor temp_projd;
  ITensor temp_conjd;


  // now do everything from above for the rest of the positions
  for (int pos = 2; pos <= N; pos++){
    // create the projectors
    projs.clear();
    auto n = new_phys_inds[pos];
    auto nn = prime(n);
    auto m = inds_to_meas[pos];
    auto mm = prime(m);

    auto [C, c] = combiner({n, m}, {"Tags", "Site,i="+std::to_string(pos)});
    // combiners.push_back(C);
    // c.addTags("Site");
    auto [CC, cc] = combiner({nn, mm}, {"Tags", "prout,Site,i="+std::to_string(pos)});
    // p_combiners.push_back(CC);
    // comb_prime_inds.push_back(cc);

    combinedSite = qudits[pos] * C;
    combinedMPS.set(pos, combinedSite);

    // leave the new physical index alone
    auto identity = ITensor(n, nn);
    for (int j = 1; j <= q; j++) {
      identity.set(n=j, nn=j, 1);
    }

    // make the projector onto each measurement result
    for (int i = 0; i < q; i++) {
      auto proj =  ITensor(m, mm);
      // set the right index to one to project onto the result
      // the +1s are because ITensors are 1-indexed
      proj.set(m(i + 1), mm(i + 1), 1);
      auto total_proj = proj * identity;
      total_proj = C * total_proj * CC;
      projs.push_back(total_proj);
    }

    std::cout <<  " - " << std::to_string(pos);
    probs.clear();
    prev_prob = 0;
    for (int i = 0; i < q; i++) {
      temp_projd = combinedMPS(pos) * projs[i];
      temp_conjd = dag(temp_projd);
      for (auto& ind : temp_conjd.inds()) {
       if (hasIndex(contracted, ind)) {
         temp_conjd.prime(ind);
       }
     }

      auto temp = temp_projd * temp_conjd;
      auto prob_scalar = temp * contracted;
      prob = norm(prob_scalar);
      probs.push_back(prob + prev_prob);
      prev_prob += prob;
    }
    result = sample(generator);

    for (int i = 0; i < q; i++) {
      if (result <= probs[i]) {
        z = z.append(std::to_string(i));

        projd = combinedMPS(pos) * projs[i];
        float meas_prob;
        if (i == 0) {
          meas_prob = probs[i];
        } else {
          meas_prob = probs[i] - probs[i - 1];
        }
        projd = projd / sqrt(meas_prob);
        combinedMPS.set(pos, projd);

        // create the projector for the "official" mps
        auto projop = ITensor(inds_to_meas[pos]);
        projop.set(inds_to_meas[pos](i + 1), 1);

        qudits[pos] = qudits[pos] * projop / sqrt(meas_prob);
        goto project;
      }
      if (i == (q - 1)) {
        println("Hasn't picked a result!!!!!! ERROR!!!");
        println("result = " + std::to_string(result));
        println("last prob = " + std::to_string(probs.back()));
        projd = combinedMPS(pos) * projs[i];
        float meas_prob;
        if (i == 0) {
          meas_prob = probs[i];
        } else {
          meas_prob = probs[i] - probs[i - 1];
        }
        projd = projd / sqrt(meas_prob);
        combinedMPS.set(pos, projd);

        // create the projector for the "official" mps
        auto projop = ITensor(inds_to_meas[pos]);
        projop.set(inds_to_meas[pos](i + 1), 1);

        qudits[pos] = qudits[pos] * projop / sqrt(meas_prob);
        goto project;
      }
    }

    // contract with itself and store
    // to store the conjugated tensor
    project: projd = combinedMPS(pos);
    conjd = dag(prime(projd, "Link"));

    contracted = projd * (contracted * conjd);
  }
  std::cout << std::endl;
  // ofile << "z=" << z << std::endl;

  return *this;
}

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

// measure by contracting everything to the right instead of orthogonalizing
// not much faster
// slightly unclear if this works-- need to check
BrickMPS& BrickMPS::measure2(std::map<int, Index> inds_to_meas, std::map<int, Index> new_phys_inds, int orth_type) {
  auto copyMPS = this->to_mps();
  for (int i = 1; i <= this->get_len(); i++) {
    qudits[i] = copyMPS(i);
  }

  // make the projectors in preparation
  // this is an Nxq 2D array, first index gives the position along the MPS
  // second index is which measurement result to project onto
  std::map<int, std::vector<ITensor>> projs;
  // std::vector<Index> comb_prime_inds;
  std::vector<ITensor> p_combiners;
  std::vector<ITensor> combiners;
  // keep track of the contracted tensors on the right
  // stored in BACKWARDS ORDER
  std::vector<ITensor> rhs_contracted;
  // record the outcome string
  string z;
  int N = this->get_len();
  println("measuring");

  // make an MPS where the physical indices are combined, into one of dimension q^2
  // this is for ease of operations and gauging
  auto sites = SiteSet(N, q*q);
  auto combinedMPS = MPS(sites);

  for (int pos = 1; pos <= N; pos++) {
    // println("pos =" + std::to_string(pos));
    auto n = new_phys_inds[pos];
    auto nn = prime(n);
    auto m = inds_to_meas[pos];
    auto mm = prime(m);

    auto [C, c] = combiner({n, m}, {"Tags", "Site,i="+std::to_string(pos)});
    combiners.push_back(C);
    // c.addTags("Site");
    auto [CC, cc] = combiner({nn, mm}, {"Tags", "prout,Site,i="+std::to_string(pos)});
    p_combiners.push_back(CC);
    // comb_prime_inds.push_back(cc);

    auto combinedSite = qudits[pos] * C;
    // println("qudits[pos]: ");
    // PrintData(qudits[pos]);
    // println("combined site");
    // PrintData(combinedSite);
    combinedMPS.set(pos, combinedSite);

    // leave the new physical index alone
    auto identity = ITensor(n, nn);
    for (int j = 1; j <= q; j++) {
      identity.set(n=j, nn=j, 1);
    }
    // make the projector onto each measurement result
    for (int i = 0; i < q; i++) {
      // println("q = " +std::to_string(q));
      auto proj =  ITensor(m, mm);
      // set the right index to one to project onto the result
      // the +1s are because ITensors are 1-indexed
      proj.set(m(i + 1), mm(i + 1), 1);
      auto total_proj = proj * identity;
      total_proj = C * total_proj * CC;
      projs[pos].push_back(total_proj);
    }
  }


  // get the probability of each measurement outcome

  // first, contract everything on the right and keep track of it all
  // again, stored in reverse order, where the rightmost contractions are FIRST
  for (int pos = N; pos >= 1; pos--) {
    auto contracted = qudits[pos] * conj(qudits[pos]);
    // PrintData(qudits[pos] * conj(qudits[pos]));
    if (pos < N) {
      contracted = contracted * rhs_contracted[N - pos - 1];
      // PrintData(contracted);
    }
    rhs_contracted.push_back(contracted);
  }

  // treat the first position differently
  // probs is a vector of length q of the probabilities of gettnig a measurement
  // result <= to i, where i is the index in the vector
  std::vector<float> probs;
  float prob;
  float prev_prob = 0;

  // PrintData(combinedMPS(3));
  // PrintData(projs[1][0]);
  // PrintData(combinedMPS(1) * projs[1][0]);
  std::cout << "pos = 1 ";
  for (int i = 0; i < q; i++) {
    prob = std::pow(norm(combinedMPS(1) * projs[1][i]), 2);
    // println("right after temporarily projecting the first qudit: ");
    // PrintData(combinedMPS(1) * projs[1][i]);
    // PrintData(projs[1][i]);
    // PrintData(combinedMPS(1) * projs[1][i]);
    // println("prob of " + std::to_string(i) + ": " + std::to_string(prob));
    probs.push_back(prob + prev_prob);
    prev_prob += prob;
  }

  // pick a measurement result

  std::random_device rand_dev;
  std::mt19937 generator(rand_dev());
  std::uniform_real_distribution<double> sample(0.0, 1.0);

  double result = sample(generator);
  // TEST WITH ALL ZEROS
  // double result = 0;
  // println("result: " + std::to_string(result));

  // store the contracted LHS tensors in contracted
  ITensor contracted;
  // temporary storage for projected tensor
  ITensor projd;
  for (int i = 0; i < q; i++) {
    // pick i=1 for TESTING
    if (result <= probs[i]) {
      z = z.append(std::to_string(i));
      projd = combinedMPS(1) * projs[1][i];
      // PrintData(projs[1][i]);

      // create the projector for the "official" mps
      auto projop = ITensor(inds_to_meas[1]);
      projop.set(inds_to_meas[1](i + 1), 1);
      qudits[1] = (qudits[1] * projop);

      // normalize by the square root of the measurement probability
      // the norm doesn't have to be 1 and probably shouldn't be
      if (i > 0) {
        // println("prob: " + std::to_string(probs[i] - probs[i - 1]));
        projd = projd / sqrt((probs[i] - probs[i - 1]));
        qudits[1] = qudits[1] / sqrt((probs[i] - probs[i - 1]));
      } else {
        // println("prob: " + std::to_string(probs[i]));
        projd = projd / sqrt(probs[i]);
        qudits[1] = qudits[1] / sqrt(probs[i]);
      }
      combinedMPS.set(1, projd);
      // println("right after projecting onto result: ");
      // PrintData(combinedMPS(1));
      // println("norm of qudits[1]: " + std::to_string(norm(qudits[1])));
      // PrintData(qudits[1]);
      break;
    }
  }

  // contract with itself and store
  // to store the conjugated tensor
  // projd = combinedMPS.ref(1);
  ITensor conjd = conj(prime(projd, "Link"));

  // println("projd inds: ");
  // for(auto& ind : projd.inds()) {
  //   println(ind);
  // }
  // println("conjd inds: ");
  // for(auto& ind : conjd.inds()) {
  //   println(ind);
  // }

  contracted = projd * conjd;
  // PrintData(combinedMPS(1));
  // PrintData(contracted);

  ITensor temp_projd;
  ITensor temp_conjd;

  // POTENTIAL ISSUE: not setting the combinedMPS to the other one
  for (int pos = 2; pos <= N; pos++){
    std::cout <<  " - " << std::to_string(pos);
    // println("pos: " + std::to_string(pos));
    probs.clear();
    // combinedMPS.position(pos);
    prev_prob = 0;
    for (int i = 0; i < q; i++) {
      temp_projd = combinedMPS(pos) * projs[pos][i];
      // println("norm of temp_projd = " + std::to_string(norm(temp_projd)));
      temp_conjd = dag(temp_projd);
      for (auto& ind : temp_conjd.inds()) {
       if (hasIndex(contracted, ind)) {
         temp_conjd.prime(ind);
       }
     }
     // println("Right after temporarily projecting: ");
    // PrintData(temp_projd);
    // PrintData(temp_conjd);

      auto prob_scalar = temp_projd * (contracted * temp_conjd);
      prob = norm(prob_scalar);
      if (pos == 11) {
        // println("prob of " + std::to_string(i) + ": " + std::to_string(prob));
      }
      // println("prob of " + std::to_string(i) + ": " + std::to_string(prob));
      // println("cumulative prob of " + std::to_string(i) + ": " + std::to_string(prev_prob + prob));
      probs.push_back(prob + prev_prob);
      prev_prob += prob;
    }
    // UNCOMMENT WHEN YOU DON'T WANT ALL ZEROS
    result = sample(generator);
    // result = 0;
    // println("last prob= " + std::to_string(probs[q-1]));

    for (int i = 0; i < q; i++) {
      if (result <= probs[i]) {
        z = z.append(std::to_string(i));

        projd = combinedMPS(pos) * projs[pos][i];
        float meas_prob;
        if (i == 0) {
          meas_prob = probs[i];
        } else {
          meas_prob = probs[i] - probs[i - 1];
        }
        projd = projd / sqrt(meas_prob);
        combinedMPS.set(pos, projd);
        // println("Right after projecting onto the result: ");
        // PrintData(combinedMPS(pos));
        // combinedMPS.normalize();
        // PrintData(combinedMPS(pos));

        // create the projector for the "official" mps
        auto projop = ITensor(inds_to_meas[pos]);
        projop.set(inds_to_meas[pos](i + 1), 1);

        qudits[pos] = qudits[pos] * projop / sqrt(meas_prob);
        // PrintData(qudits[pos]);
        goto project;
      }
      if (i == (q - 1)) {
        println("Hasn't picked a result!!!!!! ERROR!!!");
        println("result = " + std::to_string(result));
        println("last prob = " + std::to_string(probs.back()));
      }
    }

    // contract with itself and store
    // to store the conjugated tensor
    project: projd = combinedMPS(pos);
    conjd = dag(prime(projd, "Link"));

    // PrintData(projd);
    // PrintData(conjd);
    contracted = projd * (contracted * conjd);
    // PrintData(contracted);
    // println("contracted norm for pos=" + std::to_string(pos) + ": " + std::to_string(norm(contracted)));
  }
  std::cout << std::endl;
  // ofile << "z=" << z << std::endl;

  return *this;
}

// Truncate the MPS to a maximum bond dimension or error and output the entanglement data to a file
// maxerror: the maximum error after truncation and orthogonalizing if using that method
// maxdim: the maximum bond dimension in the MPS after orthogonalizing, if using that method
// method: 1 for truncating based on error, anything else for bond dimension
// ofile: the file which the entanglement data is written to
BrickMPS& BrickMPS::truncate(float maxerror, int maxdim, int method, std::ofstream& ofile) {
  auto copyMPS = this->to_mps();
  auto otherCopyMPS = this->to_mps();
  println("Max link dim: before = " + std::to_string(maxLinkDim(copyMPS)));
  if (method == 1) {
    copyMPS.orthogonalize({"Cutoff", maxerror});
  } else {
    copyMPS.orthogonalize({"MaxDim", maxdim});
  }
  copyMPS.normalize();

  // get entanglement
  entanglement_data ent = get_entanglement(copyMPS);

  ofile << "entanglement entropies" << std::endl;

  int n = this->get_len();
  for (int i = 0; i < n - 1; i++) {
    // std::cout << ent.entropies[i] << " - ";
    ofile << ent.entropies[i] << ",";
  }
  // sort the singular values
  std::sort (ent.svals[int(n/2)].begin(), ent.svals[int(n/2)].end());
  ofile << std::endl << "singular values in the middle: " << std::endl;
  for (int i = 0; i < ent.svals[n/2].size(); i++) {
    ofile << ent.svals[int(n / 2)][i] << ",";
    // ofile << ent.entropies[i] << ",";
  }
  ofile << std::endl;

  // println("norm before normalizing = " + std::to_string(norm(copyMPS)));
  // println("Max link dim without truncation: " + std::to_string(maxLinkDim(otherCopyMPS)));
  println("Max link dim: after = " + std::to_string(maxLinkDim(copyMPS)));
  for (int i = 1; i <= n; i++) {
    qudits[i] = copyMPS(i);
  }
  return *this;
}

// Truncate the MPS to a maximum truncation error and write the overlap
// between the newly truncated MPS and the old one into ofile
// Mostly used for testing to see if we truncate too much
BrickMPS& BrickMPS::truncate_err(float maxerror, std::ofstream& ofile) {
  auto copyMPS = this->to_mps();
  auto otherCopyMPS = this->to_mps();
  println("Max link dim: before = " + std::to_string(maxLinkDim(copyMPS)));
  copyMPS.orthogonalize({"Cutoff", maxerror});
  copyMPS.normalize();
  otherCopyMPS.orthogonalize();
  otherCopyMPS.normalize();

  auto overlap = norm(innerC(copyMPS, otherCopyMPS));
  ofile << "overlap" << std::endl << std::to_string(overlap) << std::endl;
  println("overlap: " + std::to_string(overlap));

  // println("norm before normalizing = " + std::to_string(norm(copyMPS)));
  copyMPS.normalize();
  // println("norm of difference after truncation: " + std::to_string(norm(copyMPS - otherCopyMPS)));
  println("Max link dim: after = " + std::to_string(maxLinkDim(copyMPS)));
  for (auto& ind : linkInds(copyMPS)) {
    // println(ind);
  }
  for (int i = 1; i <= this->get_len(); i++) {
    qudits[i] = copyMPS(i);
  }
  return *this;
}

// same as truncate_err above but orthogonalizing to a maximum dimension instead
BrickMPS& BrickMPS::truncate_maxdim(int maxdim, std::ofstream& ofile) {
  auto copyMPS = this->to_mps();
  auto otherCopyMPS = this->to_mps();
  println("Max link dim: before = " + std::to_string(maxLinkDim(copyMPS)));
  copyMPS.orthogonalize({"MaxDim", maxdim});
  copyMPS.normalize();
  otherCopyMPS.orthogonalize();
  otherCopyMPS.normalize();

  auto overlap = norm(innerC(copyMPS, otherCopyMPS));
  ofile << "overlap" << std::endl << std::to_string(overlap) << std::endl;
  println("overlap: " + std::to_string(overlap));

  // println("norm before normalizing = " + std::to_string(norm(copyMPS)));
  copyMPS.normalize();
  // println("norm of difference after truncation: " + std::to_string(norm(copyMPS - otherCopyMPS)));
  println("Max link dim: after = " + std::to_string(maxLinkDim(copyMPS)));
  for (auto& ind : linkInds(copyMPS)) {
    // println(ind);
  }
  for (int i = 1; i <= this->get_len(); i++) {
    qudits[i] = copyMPS(i);
  }
  return *this;
}

// Do one iteration of the SEBD algorithm
//  1. create a new MPS representing the new qudits
//  2. randomize it
//  3. glue it according to the architecture to the old MPS with random unitaries
//  4. measure the old physical indices
// iter: which iteration of the algorithm we're on (important for keeping track
//       of where we are in the architecture, for gluing tensors together)
// udata: the unitary gates of the circuit
// circuit_type: 1 for brickwork, 2 for cluster state, 3 for more-entangling
//               brickwork, 4 for less-entangling brickwork
// orth_type: gets passed to measure()
BrickMPS& BrickMPS::iterate(int iter, UnitaryData udata, int circuit_type, int orth_type) {
  int n = this->get_len();
  // 1. create a new MPS for the new qudits
  BrickMPS* newQudits = new BrickMPS(std::to_string(iter), n, q, true);
  auto newMPS = newQudits->to_mps();
  delete newQudits;

  // 2. randomize it by applying random unitaries between each adjacent qudit pair
  for (int pos = 1; pos < n; pos++) {
    auto T1 = newMPS(pos);
    auto T2 = newMPS(pos + 1);
    auto combined = this->apply_unitary(T1, T2, udata[pos]);
    // un-stick the adjacent tensors by SVDing them
    // get the indices for the rows of the SVD
    Index site1 = siteIndex(newMPS, pos);
    ITensor U, A;
    if (pos > 1) {
      Index bond = leftLinkIndex(newMPS, pos);
      auto [W, S, V] = svd(combined, {site1, bond}, {"Cutoff", 1E-18});
      U = W;
      A = S * V;
    } else {
      auto [W, S, V] = svd(combined, {site1}, {"Cutoff", 1E-18});
      A = S * V;
      U = W;
    }
    newMPS.set(pos, U);
    newMPS.set(pos + 1, A);
    // add the "bond=i" tag to the newly created index - optional?
    // rightLinkIndex(newMPS, pos).addTags("bond=" + pos);

  }
  newMPS.position(1);

  std::map<int, Index> inds_to_meas;
  std::map<int, Index> new_phys_inds;


  // 3. Attach the chain of new qudits to the old ones according to the architecture
  // counter keeps track of which unitary matrix we're on
  int counter = 1;
  // 1 means bricwork
  if (circuit_type == 1) {
    for (int pos = 1; pos <= n; pos++) {
      // println("applying a unitary to pos " + std::to_string(pos));
      // keep track of which indices to measure and which to keep
      for (auto& ind : inds(qudits[pos])) {
        if (hasTags(ind, "Site")) inds_to_meas[pos] = ind;
      }
      for (auto& ind : inds(newMPS(pos))) {
        if (hasTags(ind, "Site")) new_phys_inds[pos] = ind;
      }
      if (iter % 2 == 1) {
        if (((pos - 1) % 8 == 1) || ((pos - 1) % 8 == 3)) {
          qudits[pos] = this->apply_unitary(qudits[pos], newMPS(pos), udata[n + counter]);
          counter++;
        } else {
          qudits[pos] = qudits[pos] * newMPS(pos);
        }
      } else {
        if (((pos - 1) % 8 == 5) || ((pos - 1) % 8 == 7)) {
          qudits[pos] = this->apply_unitary(qudits[pos], newMPS(pos), udata[n + counter]);
          counter++;
        } else {
          qudits[pos] = qudits[pos] * newMPS(pos);
        }
      }
    }
  } else if (circuit_type == 2) {
    // circuit_type == 2 means that we're doing a cluster state (NOT BRICKWORK)
    for (int pos = 1; pos <= n; pos++) {
      // keep track of which indices to measure and which to keep
      for (auto& ind : inds(qudits[pos])) {
        if (hasTags(ind, "Site")) inds_to_meas[pos] = ind;
      }
      for (auto& ind : inds(newMPS(pos))) {
        if (hasTags(ind, "Site")) new_phys_inds[pos] = ind;
      }
      // apply the gates to EVERY qudit
      qudits[pos] = this->apply_unitary(qudits[pos], newMPS(pos), udata[n + counter]);
      counter++;
    }
  } else if (circuit_type == 3) {
    // do an in-between where instead of brickwork you do every other one (less entanglement?)
    for (int pos = 1; pos <= n; pos++) {
      // keep track of which indices to measure and which to keep
      for (auto& ind : inds(qudits[pos])) {
        if (hasTags(ind, "Site")) inds_to_meas[pos] = ind;
      }
      for (auto& ind : inds(newMPS(pos))) {
        if (hasTags(ind, "Site")) new_phys_inds[pos] = ind;
      }
      if (((pos - 1) % 8 == 1) || ((pos - 1) % 8 == 3) || ((pos - 1) % 8 == 5) || ((pos - 1) % 8 == 7)) {
        qudits[pos] = this->apply_unitary(qudits[pos], newMPS(pos), udata[n + counter]);
        counter++;
      } else {
        qudits[pos] = qudits[pos] * newMPS(pos);
      }
    }
  } else if (circuit_type == 4) {
    // do an even less entangling architecture, only HALF the connections of bw
    for (int pos = 1; pos <= n; pos++) {
      // keep track of which indices to measure and which to keep
      for (auto& ind : inds(qudits[pos])) {
        if (hasTags(ind, "Site")) inds_to_meas[pos] = ind;
      }
      for (auto& ind : inds(newMPS(pos))) {
        if (hasTags(ind, "Site")) new_phys_inds[pos] = ind;
      }
      if (iter % 2 == 1) {
        if (((pos - 1) % 8 == 1)) {
          qudits[pos] = this->apply_unitary(qudits[pos], newMPS(pos), udata[n + counter]);
          counter++;
        } else {
          qudits[pos] = qudits[pos] * newMPS(pos);
        }
      } else {
        if (((pos - 1) % 8 == 5)) {
          qudits[pos] = this->apply_unitary(qudits[pos], newMPS(pos), udata[n + counter]);
          counter++;
        } else {
          qudits[pos] = qudits[pos] * newMPS(pos);
        }
      }
    }
  }

  // combine link indices
  for (int pos = 1; pos < n; pos++) {
    auto common_links = commonInds(qudits[pos], qudits[pos + 1]);
    auto[C, c] = combiner(common_links);

    qudits[pos] = qudits[pos] * C;
    qudits[pos + 1] = qudits[pos + 1] * C;
  }


  // measure the old indices
  *this = this->measure(inds_to_meas, new_phys_inds, orth_type);

  return *this;
}
