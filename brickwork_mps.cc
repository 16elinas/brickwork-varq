#include "brickwork_mps.h"


// Create a new chain of qudits
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
    // qudits[1] = qudits[1] / norm(qudits[1]);

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
    // qudits[len] = qudits[len] / norm(qudits[len]);

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


// output tht entanglement data of the current state of the circuit
// in the form of entanglement_data which has a list of the
// total von Neumann entropies across the bonds and also
// the lists of all the singular values
// the lists are ZERO INDEXED
entanglement_data BrickMPS::get_entanglement() {
  int N = this->get_len();
  entanglement_data out;
  std::vector<std::vector<double>> svals;
  std::vector<double> entropies;

  // auto psi = this->to_mps();

  // variables for each iteration
  double ent;
  std::vector<double> bond_svals;

  // for (int b = 1; b < this->get_len(); b++) {
  //   ent = 0;
  //   bond_svals.clear();
  //
  //   // put the MPS in the right position
  //   psi.position(b);
  //
  //   //SVD this wavefunction to get the spectrum
  //   //of density-matrix eigenvalues
  //   auto l = leftLinkIndex(psi,b);
  //   auto s = siteIndex(psi,b);
  //   auto [U,S,V] = svd(psi(b),{l,s});
  //   auto u = commonIndex(U,S);
  //
  //   //Apply von Neumann formula
  //   //to the squares of the singular values
  //   println("position: " + std::to_string(b));
  //   for(auto n : range1(dim(u))) {
  //     auto Sn = elt(S,n,n);
  //     auto val = sqr(Sn);
  //     bond_svals.push_back(val);
  //     println("singular value: " + std::to_string(val));
  //     if(val > 1E-12) ent += -val*log(val);
  //   }
  //   svals.push_back(bond_svals);
  //   entropies.push_back(ent);
  // }

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

// apply a unitary between two ITensors and return the product
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
  // PrintData(U);

  // println("norm of U=" + std::to_string(norm(U)));
  // println("norm of T1=" + std::to_string(norm(T1)));
  // println("norm of T2=" + std::to_string(norm(T2)));


   auto out = T1 * T2;
   // println("norm of out=" + std::to_string(norm(out)));
   out *= U;
   out.noPrime();

   // PrintData(out);

   return out;
}


// measure the indices to be measured and orthogonalize
// inds_to_meas: a 1-indexed array of the indices to be measured
// new_phys_inds: the indices to be left alone
// for each index:
// 1. pick a measurement result, based on magnitude of projected qudit
// 2. project the qudit onto that result, for the index to be measured and normalize
// 3. contract into the next qudit to prepare
BrickMPS& BrickMPS::measure(std::map<int, Index> inds_to_meas, std::map<int, Index> new_phys_inds) {
  // make the projectors in preparation
  // this is an Nxq 2D array, first index gives the position along the MPS
  // second index is which measurement result to project onto
  std::map<int, std::vector<ITensor>> projs;
  // std::vector<Index> comb_prime_inds;
  std::vector<ITensor> p_combiners;
  std::vector<ITensor> combiners;
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

  // for (int i = 1; i <= N; i++) {
  //   PrintData(combinedMPS(i));
  //   // println("norm of combinedMPS(" + std::to_string(i) + "): " + std::to_string(norm(combinedMPS(i))));
  // }
  // combinedMPS.orthogonalize();
  // combinedMPS.normalize();
  // combinedMPS.position(1);
  // for (int i = 1; i <= N; i++) {
  //   PrintData(combinedMPS(i));
  //   // println("norm of combinedMPS(" + std::to_string(i) + "): " + std::to_string(norm(combinedMPS(i))));
  // }

  // for (auto& ind : linkInds(combinedMPS)) {
  //   println(ind);
  // }

  // get the probability of each measurement outcome

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
    if (result < probs[i]) {
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
    std::cout <<  " - " << std::to_string(pos) << std::endl;
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

      auto prob_scalar = temp_projd * contracted * temp_conjd;
      prob = norm(prob_scalar);
      // println("prob of " + std::to_string(i) + ": " + std::to_string(prob));
      // println("cumulative prob of " + std::to_string(i) + ": " + std::to_string(prev_prob + prob));
      probs.push_back(prob + prev_prob);
      prev_prob += prob;
    }
    // UNCOMMENT WHEN YOU DON'T WANT ALL ZEROS
    result = sample(generator);
    // println("last prob= " + std::to_string(probs[q-1]));

    for (int i = 0; i < q; i++) {
      if (result <= probs[i]) {
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

  return *this;
}

// Truncate the MPS to a maximum bond dimension
//TODO: this is not fully done, do it
BrickMPS& BrickMPS::truncate(int maxdim) {
  auto copyMPS = this->to_mps();
  auto otherCopyMPS = this->to_mps();
  println("Max link dim: before = " + std::to_string(maxLinkDim(copyMPS)));
  copyMPS.orthogonalize({"MaxDim", maxdim});
  // otherCopyMPS.orthogonalize();
  copyMPS.normalize();
  // println("Max link dim without truncation: " + std::to_string(maxLinkDim(otherCopyMPS)));
  println("Max link dim: after = " + std::to_string(maxLinkDim(copyMPS)));
  // println("difference in norms after truncation: " + std::to_string(norm(this->to_mps() - copyMPS)));
  for (int i = 1; i <= this->get_len(); i++) {
    qudits[i] = copyMPS(i);
  }
  return *this;
}

// Truncate the MPS to a maximum truncation error
BrickMPS& BrickMPS::truncate(float maxerror) {
  auto copyMPS = this->to_mps();
  // auto origMPS = copyMPS;
  // origMPS.orthogonalize();
  // println("Norm before forceful normalization: " + std::to_string(norm(origMPS)));
  println("Max link dim: before = " + std::to_string(maxLinkDim(copyMPS)));
  copyMPS.orthogonalize({"Cutoff", maxerror});
  // origMPS.normalize();
  copyMPS.normalize();
  // println("norm of difference after truncation: " + std::to_string(norm(origMPS - copyMPS)));
  println("Max link dim: after = " + std::to_string(maxLinkDim(copyMPS)));
  for (int i = 1; i <= this->get_len(); i++) {
    qudits[i] = copyMPS(i);
  }
  return *this;
}

// do one iteration of the SEBD algorithm
// 1. create a new MPS representing the new qudits
// 2. randomize it
// 3. glue it to the old MPS with random unitaries according to the architecture
// 4. measure the old physical indices
BrickMPS& BrickMPS::iterate(int iter, UnitaryData udata) {
  int n = this->get_len();
  // 1. create a new MPS for the new qudits
  BrickMPS newQudits = BrickMPS(std::to_string(iter), n, q, true);
  auto newMPS = newQudits.to_mps();


  // print all the qudit indices
  // for (int pos = 1; pos <= n; pos++) {
  //   println("----qudits[" + std::to_string(pos) + "]:");
  //   for (auto& index : inds(qudits[pos])) {
  //     println(index);
  //   }
  //   // PrintData(qudits[pos]);
  //   println("norm=" + std::to_string(norm(qudits[pos])));
  // }

  // 2. randomize it by applying random unitaries between each adjacent qudit pair
  for (int pos = 1; pos < n; pos++) {
    auto T1 = newMPS(pos);
    auto T2 = newMPS(pos + 1);
    auto combined = this->apply_unitary(T1, T2, udata[pos]);
    // println("norm of post-unitary: " + std::to_string(norm(combined)));
    // PrintData(combined);
    // un-stick them by SVDing it
    // get the indices for the rows of the SVD
    Index site1 = siteIndex(newMPS, pos);
    ITensor U, A;
    if (pos > 1) {
      Index bond = leftLinkIndex(newMPS, pos);
      auto [W, S, V] = svd(combined, {site1, bond}, {"Cutoff", 1E-15});
      U = W;
      A = S * V;
    } else {
      auto [W, S, V] = svd(combined, {site1}, {"Cutoff", 1E-15});
      A = S * V;
      U = W;
    }
    newMPS.set(pos, U);
    // println("norm of post-unitary U: " + std::to_string(norm(combined)));

    newMPS.set(pos + 1, A);
    // println("norm of post-unitary A: " + std::to_string(norm(combined)));

    // add the "bond=i" tag to the newly created index
    // rightLinkIndex(newMPS, pos).addTags("bond=" + pos);

  }
  newMPS.orthogonalize();
  // println(norm(newMPS));
  // newMPS.normalize();


  // for (int ipos = 1; ipos <= n; ipos++) {
  //   println("----newMPS(" + std::to_string(ipos) + "):");
  //   for (auto& index : inds(newMPS(ipos))) {
  //     println(index);
  //   }
  //   println("norm=" + std::to_string(norm(newMPS(ipos))));
  // }



  std::map<int, Index> inds_to_meas;
  std::map<int, Index> new_phys_inds;


  // 3. Attach the chain of new qudits to the old ones according to the architecture
  // counter keeps track of which unitary matrix we're on
  int counter = 1;
  for (int pos = 1; pos <= n; pos++) {
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

  // combine link indices
  for (int pos = 1; pos < n; pos++) {
    auto common_links = commonInds(qudits[pos], qudits[pos + 1]);
    auto[C, c] = combiner(common_links);

    qudits[pos] = qudits[pos] * C;
    qudits[pos + 1] = qudits[pos + 1] * C;
  }



  // // print all the qudit indices
  // for (int pos = 1; pos <= n; pos++) {
  //   println("----qudits[" + std::to_string(pos) + "]:");
  //   for (auto& index : inds(qudits[pos])) {
  //     println(index);
  //   }
  //   // PrintData(qudits[pos]);
  //   println("norm=" + std::to_string(norm(qudits[pos])));
  // }

  // measure the indices
  *this = this->measure(inds_to_meas, new_phys_inds);

  // truncate the bonds


  return *this;
}
