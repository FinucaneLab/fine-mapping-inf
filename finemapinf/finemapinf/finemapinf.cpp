#include "finemapinf.hpp"
#include <math.h>
#include <stdlib.h>
#include <limits>
#include <vector>
#include <queue>
#include <iostream>
#include <chrono>

// Use external Eigen library For computing matrix inverse and determinant
#include <Eigen/Dense>
#include <Eigen/LU>

// Time individual components of the computation, for code optimization
// #define TIMEIT

using namespace FINEMAPINF;
using std::chrono::steady_clock;
using std::chrono::duration_cast;
using std::chrono::microseconds;

#ifdef TIMEIT
steady_clock::time_point START;
steady_clock::time_point END;
microseconds SSS_T;
microseconds MODEL_T;
microseconds PAIR_T;
microseconds XOMEGA_T;
microseconds VDV_T;
microseconds LOGQ_T;
microseconds MOM_T;
microseconds PIP_T;
#endif

// Add two values in log-space
inline double logaddexp(double a, double b) {
  return (a>b?a:b)+log(1+exp(-fabs(b-a)));
}

// For printing models to std::out
std::ostream& operator<<(std::ostream& os, const Model& m) {
  os << "(";
  for (unsigned i = 0; i < m.L; ++i) {
    os << m[i];
    if (i+1 != m.L) os << ",";
  }
  os << ")";
  return os;
}

// // For printing hash tables to std::out
// std::ostream& operator<<(std::ostream& os, const ModelHash& tab) {
//   for (unsigned i = 0; i < tab.used.size(); ++i) {
//     os << "  hash " << tab.used[i] << std::endl;
//     Model* m = tab.keys[tab.used[i]];
//     while (m != NULL) {
//       os << "    model " << *m << std::endl;
//       os << "    logq " << m->logq << std::endl;
//       m = m->next;
//     }
//   }
//   return os;
// }
// 
// std::ostream& operator<<(std::ostream& os, const PairHash& tab) {
//   for (unsigned i = 0; i < tab.used.size(); ++i) {
//     os << "  hash " << tab.used[i] << std::endl;
//     Pair* p = tab.keys[tab.used[i]];
//     while (p != NULL) {
//       os << "    pair " << p->i << " " << p->j << std::endl;
//       os << "    xOmegax " << p->xOmegax << std::endl;
//     }
//   }
//   return os;
// }

FINEMAP::FINEMAP(unsigned _n, unsigned _p, double _sigmasq, double _tausq,
        const double* _pi0, const double* _S, unsigned _nssq,
        unsigned _Lmax, double _meansq, const double* _z,
        const double* _V, const double* _Dsq, int _verbose,
        unsigned _models_nbins, unsigned _pairs_nbins,
        double* _PIP, double* _beta, double* _se, double* _alpha) :
    n(_n),p(_p),Lmax(_Lmax),nssq(_nssq),sigmasq(_sigmasq),tausq(_tausq),
    logqsum(NEGINF),
    verbose(_verbose),V(_V),Dsq(_Dsq),S(_S),PIP(_PIP),beta(_beta),se(_se),
    alpha(_alpha), pairs(_pairs_nbins),models(_models_nbins) {
  ysq = _n*_meansq;
  Xysq = 0;
  Omega = new double[p];
  Xy = new double[p];
  VXy = new double[p];
  XXXy = new double[p];
  xOmegay = new double[p];
  logodds = new double[p];
  logpi0c = NEGINF;
  double sqrtn = sqrt(n);
  // Precompute log(pi0) and log(1-pi0)
  // Precompute X'y, ||X'y||_2^2, and A-matrix for method-of-moments
  MoMA[0] = n;
  MoMA[1] = 0;
  MoMA[3] = 0;
  for (unsigned i = 0; i < p; ++i) {
    logpi0c = logaddexp(logpi0c,log(1-_pi0[i]));
    logodds[i] = (_pi0[i] > 0 ? log(_pi0[i]/(1-_pi0[i])) : NEGINF);
    Omega[i] = 1.0/(tausq*Dsq[i]+sigmasq);
    xOmegay[i] = NaN;
    Xy[i] = sqrtn*_z[i];
    Xysq += _z[i]*_z[i]*n;
    MoMA[1] += Dsq[i];
    MoMA[3] += Dsq[i]*Dsq[i];
  }
  MoMA[2] = MoMA[1];
  // Accumulate unnormalized values in b-vector for method-of-moments
  MoMb[0] = 0;
  MoMb[1] = 0;
  // Precompute V'X'y
  for (unsigned i = 0; i < p; ++i) {
    VXy[i] = 0;
    for (unsigned j = 0; j < p; ++j)
      VXy[i] += V[j*p+i]*Xy[j];
  }
  // Precompute (X'X)(X'y)
  for (unsigned i = 0; i < p; ++i) {
    XXXy[i] = 0;
    for (unsigned j = 0; j < p; ++j)
      XXXy[i] += V[i*p+j]*Dsq[j]*VXy[j];
  }
  // Initialize temporary storage
  logpissq = new double[_nssq];
  mus = new double[_Lmax*_nssq];
  omegainvs = new double[_Lmax*_Lmax*_nssq];
  Eb = new double[_Lmax];
  Ebb = new double[_Lmax*_Lmax];
  nbrs = new Model[p*(_Lmax+1)];
  for (unsigned i = 0; i < p*(_Lmax+1); ++i)
    nbrs[i].inds = new unsigned[_Lmax];
  last_logqsum = new double[p];
  // Initialize outputs
  for (unsigned i = 0; i < p; ++i) {
    last_logqsum[i] = NEGINF;
    PIP[i] = 0;
    beta[i] = 0;
    se[i] = 0;
  }
}

FINEMAP::~FINEMAP() {
  delete[] Omega;
  delete[] Xy;
  delete[] VXy;
  delete[] XXXy;
  delete[] logpissq;
  delete[] mus;
  delete[] omegainvs;
  delete[] logodds;
  delete[] Eb;
  delete[] Ebb;
  delete[] xOmegay;
  for (unsigned i = 0; i < p*(Lmax+1); ++i)
    delete[] nbrs[i].inds;
  delete[] nbrs;
  delete[] last_logqsum;
}

double FINEMAP::compute_xOmegax(unsigned i, unsigned j) {
#ifdef TIMEIT
  START = steady_clock::now();
#endif
  double val = 0;
  for (unsigned k = 0; k < p; ++k)
    val += V[i*p+k]*V[j*p+k]*Dsq[k]*Omega[k];
#ifdef TIMEIT
  END = steady_clock::now();
  XOMEGA_T += duration_cast<microseconds>(END - START);
#endif
  return val;
}

double FINEMAP::compute_xOmegay(unsigned j) {
#ifdef TIMEIT
  START = steady_clock::now();
#endif
  double val = 0;
  for (unsigned k = 0; k < p; ++k)
    val += V[j*p+k]*VXy[k]*Omega[k];
#ifdef TIMEIT
  END = steady_clock::now();
  XOMEGA_T += duration_cast<microseconds>(END - START);
#endif
  return val;
}

double FINEMAP::compute_vD2v(unsigned i, unsigned j) {
#ifdef TIMEIT
  START = steady_clock::now();
#endif
  double val = 0;
  for (unsigned k = 0; k < p; ++k)
    val += V[i*p+k]*V[j*p+k]*Dsq[k];
#ifdef TIMEIT
  END = steady_clock::now();
  VDV_T += duration_cast<microseconds>(END - START);
#endif
  return val;
}

double FINEMAP::compute_vD4v(unsigned i, unsigned j) {
#ifdef TIMEIT
  START = steady_clock::now();
#endif
  double val = 0;
  for (unsigned k = 0; k < p; ++k)
    val += V[i*p+k]*V[j*p+k]*Dsq[k]*Dsq[k];
#ifdef TIMEIT
  END = steady_clock::now();
  VDV_T += duration_cast<microseconds>(END - START);
#endif
  return val;
}

double FINEMAP::compute_logq(const Model& Gamma, bool update_MoM) {
  // Compute or look up x_i'Omega x_j, x_i'Omega y, v_i'D^2v_j, v_i'D^4v_j
  unsigned L = Gamma.L;
  Pair** ptab = new Pair*[L*L];
  for (unsigned i = 0; i < L; ++i) {
    if (std::isnan(xOmegay[Gamma[i]]))
      xOmegay[Gamma[i]] = compute_xOmegay(Gamma[i]);
    for (unsigned j = i; j < L; ++j) {
#ifdef TIMEIT
      START = steady_clock::now();
#endif
      ptab[i*L+j] = pairs.find(Gamma[i],Gamma[j]);
#ifdef TIMEIT
      END = steady_clock::now();
      PAIR_T += duration_cast<microseconds>(END - START);
#endif
      if (ptab[i*L+j] == NULL) {
        double v1 = compute_xOmegax(Gamma[i],Gamma[j]);
        // Omit computation of v_i'D^2v_j, v_i'D^4v_j if not update_MoM
        double v2 = (update_MoM ? compute_vD2v(Gamma[i],Gamma[j]) : NaN);
        double v3 = (update_MoM ? compute_vD4v(Gamma[i],Gamma[j]) : NaN);
#ifdef TIMEIT
        START = steady_clock::now();
#endif
        ptab[i*L+j] = pairs.insert(Gamma[i],Gamma[j],v1,v2,v3);
#ifdef TIMEIT
        END = steady_clock::now();
        PAIR_T += duration_cast<microseconds>(END - START);
#endif
      }
    }
  }
  // Compute model log(q) and posterior first and second moment of beta
#ifdef TIMEIT
  START = steady_clock::now();
#endif
  Eigen::MatrixXd omega(L,L);
  Eigen::MatrixXd RHS = Eigen::MatrixXd::Zero(L,L+1);
  for (unsigned i = 0; i < L; ++i) {
    RHS(i,0) = xOmegay[Gamma[i]];
    RHS(i,i+1) = 1;
    for (unsigned j = i; j < L; ++j) {
      omega(i,j) = ptab[i*L+j]->xOmegax;
      omega(j,i) = ptab[i*L+j]->xOmegax;
    }
  }
  double logpisum = NEGINF;
  Eigen::MatrixXd AinvB(L,L+1);
  for (unsigned i = 0; i < nssq; ++i) {
    logpissq[i] = 0;
    if (L > 0) {
      // Set omega = X_Gamma' Omega X_Gamma + diag(1/S[i])
      double delta = 1/S[i];
      if (i > 0) delta -= 1/S[i-1];
      for (unsigned j = 0; j < L; ++j) omega(j,j) += delta;
      // Simultaneously compute omega^{-1} (X_Gamma' Omega y) and omega^{-1}
      AinvB = omega.ldlt().solve(RHS);
      // Save mu, omega^{-1}, and unnormalized log(pi(s^2 | Gamma))
      for (unsigned j = 0; j < L; ++j) {
        mus[j*nssq+i] = AinvB(j,0);
        for (unsigned k = 0; k < L; ++k)
          omegainvs[j*nssq*L+k*nssq+i] = AinvB(j,k+1);
        logpissq[i] += 0.5*AinvB(j,0)*RHS(j,0);
      }
      logpissq[i] -= 0.5*(log(omega.determinant())+L*log(S[i]));
    }
    logpisum = logaddexp(logpisum, logpissq[i]);
  }
  // Value of log(q)
  double logqval = logpi0c+logpisum;
  for (unsigned i = 0; i < L; ++i)
    logqval += logodds[Gamma[i]];
  // Update logqsum and total posterior moments of b
  double prev_logqsum = logqsum;
  logqsum = logaddexp(logqsum,logqval);
  for (unsigned i = 0; i < nssq; ++i) {
    // Use temporarily to save pi(s^2 | Gamma) with normalization
    logpissq[i] = exp(logpissq[i]-logpisum);
  }
  for (unsigned i = 0; i < L; ++i) {
    Eb[i] = 0;
    for (unsigned k = 0; k < nssq; ++k)
      Eb[i] += logpissq[k]*mus[i*nssq+k];
    for (unsigned j = i; j < L; ++j) {
      Ebb[i*L+j] = 0;
      for (unsigned k = 0; k < nssq; ++k) {
        Ebb[i*L+j] += logpissq[k]*mus[i*nssq+k]*mus[j*nssq+k];
        Ebb[i*L+j] += logpissq[k]*omegainvs[i*nssq*L+j*nssq+k];
      }
      if (!update_MoM) break; // Only update diagonal j=i if not update_MoM
    }
  }
#ifdef TIMEIT
  END = steady_clock::now();
  LOGQ_T += duration_cast<microseconds>(END - START);
#endif
  // Update b-vector for MoM
#ifdef TIMEIT
  START = steady_clock::now();
#endif
  double wt = exp(logqval-logqsum);
  double rewt = exp(prev_logqsum-logqsum);
  if (update_MoM) {
    MoMb[0] *= rewt;
    MoMb[1] *= rewt;
    for (unsigned i = 0; i < L; ++i) {
      MoMb[0] -= 2*wt*Eb[i]*Xy[Gamma[i]];
      MoMb[1] -= 2*wt*Eb[i]*XXXy[Gamma[i]];
      MoMb[0] += wt*Ebb[i*L+i]*ptab[i*L+i]->vD2v;
      MoMb[1] += wt*Ebb[i*L+i]*ptab[i*L+i]->vD4v;
      for (unsigned j = i+1; j < L; ++j) {
        MoMb[0] += 2*wt*Ebb[i*L+j]*ptab[i*L+j]->vD2v;
        MoMb[1] += 2*wt*Ebb[i*L+j]*ptab[i*L+j]->vD4v;
      }
    }
  }
#ifdef TIMEIT
  END = steady_clock::now();
  MOM_T += duration_cast<microseconds>(END - START);
#endif
  // Update PIPs and posterior means and variances of b
#ifdef TIMEIT
  START = steady_clock::now();
#endif
  for (unsigned l = 0; l < L; ++l) {
    rewt = exp(last_logqsum[Gamma[l]]-logqsum);
    PIP[Gamma[l]] *= rewt;
    beta[Gamma[l]] *= rewt;
    se[Gamma[l]] *= rewt;
    PIP[Gamma[l]] += wt;
    beta[Gamma[l]] += Eb[l]*wt;
    se[Gamma[l]] += (Ebb[l*L+l]-Eb[l]*Eb[l])*wt;
    last_logqsum[Gamma[l]] = logqsum;
  }
#ifdef TIMEIT
  END = steady_clock::now();
  PIP_T += duration_cast<microseconds>(END - START);
#endif
  //std::cout << "logq: " << Gamma << " " << logqval << std::endl;
  delete[] ptab;
  return logqval; // Return value of log(q)
}

void FINEMAP::update_variance(bool update_MoM) {
  if (verbose >= 1)
    std::cout << " Computing update for (sigma^2,tau^2)" << std::endl;
  // Finalize b-vector for method-of-moments
#ifdef TIMEIT
  START = steady_clock::now();
#endif
  MoMb[0] = MoMb[0] + ysq;
  MoMb[1] = MoMb[1] + Xysq;
  // Explicitly compute (sigma^2,tau^2) = A^{-1}b
  double det = MoMA[0]*MoMA[3]-MoMA[1]*MoMA[2];
  sigmasq = (MoMA[3]*MoMb[0]-MoMA[1]*MoMb[1])/det;
  tausq = (-MoMA[2]*MoMb[0]+MoMA[0]*MoMb[1])/det;
  if (sigmasq < 0 || tausq < 0) {
    // Update only sigma^2, and set tau^2 = 0
    sigmasq = MoMb[0]/MoMA[0];
    tausq = 0;
  }
#ifdef TIMEIT
  END = steady_clock::now();
  MOM_T += duration_cast<microseconds>(END - START);
#endif
  if (verbose >= 1)
    std::cout << " Update (sigma^2,tau^2) = ("
        << sigmasq << "," << tausq << ")" << std::endl;
  // Recompute xOmegax and xOmegay
  if (verbose >= 1)
    std::cout << " Re-computing probabilities of "
        << models.size() << " explored models" << std::endl;
  for (unsigned i = 0; i < p; ++i)
    Omega[i] = 1.0/(tausq*Dsq[i]+sigmasq);
  for (unsigned i = 0; i < p; ++i) {
    if (!std::isnan(xOmegay[i])) xOmegay[i] = compute_xOmegay(i);
  }
  for (unsigned i = 0; i < pairs.used.size(); ++i) {
    Pair* pair = pairs.keys[pairs.used[i]];
    while (pair != NULL) {
      pair->xOmegax = compute_xOmegax(pair->i,pair->j);
      pair = pair->next;
    }
  }
  // Reset logqsum, MoMb, PIP, beta, se
  // Recompute logq values for existing models, and update these quantities 
  logqsum = NEGINF;
  MoMb[0] = 0;
  MoMb[1] = 0;
  for (unsigned i = 0; i < p; ++i) {
    last_logqsum[i] = NEGINF;
    PIP[i] = 0;
    beta[i] = 0;
    se[i] = 0;
  }
  for (unsigned i = 0; i < models.used.size(); ++i) {
    Model* m = models.keys[models.used[i]];
    while (m != NULL) {
      m->logq = compute_logq(*m, update_MoM);
      m = m->next;
    }
  }
}

void FINEMAP::finemap(const int* sched_sss, unsigned n_conv_sss,
    double prob_tol_sss, unsigned seed) {
  steady_clock::time_point time_begin = steady_clock::now();
  unsigned total_iters = 0;
  // Seed RNG for SSS
  srand(seed);
  // Initialize SSS with empty model
  Model Gamma;
  Gamma.inds = new unsigned[Lmax];
  for (unsigned epoch = 0; ; ++epoch) {
    if (sched_sss[epoch] < 0) break; // Detect sentinel values of -1
    if (verbose >= 1) {
      std::cout << "Starting epoch " << (epoch+1) << std::endl;
      std::cout << " Performing SSS" << std::endl;
    }
    // Update statistics for MoM if not the last epoch
    bool update_MoM = (sched_sss[epoch+1] >= 0);
    // Save total log(q) in each iteration
    std::vector<double> totlogq;
    totlogq.reserve(sched_sss[epoch]);
    for (int it = 0; it < sched_sss[epoch]; ++it) {
      ++total_iters;
#ifdef TIMEIT
      START = steady_clock::now();
#endif
      // Compute all neighboring models
      unsigned idx = 0;
      unsigned pos = 0;
      for (unsigned i = 0; i < p; ++i) {
        if (idx == Gamma.L || i < Gamma[idx]) { // i not in Gamma
          if (Gamma.L < Lmax) {
            // Add i to Gamma
            Model& Gammap = nbrs[pos];
            Gammap.L = Gamma.L+1;
            for (unsigned l = 0; l < idx; ++l) Gammap[l] = Gamma[l];
            Gammap[idx] = i;
            for (unsigned l = idx; l < Gamma.L; ++l) Gammap[l+1] = Gamma[l];
            ++pos;
          }
          for (unsigned l = 0; l < Gamma.L; ++l) {
            // Swap i for Gamma[l] in current model
            Model& Gammap = nbrs[pos];
            Gammap.L = Gamma.L;
            if (l < idx) {
              for (unsigned k = 0; k < l; ++k) Gammap[k] = Gamma[k];
              for (unsigned k = l+1; k < idx; ++k) Gammap[k-1] = Gamma[k];
              Gammap[idx-1] = i;
              for (unsigned k = idx; k < Gamma.L; ++k) Gammap[k] = Gamma[k];
            } else {
              for (unsigned k = 0; k < idx; ++k) Gammap[k] = Gamma[k];
              Gammap[idx] = i;
              for (unsigned k = idx; k < l; ++k) Gammap[k+1] = Gamma[k];
              for (unsigned k = l+1; k < Gamma.L; ++k) Gammap[k] = Gamma[k];
            }
            ++pos;
          }
        } else if (idx < Gamma.L && i == Gamma[idx]) {
          // Delete i from current model
          Model& Gammap = nbrs[pos];
          Gammap.L = Gamma.L-1;
          for (unsigned l = 0; l < idx; ++l) Gammap[l] = Gamma[l];
          for (unsigned l = idx+1; l < Gamma.L; ++l) Gammap[l-1] = Gamma[l];
          ++pos;
          ++idx;
        }
      }
#ifdef TIMEIT
      END = steady_clock::now();
      SSS_T += duration_cast<microseconds>(END - START);
#endif
      // Sample new model from neighbors
      double logqsumnbrs = NEGINF;
      for (unsigned i = 0; i < pos; ++i) {
#ifdef TIMEIT
        START = steady_clock::now();
#endif
        Model* m = models.find(nbrs[i]);
#ifdef TIMEIT
        END = steady_clock::now();
        MODEL_T += duration_cast<microseconds>(END - START);
#endif
        if (m != NULL) nbrs[i].logq = m->logq;
        else {
          nbrs[i].logq = compute_logq(nbrs[i], update_MoM);
#ifdef TIMEIT
          START = steady_clock::now();
#endif
          models.insert(nbrs[i],nbrs[i].logq);
#ifdef TIMEIT
          END = steady_clock::now();
          MODEL_T += duration_cast<microseconds>(END - START);
#endif
        }
        logqsumnbrs = logaddexp(logqsumnbrs,nbrs[i].logq);
      }
#ifdef TIMEIT
      START = steady_clock::now();
#endif
      double runif = ((double)rand() / RAND_MAX);
      double nbr_tot = 0;
      for (unsigned i = 0; i < pos; ++i) {
        nbr_tot += exp(nbrs[i].logq-logqsumnbrs);
        if (nbr_tot > runif) {
          // Set Gamma to nbrs[i]
          Gamma.L = nbrs[i].L;
          for (unsigned l = 0; l < Gamma.L; ++l) Gamma[l] = nbrs[i][l];
          Gamma.logq = nbrs[i].logq;
          break;
        }
      }
#ifdef TIMEIT
      END = steady_clock::now();
      SSS_T += duration_cast<microseconds>(END - START);
#endif
      // Save total model log(q), and assess convergence
      totlogq.push_back(logqsum);
      if (totlogq.size() > n_conv_sss) {
        double added = 1-exp(totlogq[totlogq.size()-n_conv_sss-1]
                -totlogq[totlogq.size()-1]);
        if (verbose >= 2) {
          std::cout << "  Iteration " << (it+1) << ", model = " << Gamma;
          std::cout << ", probability = "
              << exp(Gamma.logq-logqsum)
              << ", added probability in last " << n_conv_sss
              << " iterations: " << added << std::endl;
        }
        if (added < prob_tol_sss) {
          if (verbose >= 1) std::cout << " Epoch " << (epoch+1)
              << " converged in " << (it+1) << " SSS iterations"
              << std::endl;
          break;
        }
      } else if (verbose >= 2) {
        std::cout << "  Iteration " << (it+1) << ", model = " << Gamma;
        std::cout << ", probability = "
            << exp(Gamma.logq-logqsum) << std::endl;
      }
      //std::cout << "  models size " << models.size() << std::endl;
      //std::cout << "  pairs size " << pairs.size() << std::endl;
      //std::cout << "models table:" << std::endl << models << std::endl;
      //std::cout << "pairs table:" << std::endl << pairs << std::endl;
    }
    // If not the last epoch, update (sigma^2,tau^2)
    if (sched_sss[epoch+1] >= 0) {
      // If not the second-to-last epoch, also update statistics for the
      // next update of (sigma^2,tau^2)
      update_MoM = (sched_sss[epoch+2] >= 0);
      update_variance(update_MoM);
    }
  }
  delete[] Gamma.inds;
#ifdef TIMEIT
  START = steady_clock::now();
#endif
  // Compute posterior mean of alpha
  for (unsigned i = 0; i < p; ++i)
    for (unsigned j = 0; j < p; ++j)
      VXy[i] -= Dsq[i]*V[j*p+i]*beta[j]; // Residualize D^2V'beta from V'X'y
  for (unsigned j = 0; j < p; ++j) {
    alpha[j] = 0;
    for (unsigned k = 0; k < p; ++k)
      alpha[j] += tausq*V[j*p+k]*VXy[k]*Omega[k];
  }
  // Condition posterior mean and std on inclusion
  for (unsigned i = 0; i < p; ++i) {
    beta[i] = beta[i] / PIP[i];
    se[i] = sqrt(se[i] / PIP[i]);
  }
#ifdef TIMEIT
  END = steady_clock::now();
  PIP_T += duration_cast<microseconds>(END - START);
#endif
  steady_clock::time_point time_end = steady_clock::now();
  microseconds time_elapsed =
      duration_cast<microseconds>(time_end-time_begin);
  if (verbose >= 1) {
    std::cout << "DONE" << std::endl;
    std::cout << "Elapsed time " << double(time_elapsed.count())/1e6 << " seconds" << std::endl;
    std::cout << models.size() << " total explored models" << std::endl;
    std::cout << pairs.size() << " total explored pairs" << std::endl;
    std::cout << total_iters << " total iterations" << std::endl;
#ifdef TIMEIT
    std::cout << "SSS_T " << double(SSS_T.count())/1e6 << std::endl;
    std::cout << "MODEL_T " << double(MODEL_T.count())/1e6 << std::endl;
    std::cout << "PAIR_T " << double(PAIR_T.count())/1e6 << std::endl;
    std::cout << "XOMEGA_T " << double(XOMEGA_T.count())/1e6 << std::endl;
    std::cout << "VDV_T " << double(VDV_T.count())/1e6 << std::endl;
    std::cout << "LOGQ_T " << double(LOGQ_T.count())/1e6 << std::endl;
    std::cout << "MOM_T " << double(MOM_T.count())/1e6 << std::endl;
    std::cout << "PIP_T " << double(PIP_T.count())/1e6 << std::endl;
#endif
  }
}

void FINEMAP::getmodels(unsigned nmodels, ModelQueue& q) {
  for (unsigned i = 0; i < models.used.size(); ++i) {
    Model* m = models.keys[models.used[i]];
    while (m != NULL) {
      if (q.size() < nmodels)
        q.push(m);
      else if (q.top()->logq < m->logq) {
        q.pop();
        q.push(m);
      }
      m = m->next;
    }
  }
}
