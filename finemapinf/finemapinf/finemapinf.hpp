#ifndef __finemap_hpp__
#define __finemap_hpp__

#include "hashtables.hpp"
#include <queue>

namespace FINEMAPINF {
 
  // FINEMAP object contains the input data, maintains the output arrays
  // and hash tables, and performs SSS and method-of-moments updates
  class FINEMAP {
    public:
      // All data inputs and outputs are passed to the constructor
      FINEMAP(unsigned _n, unsigned _p, double _sigmasq, double _tausq,
              const double* _pi0, const double* _S, unsigned _nssq,
              unsigned _Lmax,
              double _meansq, const double* _z, const double* _V,
              const double* _Dsq, int _verbose,
              unsigned _models_nbins, unsigned _pairs_nbins,
              double* _PIP, double* _beta, double* _se, double* _alpha);

      ~FINEMAP();

      // Compute log(q) (unnormalized posterior probability) for model Gamma
      // Also updates PIP, beta, and se using this model
      // Also updates MoM b-vector using this model, if update_MoM = true
      double compute_logq(const Model& Gamma, bool update_MoM);

      // Compute x_i'Omega x_j, x_i'Omega y, v_i'D^2 v_j, v_i'D^4 v_j
      double compute_xOmegax(unsigned i, unsigned j);
      double compute_xOmegay(unsigned j);
      double compute_vD2v(unsigned i, unsigned j);
      double compute_vD4v(unsigned i, unsigned j);

      // Perform method-of-moments update for (sigma^2,tau^2)
      // Recomputes log(q) and PIP, beta, se values for existing models
      // Also recomputes MoM b-vector for next update, if update_MoM = true
      void update_variance(bool update_MoM);

      unsigned n;     // Number samples
      unsigned p;     // Number SNPs
      unsigned Lmax;  // Max number of causal variants
      unsigned nssq;  // Number of prior effect-size variance values
      double sigmasq; // sigma^2 variance parameter
      double tausq;   // tau^2 variance parameter
      double* logodds;  // log(pi0/(1-pi0))
      double logpi0c; // sum of log(1-pi0) across all SNPs
      double logqsum; // log-sum of q for all explored models
      double ysq;     // ||y||_2^2
      double Xysq;    // ||X'y||_2^2
      int verbose;    

      double MoMA[4]; // 2x2 matrix A; method-of-moments solves A^{-1}b
      double MoMb[2]; // 2x1 vector b; method-of-moments solves A^{-1}b

      /* Memory externally managed */
      const double* V;    // From eigendecomposition VD^2V' = n*LD
      const double* Dsq;  // From eigendecomposition VD^2V' = n*LD
      const double* S;    // Prior effect-size variances
      double* PIP;        // Output PIPs
      double* beta;       // Output posterior means
      double* se;         // Output posterior stds
      double* alpha;      // Output posterior mean of alpha

      /* Memory internally managed */
      double* Omega;      // diag(1/(tau^2 D^2+sigma^2 I))
      double* Xy;         // X'y = sqrt(n) * z-scores
      double* VXy;        // V'X'y
      double* XXXy;       // (X'X)(X'y)
      double* logpissq;   // Temporary storage
      double* mus;        // Temporary storage
      double* omegainvs;  // Temporary storage
      double* Eb;         // Temporary storage
      double* Ebb;        // Temporary storage
      Model* nbrs;        // Temporary storage
      double* last_logqsum; // Temporary storage
     
      // Hash tables
      double* xOmegay;    // Values of x_i'Omega y
      PairHash pairs;     // Values of x_i'Omega x_j, v_i'D^2v_j, v_i'D^4v_j
      ModelHash models;   // Values of log(q)

      // Main finemapping function
      // Parameters of SSS and random seed for SSS are input here
      void finemap(const int* sched_sss, unsigned n_conv_sss,
              double prob_tol_sss, unsigned seed);

      // Filter to list of most probable causal models
      typedef std::priority_queue<Model*,std::vector<Model*>,ModelCompare>
          ModelQueue;
      void getmodels(unsigned nmodels, ModelQueue& q);
  };
}

#endif
