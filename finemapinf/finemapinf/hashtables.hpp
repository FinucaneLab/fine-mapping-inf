#ifndef __hashtables_hpp__
#define __hashtables_hpp__

#include <vector>
#include <limits>
#include <cstddef>

namespace FINEMAPINF {

  const double NEGINF = -std::numeric_limits<double>::max();
  const double NaN = std::numeric_limits<double>::quiet_NaN();

  // Model object, containing tuple of indices, log(q) value, and pointer
  // to implement linked list
  //
  // (Memory for inds is managed explicitly, as default memory management
  // using std::vector is too high)
  struct Model {

    // Memory for inds should be allocated and freed externally
    Model() : L(0), inds(NULL), logq(NaN), next(NULL) { }

    unsigned& operator[](unsigned l) { return inds[l]; }
    const unsigned& operator[](unsigned l) const { return inds[l]; }

    // Compare two models based on tuple of indices
    bool operator==(const Model& m) {
      if (L != m.L) return false;
      for (unsigned l = 0; l < m.L; ++l) {
        if (inds[l] != m.inds[l]) return false;
      }
      return true;
    }

    unsigned L;      // number of indices
    unsigned* inds;  // array of indices; always assumed sorted
    double logq;     // log(q) value; NaN if not yet computed
    Model* next;     // next Model, if used in linked list
  };

  struct ModelCompare {
    bool operator() (const Model* lhs, const Model* rhs) const {
      return lhs->logq > rhs->logq;
    }
  };

  // Explicit implementation of hash table for visited models
  // (Memory overhead using std::unordered_map is too high)
  struct ModelHash {
    ModelHash(unsigned nbins);
    ~ModelHash();

    unsigned hash(const Model& m);
    unsigned size();

    // return pointer to matching Model in the table, or NULL if not found
    Model* find(const Model& key);

    // add new Model with log(q) = val; assumes not yet in table
    Model* insert(const Model& key, double val);
    
    Model** keys;                // array of _nbins linked lists
    std::vector<unsigned> used;  // which hash keys have been used
    unsigned _nbins;             // number of buckets of hash table
    unsigned _size;              // number of models
  };

  // Object representing pair of indices (i,j) with values x_i' Omega x_j,
  // v_i' D^2 v_j, and v_i' D^4 v_j, and pointer to implement linked list
  struct Pair {
    Pair() : i(0), j(0), xOmegax(NaN), vD2v(NaN), vD4v(NaN), next(NULL) { }
    bool operator==(const Pair& p) { return (i==p.i&&j==p.j); }
    
    unsigned i;      // first index
    unsigned j;      // second index
    double xOmegax;  // x_i' Omega x_j value; NaN if not yet computed
    double vD2v;     // v_i' D^2 v_j value; NaN if not yet computed
    double vD4v;     // v_i' D^4 v_j value; NaN if not yet computed
    Pair* next;      // next Pair, if used in linked list
  };

  // Explicit implementation of hash table for visited pairs
  // (Memory overhead using std::unordered_map is too high)
  class PairHash {
    public:
      PairHash(unsigned nbins);
      ~PairHash();

      unsigned hash(unsigned i, unsigned j);
      unsigned size();

      // return pointer to matching Pair in the table, or NULL if not found
      Pair* find(unsigned i, unsigned j);

      // add new Pair with (xOmegax,vD2v,vD4v)=(v1,v2,v3);
      // assumes not yet in table
      Pair* insert(unsigned i, unsigned j, double v1, double v2, double v3);

      Pair** keys;                 // array of _nbins linked lists
      std::vector<unsigned> used;  // which hash keys have been used
      unsigned _nbins;             // number of buckets of hash table
      unsigned _size;              // number of pairs
  };
}

#endif
