#include "hashtables.hpp"
#include <iostream>

using namespace FINEMAPINF;

ModelHash::ModelHash(unsigned nbins) : _nbins(nbins), _size(0) {
  keys = new Model*[_nbins];
  for (unsigned i = 0; i < _nbins; ++i) keys[i] = NULL;
  used.reserve(_nbins); // Reserve memory for maximum possible size
}

ModelHash::~ModelHash() {
  // Delete all models and their indices in the table
  for (unsigned i = 0; i < used.size(); ++i) {
    Model* m = keys[used[i]];
    while (m != NULL) {
      Model* next = m->next;
      delete[] m->inds;
      delete m;
      m = next;
    }
  }
  delete[] keys;
}

unsigned ModelHash::hash(const Model& m) {
  // Hash function from array of unsigned to unsigned in [0,_nbins-1]
  unsigned tot = m.L;
  for (unsigned i = 0; i < m.L; ++i) {
    tot ^= m[i] + 0x9e3779b9 + (tot << 6) + (tot >> 2);
  }
  return tot % _nbins;
}

unsigned ModelHash::size() {
  return _size;
}

Model* ModelHash::find(const Model& key) {
  unsigned h = hash(key);
  Model* m = keys[h];
  while (m != NULL) {
    if (*m == key) return m;
    m = m->next;
  }
  return NULL;
}

Model* ModelHash::insert(const Model& key, double val) {
  unsigned h = hash(key);
  // Create copy of Model
  Model* m = new Model();
  m->L = key.L;
  if (m->L > 0) { m->inds = new unsigned[m->L]; }
  for (unsigned l = 0; l < key.L; ++l) (*m)[l] = key[l];
  m->logq = val;
  // Add to the start of linked list for h
  m->next = keys[h];
  if (m->next == NULL) used.push_back(h);
  keys[h] = m;
  ++_size;
  return m;
}

PairHash::PairHash(unsigned nbins) : _nbins(nbins), _size(0) {
  keys = new Pair*[_nbins];
  for (unsigned i = 0; i < _nbins; ++i) keys[i] = NULL;
  used.reserve(_nbins); // Reserve memory for maximum possible size
}

PairHash::~PairHash() {
  // Delete all pairs in the table
  for (unsigned i = 0; i < used.size(); ++i) {
    Pair* p = keys[used[i]];
    while (p != NULL) {
      Pair* next = p->next;
      delete p;
      p = next;
    }
  }
  delete[] keys;
}

unsigned PairHash::hash(unsigned i, unsigned j) {
  // Hash function from pair of unsigned to unsigned in [0,_nbins-1]
  unsigned tot = 2;
  tot ^= i + 0x9e3779b9 + (tot << 6) + (tot >> 2);
  tot ^= j + 0x9e3779b9 + (tot << 6) + (tot >> 2);
  return tot % _nbins;
}

unsigned PairHash::size() {
  return _size;
}

Pair* PairHash::find(unsigned i, unsigned j) {
  unsigned h = hash(i,j);
  Pair* p = keys[h];
  while (p != NULL) {
    if (p->i == i && p->j == j) return p;
    p = p->next;
  }
  return NULL;
}

Pair* PairHash::insert(unsigned i, unsigned j, double v1, double v2,
    double v3) {
  unsigned h = hash(i,j);
  // Create copy of Pair
  Pair* p = new Pair();
  p->i = i;
  p->j = j;
  p->xOmegax = v1;
  p->vD2v = v2;
  p->vD4v = v3;
  // Add to the start of linked list for h
  p->next = keys[h];
  if (p->next == NULL) used.push_back(h);
  keys[h] = p;
  ++_size;
  return p;
}
