#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#include <Python.h>
#include <numpy/arrayobject.h>
#include <math.h>
#include "finemapinf.hpp"

PyObject* py_finemap(PyObject* self, PyObject* args) {
  // Inputs
  int n;
  double meansq;
  PyArrayObject* py_z = NULL;
  PyArrayObject* py_V = NULL;
  PyArrayObject* py_Dsq = NULL;
  int Lmax;
  double sigmasq;
  double tausq;
  PyArrayObject* py_pi0 = NULL;
  PyArrayObject* py_S = NULL;
  PyArrayObject* py_sched_sss = NULL;
  int n_conv_sss;
  double prob_tol_sss;
  int nmodels;
  int seed;
  int models_nbins;
  int pairs_nbins;
  int verbose;
  // Outputs
  PyArrayObject* py_PIP = NULL;
  PyArrayObject* py_beta = NULL;
  PyArrayObject* py_se = NULL;
  PyArrayObject* py_alpha = NULL;
  if (!PyArg_ParseTuple(args, "idOOOiddOOOidiiiiiOOOO",
              &n, &meansq, &py_z, &py_V, &py_Dsq, &Lmax, &sigmasq, &tausq,
              &py_pi0, &py_S, &py_sched_sss, &n_conv_sss, &prob_tol_sss,
              &nmodels, &seed, &verbose, &models_nbins, &pairs_nbins, 
              &py_PIP, &py_beta, &py_se, &py_alpha))
    return NULL;
  // Dimensions
  int p = (int)PyArray_DIM(py_z,0);
  int nssq = (int)PyArray_DIM(py_S,0);
  // Cast to C arrays
  double* z = static_cast<double*>(PyArray_DATA(py_z));
  double* V = static_cast<double*>(PyArray_DATA(py_V));
  double* Dsq = static_cast<double*>(PyArray_DATA(py_Dsq));
  double* pi0 = static_cast<double*>(PyArray_DATA(py_pi0));
  double* S = static_cast<double*>(PyArray_DATA(py_S));
  int* sched_sss = static_cast<int*>(PyArray_DATA(py_sched_sss));
  double* PIP = static_cast<double*>(PyArray_DATA(py_PIP));
  double* beta = static_cast<double*>(PyArray_DATA(py_beta));
  double* se = static_cast<double*>(PyArray_DATA(py_se));
  double* alpha = static_cast<double*>(PyArray_DATA(py_alpha));
  // Run finemap
  FINEMAPINF::FINEMAP FM(n,p,sigmasq,tausq,pi0,S,nssq,(unsigned)Lmax,meansq,
          z,V,Dsq,verbose,models_nbins,pairs_nbins,PIP,beta,se,alpha);
  FM.finemap(sched_sss,n_conv_sss,prob_tol_sss,seed);
  FINEMAPINF::FINEMAP::ModelQueue q;
  FM.getmodels(nmodels,q);
  PyObject* list = PyList_New(0);
  while (q.size() > 0) {
    FINEMAPINF::Model* m = q.top();
    q.pop();
    PyObject* model = PyList_New(m->L);
    for (unsigned i = 0; i < m->L; ++i)
      PyList_SET_ITEM(model, i, PyLong_FromLong((long)m->inds[i]));
    double prob = exp(m->logq - FM.logqsum);
    PyObject* tup = PyTuple_Pack(2,PyFloat_FromDouble(prob),model);
    PyList_Append(list, tup);
  }
  PyList_Reverse(list);
  // Return sigma^2, tau^2, and list of models; remaining outputs are
  // saved in PIP, beta, se
  return Py_BuildValue("ddO",FM.sigmasq,FM.tausq,list);
}

static PyMethodDef methods[] = {
  {"finemap", py_finemap, METH_VARARGS, "Run FINEMAP-inf"},
  {NULL, NULL, 0, NULL}
};

static struct PyModuleDef cMod = {
  PyModuleDef_HEAD_INIT,
  "_c_funcs",
  "",
  -1,
  methods
};

PyMODINIT_FUNC PyInit__c_funcs(void) {
  import_array();
  return PyModule_Create(&cMod);
}
