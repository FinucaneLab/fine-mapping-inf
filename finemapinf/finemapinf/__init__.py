import numpy as np
import scipy.linalg
import _c_funcs

def finemap(z,meansq,n,L,LD=None,V=None,Dsq=None,
    sigmasq=1,tausq=0,pi0=None,S=[0.05**2],
    sched_sss=[100,100,10000],n_conv_sss=100,nmodels=1000,prob_tol_sss=1e-3,
    seed=21,models_nbins=None,pairs_nbins=None,verbose=1):
  '''FINEMAP with infinitesimal random effects

  z -- vector of z-scores (equal to X'y/sqrt(n))
  meansq -- average squared magnitude of y (equal to y'y/n)
  n -- sample size
  L -- maximum number of causal effects

  LD -- LD matrix (equal to X'X/n)
  V -- precomputed p x p matrix of eigenvectors of X'X
  Dsq -- precomputed length-p vector of eigenvalues of X'X
  (Must provide either LD or the pair (V,Dsq))

  sigmasq -- initial value for sigma^2
  tausq -- initial value for tau^2
  pi0 -- length-p vector of prior causal probabilities
       Default: 1/p for each SNP
  S -- grid of effect size variances s^2; prior for s^2 is uniform on this grid
       Default: single value s^2 = 0.05^2

  sched_sss -- maximum number of iterations of SSS to perform in each epoch
       The parameters (sigma^2,tau^2) are updated between consecutive epochs
       Default: 3 epochs of 100, 100, 10000 iterations, with two updates of
                (sigma^2,tau^2) in between
  n_conv_sss, prob_tol_sss -- convergence criteria within each epoch
       If the total probability of newly explored models in the last
       n_conv_sss iterations of SSS is < prob_tol_sss, terminate this epoch
  nmodels -- number of most probable causal models to return

  seed -- random seed for SSS
  models_nbins -- Number of bins in hash table of explored models
       Default: min(p*1000,1e8)
  pairs_nbins -- Number of bins in hash table of explored pairs (i,j)
       Default: min(p*100,1e7)

  verbose -- 0 for no output, 1 for basic output, 2 for per-iteration output

  Returns: Dictionary with keys
    PIP -- length-p vector of PIPs
    beta -- length-p vector of posterior mean of beta conditional on causal
            i.e. E[beta_j | X, y, gamma_j = 1]
    se -- length-p vector of posterior std.dev. of beta conditional on causal
            i.e. sqrt(Var[beta_j | X, y, gamma_j = 1])
    sigmasq -- final value of sigma^2
    tausq -- final value of tau^2
    alpha -- length-p vector of posterior mean of infinitesimal effects
    models -- list of most probable causal models and log-posterior values
  '''
  # Compute eigendecomposition of LD
  p = len(z)
  if (V is None or Dsq is None) and LD is None:
    raise RuntimeError('Missing LD')
  elif V is None or Dsq is None:
    if LD.shape != (p,p):
      raise RuntimeError('Incompatible z-score and LD dimensions')
    eigvals,V = scipy.linalg.eigh(LD)
    Dsq = np.maximum(n*eigvals,0)
  else:
    if V.shape != (p,p) or Dsq.shape != (p,):
      raise RuntimeError('Incompatible z-score and V,Dsq dimensions')
    Dsq = np.maximum(Dsq,0)
  # Initialize prior parameters
  if pi0 is None:
    pi0 = np.ones(p,dtype=np.double,order='C')*1.0/p
  else:
    pi0 = np.array(pi0,dtype=np.double,order='C')
  S = np.array(S,dtype=np.double,order='C')
  # Initialize SSS parameters
  sched_sss = list(sched_sss) + [-1] # Add sentinel
  sched_sss = np.array(sched_sss,dtype=np.intc,order='C')
  if models_nbins is None: models_nbins = min(p*1000, 1e8)
  if pairs_nbins is None: pairs_nbins = min(p*100,1e7)
  # Ensure that input numpy input arrays are double-precision and C-contiguous
  z = z.astype(np.double,'C',copy=False)
  V = V.astype(np.double,'C',copy=False)
  Dsq = Dsq.astype(np.double,'C',copy=False)
  # Initialize output arrays
  PIP = np.zeros(p,dtype=np.double,order='C')
  beta = np.zeros(p,dtype=np.double,order='C')
  se = np.zeros(p,dtype=np.double,order='C')
  alpha = np.zeros(p,dtype=np.double,order='C')
  # Run FINEMAP from C++
  (sigmasq,tausq,models) = _c_funcs.finemap(n,meansq,z,V,Dsq,L,sigmasq,tausq,
      pi0,S,sched_sss,n_conv_sss,prob_tol_sss,nmodels,seed,verbose,models_nbins,
      pairs_nbins,PIP,beta,se,alpha)
  return {'PIP':PIP, 'beta':beta, 'se':se, 'sigmasq':sigmasq, 'tausq':tausq,
    'alpha': alpha, 'models': models}

