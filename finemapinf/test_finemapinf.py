import numpy as np
from finemapinf import finemap

verbose = 1

np.random.seed(21)
n = 5000
p = 500
MAF = 0.1
Ltrue = 5
ssq = 0.01
sigmasq = 1
X = np.random.binomial(2,MAF,size=(n,p))
X = X - X.mean(axis=0)
X = X / X.std(axis=0)
LD = (X.T.dot(X))/n
b = np.zeros(p)
inds = np.random.choice(p,size=Ltrue,replace=False)
b[inds] = np.random.normal(size=Ltrue) * np.sqrt(ssq)
order = np.argsort(inds)
print('True effects: %s' % str(inds[order]))
print('Effect sizes: %s\n' % str(b[inds[order]]))

# Generate data without infinitesimal effects
print('###### Simulating without infinitesimal effects ######')
effects = X.dot(b)
y = effects + np.random.normal(size=n) * np.sqrt(sigmasq)
print('Total fraction of variance explained by SNPs: %f\n' % (np.var(effects)/np.var(y)))
z = (X.T.dot(y))/np.sqrt(n)
meansq = y.T.dot(y)/n

print('FINEMAP-inf with fixed sigma^2=<y^2> and tau^2=0')
L = 5
output = finemap(z,meansq,n,L,LD=LD,verbose=verbose,
    sigmasq=np.mean(y**2),tausq=0,sched_sss=[10000],nmodels=10)
PIP = output['PIP']
beta = output['beta']
print('SNPs with PIP > 0.1: ' + str(np.nonzero(PIP > 0.1)[0]))
print('Corresponding PIPs: %s' % str(np.round(PIP[np.nonzero(PIP > 0.1)],2)))
print('Corresponding effect sizes: %s' % str(np.round(beta[np.nonzero(PIP > 0.1)],5)))
print('Estimated (sigma^2,tau^2): (%f,%f)\n' % (output['sigmasq'],output['tausq']))

print('FINEMAP-inf with estimated sigma^2 and tau^2')
L = 5
output = finemap(z,meansq,n,L,LD=LD,verbose=verbose,nmodels=10)
PIP = output['PIP']
beta = output['beta']
print('SNPs with PIP > 0.1: ' + str(np.nonzero(PIP > 0.1)[0]))
print('Corresponding PIPs: %s' % str(np.round(PIP[np.nonzero(PIP > 0.1)],2)))
print('Corresponding effect sizes: %s' % str(np.round(beta[np.nonzero(PIP > 0.1)],5)))
print('Estimated (sigma^2,tau^2): (%f,%f)\n' % (output['sigmasq'],output['tausq']))

# Generate data with strong infinitesimal effects
print('###### Simulating with strong infinitesimal effects, tau^2 = 1e-3 ######')
tausq = 1e-3
effects = X.dot(b) + X.dot(np.random.normal(size=p) * np.sqrt(tausq))
y = effects + np.random.normal(size=n) * np.sqrt(sigmasq)
print('Total fraction of variance explained by SNPs: %f\n' % (np.var(effects)/np.var(y)))
LD = (X.T.dot(X))/n
z = (X.T.dot(y))/np.sqrt(n)
meansq = y.T.dot(y)/n

print('FINEMAP-inf with fixed sigma^2=<y^2> and tau^2=0')
L = 5
output = finemap(z,meansq,n,L,LD=LD,verbose=verbose,
    sigmasq=np.mean(y**2),tausq=0,sched_sss=[10000],nmodels=10)
PIP = output['PIP']
beta = output['beta']
print('SNPs with PIP > 0.1: ' + str(np.nonzero(PIP > 0.1)[0]))
print('Corresponding PIPs: %s' % str(np.round(PIP[np.nonzero(PIP > 0.1)],2)))
print('Corresponding effect sizes: %s' % str(np.round(beta[np.nonzero(PIP > 0.1)],5)))
print('Estimated (sigma^2,tau^2): (%f,%f)\n' % (output['sigmasq'],output['tausq']))

print('FINEMAP-inf with estimated sigma^2 and tau^2')
L = 5
output = finemap(z,meansq,n,L,LD=LD,verbose=verbose,nmodels=10)
PIP = output['PIP']
beta = output['beta']
print('SNPs with PIP > 0.1: ' + str(np.nonzero(PIP > 0.1)[0]))
print('Corresponding PIPs: %s' % str(np.round(PIP[np.nonzero(PIP > 0.1)],2)))
print('Corresponding effect sizes: %s' % str(np.round(beta[np.nonzero(PIP > 0.1)],5)))
print('Estimated (sigma^2,tau^2): (%f,%f)\n' % (output['sigmasq'],output['tausq']))

# Generate data with weak infinitesimal effects
print('###### Simulating with weak infinitesimal effects, tau^2 = 1e-4 ######')
tausq = 1e-4
effects = X.dot(b) + X.dot(np.random.normal(size=p) * np.sqrt(tausq))
y = effects + np.random.normal(size=n) * np.sqrt(sigmasq)
print('Total fraction of variance explained by SNPs: %f\n' % (np.var(effects)/np.var(y)))
z = (X.T.dot(y))/np.sqrt(n)
meansq = y.T.dot(y)/n

print('FINEMAP-inf with fixed sigma^2=<y^2> and tau^2=0')
L = 5
output = finemap(z,meansq,n,L,LD=LD,verbose=verbose,
    sigmasq=np.mean(y**2),tausq=0,sched_sss=[10000],nmodels=10)
PIP = output['PIP']
beta = output['beta']
print('SNPs with PIP > 0.1: ' + str(np.nonzero(PIP > 0.1)[0]))
print('Corresponding PIPs: %s' % str(np.round(PIP[np.nonzero(PIP > 0.1)],2)))
print('Corresponding effect sizes: %s' % str(np.round(beta[np.nonzero(PIP > 0.1)],5)))
print('Estimated (sigma^2,tau^2): (%f,%f)\n' % (output['sigmasq'],output['tausq']))

print('FINEMAP-inf with estimated sigma^2 and tau^2')
L = 5
output = finemap(z,meansq,n,L,LD=LD,verbose=verbose,nmodels=10)
PIP = output['PIP']
beta = output['beta']
print('SNPs with PIP > 0.1: ' + str(np.nonzero(PIP > 0.1)[0]))
print('Corresponding PIPs: %s' % str(np.round(PIP[np.nonzero(PIP > 0.1)],2)))
print('Corresponding effect sizes: %s' % str(np.round(beta[np.nonzero(PIP > 0.1)],5)))
print('Estimated (sigma^2,tau^2): (%f,%f)\n' % (output['sigmasq'],output['tausq']))

