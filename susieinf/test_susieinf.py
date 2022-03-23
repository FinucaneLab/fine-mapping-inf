import numpy as np
import susieinf

verbose = False

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

print('SuSiE with MLE-estimated sigma^2 and tau^2=0')
L = 5
output = susieinf.susie(z,meansq,n,L,LD=LD,est_tausq=False,tausq=0,method='MLE',verbose=verbose)
PIP = 1-(1-output['PIP']).prod(axis=1)
cred = susieinf.cred(output['PIP'],LD=LD)
print('SNPs with PIP > 0.1: ' + str(np.nonzero(PIP > 0.1)[0]))
print('Corresponding PIPs: %s' % str(np.round(PIP[np.nonzero(PIP > 0.1)],2)))
print('Credible sets: ' + str(cred))
print('Estimated sigma^2: %f\n' % output['sigmasq'])

print('SuSiE with MLE-estimated sigma^2 and tau^2')
L = 5
output = susieinf.susie(z,meansq,n,L,LD=LD,method='MLE',verbose=verbose)
PIP = 1-(1-output['PIP']).prod(axis=1)
cred = susieinf.cred(output['PIP'],LD=LD)
print('SNPs with PIP > 0.1: ' + str(np.nonzero(PIP > 0.1)[0]))
print('Corresponding PIPs: %s' % str(np.round(PIP[np.nonzero(PIP > 0.1)],2)))
print('Credible sets: ' + str(cred))
print('Estimated (sigma^2,tau^2): (%f,%f)\n' % (output['sigmasq'],output['tausq']))

print('SuSiE with moments-estimated sigma^2 and tau^2')
L = 5
output = susieinf.susie(z,meansq,n,L,LD=LD,method='moments',verbose=verbose)
PIP = 1-(1-output['PIP']).prod(axis=1)
cred = susieinf.cred(output['PIP'],LD=LD)
print('SNPs with PIP > 0.1: ' + str(np.nonzero(PIP > 0.1)[0]))
print('Corresponding PIPs: %s' % str(np.round(PIP[np.nonzero(PIP > 0.1)],2)))
print('Credible sets: ' + str(cred))
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

print('SuSiE with MLE-estimated sigma^2 and tau^2=0')
L = 5
output = susieinf.susie(z,meansq,n,L,LD=LD,est_tausq=False,tausq=0,method='MLE',verbose=verbose)
PIP = 1-(1-output['PIP']).prod(axis=1)
cred = susieinf.cred(output['PIP'],LD=LD)
print('SNPs with PIP > 0.1: ' + str(np.nonzero(PIP > 0.1)[0]))
print('Corresponding PIPs: %s' % str(np.round(PIP[np.nonzero(PIP > 0.1)],2)))
print('Credible sets: ' + str(cred))
print('Estimated sigma^2: %f\n' % output['sigmasq'])

print('SuSiE with MLE-estimated sigma^2 and tau^2')
L = 5
output = susieinf.susie(z,meansq,n,L,LD=LD,method='MLE',verbose=verbose)
PIP = 1-(1-output['PIP']).prod(axis=1)
cred = susieinf.cred(output['PIP'],LD=LD)
print('SNPs with PIP > 0.1: ' + str(np.nonzero(PIP > 0.1)[0]))
print('Corresponding PIPs: %s' % str(np.round(PIP[np.nonzero(PIP > 0.1)],2)))
print('Credible sets: ' + str(cred))
print('Estimated (sigma^2,tau^2): (%f,%f)\n' % (output['sigmasq'],output['tausq']))

print('SuSiE with moments-estimated sigma^2 and tau^2')
L = 5
output = susieinf.susie(z,meansq,n,L,LD=LD,method='moments',verbose=verbose)
PIP = 1-(1-output['PIP']).prod(axis=1)
cred = susieinf.cred(output['PIP'],LD=LD)
print('SNPs with PIP > 0.1: ' + str(np.nonzero(PIP > 0.1)[0]))
print('Corresponding PIPs: %s' % str(np.round(PIP[np.nonzero(PIP > 0.1)],2)))
print('Credible sets: ' + str(cred))
print('Estimated (sigma^2,tau^2): (%f,%f)\n' % (output['sigmasq'],output['tausq']))

# Generate data with weak infinitesimal effects
print('###### Simulating with weak infinitesimal effects, tau^2 = 1e-4 ######')
tausq = 1e-4
effects = X.dot(b) + X.dot(np.random.normal(size=p) * np.sqrt(tausq))
y = effects + np.random.normal(size=n) * np.sqrt(sigmasq)
print('Total fraction of variance explained by SNPs: %f\n' % (np.var(effects)/np.var(y)))
z = (X.T.dot(y))/np.sqrt(n)
meansq = y.T.dot(y)/n

print('SuSiE with MLE-estimated sigma^2 and tau^2=0')
L = 5
output = susieinf.susie(z,meansq,n,L,LD=LD,est_tausq=False,tausq=0,method='MLE',verbose=verbose)
PIP = 1-(1-output['PIP']).prod(axis=1)
cred = susieinf.cred(output['PIP'],LD=LD)
print('SNPs with PIP > 0.1: ' + str(np.nonzero(PIP > 0.1)[0]))
print('Corresponding PIPs: %s' % str(np.round(PIP[np.nonzero(PIP > 0.1)],2)))
print('Credible sets: ' + str(cred))
print('Estimated sigma^2: %f\n' % output['sigmasq'])

print('SuSiE with MLE-estimated sigma^2 and tau^2')
L = 5
output = susieinf.susie(z,meansq,n,L,LD=LD,method='MLE',verbose=verbose)
PIP = 1-(1-output['PIP']).prod(axis=1)
cred = susieinf.cred(output['PIP'],LD=LD)
print('SNPs with PIP > 0.1: ' + str(np.nonzero(PIP > 0.1)[0]))
print('Corresponding PIPs: %s' % str(np.round(PIP[np.nonzero(PIP > 0.1)],2)))
print('Credible sets: ' + str(cred))
print('Estimated (sigma^2,tau^2): (%f,%f)\n' % (output['sigmasq'],output['tausq']))

print('SuSiE with moments-estimated sigma^2 and tau^2')
L = 5
output = susieinf.susie(z,meansq,n,L,LD=LD,method='moments',verbose=verbose)
PIP = 1-(1-output['PIP']).prod(axis=1)
cred = susieinf.cred(output['PIP'],LD=LD)
print('SNPs with PIP > 0.1: ' + str(np.nonzero(PIP > 0.1)[0]))
print('Corresponding PIPs: %s' % str(np.round(PIP[np.nonzero(PIP > 0.1)],2)))
print('Credible sets: ' + str(cred))
print('Estimated (sigma^2,tau^2): (%f,%f)\n' % (output['sigmasq'],output['tausq']))

