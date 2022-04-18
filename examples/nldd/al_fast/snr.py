# Computing rotation amplitude constraints

import numpy as np, matplotlib as mpl
from matplotlib.pyplot import *

# model
l   = np.arange(2,3001,1)
caa = 2*np.pi*1e-5/(l**2+l)

# sensitivity in uK-ac
s = np.arange(1.,11.,2.)

# resolution in ac
t = np.array([30.,25.,15.,3.,])
a = 0.0
#naa = np.loadtxt('dat/plk_a0.0.dat').T[1]
#plk = 1./np.sqrt(np.sum(.7*(l+.5)*caa**2/(naa[1:])**2))

lmin = 4
mpl.rcParams.update({'font.size': 15})
for j, tj in enumerate(t):
  xlabel(r'Sensitivity [$\mu$K-arcmin]')
  ylabel(r'$\sigma(A_{CB})$')
  yscale('log')
  ylim(.00005,1.)
  grid(True)
  snr = np.zeros(len(s))
  for i, si in enumerate(s):
    naa = np.loadtxt('dat/'+str(tj)+'_'+str(si)+'_a'+str(a)+'.dat').T[1]
    #snr[i] = 1./np.sqrt(np.sum(.01*(l+.5)*caa**2/(caa+naa)**2))
    snr[i] = 1./np.sqrt(np.sum(.01*(l[lmin:]+.5)*caa[lmin:]**2/(naa[lmiin:])**2))
  plot(s,snr,'r-',label='1%')
  plot(s,snr/np.sqrt(10.),'g--',label='10%')
  plot(s,snr/np.sqrt(70.),'b:',label='70%')
  figtext(.2,.75,r'$\theta_{\rm FWHM}=$'+str(tj)+' arcmin')
  #axhline(plk,ls='--',color='k',label='Planck')
  legend(loc=4,frameon=False)
  savefig('fig_snr_fwhm'+str(tj)+'_a'+str(a)+'.png',bbox_inches='tight')
  show()
  clf()

