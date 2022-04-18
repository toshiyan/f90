# Computing rotation noise spectrum

import numpy as np, matplotlib as mpl
from matplotlib.pyplot import *

# sensitivity in uK-ac
s = np.arange(1.,11.,2.)

# resolution in ac
t = np.array([30.,25.,15.,3.,])
c = ['r','g','b','m','k']
a = 0.0

mpl.rcParams.update({'font.size': 15})
for j, tj in enumerate(t):
  xlabel(r'$L$')
  ylabel(r'$L(N_L^{\alpha\alpha}/2\pi)^{1/2}$ [deg$^2$]')
  xlim(2,1000)
  xscale('log')
  yscale('log')
  ylim(.001,100.)
  grid(True)
  snr = np.zeros(len(s))
  for i, si in enumerate(s):
    L, naa = np.loadtxt('dat/'+str(tj)+'_'+str(si)+'_a'+str(a)+'.dat').T
    plot(L,L*np.sqrt(naa/(2*np.pi))*(180/np.pi),'-',color=c[i],label=str(si)+' $\mu$K\'')
  figtext(.2,.85,r'$\theta_{\rm FWHM}=$'+str(tj)+' arcmin')
  #L, naa = np.loadtxt('dat/plk_a0.0.dat').T
  #plot(L,L*np.sqrt(naa/(2*np.pi))*(180/np.pi),'k--',label='Planck')
  legend(loc=4,frameon=False)
  savefig('fig_naa_fwhm'+str(tj)+'_a'+str(a)+'.png',bbox_inches='tight')
  show()
  clf()

