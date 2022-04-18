# Computing polarization rotation amplitude constraints as a function of cmb noise (uK')

import os, numpy as np, matplotlib as mpl
from matplotlib.pyplot import *

def add(text,f,ini=False):
  if ini:     os.system('echo "'+text+'" > '+f)
  if not ini: os.system('echo "'+text+'" >> '+f)

def run(tag,lmax,fwhm,sP,alpha):
  p = 'params_'+tag+'.ini'
  add('ucl = ../../dat/fid_P13.dat',p,ini=True)
  add('lcl = ../../dat/lensedfid_P13.dat',p)
  add('oL = 2, 300',p)
  add('rL = 150, '+str(lmax),p)
  add('CV = F',p)
  add('fwhm = '+str(fwhm),p)
  add('sP = '+str(sP),p)
  add('alpha = 1.',p)
  os.system('./exe '+p)
  os.system('mv al.dat dat/al_'+tag+'.dat')

def loop(t,s):
  for j, tj in enumerate(t):
    lmax = np.int(500*30/tj)
    for i, si in enumerate(s):
      print tj, si
      run(str(tj)+'_'+str(si),lmax,tj,si)

# model
l   = np.arange(2,301,1)
caa = 2*np.pi*1e-4/(l**2+l)

# sensitivity in uK-ac
s = np.arange(1.,11.,2.)

# resolution in ac
t = np.arange(30.,20.,-5.)

loop(t,s)

mpl.rcParams.update({'font.size': 15})
xlabel(r'Sensitivity [$\mu$K-arcmin]')
#ylabel(r'$\sigma(A_{CB})$')
ylabel(r'SNR')
yscale('log')
grid(True)
c = ['r','g','b','m']
for j, tj in enumerate(t):
  snr = np.zeros(len(s))
  for i, si in enumerate(s):
    naa = np.loadtxt('dat/al_'+str(tj)+'_'+str(si)+'.dat').T[1]
    snr[i] = np.sqrt(np.sum(.01*(l+.5)*caa**2/(caa+naa)**2))
  plot(s,snr,'-',color=c[j],label=r'$\theta_{\rm FWHM}=$'+str(tj)+' arcmin')
  plot(s,snr*np.sqrt(5.),'--',color=c[j])
  plot(s,snr*np.sqrt(10.),':',color=c[j])
figtext(.2,.85,'solid: 1%')
figtext(.2,.8,' dashed: 5%')
figtext(.2,.75,'dotted: 10%')
legend(loc=0,frameon=False)
savefig('test.png',bbox_inches='tight')
show()

