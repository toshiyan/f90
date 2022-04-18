# Computing rotation reconstruction noise

import os, numpy as np

def add(text,f,ini=False):
  if ini:     os.system('echo "'+text+'" > '+f)
  if not ini: os.system('echo "'+text+'" >> '+f)

def run(tag,fwhm,sP,alpha):
  p = 'params_'+tag+'.ini'
  add('ucl = ../../dat/fid_P15.dat',p,ini=True)
  add('lcl = ../../dat/lensedfid_P15.dat',p)
  add('oL = 2, 3000',p)
  add('rL = 150, 3000',p)
  add('CV = F',p)
  add('fwhm = '+str(fwhm),p)
  add('sP = '+str(sP),p)
  add('alpha = '+str(alpha),p)
  os.system('./exe '+p)
  os.system('mv Aaa.dat dat/'+tag+'.dat')

def loop(t,s,a):
  for j, tj in enumerate(t):
    for i, si in enumerate(s):
      print tj, si
      run(str(tj)+'_'+str(si)+'_a'+str(a),tj,si,a)

# sensitivity in uK-ac
s = np.arange(1.,11.,2.)

# resolution in ac
t = np.array([30.,25.,15.,3.,])
a = .9

loop(t,s,a)

