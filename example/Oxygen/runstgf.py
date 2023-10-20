'''
This script smooths the s-matrix and plot figures
Rui Jin
'''

import numpy as np
import matplotlib.pyplot as plt
import os
from subprocess import Popen, PIPE
from matplotlib.ticker import (MultipleLocator, AutoMinorLocator)

smooth = '/home/ruijin/RR-eigen-dev/s-matrix-prepare/smooth.x'
stgf = '/home/ruijin/RR-eigen-dev/RR-eigen/stgf.x'
stgbf = '/home/ruijin/RR-eigen-dev/RR-eigen/stgbf.x'
prebf = '/home/ruijin/RR-eigen-dev/RR-eigen/prebf.x'
colorlist = ["lime","orange","violet",\
       "green","cyan","red",\
       "blue","gray", "maroon", "chartreuse","cornflowerblue"]

def ordinal(i) :
   if i == 1 : return r'$1^{st}$'
   elif i == 2 : return r'$2^{nd}$'
   elif i == 3 : return r'$3^{rd}$'
   else : return '{}'.format(i)+r'$^{th}$'


def get_channel(path) :
  files = filter(os.path.isdir, [os.path.join(path, f) for f in os.listdir(path)] )
  res = []
  for ifile in files :
    if 'bk-' in ifile : continue
    res.append(ifile)
  return res
def symb_jpi(sym) :
  if 'o' in sym :
    j = sym.split('o')[0]
    pi = 1
  elif 'e' in sym : 
    j = sym.split('e')[0]
    pi = 0
  elif '+' in sym : 
    j = sym.split('+')[0]
    pi = 0
  elif '-' in sym : 
    j = sym.split('-')[0]
    pi = 1
  else : 
    print('check the symmetry block symbol', sym)
    print('it should be j+, j-, jo, je')
    exit()
  return float(j), pi

sampleinput = 'stgf.sample'
sym_name = {'3P2e':'2e', '3P1e':'1e', '3P0e':'0e', '1D2e':'2e', '1S0e':'0e'}
def make_stgf(E0, E1, dE, chanorder, sym, sym_f, nch, run = True, trans = True) :
  j,pi = symb_jpi(sym_f)
  with open(sampleinput, 'r') as sample :
    lines = sample.readlines()
  if not os.path.isdir(sym+"/{}".format(nch)) : os.makedirs(sym+"/{}".format(nch))
  finput = sym+"/{}".format(nch)+'/stgf.inp'
  fo = open(finput, 'w')
  nE = int((E1-E0)/dE) + 1
  lenth = len(lines)
  lines[lenth-2] = ' 0 {} {} 3 0.34\n'.format(int(2.0 * j), pi)
  for i, line in enumerate(lines) :
    if 'E0' in line :
      lines[i] = 'E0={}\n'.format(E0)
    elif 'MXE' in line : 
      lines[i] = 'MXE={}\n'.format(nE)
    elif 'EINCR' in line :
      lines[i] = 'EINCR={}\n'.format(dE)
    fo.write(lines[i])
  print(finput+' modified')#lines)
  fo.close()
  fopen = sym+"/{}".format(nch)+'/openchan'
  fo = open(fopen, 'w')
  fo.write("{} ".format(len(chanorder)))
  for ic in chanorder : fo.write('{} '.format(ic))
  fo.write('\n')
  fo.close()
  print(fopen+' modified')
  if run :
    cwd1 = os.getcwd()
    os.system('cp -v '+finput+' '+sym)
    os.system('cp -v '+fopen+' '+sym)
    os.chdir('./'+sym)
    print(os.getcwd())
    print(stgf+' < stgf.inp > stgf.log 2>&1')
    os.system(stgf+' < stgf.inp > stgf.log 2>&1')
    if trans :
      print(prebf+' < prebf.inp > prebf.log 2>&1')
      os.system(prebf+' < prebf.inp > prebf.log 2>&1')
      print(stgbf+' < stgbf.inp > stgbf.log 2>&1')
      os.system(stgbf+' < stgbf.inp > stgbf.log 2>&1')

    os.system('mv -v ESUM* routf eigenqd MIU-f.OUT OMEGA RK.OUT sigpw.dat stgf.log {}/'.format(nch))

    if trans :
      os.system('mv -v routbf Da*  stgbf.log prebf.log  {}/'.format(nch))
    os.chdir(cwd1)
    print(os.getcwd())

Z = 1.0

# ['0o', '0e', '1o', '1e', '2o','2e'] 
open_ctrl={'0o':{'2ch':{'E': [-0.27,0.05,0.001], 'order' : [1,5]},
                 '5ch':{'E':[0.0, 0.76, 0.001], 'order': [1,5,2,3,4]},
                 '7ch':{'E':[0.75,1.0, 0.001], 'order': [1,5,2,3,4,7,8]},
                 '8ch':{'E':[0.98,1.11,0.001], 'order': [1,5,2,3,4,7,8, 6]},
                 '12ch':{'E':[1.09,1.50,0.001], 'order':[1,5,2,3,4,7,8, 6, 10, 11, 12, 9]},
                 '13ch':{'E':[1.50,2.0,0.001], 'order':[1,5,2,3,4,7,8, 6, 10, 11, 12, 9, 13]}
              },
         '0e':{'1ch':{'E':[-0.5,-0.1,0.001], 'order':[1]},
               '4ch':{'E':[-0.25,0.0,0.001], 'order':[1, 3, 5, 4]},
               '5ch':{'E':[-0.1,0.2,0.001], 'order':[1, 3, 5, 4, 2]},
               '6ch':{'E':[0.1,0.75,0.001], 'order':[1, 3, 5, 4, 2, 8]},
               '9ch':{'E':[0.7,1.3,0.001], 'order': [1, 3, 5, 4, 2, 8,6, 7, 11]},
               '13ch':{'E':[1.29,2.0,0.001], 'order': [1, 3, 5, 4, 2, 8, 6, 7, 11, 9, 10, 12, 13]}
              },
         '1o':{
               '1ch':{'E':[-0.5,-0.1,0.001], 'order':[1]},
               #'4ch':{'E':[-0.25,0.0,0.001], 'order':[1,]},
               '6ch':{'E':[-0.15,0.01,0.001],'order':[1,2,3,7,10,13]},
              '10ch':{'E':[0.05,0.200,0.001],'order':[1,2,3,7,10,13,4,5,8,9]},
              #'13ch':{'E':[0.16,0.27,0.001], 'order':[1,2,3,7,10,13,4,5,8,9,11,12,14]},  'not useful'
              '14ch':{'E':[0.16,0.82,0.001], 'order':[1,2,3,7,10,13,4,5,8,9,11,12,14,6]}, 
              '19ch':{'E':[0.50,1.00,0.001], 'order':[1,2,3,7,10,13,4,5,8,9,11,12,14,6,15,18,19,21,22]}, 
              '22ch':{'E':[0.90,1.20,0.001], 'order':[1,2,3,7,10,13,4,5,8,9,11,12,14,6,15,18,19,21,22,16,17,20]}, 
              '25ch':{'E':[1.10,1.45,0.001], 'order':[1,2,3,7,10,13,4,5,8,9,11,12,14,6,15,18,19,21,22,16,17,20,26,23,27]}, 
              '28ch':{'E':[1.35,1.55,0.001], 'order':[1,2,3,7,10,13,4,5,8,9,11,12,14,6,15,18,19,21,22,16,17,20,26,23,27,24,25,28]},
              '34ch':{'E':[1.55,2.3,0.001], 'order':[1,2,3,7,10,13,4,5,8,9,11,12,14,6,15,18,19,21,22,16,17,20,26,23,27,24,25,28, 29, 30, 31, 32, 33, 34]} 
              },
         '1e':{'2ch':{'E':[-0.5,-0.1,0.001],  'order':[1,2]},
               #'6ch':{'E':[-0.15, 0.02,0.001],'order':[1,2,3,4,7,8]},
              '10ch':{'E':[-0.1,0.1,0.001],   'order':[1,2,3,4,7,8,10,13,11,14]},
              '13ch':{'E':[0.0, 0.25,0.001],  'order':[1,2,3,4,7,8,10,13,11,14,5,6,9]},
              '14ch':{'E':[0.2, 0.70,0.001],  'order':[1,2,3,4,7,8,10,13,11,14,5,6,9,12]},
              '21ch':{'E':[0.6, 0.95,0.001],  'order':[1,2,3,4,7,8,10,13,11,14,5,6,9,12,18,21,15,16,19,20,22]},
              '28ch':{'E':[0.9, 1.4,0.001],   'order':[1,2,3,4,7,8,10,13,11,14,5,6,9,12,18,21,15,16,19,20,22,17,26,23,24,27,28,29]},
              '31ch':{'E':[1.4, 1.6,0.001],   'order':[1,2,3,4,7,8,10,13,11,14,5,6,9,12,18,21,15,16,19,20,22,17,26,23,24,27,28,29,25,31,34]}
             # '35ch':{'E':[1.4, 1.6,0.001],   'order':[1,2,3,4,7,8,10,13,11,14,5,6,9,12,18,21,15,16,19,20,22,17,26,23,24,27,28,29,25,31,34,]}
              },
         '2o':{
               '1ch':{'E':[-0.5,-0.2,0.001],'order':[1]},
               '5ch':{'E':[-0.3, 0.0,0.001],'order':[1,2,3,5,10]},
              '11ch':{'E':[-0.15, 0.15,0.001],'order':[1,2,3,5,10,4,14,6,7,11,12]},
              '15ch':{'E':[ 0.11, 0.17,0.001],'order':[1,2,3,5,10,4,14,6,7,11,12,15,16,18,19]},
              '18ch':{'E':[ 0.15, 0.30,0.001],'order':[1,2,3,5,10,4,14,6,7,11,12,15,16,18,19,8,9,13]},
              '19ch':{'E':[ 0.25, 0.72,0.001],'order':[1,2,3,5,10,4,14,6,7,11,12,15,16,18,19,8,9,13,17]},
              '24ch':{'E':[ 0.70, 1.00,0.001],'order':[1,2,3,5,10,4,14,6,7,11,12,15,16,18,19,8,9,13,17,20,21,24,25,28]},
              #'29ch':{'E':[ 0.80, 1.20,0.001], 'order':[1,2,3,5,10,4,14,6,7,11,12,15,16,18,19,8,9,13,17,20,21,24,25,28,22,23,26,27,29]},
              '33ch':{'E':[ 0.98, 1.40,0.001], 'order':[1,2,3,5,10,4,14,6,7,11,12,15,16,18,19,8,9,13,17,20,21,24,25,28,22,23,26,27,29,30,34,31,35]},
              '37ch':{'E':[ 1.38, 1.50,0.001], 'order':[1,2,3,5,10,4,14,6,7,11,12,15,16,18,19,8,9,13,17,20,21,24,25,28,22,23,26,27,29,30,34,31,35,32,33,36,37]},
              },
         '2e':{
              '2ch':{'E':[-0.5, -0.1,0.001],'order':[1,2]},
             '11ch':{'E':[-0.2,  0.12,0.001],'order':[1,2,3,4,5,6,9,10,13,14,17]},
             '15ch':{'E':[ 0.1,  0.3,0.001],'order': [1,2,3,4,5,6,9,10,13,14,17,7,8,11,12]},
             #'18ch':{'E':[ 0.22, 0.67,0.001],'order':[1,2,3,4,5,6,9,10,13,14,17,7,8,11,12,15,16,18]},
             '20ch':{'E':[ 0.22, 0.90,0.001],'order':[1,2,3,4,5,6,9,10,13,14,17,7,8,11,12,15,16,18,19,24]},
             '26ch':{'E':[ 0.80, 1.10,0.001],'order':[1,2,3,4,5,6,9,10,13,14,17,7,8,11,12,15,16,18,19,24,20,21,25,26,28,29]},
             '31ch':{'E':[ 1.00, 1.35,0.001],'order':[1,2,3,4,5,6,9,10,13,14,17,7,8,11,12,15,16,18,19,24,20,21,25,26,28,29,23,22,27,30,35]},
             '35ch':{'E':[ 1.30, 1.50,0.001],'order':[1,2,3,4,5,6,9,10,13,14,17,7,8,11,12,15,16,18,19,24,20,21,25,26,28,29,23,22,27,30,35,31,32,36,37]},
              }
        }
calc = 'ctrl_list' # or update 

transit = ['1S0e-1o', '1D2e-1o', "3P2e-1o", '3P1e-1o', '3P0e-1o']
sel_ch = ['34ch']
#transit = ['0e-1o', '2e-1o']
#print('finall symm-blocks',symms)
for iax, isym in enumerate(transit) : #symms) :
  symm_f = isym.split('-')[1]
  if (symm_f not in open_ctrl) : continue
  elif (calc == 'update') :
    channels = get_channel('./'+symm_f)
  elif (calc == 'ctrl_list') :
    channels = ["./{}/{}".format(symm_f,i) for i in open_ctrl[symm_f].keys()] 
     
  maxchan = max([int(ich.split('/')[2].split('ch')[0]) for ich in channels])
  #ch_name = [int(ich.split('/')[2].split('ch')[0]) for ich in channels]
  for isec, ich in enumerate(channels) : 
    print(open_ctrl[symm_f])
    ch_name = ich.split('/')[2]
    print(ch_name)
    if ch_name not in open_ctrl[symm_f] : continue
    [E0,E1,dE] = open_ctrl[symm_f][ch_name]['E']
    chanorder = open_ctrl[symm_f][ch_name]['order']
    print(symm_f, ich, ch_name, E0,E1,dE)
    if ch_name in sel_ch :
      make_stgf(E0, E1, dE, chanorder, isym, symm_f, ch_name, run = True, trans = True)
