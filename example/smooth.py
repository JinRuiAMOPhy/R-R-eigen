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
colorlist = ["lime","orange","violet",\
       "green","cyan","red",\
       "blue","gray", "maroon", "chartreuse","cornflowerblue"]
#colorlist = ['#aa489d', '#ad65a7', '#ae7fb2', '#ad97bd', '#a8b0c9', '#9fc9d6', '#8ce2e4', '#00ffff', '#ffd980', '#ffb176', '#f78b6b', '#e96560', '#d44054', '#b81b48', '#93003a']

figsize = (12,7.24*4) # (3.4,2.05)
lab_size = 18; font_size = 12; ticks_size = 16
xrel = 0.07; yrel = 0.8
def ordinal(i) :
   if i == 1 : return r'$1^{st}$'
   elif i == 2 : return r'$2^{nd}$'
   elif i == 3 : return r'$3^{rd}$'
   else : return '{}'.format(i)+r'$^{th}$'

def tidy_axis(**iax) :
   ax = iax['ax']
   ax.tick_params(axis='both', which='both', labelsize= ticks_size)
   xlabel = iax['xlab']; ylabel = iax['ylab']
   if xlabel != '' : ax.set_xlabel(xlabel, fontsize = lab_size)
   ax.set_ylabel(ylabel, fontsize = lab_size)
   #ax.xaxis.set_major_locator(MultipleLocator(5))
   ax.xaxis.set_major_formatter('{x:.3f}')
   ax.yaxis.set_minor_locator(AutoMinorLocator())
   ax.xaxis.set_minor_locator(AutoMinorLocator())
   xmin = iax['xmin']; xmax = iax['xmax']
   ax.set_xlim(xmin, xmax)

def tidy_text(text, **iax) :
   
   xmin, xmax = iax['ax'].get_xlim(); xrel = iax['xrel']
   ymin, ymax = iax['ax'].get_ylim(); yrel = iax['yrel']
   #print(xmin, xmax, ymin, ymax)
   x = xmin + xrel * (xmax - xmin)
   y = ymin + yrel * (ymax - ymin)
   ax.text(x, y, text, fontsize = font_size)




def get_channel(path) :
  files = filter(os.path.isdir, [os.path.join(path, f) for f in os.listdir(path)] )
  res = []
  for ifile in files :
    if 'bk-' in ifile : continue
    res.append(ifile)
  return res


colorscheme = 'chan' #'sec'
action = 'plot'#'; smooth:1o/all' #'6ch,10ch,14ch,19ch,22ch,25ch,28ch' #'smooth;plot'
Z = 1.0
symms = ['0o','0e', '1o', '1e', '2o','2e'] 
update_sym = ['0o', '0e' , '1o', '1e', '2o','2e']
fig = plt.figure(figsize= figsize )
gs = fig.add_gridspec(len(symms),1, hspace=0.3)
axs =  gs.subplots(sharex = False)

Rydbergcminv = 109737.3177
IP = np.array([109837.02, 136647.57, 136667.59, 150305.03, 150307.02,
      229674.23, 229837.45, 229919.88, 275825.48, 275833.52,
      305547.49, 322430.84, 322599.27]) /Rydbergcminv * 13.6058

for iax, isym in enumerate(symms) :
  channels = get_channel('./'+isym)
  print(channels )#, maxchan)
  maxchan = max([int(ich.split('/')[2].split('ch')[0]) for ich in channels])
  ax = axs[iax]
  for isec, ich in enumerate(channels) : 
    if 'smooth' in action and \
         ((ich.split('/')[2] in action ) or 'all' in action) and isym in action :
      if not os.path.isfile('./'+isym+'/smooth.inp') : 
        print('./'+isym+'/smooth.inp is missing')
        exit()
      else : os.system('cp ./'+isym+'/smooth.inp '+ich+'/.')
      cmdl = smooth+' -p '+ich+'/'
      print(cmdl)
      process = Popen(cmdl.split(), stdout=PIPE)
      (output, err) = process.communicate()
      exit_code = process.wait()
      if err is not None : print(err.decode("utf-8"))
      else : print ("*"*10+'succeed'+'*'*10)
      log = open(ich+'/smooth.log', 'w')
      if output is not None : log.write("{:s}".format(output.decode("utf-8")))
      if err is not None : log.write("{:s}".format(err.decode("utf-8")))
      log.close()

    # read in s-matrix even when error happens.
    smat = np.loadtxt(ich+'/miuang.out', unpack=False, comments='#')
    with open(ich+'/debug.log', 'r') as fl :
      lines = fl.readlines()
    eigenchan = lines[-1].split()[1:]
    nch = int(ich.split('/')[2].split('ch')[0])
    for i in range(2, nch+2) :
      if colorscheme == 'sec' : ic = isec % len(colorlist)
      elif  colorscheme == 'chan' : ic = (i-2) % len(colorlist)
      if nch == maxchan :
        #ax.plot(smat[:,1] * Z**2., smat[:,i], 
        ax.plot(IP[0] + 13.6058 * smat[:,1] * Z**2., smat[:,i], 
              color = colorlist[ic], label = r'{}'.format(eigenchan[i-2])) #ordinal(i-1))
        #ax.set_xticks(smat[:,1] * Z**2.)
      else :
        ax.plot(IP[0] + 13.6058 *smat[:,1] * Z**2., smat[:,i], 
              color = colorlist[ic])
        #ax.set_xticks(smat[:,1] * Z**2.)
    xmin = IP[0] + 13.6058 * smat[0,1] * Z**2.; xmax = IP[0] + 13.6058 * smat[-1,1] * Z**2. ; 
    xcenter = 0.5 * xmin + 0.5 * xmax 
    ax.plot([xmin,xmin],[1.0,1.1], color = colorlist[isec])
    ax.plot([xmax,xmax],[1.0,1.1], '--',color = colorlist[isec])
    ax.arrow(xmin,1.0,xmax-xmin,0.0, color = colorlist[isec] )
    ax.text(xcenter, 1.05,'{} ch'.format(nch), color = colorlist[isec])
    if nch == maxchan and colorscheme == 'chan' : ax.legend(loc='upper left', fontsize= font_size) #, ncol=3)
    if isec == 0 :
      ax_dict = {}; ax_dict['xrel'] = xrel; ax_dict['yrel'] = yrel
      ax_dict['xmin'] = IP[0] + 13.6058 * (-0.5); ax_dict['xmax'] = IP[0] + 13.6058 * 1.5
      #ax.set_ylim(-1., 1.)
      ax.tick_params(axis='both', which='both', labelsize= ticks_size)
      ax_dict['ax'] = ax; ax_dict['ylab'] = r'$\mu_\alpha$'; ax_dict['xlab'] = ''
      ax_dict['xlab'] = 'Excitation energy (eV)'
      tidy_axis(**ax_dict)
      text = '({})'.format(chr(ord('a') + iax)) +'symmetry :' +isym; 
      tidy_text(text, **ax_dict)
plt.savefig('miu.pdf', format = 'pdf')

