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

figsize = (12,7.24*1) # (3.4,2.05)
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
   if 'xmin' in iax :
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


Z = 1.0
#symms = ['0o','0e', '1o', '1e', '2o','2e'] 
#update_sym = ['0o', '0e' , '1o', '1e', '2o','2e']
#fig = plt.figure(figsize= figsize )
#gs = fig.add_gridspec(len(symms),1, hspace=0.3)
#axs =  gs.subplots(sharex = False)
fig, ax = plt.subplots(figsize= figsize )
Rydbergcminv = 109737.3177
IP = np.array([109837.02, 136647.57, 136667.59, 150305.03, 150307.02,
      229674.23, 229837.45, 229919.88, 275825.48, 275833.52,
      305547.49, 322430.84, 322599.27]) /Rydbergcminv * 13.6058
IPn = [r'${}^4S_{3/2}^o$', r'${}^2D_{5/2}^o$', r'${}^2D_{3/2}^o$', r'${}^2P_{3/2}^o$',r'${}^2P_{1/2}^o$',
       r'${}^4P_{5/2}^e$', r'${}^4P_{3/2}^e$', r'${}^4P_{1/2}^e$', r'${}^2D_{5/2}^e$',r'${}^2D_{3/2}^e$',
       r'${}^2S_{1/2}^e$', r'${}^2P_{3/2}^e$', r'${}^2P_{1/2}^e$']
# read in s-matrix even when error happens.
#channels = ['6ch','10ch', '14ch', '19ch']; nseg = len(channels)
channels = ['6ch','10ch','14ch','14ch-ok', '14ch-ok+', '14ch+', '19ch', '19ch+', '19ch++', '22ch', '25ch']; nseg = len(channels)
#channels = ['6ch','10ch', '14ch+', '19ch', '19ch+', '19ch++', '25ch']; nseg = len(channels)
#indx = [2, 2, 6, 6, 6, 6, 9]
# 6ch -> ipx = 2; 10ch -> ipx =2; 19ch -> ipx = 6
indx = [2, 2, 2, 4, 4, 6, 6, 6, 6, 9, 9]
for isec, ich in enumerate(channels) :
  os.system("/home/ruijin/RR-eigen-dev/mqdt-util/mqdt-post/postmqdt.exe -p {}".format(ich))
  theo = np.loadtxt(ich+'/theo-plot.out', unpack=False, comments='#', dtype = float)
  if os.path.isfile(ich+'/nux_nuy.out') : yes_dis = True
  else : yes_dis = False
  if yes_dis : nxny = np.loadtxt(ich+'/nux_nuy.out', unpack=False, comments='#', dtype = float)
  tau = np.loadtxt(ich+'/tau2-plot.out', unpack=False, comments='#', dtype = float)
  with open(ich+'/debug.log', 'r') as fl :
    lines = fl.readlines()
  eigenchan = lines[-1].split()[1:]
  #nch = int(ich.split('/')[2].split('ch')[0])
  nch = int(ich.split('ch')[0])
  print(ich,'theo')
  if yes_dis : 
    En2 = IP[indx[isec]] - 1.0 / nxny[:,0]**2 * 13.6058
    #print(En2[0],'<Eph<',En2[-1])
    #ax.plot(nxny[:,0], nxny[:,1], color = 'gray', label = 'Energy relation')
    ax.plot(En2, nxny[:,1], color = 'gray', label = 'Energy relation')
  for i in range(1,nch+1) :
    ic = i % len(colorlist)
    En0 = IP[indx[isec]] - 1.0 / tau[:,0]**2 * 13.6058
    if i == 1 :print(En0[0],'<Eph<',En0[-1])
    ax.plot(En0, tau[:,i] % 1.0, color = colorlist[ic])
    #ax.plot(tau[:,0], tau[:,i], color = colorlist[ic])
    sizemat = theo.shape; nd = len(sizemat); 
    if isec == nseg and nd > 1 :
      En = IP[indx[isec]] - 1.0 / theo[:,0]**2 * 13.6058
      #ax.plot(smat[:,1] * Z**2., smat[:,i], 
      #ax.scatter(theo[:,0], theo[:,i], alpha = 0.8,
      ax.scatter(En, theo[:,i], alpha = 1.,
            color = colorlist[ic], s= 15, label = r'{}'.format(eigenchan[i-1])) #ordinal(i-1))
      #ax.set_xticks(smat[:,1] * Z**2.)
    elif nd > 1:
      En = IP[indx[isec]] - 1.0 / theo[:,0]**2 * 13.6058
      #ax.scatter(theo[:,0], theo[:,i], alpha = 0.8, s= 10,
      ax.scatter(En, theo[:,i], alpha = 1., s= 15,
            color = colorlist[ic])
    elif nd == 1 and sizemat[0]> 0:
      En = IP[indx[isec]] - 1.0 / theo[0]**2 * 13.6058
      ax.scatter([En], [theo[i]], alpha = 1., s= 15,
            color = colorlist[ic])
# channel labels
  xmin = En0[0]; xmax = En0[-1] 
  xcenter = 0.5 * xmin + 0.5 * xmax 
  ax.plot([xmin,xmin],[1.0,1.1], color = colorlist[isec])
  ax.plot([xmax,xmax],[1.0,1.1], '--',color = colorlist[isec], linewidth = 3)
  ax.arrow(xmin,1.0,xmax-xmin,0.0, color = colorlist[isec] )
  ax.text(xcenter, 1.05,'{} ch'.format(ich), color = colorlist[isec])
# 
  if isec == nseg : ax.legend(loc='upper left', fontsize= font_size) #, ncol=3)

# IPs and labels
for i, ip in enumerate(IP) :
  ic = i % len(colorlist)
  ax.plot([ip,ip],[-0.1,0.0], color = colorlist[ic])
  ax.text(ip, -0.15,'{}'.format(IPn[i]), color = colorlist[ic])
# 
ax_dict = {}; ax_dict['xrel'] = xrel; ax_dict['yrel'] = yrel
#ax_dict['xmin'] = -0.5; ax_dict['xmax'] = 1.5
#ax.set_ylim(-1., 1.)
ax.tick_params(axis='both', which='both', labelsize= ticks_size)
ax_dict['ax'] = ax; ax_dict['ylab'] = r'$\tau_\rho$'; ax_dict['xlab'] = ''
ax_dict['xlab'] = 'Excitation energy (eV)'
tidy_axis(**ax_dict)
#text = '({})'.format(chr(ord('a') + iax)) ; 
#tidy_text(text, **ax_dict)
plt.savefig('mqdt.pdf', format = 'pdf')

