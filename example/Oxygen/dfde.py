'''
This script smooths the s-matrix and plot figures
Rui Jin
'''

import numpy as np
import matplotlib.pyplot as plt
import os
import numba
from numba import jit
from subprocess import Popen, PIPE
from matplotlib.ticker import (MultipleLocator, AutoMinorLocator)

postmqdt = '/home/ruijin/RR-eigen-dev/mqdt-util/mqdt-post/postmqdt.exe'
colorlist = ["lime","orange","violet",\
       "green","cyan","red",\
       "blue","gray", "maroon", "chartreuse","cornflowerblue"]
tr_cl_map = {'3P0e->1o':0, '3P1e->1o':1, '3P2e->1o':3, '1D2e->1o':4, '1S0e->1o':5}
#colorlist = ['#aa489d', '#ad65a7', '#ae7fb2', '#ad97bd', '#a8b0c9', '#9fc9d6', '#8ce2e4', '#00ffff', '#ffd980', '#ffb176', '#f78b6b', '#e96560', '#d44054', '#b81b48', '#93003a']

global indx

Rydbergcminv = 109737.3177
IP = np.array([109837.02, 136647.57, 136667.59, 150305.03, 150307.02,
      229674.23, 229837.45, 229919.88, 275825.48, 275833.52,
      305547.49, 322430.84, 322599.27]) /Rydbergcminv * 13.6058
IPn = [r'${}^4S_{3/2}^o$', r'${}^2D_{5/2}^o$', r'${}^2D_{3/2}^o$', r'${}^2P_{3/2}^o$',r'${}^2P_{1/2}^o$',
       r'${}^4P_{5/2}^e$', r'${}^4P_{3/2}^e$', r'${}^4P_{1/2}^e$', r'${}^2D_{5/2}^e$',r'${}^2D_{3/2}^e$',
       r'${}^2S_{1/2}^e$', r'${}^2P_{3/2}^e$', r'${}^2P_{1/2}^e$']
sym_latex = {'3P2e':'${}^{3}P^e_2$', '3P1e':'${}^{3}P^e_1$', '3P0e':'${}^{3}P^e_0$',
             '1D2e':'${}^{1}D^e_2$', '1S0e':'${}^{1}S^e_0$'}
lab_size = 18; font_size = 12; ticks_size = 16
xrel = 0.07; yrel = 0.8

xmin = 22; xmax = 40
xmin = 32.5; xmax = 35
#xmin = 25.5; xmax = 26.5

colorscheme = 'chan' #'sec'
#action = 'plot, post:all' #'6ch,10ch,14ch,19ch,22ch,25ch,28ch' #'smooth;plot'
#action = 'plot, post:19ch, 19ch+, 19ch++, 28ch, 28ch+, 28ch-, 34ch, 34ch+, 34ch+_fine1' 
#action = 'plot, post:25ch,28ch,28ch+,34ch;' #'6ch,10ch,14ch,19ch,22ch,25ch,28ch' #'smooth;plot'
#action = 'plot; post:all; fold:1o/all' #'6ch,10ch,14ch,19ch,22ch,25ch,28ch' #'smooth;plot'
#action = 'plot; fold' #, :1o/all' #'6ch,10ch,14ch,19ch,22ch,25ch,28ch' #'smooth;plot'
action = 'plot; post:34ch+++, 34ch-2fast, 34ch-3fast, 34ch-4fast, 34ch-5fast' #, :1o/all' #'6ch,10ch,14ch,19ch,22ch,25ch,28ch' #'smooth;plot'
#action = 'plot; post:34ch+++, 34ch-2fast' #, :1o/all' #'6ch,10ch,14ch,19ch,22ch,25ch,28ch' #'smooth;plot'
action = 'plot; fold-reuse'
Z = 1.0
#symms = ['0e-1o', '1e-1o', '2e-1o']  # for 3P_{2,1,0}
#Egrnd ={'2e':-1.0009086, '1e':-0.99946638,'0e':-0.99884023} this manifold is for 3P_{2,1,0}
#symms = ['0e-1o', '2e-1o'] 
#Egrnd ={'2e':-0.85631,'0e':-0.69296793} # These are for 1D2e and 1S0e 
symms = ['3P2e-1o', '3P1e-1o', '3P0e-1o', '1D2e-1o', '1S0e-1o']
#symms = ['3P2e-1o']#, '3P1e-1o', '3P0e-1o', '1D2e-1o', '1S0e-1o']
Egrnd = {'3P2e':-1.0009086, '3P1e':-0.99946638, '3P0e':-0.99884023, '1S0e':-0.69296793, '1D2e':-0.85631}
#segments = {'1o':['1ch', '6ch', '10ch', '14ch', '14ch-ok', '14ch-ok+', '14ch+', '19ch', '19ch+', '19ch++', '22ch', '25ch']}
segments = {'1o':['6ch', '10ch',   '14ch', '14ch-ok', '14ch-ok+', '14ch+', '19ch', 
                  '19ch+','19ch++','22ch', '25ch', '28ch-', '28ch', '28ch+', '34ch-1',
                  '34ch-2fast', '34ch-3fast', '34ch-4fast', '34ch-5fast',
                  '34ch+++', '34ch+', '34ch+_fine1']}
directory_mqdt = {'1o':'1o-mqdt'}
#segments = {'1o':[ '10ch']} 
#indx = [2, 2, 2, 2, 4, 4, 6, 6, 6, 6, 9, 9]
indx =  np.array([2,     2,        2,      4,         4,          6,       6, 
                   6,     6,        9,      9,      9,       9,      11,     11,       
                   11,           11,           11,          11,     
                   11,        12,       12] ) - 1 # fortran index start from 1#!!!!

reset_xrange = True, 32.5, 34.5



global use_erange
use_erange = {'6ch':'use', '10ch':'use',   '14ch':'use', '14ch-ok':'use', '14ch-ok+':'use', 
              '14ch+':'use', '19ch':'block:[1.54624:1.54877]', '19ch+':'use','19ch++':'use','22ch':'use', '25ch':'use',
              #'28ch-':'use', '28ch':'insert:28ch-;', '28ch+':'use', '34ch':'use', '34ch+':'insert:34ch+_fine1', 
              '28ch-':'use', '28ch':'insert:28ch-;', '28ch+':'use', '34ch-1':'use','34ch-2':'use', 
              '34ch-2fast':'use',  '34ch-3fast':'use',  '34ch-4fast':'use',  '34ch-5fast':'use',
              '34ch+++':'use', '34ch+':'use', 
              '34ch+_fine1':'use'
             }



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
   xmin, xmax = iax['ax'].get_xlim()
   if 'xmin' in iax :
     xmin = iax['xmin']; 
   if 'xmax' in iax :
     xmax = iax['xmax']
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


class MQDT :
  def __init__(self, path) :
    channels = get_channel(path)
 
#sigma FHWM  
def gauss(x, mu, sigma) :
   y = np.exp( -(x-mu)**2 / (sigma**2/4/np.log(2)) )
   y = y / sigma / np.sqrt(np.pi/4/np.log(2))
   return y
def fold(xs, ys, x, sigma) :
  # x is uniform grid
  y = np.zeros_like(x, dtype = np.float64)
  for i, xi in enumerate(x) :
    y[i] = np.trapz(ys * gauss(xs,xi,sigma), x = xs)
  return y  
  

def plot_dfde(ax, dir_i, isec, seg, nseg, sym_i, sym_f, action, legend, width = 0.05) : # width FHWM in eV
  sig = np.loadtxt(dir_i+'/{}-{}-{}.os'.format(sym_i, sym_f,seg), unpack=False, comments='#', dtype = np.float64)
# not that the energy in os.out is I_x - Z^2/nx**2 (with I1 set to zero)
  xs = sig[:,0] * 13.6058 
  print('dfde', sym_i+ '->' +sym_f,xs[0],xs[-1])
  ys = sig[:,1]
  #xn = np.linspace(xs[0]+width,xs[-1]-width,num = 1000) 
  xn = np.linspace(xs[0],xs[-1],num = 1000) 
  cl = tr_cl_map[sym_i+'->'+sym_f]
  if 'fold' in action :
    #print(width,'width')
    fn = dir_i+'/{}-{}-{}-fold-{:.2e}eV.os'.format(sym_i, sym_f,seg, width) 
    if os.path.isfile(fn) and 'fold-reuse' in action: 
      xp, yp = np.loadtxt(fn, unpack = True, comments='#', dtype = np.float64)
    else :
      xp = xn; yp = fold(xs, ys, xn, width)
      sig_fold = np.vstack((xp,yp)).transpose()
      np.savetxt(fn, sig_fold, fmt='%.18e', delimiter=' ', newline='\n', header='', footer='', comments='# ', encoding=None)
    ax.plot(xs, ys, ':', color = colorlist[cl], alpha = 0.1)
    #exit()
  else :
    xp = xs; yp = ys
  if legend:
    ax.plot(xp, yp, color = colorlist[cl], label = sym_latex[sym_i] +'->'+sym_f)
  else :
    ax.plot(xp, yp, color = colorlist[cl])
  return ax

global states_in_chan 
states_in_chan = np.zeros(34) # max eigenchan for final sym 1o

# for Oxgen only
least_pqn = {'s':3,'p':3,'d':3,'f':4,'g':5,'h':6}
def get_pqn(eigenchan, ich) :
  terms = eigenchan.split('}')
  orb = terms[2].split('_')[0]
  n = least_pqn[orb] + states_in_chan[ich]
  states_in_chan[ich] += 1
  terms[2] = str(int(n)) + terms[2]
  #print('}'.join(terms), states_in_chan[ich])
  return '}'.join(terms) , states_in_chan[ich]

def sel_tau_2_plot(plot_erange, tau0, sym_f) :
  if plot_erange is None : tau = tau0
  elif 'use' in plot_erange :  tau = tau0
  elif 'block' in plot_erange : 
    s = plot_erange
    s__ = s[s.find("[")+1:s.find("]")] 
    id_ = np.array([True] * len(tau0))
    print('block',s__)
    for __ in s__.split(';') :
      [es,ee] = [float(i) for i in __.split(':')]
      id_ = id_ * (tau0[:,0] >= es) * (tau0[:,0] <= ee) 
    id_ = [ not i for i in id_ ]
    tau = tau0[id_]
  elif 'insert' in plot_erange : 
    final_mqdt_dir = directory_mqdt[sym_f]
    s__ = plot_erange.split(':')[1].split(';') 
    id_ = np.array([True] * len(tau0))
    for s_ in s__ :
      if s_ == '' : continue
      s = final_mqdt_dir +'/'+s_ + '/debug-mqdt.out' 
      print('search file',s)
      __ = ''
      with open(s, 'r') as FL : 
        while True : 
          line = FL.readline()
          if '# E range in terms of nu' in line : 
            xx = line.split(')')[1]
            if '(' in xx : __ = xx[:xx.find('(')] 
            else : __ = xx
            break
      print('insert',xx)
      [es,ee] = [float(i) for i in __.split(':')]
      id_ = id_ * (tau0[:,0] >= es) * (tau0[:,0] <= ee) 
    id_ = [ not i for i in id_ ]
    tau = tau0[id_]
  return tau

def plot_tau(ax, dir_i, isec, seg, nseg, Eg, do_plot = True, plot_erange = None, sym_f = None) :
  theo = np.loadtxt(dir_i+'/theo-plot.out', unpack=False, comments='#', dtype = float)
  if os.path.isfile(dir_i+'/nux_nuy.out') : yes_dis = True
  else : yes_dis = False
  if yes_dis : nxny = np.loadtxt(dir_i+'/nux_nuy.out', unpack=False, comments='#', dtype = float)
  tau0 = np.loadtxt(dir_i+'/tau2-plot.out', unpack=False, comments='#', dtype = float)
  tau = sel_tau_2_plot(plot_erange, tau0, sym_f)
  print(seg, tau0.shape, tau.shape)
  with open(dir_i+'/debug.log', 'r') as fl :
    lines = fl.readlines()
  eigenchan = lines[-1].split()[1:]
  #nch = int(ich.split('/')[2].split('ch')[0])
  nch = int(seg.split('ch')[0])
  print(seg,'theo',dir_i, theo.shape)
  elev = np.empty(0, dtype = float)
  notes = np.empty(0, dtype = str)
  notes__ = np.empty(0, dtype = str)
  max_n = np.empty(0, dtype = int)


  if yes_dis : 
    En2 = IP[indx[isec]] - 1.0 / nxny[:,0]**2 * 13.6058 - Eg * 13.6058
    if do_plot : ax.plot(En2, nxny[:,1], color = 'gray', label = 'Energy relation')
  for i in range(1,nch+1) :
    ic = i % len(colorlist)
    En0 = IP[indx[isec]] - 1.0 / tau[:,0]**2 * 13.6058 - Eg * 13.6058
    if i == 1 :print(En0[0],'<Eph<',En0[-1])
    if do_plot : ax.plot(En0, tau[:,i] % 1.0, color = colorlist[ic])
    sizemat = theo.shape; nd = len(sizemat); 
    if isec == nseg and nd > 1 :
      En = IP[indx[isec]] - 1.0 / theo[:,0]**2 * 13.6058 - Eg * 13.6058
      if do_plot : ax.scatter(En, theo[:,i], alpha = 1.,
            color = colorlist[ic], s= 15, label = r'{}'.format(eigenchan[i-1])) #ordinal(i-1))
      elev = np.concatenate((elev, [En[j] for j, v in enumerate(theo[:,i]) if not np.isnan(v)] ))
      notes = np.concatenate((notes, [eigenchan[i-1] for v in theo[:,i] if not np.isnan(v)] ))
      notes__ = np.concatenate((notes__, [get_pqn(eigenchan[i-1], i-1)[0] for v in theo[:,i] if not np.isnan(v)] ))
      max_n = np.concatenate((max_n, [get_pqn(eigenchan[i-1], i-1)[1] for v in theo[:,i] if not np.isnan(v)] ))
    elif nd > 1:
      En = IP[indx[isec]] - 1.0 / theo[:,0]**2 * 13.6058 - Eg * 13.6058
      if do_plot : ax.scatter(En, theo[:,i], alpha = 1., s= 15,
            color = colorlist[ic])
      elev = np.concatenate((elev, [En[j] for j, v in enumerate(theo[:,i]) if not np.isnan(v)] ))
      notes = np.concatenate((notes, [eigenchan[i-1] for v in theo[:,i] if not np.isnan(v)] ))
      notes__ = np.concatenate((notes__, [get_pqn(eigenchan[i-1], i-1)[0] for v in theo[:,i] if not np.isnan(v)] ))
      max_n = np.concatenate((max_n, [get_pqn(eigenchan[i-1], i-1)[1] for v in theo[:,i] if not np.isnan(v)] ))
    elif nd == 1 and sizemat[0]> 0:
      En = IP[indx[isec]] - 1.0 / theo[0]**2 * 13.6058 - Eg * 13.6058
      if do_plot : ax.scatter([En], [theo[i]], alpha = 1., s= 15,
            color = colorlist[ic])
      elev = np.concatenate((elev, [En]))
      notes = np.concatenate((notes, [eigenchan[i-1]] ))
      notes__ = np.concatenate((notes__, [get_pqn(eigenchan[i-1], i-1)[0]] ))
      max_n = np.concatenate((max_n, [get_pqn(eigenchan[i-1], i-1)[1]] ))
    
# channel labels
  xmin = En0[0]; xmax = En0[-1] 
  xcenter = 0.5 * xmin + 0.5 * xmax 
  if do_plot :
    ax.plot([xmin,xmin],[1.05,1.1], color = colorlist[isec%len(colorlist)])
    ax.plot([xmax,xmax],[1.05,1.1], '--',color = colorlist[isec%len(colorlist)], linewidth = 3)
    ax.arrow(xmin,1.05,xmax-xmin,0.0, color = colorlist[isec%len(colorlist)] )
    ax.text(xcenter, 1.07,'{} ch'.format(nch), color = colorlist[isec%len(colorlist)])
  return ax, elev, notes, notes__, max_n
# 
   # ax[0].legend(loc='upper left', fontsize= font_size) #, ncol=3)




sel_mqdt_plot = 0

symm_f_list = {}



for isym in symms :
  if '-' in isym : sym_f = isym.split('-')[1]
  else : exit()
  if '-' in isym : sym_i = isym.split('-')[0]
  else : exit()
  if not sym_f in symm_f_list : symm_f_list[sym_f] = [sym_i]
  else : symm_f_list[sym_f].append(sym_i)
print(symm_f_list,len(symm_f_list))

resonances = np.empty(0, dtype = float)
trans_note = np.empty(0, dtype = str)
pqn        = np.empty(0, dtype = str)
max_n_pq   = np.empty(0, dtype = int)

for sym_f in symm_f_list.keys() :
  figsize = (16,12.24) # (3.4,2.05)
  fig = plt.figure(figsize= figsize )
  gs = fig.add_gridspec(1+len(symm_f_list[sym_f]),1, hspace=0.1) #(a) tau, (b) dfde
  axs =  gs.subplots(sharex = True)
  figname = 'mqdt-{}'.format(sym_f)
  if 'fold' in action :
    figname +=  '-fold'
  if reset_xrange[0] : 
    figname += '-range-{:.3f}-{:.3f}'.format(reset_xrange[1], reset_xrange[2])
  else :
    figname += 'full-range'
  figname += '.pdf'
  istate = 0; n_istate = len(symm_f_list[sym_f])
  for  n_sym_i, sym_i in enumerate(symm_f_list[sym_f]) :
    states_in_chan = np.zeros(34) # max eigenchan for final sym 1o
    print(sym_i,'->',sym_f)
    segs = segments[sym_f] #get_channel('./{}-mqdt'.format(sym_f)) 
    nseg = len(segs)
    maxchan = max([int(seg.split('ch')[0]) for seg in segs])
    print(segs);# exit()
    for isec, seg in enumerate(segs) :
      #if 'post' in action and (seg == '22ch'):
      if 'post' in action  and (seg in action or 'all' in action) :
        dm = '{}-{}/{}ch/DalphaL-smooth.1.out'.format(sym_i,sym_f,seg.split('ch')[0])
        cmd = 'cp -v {} {}-mqdt/{}/dalfa.in'.format(dm, sym_f, seg)
        os.system(cmd)
        cmd = postmqdt+' -p {}-mqdt/{}/ -eg {}'.format(sym_f, seg, Egrnd[sym_i])
        print(cmd)
        process = Popen(cmd.split(), stdout=PIPE)
        (output, err) = process.communicate()
        exit_code = process.wait()
        if err is not None : print("*"*10+'something wrong'+'*'*10)#err.decode("utf-8"))
        else : print ("*"*10+'succeed'+'*'*10)
        log = open(sym_f+'-mqdt/'+seg+'/post-dfde.log', 'w')
        if output is not None : log.write("{:s}".format(output.decode("utf-8")))
        if err is not None : log.write("{:s}".format(err.decode("utf-8")))
        log.close()
        dir_i = '{}-mqdt/{}/'.format(sym_f, seg); dir_f = '{}-mqdt/{}/{}-{}-{}.os'.format(sym_f, seg, sym_i,sym_f, seg)
        cmd = 'cp -v {}/os.out {}'.format(dir_i, dir_f)
        os.system(cmd)
      # end if post 
      dir_i = '{}-mqdt/{}/'.format(sym_f, seg); dir_f = '{}-mqdt/{}/{}-{}-{}.os'.format(sym_f, seg, sym_i,sym_f, seg)
      if isec == nseg - 1 : legend = True
      else : legend = False
      print('isec',isec,seg,'nseg',nseg,'istate', istate, 'n_istate', n_istate - 1)
      if istate == sel_mqdt_plot:  # only plot the final channel QDT onece for different initial transitions.
        print('plot:',seg)
        do_plot =  True
      else : do_plot = False
      plot_erange = use_erange[seg]
      axs[0], elev, notes, notes__, max_n = plot_tau(axs[0], dir_i, isec, seg, nseg, Egrnd[sym_i]- Egrnd['3P2e'], 
           do_plot, plot_erange, sym_f) 
      #print('elev', elev)
      if reset_xrange[0] : axs[0].set_xlim(reset_xrange[1], reset_xrange[2])
      resonances = np.concatenate((resonances, elev)) # in eV
      notes = [sym_latex[sym_i]+'$\\to$'+lab for lab in notes]
      trans_note = np.concatenate((trans_note, notes))
      pqn = np.concatenate((pqn, notes__))
      max_n_pq = np.concatenate((max_n_pq, max_n))
      axs[n_sym_i+1] = plot_dfde(axs[n_sym_i+1], dir_i, isec, seg, nseg, sym_i, sym_f, action, legend)
      if reset_xrange[0] : axs[n_sym_i+1].set_xlim(reset_xrange[1], reset_xrange[2])
    istate += 1

  # IPs and labels
  for i, ip in enumerate(IP) :
    ic = i % len(colorlist)
    axs[0].plot([ip,ip],[-0.1,0.0], color = colorlist[ic])
    axs[0].text(ip, -0.15,'{}'.format(IPn[i]), color = colorlist[ic])
  #
 
  ax_dict = {}; ax_dict['xrel'] = xrel; ax_dict['yrel'] = yrel
  #ax_dict['xmin'] = 20.0
  #ax_dict['xmin'] = -0.5; ax_dict['xmax'] = 1.5
  axs[0].tick_params(axis='both', which='both', labelsize= ticks_size)
  ax_dict['ax'] = axs[0]; ax_dict['ylab'] = r'Phase shift $\tau_\rho$'; ax_dict['xlab'] = ''
  ax_dict['xlab'] = '' 
  tidy_axis(**ax_dict)

  ax_dict = {}; ax_dict['xrel'] = xrel; ax_dict['yrel'] = yrel
  #ax_dict['xmin'] = 20.0
  for iax in range(1, 1+len(symm_f_list[sym_f])) :
    axs[iax].set_ylim(0.0, 2.e2)
    axs[iax].set_xlim(xmin, xmax)
    axs[iax].tick_params(axis='both', which='both', labelsize= ticks_size)
    ax_dict['ax'] = axs[iax]; 
    if iax == int(len(symm_f_list[sym_f]) / 2) :
      ax_dict['ylab'] = 'Oscilator strength (Density) (arb. unit)'; 
    else : ax_dict['ylab'] = ''
    if iax == len(symm_f_list[sym_f]) : ax_dict['xlab'] = 'Photon energy (eV)' 
    else : ax_dict['xlab'] = ''
    tidy_axis(**ax_dict)
    axs[iax].legend(loc='upper right')
    #axs[iax].set_yscale('log')
    #text = '({})'.format(chr(ord('a') + iax)) ; 
    #tidy_text(text, **ax_dict)
  plt.savefig(figname, format = 'pdf')
  
trans_latex = open('trans.tex', 'w')

lev_sort = np.argsort(resonances)
new_notes = trans_note[lev_sort]
new_lev = resonances[lev_sort]
new_pqn = pqn[lev_sort]
new_n = max_n_pq[lev_sort]
for e, lab, n, ns in zip(new_lev, new_notes, new_pqn, new_n) :
  if e >= xmin and e <= xmax and ns < 10 :
    trans_latex.write("{}, {}, {} \\newline\n".format(e, lab, n))
    
