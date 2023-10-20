#  Z  A   mass 
import os
import copy
import numpy as np
from subprocess import Popen, PIPE

import sys
sys.path.append('/home/ruijin/')
from mypylib import grasp_mcscf as grasp

dbg = {'cfg':False, 'wave':False, 'nuc':False,
       'angular':False, 'scf':False}

def main() :
  Reigen = {'idiag':0, 'ipolph': 2, 'n_phy_targ':4, 'sym_blocks':[grasp.symmetry(0,1), grasp.symmetry(1,-1)]}
  mod1 = grasp.calc_setup('mod1', clear = True, clear_what = ['All'], **Reigen)
  mod1.Clight(137.035999139000)
  dbg_control = dbg # defined at the head of file
  mod1.Debug(dbg_control)
  mod1.Grid({'RHN':2.500000000000000E-007,
             'H'  :5.000000000000000E-002,
             'HP' :0.000000000000000E-001,
             'N'  :590})
  mod1.NucInfo(grasp.ElemTab['O'], NucDetail = {"NucSpin":1.0, "Dipole":1.0, 
                                           "Qpole":1.0,  
                                           "RNucRMS":0.1, "RNucSkin":0.05})
  mod1.CSFBasis({#'InactShell':'1s2 ', 'J':1, 
                  'CloseCore':'He', 'JLow':1, 'JHigh':5,  # 2J
                  'config':['orbs:(n1,n2);ref:2s2_2p3',
                            'orbs:(n1,n2);ref:2s1_2p4',
                            'orbs:(n1,n2);ref:2p5']
                 })
  #mod1.WaveFunc(Control = {'type':'n1,n2:TF;3sp,3d-:H;3d+:./rwfn.out'}) 
  mod1.WaveFunc(Control = {'type':'n1,n2:TF'})
  mod1.Angular(Control ={'kind':'parallel', 'interact':'full'})
  mod1.SCF(Control = {#'QED':{},
                      'maxiter':100,
                      'spec_orb':['*'],
                      'update_orb':['*'],
                      'accuracy':'1.e-8',
                      #'block' : {'weight':{'kind':'user', 'ratio':[0.7,0.3]}, 
                      'block' : {'weight':{'kind':'standard'},
                                 'level': {'3/2-':[1],'5/2-':[]}, # opt. 1s,2s,2p
                      'order' : 'update' # or scc -- self consistency connected
                    # 'converg':{}   
                                }
                     })
  mod1.CI(Control = {} 
         )
  #mod1.execute('Nuc,CSF,split,Wave,Angular,SCF,CI') 
  # not every module need to be executed.
  #mod1.execute('Nuc,csf,Wave,Angular, scf') 
  
  del mod1

  Reigen = {'idiag':0, 'ipolph': 2, 'n_phy_targ':9, 'sym_blocks':[grasp.symmetry(0,1), grasp.symmetry(1,-1)]}
  mod2 = grasp.calc_setup('mod2', clear = True, clear_what = ['All'], **Reigen)
  mod2.NucInfo(grasp.ElemTab['O'], NucDetail = {"NucSpin":1.0, "Dipole":1.0, 
                                          "Qpole":1.0,  
                                          "RNucRMS":0.1, "RNucSkin":0.05})
  mod2.CSFBasis({#'InactShell':'1s2 ', 'J':1, 
                  'CloseCore':'He', 'JLow':1, 'JHigh':5,  # 2J
                  'config':['orbs:(n1,n2,3s);ref:2s2_2p3;Excite:n2->n3(S,D)',
                            'orbs:(n1,n2,3s);ref:2s1_2p4;Excite:n2->n3(S,D)',
                            'orbs:(n1,n2,3s);ref:2p5;Excite:n2->n3(S,D)'] 
                 })
  #mod1.WaveFunc(Control = {'type':'n1,n2:TF;3sp,3d-:H;3d+:./rwfn.inp'}) 
  mod2.WaveFunc(Control = {'type':'n3:TF;n1,n2:mod1/rwfn.out'})
  mod2.Angular(Control ={'kind':'parallel', 'interact':'full'})
  mod2.SCF(Control = {#'QED':{},
                      'maxiter':100,
                      'spec_orb':['*'],
                      'update_orb':['3*'],
                      'accuracy':'1.e-8',
              # 3s, 3p is only slightly mixed in first few levels, one can manually add more weight.
              # to achieve SCF convergence
                      #'block' : {'weight':{'kind':'user', 'ratio':[0.3,0.2,0.3,0.2]}, 
                      'block' : {'weight':{'kind':'standard'}, 
                                 'level': {'1/2+':[1],   # opt. 3s
                                           '1/2-':[1]},  # opt. 1s,2s,2p
                                           #'3/2-':[1], '5/2-':[],
                                           #'5/2+':[1,2,3]}, 
                      'order' : 'update' # or scc -- self consistency connected
                    # 'converg':{}   
                                }
                     })
  mod2.CI(Control = {} 
         )

  del mod2

  Reigen = {'idiag':0, 'ipolph': 2, 'n_phy_targ':9, 'sym_blocks':[grasp.symmetry(0,1), grasp.symmetry(1,-1)]}
  mod3 = grasp.calc_setup('mod3', clear = True, clear_what = ['All'], **Reigen)
  mod3.NucInfo(grasp.ElemTab['O'], NucDetail = {"NucSpin":1.0, "Dipole":1.0, 
                                          "Qpole":1.0,  
                                          "RNucRMS":0.1, "RNucSkin":0.05})
  mod3.CSFBasis({#'InactShell':'1s2 ', 'J':1, 
                  'CloseCore':'He', 'JLow':1, 'JHigh':5,  # 2J
                  'config':['orbs:(n1,n2,3spd);ref:2s2_2p3;Excite:n2->n3(S,D)',
                            'orbs:(n1,n2,3spd);ref:2s1_2p4;Excite:n2->n3(S,D)',
                            'orbs:(n1,n2,3spd);ref:2p5;Excite:n2->n3(S,D)'] 
                 })
  #mod1.WaveFunc(Control = {'type':'n1,n2:TF;3sp,3d-:H;3d+:./rwfn.inp'}) 
  mod3.WaveFunc(Control = {'type':'3pd:TF;n1,n2,3s:mod2/rwfn.out'})
  mod3.Angular(Control ={'kind':'parallel', 'interact':'full'})
  mod3.SCF(Control = {#'QED':{},
                      'maxiter':100,
                      'spec_orb':['*'],
                      'update_orb':['3p*', '3d*'],
                      'accuracy':'1.e-8',
              # 3s, 3p is only slightly mixed in first few levels, one can manually add more weight.
              # to achieve SCF convergence
                      #'block' : {'weight':{'kind':'user', 'ratio':[0.3,0.2,0.3,0.2]}, 
                      'block' : {'weight':{'kind':'standard'}, 
                                 'level': {'1/2+':[1],   # opt. 3s
                                           '1/2-':[1],   # opt. 1s,2s,2p,3p
                                           '3/2+':[1]},  # opt. 3d
                                           #'3/2-':[1], '5/2-':[],
                                           #'5/2+':[1,2,3]}, 
                      'order' : 'update' # or scc -- self consistency connected
                    # 'converg':{}   
                                }
                     })
  mod3.CI(Control = {} 
         )

  del mod3


  Reigen = {'idiag':0, 'ipolph': 2, 'n_phy_targ':9, 'sym_blocks':[grasp.symmetry(0,1), grasp.symmetry(1,-1)]}
  mod4 = grasp.calc_setup('mod4', clear = True, clear_what = ['All'], **Reigen)
  mod4.NucInfo(grasp.ElemTab['O'], NucDetail = {"NucSpin":1.0, "Dipole":1.0, 
                                          "Qpole":1.0,  
                                          "RNucRMS":0.1, "RNucSkin":0.05})
  mod4.CSFBasis({#'InactShell':'1s2 ', 'J':1, 
                  'CloseCore':'He', 'JLow':1, 'JHigh':5,  # 2J
                  'config':['orbs:(n1,n2,3spd,n4);ref:2s2_2p3;Excite:n2->n3,n4(S,D)',
                            'orbs:(n1,n2,3spd,n4);ref:2s1_2p4;Excite:n2->n3,n4(S,D)',
                            'orbs:(n1,n2,3spd,n4);ref:2p5;Excite:n2->n3,n4(S,D)'] 
                 })
  #mod1.WaveFunc(Control = {'type':'n1,n2:TF;3sp,3d-:H;3d+:./rwfn.inp'}) 
  mod4.WaveFunc(Control = {'type':'n4:TF;n1,n2,n3:mod3/rwfn.out'})
  mod4.Angular(Control ={'kind':'parallel', 'interact':'full'})
  mod4.SCF(Control = {#'QED':{},
                      'maxiter':100,
                      'spec_orb':[],
                      'update_orb':['4*'],
                      'accuracy':'1.e-8',
              # 3s, 3p is only slightly mixed in first few levels, one can manually add more weight.
              # to achieve SCF convergence
                      #'block' : {'weight':{'kind':'user', 'ratio':[0.3,0.2,0.3,0.2]}, 
                      'block' : {'weight':{'kind':'standard'}, 
                                 'level': {'1/2+':[1],   # opt. 3s
                                           '1/2-':[1],   # opt. 1s,2s,2p,3p
                                           '3/2+':[1]},  # opt. 3d
                                           #'3/2-':[1], '5/2-':[],
                                           #'5/2+':[1,2,3]}, 
                      'order' : 'update' # or scc -- self consistency connected
                    # 'converg':{}   
                                }
                     })
  mod4.CI(Control = {} 
         )

  del mod4


  exit()

  Reigen = {'idiag':0, 'ipolph': 2, 'n_phy_targ':9, 'sym_blocks':[grasp.symmetry(0,1), grasp.symmetry(1,-1)]}
  mod3 = grasp.calc_setup('mod3', clear = True, clear_what = ['All'], **Reigen)
  mod3.NucInfo(grasp.ElemTab['O'], NucDetail = {"NucSpin":1.0, "Dipole":1.0, 
                                          "Qpole":1.0,  
                                          "RNucRMS":0.1, "RNucSkin":0.05})
  mod3.CSFBasis({#'InactShell':'1s2 ', 'J':1, 
                  'CloseCore':'He', 'JLow':1, 'JHigh':7,  # 2J
                  'config':['orbs:(n1,n2,n3);ref:2s2_2p3;Excite:n2->n3(S,D)',
                            'orbs:(n1,n2,n3);ref:2s1_2p4;Excite:n2->n3(S,D)',
                            'orbs:(n1,n2,n3);ref:2p5;Excite:n2->n3(S,D)'] 
                 })
  mod3.WaveFunc(Control = {'type':'3d:TF;n1,n2,3sp:mod2/rwfn.out'})
  mod3.Angular(Control ={'kind':'parallel', 'interact':'full'})
  mod3.SCF(Control = {#'QED':{},
                      'maxiter':100,
                      'spec_orb':[],
                      'update_orb':['3d*'],
                      'accuracy':'1.e-8',
                      #'block' : {'weight':{'kind':'user', 'ratio':[0.3,0.2,0.3,0.2]}, 
                      'block' : {'weight':{'kind':'standard'},
                                 'level': {'1/2+':[1,2],   #1/2 to opt. 3s
                                           '1/2-':[1,2],  # 1/2- to opt. 3p
                                           '3/2-':[1], '5/2-':[],
                                           '7/2+':[1,2]},
                      'order' : 'update' # or scc -- self consistency connected
                    # 'converg':{}   
                                }
                     })
  mod3.CI(Control = {} 
         )
  #mod1.execute('Nuc,CSF,split,Wave,Angular,SCF,CI,Rmat') 
  # not every module need to be executed.
  #mod2.execute('Nuc,csf,Wave,Angular, scf') #,split,Wave,Angular,SCF,CI,Rmat') 
#
# main 
if __name__ == '__main__' :
  main()
