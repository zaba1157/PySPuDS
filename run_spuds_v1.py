# -*- coding: utf-8 -*-
"""
Created on Wed Sep  4 16:50:58 2019

@author: zaba1157
"""

import os
from shutil import move
import warnings
import subprocess
import re
from pymatgen.io.cif import CifBlock,CifFile,CifParser
from collections import OrderedDict
from monty.io import zopen
from pymatgen.io.ase import AseAtomsAdaptor
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from pymatgen.transformations.advanced_transformations import OrderDisorderedStructureTransformation
from ase.io import read
from pymatgen.io.vasp.inputs import Poscar
###
import numpy as np
import pandas as pd
from pymatgen.core.periodic_table import Specie, Element
#from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from pymatgen.symmetry.structure import SymmetrizedStructure
#import warnings
#from functools import lru_cache
###

OD = OrderDisorderedStructureTransformation()
opj = os.path.join





class SPuDS:
    '''
    Inputs:
        A: dict of A cations and oxidation states, {Ce:4,La:3} or {Gd:3}
        B: dict of B cations and oxidation states, {Al:3,Mn:4} or {Zr:4}
        X: dict of X anions and oxidation states, {O:2} #assumes (-) oxi_state
        
    Allowed menu items for > ternary perovskites:
    
    3. AA'B2X6
    6. A2BB'X6
    7. AA'BB'X6
    
    
    '''
    def __init__(self, As, Bs, X, store_dir = 'SPuDS_output'):
        self.As = As
        self.Bs = Bs
        self.X = X
        self.Xel = [el for el in X.keys()][0]
        els = [el for el in As.keys()]
        els.extend([el for el in Bs.keys()])
        els.extend([el for el in X.keys()])
        self.ellist = els
        self.store_dir = store_dir
        if os.path.exists(self.store_dir) == False:
            os.mkdir(self.store_dir)
        self.allowed_menus = [3,6,7]
        self.n = len(As.keys()) + len(Bs.keys())
        self.pwd = os.getcwd()

        oxi_states = []
        for k,v in As.items():
            oxi_states.append(v)
        for k,v in Bs.items():
            oxi_states.append(v)
        for k,v in X.items():
            oxi_states.append(-1*v)
            
        self.oxilist = oxi_states
    
        if len(X.keys()) > 1:
            raise Exception('Only single anion supported')
        
        if len(As.keys()) == 2 and len(Bs.keys()) == 2:
            self.menu = 7
            self.allowed_tilts =[3,10,14,20,21,22,23]
            els = [el for el in As.keys()]
            els.extend([el for el in Bs.keys()])
            els.extend([el for el in X.keys()])
            form = "".join(els)
            self.formula = form+'6'
            
        elif len(As.keys()) == 1 and len(Bs.keys()) == 2:
            self.menu = 6
            self.allowed_tilts =[3,10,14,20,21,22,23]
            els = [el for el in As.keys()]
            els[0] = els[0]+'2'
            els.extend([el for el in Bs.keys()])
            els.extend([el for el in X.keys()])
            form = "".join(els)
            self.formula = form+'6' 
            
        elif len(As.keys()) == 2 and len(Bs.keys()) == 1:
            self.menu = 3
            self.allowed_tilts =[10,17,20,21,22,23]            
            els = [el for el in As.keys()]            
            els.extend([el for el in Bs.keys()])
            els[2] = els[2]+'2'
            els.extend([el for el in X.keys()])
            form = "".join(els)
            self.formula = form+'6'
            
        elif len(As.keys()) == 1 and len(Bs.keys()) == 1:
            self.menu = 1
            self.allowed_tilts =[3,5,10,13,14,16,17,20,21,22,23]            
            els = [el for el in As.keys()]            
            els.extend([el for el in Bs.keys()])
            els.extend([el for el in X.keys()])
            form = "".join(els)
            self.formula = form+'3'            
            
            
        else:
            raise Exception('Model not supported')
        
    
    def rocksalt_B_order(self):
        '''
        Double perovskites 1:1 M-site
        No J-T distortions
        Temp only @ 298K (dR/dT = 0.0)
        '''
        self.allowed_tilts =[3,10,14,20,21,22,23]
        self.JT = 0
        self.JT_params = [-0.1,-0.04,+0.18]
        self.T = 298
        self.dRdT = [0.000] * self.n
        
    def gen_default_input(self,tilt):
        '''
        
        No J-T distortions
        Temp only @ 298K (dR/dT = 0.0)
        
        Inputs: tilt -- glazer tilt number (int)
        
        '''
        self.tilt = tilt
        #SPuDS.rocksalt_B_order(self)
        
        self.JT = 0
        self.JT_params = ['-0.1','-0.04','+0.18']
        self.T = 298
        self.dRdT = [0] * self.n
        #if tilt not in self.allowed_tilts:
        #    raise Exception('Tilt not in allowed tilts')
        writelines = []
        writelines.append(self.menu)
        writelines.append(tilt)
        writelines.append(self.JT)
        writelines.append(0)
        for j in self.JT_params:
            writelines.append(j)
        writelines.append(self.formula)
        for oxi in self.oxilist:
            if oxi < 0:
                oxi = oxi*-1
            writelines.append(oxi)
            
        writelines.append(self.T)
        for drdt in self.dRdT:
            writelines.append(drdt)
        
        with open('input.txt','w') as f:
            for line in writelines:
                f.write(str(line)+"\n")

    def gen_JT_input(self,tilt):
        '''
        
        No J-T distortions
        Temp only @ 298K (dR/dT = 0.0)
        
        Inputs: tilt -- glazer tilt number (int)
        
        '''
        self.tilt = tilt
        #SPuDS.rocksalt_B_order(self)
        
        self.JT = 1
        self.JT_params = ['-0.1','-0.04','+0.18']
        self.T = 298
        self.dRdT = [0] * self.n
        #if tilt not in self.allowed_tilts:
        #    raise Exception('Tilt not in allowed tilts')
        writelines = []
        writelines.append(self.menu)
        writelines.append(tilt)
        writelines.append(self.JT)
        writelines.append(0)
        for j in self.JT_params:
            writelines.append(j)
        writelines.append(self.formula)
        for oxi in self.oxilist:
            if oxi < 0:
                oxi = oxi*-1
            writelines.append(oxi)
            
        writelines.append(self.T)
        for drdt in self.dRdT:
            writelines.append(drdt)
        
        with open('input.txt','w') as f:
            for line in writelines:
                f.write(str(line)+"\n")                
                
    def run(self):
        '''
        Runs the input.txt file
        '''
        if self.menu == 6 or self.menu == 7 and self.tilt == 23:
            inputs = "\n".join(["0","0","1"]) #optimize GII
        else:
            inputs = "\n".join(["0","0"])
            
        result = subprocess.run('spuds', input=inputs, capture_output=True, text=True)
        print(result.stdout)
        
    def make_cif(self):
        '''
        Generates a .cif file from the SPuDS output.txt
        
        '''
        cif_file = str(self.tilt)+'_'+self.formula+'.cif'
        self.cif_file = cif_file
        
        def get_lat(tilt,a_len,b_len,c_len,alp,bet,gam):
            if tilt == 3:
                a = a_len
                b = a_len
                c = a_len
                alpha = 90
                beta = 90
                gamma = 90
            if tilt == 5:
                a = a_len
                b = a_len
                c = c_len
                alpha = 90
                beta = 90
                gamma = 90                
            if tilt == 10:
                a = a_len
                b = b_len
                c = c_len
                alpha = 90
                beta = 90
                gamma = 90
            if tilt == 17:
                a = a_len
                b = b_len
                c = c_len
                alpha = 90
                beta = 90
                gamma = 90                
            if tilt == 20:
                a = a_len
                b = b_len
                c = c_len
                alpha = 90
                beta = 90
                gamma = 90
            if tilt == 21:
                a = a_len
                b = a_len
                c = c_len
                alpha = 90
                beta = 90
                gamma = 90
            if tilt == 22:
                a = a_len
                b = a_len
                c = c_len
                alpha = 90
                beta = 90
                gamma = 90
            if tilt == 23:
                a = a_len
                b = a_len
                c = a_len
                alpha = 90
                beta = 90
                gamma = 90
            return a,b,c,alpha,beta,gamma
        
        
        def get_lat_67(tilt,a_len,b_len,c_len,alp,bet,gam):
            if tilt == 3:
                a = a_len
                b = a_len
                c = a_len
                alpha = 90
                beta = 90
                gamma = 90
            if tilt == 10:
                a = a_len
                b = b_len
                c = c_len
                alpha = bet
                beta = bet
                gamma = bet
            if tilt == 14:
                a = a_len
                b = a_len
                c = a_len
                alpha = alp
                beta = alp
                gamma = alp
            if tilt == 20:
                a = a_len
                b = b_len
                c = c_len
                alpha = bet
                beta = bet
                gamma = bet
            if tilt == 21:
                a = a_len
                b = a_len
                c = c_len
                alpha = 90
                beta = 90
                gamma = 90
            if tilt == 22:
                a = a_len
                b = a_len
                c = c_len
                alpha = 90
                beta = 90
                gamma = 90
            if tilt == 23:
                a = a_len
                b = a_len
                c = a_len
                alpha = 90
                beta = 90
                gamma = 90
            return a,b,c,alpha,beta,gamma
    
        def get_lat_1(tilt,a_len,b_len,c_len,alp,bet,gam):
            if tilt == 3:
                a = a_len
                b = a_len
                c = a_len
                alpha = 90
                beta = 90
                gamma = 90
            if tilt == 5:
                a = a_len
                b = a_len
                c = c_len
                alpha = 90
                beta = 90
                gamma = 90
            if tilt == 10:
                a = a_len
                b = b_len
                c = c_len
                alpha = 90
                beta = 90
                gamma = 90  
            if tilt == 13:
                a = a_len
                b = b_len
                c = c_len
                alpha = 90
                beta = bet
                gamma = 90        
            if tilt == 14:
                a = a_len
                b = a_len
                c = c_len
                alpha = 90
                beta = 90
                gamma = 120
            if tilt == 16:
                a = a_len
                b = a_len
                c = c_len
                alpha = 90
                beta = 90
                gamma = 90 
            if tilt == 17:
                a = a_len
                b = b_len
                c = c_len
                alpha = 90
                beta = 90
                gamma = 90         
            if tilt == 20:
                a = a_len
                b = b_len
                c = c_len
                alpha = 90
                beta = 90
                gamma = 90 
            if tilt == 21:
                a = a_len
                b = a_len
                c = c_len
                alpha = 90
                beta = 90
                gamma = 90
            if tilt == 22:
                a = a_len
                b = a_len
                c = c_len
                alpha = 90
                beta = 90
                gamma = 90
            if tilt == 23:
                a = a_len
                b = a_len
                c = a_len
                alpha = 90
                beta = 90
                gamma = 90
            return a,b,c,alpha,beta,gamma
         
    
        a_lat,b_lat,c_lat,alp,bet,gam = None,None,None,None,None,None
        struct_data = []
        #struct_data2 = []
        with open('output.txt','r') as f:
            linecount = 0
            count=0
            for line in f:
                linelist = line.split()
                if linelist == ['**********************************************************************']:
                    count+=1
                        
                if count == 2:
                    if len(linelist) == 4:
                        if linelist[0] == 'Space' and linelist[1] == 'group':
                            key = " ".join([linelist[2],linelist[3]])
                    if len(linelist) == 3:
                        if linelist[0] == 'a':
                            a_lat = linelist[2]
                        if linelist[0] == 'b':
                            b_lat = linelist[2]
                        if linelist[0] == 'c':
                            c_lat = linelist[2]
                            
                        if linelist[0] == 'alpha':
                            alp = linelist[2]
                        if linelist[0] == 'beta':
                            bet = linelist[2]
                        if linelist[0] == 'gamma':
                            gam = linelist[2]                    
                    if len(linelist) > 3:
                        for el in self.ellist:
                            if linelist[0] == el:
                                
                                #if el == 'O':
                                if el == self.Xel:
                                    occup = 1.00
                                else:
                                    occup = linelist[5]
                                    
                                multsite = re.split('(\d+)',linelist[1])
                                mult = multsite[1]
                                wycoff = multsite[2]
                                x = linelist[2]
                                y = linelist[3]
                                z = linelist[4]
                                struct_data.append([el,mult,wycoff,x,y,z,occup])
                            
                            
                linecount+=1
                
                if self.menu == 6 or self.menu == 7:
                    a,b,c,alpha,beta,gamma = get_lat_67(self.tilt,a_lat,b_lat,c_lat,alp,bet,gam)
                elif self.menu == 1:
                    a,b,c,alpha,beta,gamma = get_lat_1(self.tilt,a_lat,b_lat,c_lat,alp,bet,gam)
                else:
                    a,b,c,alpha,beta,gamma = get_lat(self.tilt,a_lat,b_lat,c_lat,alp,bet,gam)
            
        header = 'SPuDS'
        
        loops = [['_space_group_symop_operation_xyz'],['_atom_type_symbol','_atom_type_oxidation_number'],
                 ['_atom_site_label','_atom_site_type_symbol','_atom_site_symmetry_multiplicity',
                  '_atom_site_Wycoff_symbol','_atom_site_fract_x','_atom_site_fract_y','_atom_site_fract_z',
                  '_atom_site_occupancy']]
        
        
        space_group_name_HM_and_number = {'Pn-3 (#201)': {'name': 'P n -3', 'number': 201},
                                           'P2(1)/n (#14)': {'name': 'P 21/n','number': 14},
                                           'R-3 (#148)': {'name': 'R -3','number': 148},
                                           'I2/m (#12)': {'name': 'I 2/m','number': 12},
                                           'P4/mnc (#128)': {'name': 'P 4/mnc','number': 128},
                                           'I4/m (#87)': {'name': 'I 4/m','number': 87},
                                           'Fm-3m (#225)': {'name': 'F m -3 m','number': 225},
                                           'P4(2)/nmc #137':{'name': 'P 42/n m c','number': 137},
                                           'Im-3 (#204)': {'name': 'I m -3', 'number': 204},
                                           'Pnma (#62)': {'name': 'P n m a', 'number': 62},
                                           'Pbnm (#62)': {'name': 'P b n m', 'number': 62},
                                           'Cmcm (#63)': {'name': 'C m c m', 'number': 63},
                                           'Imma (#74)': {'name': 'I m m a', 'number': 74},
                                           'P4/mbm (#127)': {'name': 'P 4/m b m', 'number': 127},
                                           'I4/mcm (#140)': {'name': 'I 4/m c m', 'number': 140},
                                           'Pm-3m (#221)': {'name': 'P m -3 m', 'number': 221},
                                            'P4(2)/nmc (#137)':{'name':'P 42/n m c','number':137}, #not right?
                                            'C2/c (#15)':{'name':'C 2/c','number':15},
                                            'R-3c (#167)':{'name':'R -3 c','number':167},
                                            'I4/mmm (#139)':{'name':'I 4/m m m','number':139}}
        space_group_sym_cell = {'Im-3 (#204)':'cubic',
                                'P4(2)/nmc #137':'tetragonal',
                                'Pnma (#62)':'orthorhombic',
                                'C2/c (#15)':'monoclinic',
                                'R-3c (#167)':'trigonal',
                                'I4/mmm (#139)':'tetragonal',
                                'Cmcm (#63)':'orthorhombic',
                                'Imma (#74)':'orthorhombic',
                                'P4/mbm (#127)':'tetragonal',
                                'I4/mcm (#140)':'tetragonal',
                                'Pm-3m (#221)':'cubic'
                                }                
                                            
        space_group_equiv_pos = {'Im-3 (#204)':['x,y,z','z,x,y','y,z,x','-y,-z,x','z,-x,-y',
                                                '-y,z,-x','-z,-x,y','-z,x,-y','y,-z,-x','-x,-y,z',
                                                'x,-y,-z','-x,y,-z','-x,-y,-z','-z,-x,-y','-y,-z,-x',
                                                'y,z,-x','-z,x,y','y,-z,x','z,x,-y','z,-x,y','-y,z,x',
                                                'x,y,-z','-x,y,z','x,-y,z','1/2+x,1/2+y,1/2+z','1/2+z,1/2+x,1/2+y',
                                                '1/2+y,1/2+z,1/2+x','1/2-y,1/2-z,1/2+x','1/2+z,1/2-x,1/2-y',
                                                '1/2-y,1/2+z,1/2-x','1/2-z,1/2-x,1/2+y','1/2-z,1/2+x,1/2-y',
                                                '1/2+y,1/2-z,1/2-x','1/2-x,1/2-y,1/2+z','1/2+x,1/2-y,1/2-z',
                                                '1/2-x,1/2+y,1/2-z','1/2-x,1/2-y,1/2-z','1/2-z,1/2-x,1/2-y',
                                                '1/2-y,1/2-z,1/2-x','1/2+y,1/2+z,1/2-x','1/2-z,1/2+x,1/2+y',
                                                '1/2+y,1/2-z,1/2+x','1/2+z,1/2+x,1/2-y','1/2+z,1/2-x,1/2+y',
                                                '1/2-y,1/2+z,1/2+x','1/2+x,1/2+y,1/2-z','1/2-x,1/2+y,1/2+z',
                                                '1/2+x,1/2-y,1/2+z'],
                                'P4(2)/nmc (#137)':['x,y,z','1/2-y,x,1/2+z','1/2-x,1/2-y,z','y,1/2-x,1/2+z',
                                                    '1/2+x,-y,-z','-x,1/2+y,-z','1/2+y,1/2+x,1/2-z','-y,-x,1/2-z',
                                                    '-x,-y,-z','1/2+y,-x,1/2-z','1/2+x,1/2+y,-z','-y,1/2+x,1/2-z',
                                                    '1/2-x,y,z','x,1/2-y,z','1/2-y,1/2-x,1/2+z','y,x,1/2+z'] ,
                                'I4/mmm (#139)':['x,y,z','-y,x,z','-x,-y,z','y,-x,z','x,-y,-z','-x,y,-z',
                                                'y,x,-z','-y,-x,-z','-x,-y,-z','y,-x,-z','x,y,-z','-y,x,-z',
                                                '-x,y,z','x,-y,z','-y,-x,z','y,x,z','1/2+x,1/2+y,1/2+z',
                                                '1/2-y,1/2+x,1/2+z','1/2-x,1/2-y,1/2+z','1/2+y,1/2-x,1/2+z',
                                                '1/2+x,1/2-y,1/2-z','1/2-x,1/2+y,1/2-z','1/2+y,1/2+x,1/2-z',
                                                '1/2-y,1/2-x,1/2-z','1/2-x,1/2-y,1/2-z','1/2+y,1/2-x,1/2-z',
                                                '1/2+x,1/2+y,1/2-z','1/2-y,1/2+x,1/2-z','1/2-x,1/2+y,1/2+z',
                                                '1/2+x,1/2-y,1/2+z','1/2-y,1/2-x,1/2+z','1/2+y,1/2+x,1/2+z'],
                                'Cmcm (#63)':['x,y,z','-x,-y,1/2+z','x,-y,-z','-x,y,1/2-z','-x,-y,-z',
                                            'x,y,1/2-z','-x,y,z','x,-y,1/2+z','1/2+x,1/2+y,z','1/2-x,1/2-y,1/2+z',
                                            '1/2+x,1/2-y,-z','1/2-x,1/2+y,1/2-z','1/2-x,1/2-y,-z','1/2+x,1/2+y,1/2-z',
                                            '1/2-x,1/2+y,z','1/2+x,1/2-y,1/2+z']
                                                               
                                }                                     
        
        
        space_group_symops = {'Pn-3 (#201)':['x,y,z','-x+1/2,-y+1/2,z', '-x+1/2,y,-z+1/2', 'x,-y+1/2,-z+1/2',
                                             'z,x,y', 'z,-x+1/2,-y+1/2', '-z+1/2,-x+1/2,y',	'-z+1/2,x,-y+1/2',
                                             'y,z,x', '-y+1/2,z,-x+1/2', 'y,-z+1/2,-x+1/2', '-y+1/2,-z+1/2,x',
                                             '-x,-y,-z', 'x+1/2,y+1/2,-z', 'x+1/2,-y,z+1/2', '-x,y+1/2,z+1/2',
                                             '-z,-x,-y', '-z,x+1/2,y+1/2', 'z+1/2,x+1/2,-y', 'z+1/2,-x,y+1/2',
                                             '-y,-z,-x',	'y+1/2,-z,x+1/2', '-y,z+1/2,x+1/2', 'y+1/2,z+1/2,-x'],
                                
                              'P2(1)/n (#14)':['x,y,z', '-x+1/2,y+1/2,-z+1/2', '-x,-y,-z', 'x+1/2,-y+1/2,z+1/2'],  #probably incorrect
                              
                              'R-3 (#148)':['x,y,z', 'z,x,y', 'y,z,x', '-x,-y,-z', '-z,-x,-y', '-y,-z,-x'],
                              
                              'I2/m (#12)':['x, y, z','-x, -y, -z','-x, y, -z','x, -y, z','x+1/2, y+1/2, z+1/2',
                                            '-x+1/2, -y+1/2, -z+1/2','-x+1/2, y+1/2, -z+1/2','x+1/2, -y+1/2, z+1/2'],
                              
                              'P4/mnc (#128)':['x,y,z','-x,-y,z', '-y,x,z', 'y,-x,z', '-x+1/2,y+1/2,-z+1/2', 'x+1/2,-y+1/2,-z+1/2', 	
                                               'y+1/2,x+1/2,-z+1/2', '-y+1/2,-x+1/2,-z+1/2', '-x,-y,-z', 'x,y,-z', 'y,-x,-z', '-y,x,-z',
                                               'x+1/2,-y+1/2,z+1/2', '-x+1/2,y+1/2,z+1/2', '-y+1/2,-x+1/2,z+1/2', 'y+1/2,x+1/2,z+1/2'],
                              
                              'I4/m (#87)':['x, y, z','-x, -y, -z','-x, -y, z','x, y, -z','-y, x, z','y, -x, -z','y, -x, z','-y, x, -z',
                                            'x+1/2, y+1/2, z+1/2','-x+1/2, -y+1/2, -z+1/2','-x+1/2, -y+1/2, z+1/2','x+1/2, y+1/2, -z+1/2',
                                            '-y+1/2, x+1/2, z+1/2','y+1/2, -x+1/2, -z+1/2','y+1/2, -x+1/2, z+1/2','-y+1/2, x+1/2, -z+1/2'],
                              
                              'Fm-3m (#225)':['x, y, z', '-x, -y, -z','-x, -y, z','x, y, -z','-x, y, -z','x, -y, z','x, -y, -z',
                                              '-x, y, z','z, x, y','-z, -x, -y','z, -x, -y','-z, x, y','-z, -x, y','z, x, -y',
                                              '-z, x, -y','z, -x, y','y, z, x','-y, -z, -x','-y, z, -x','y, -z, x','y, -z, -x',
                                              '-y, z, x','-y, -z, x','y, z, -x','y, x, -z','-y, -x, z','-y, -x, -z','y, x, z',
                                              'y, -x, z','-y, x, -z','-y, x, z','y, -x, -z','x, z, -y','-x, -z, y','-x, z, y',
                                              'x, -z, -y','-x, -z, -y','x, z, y','x, -z, y','-x, z, -y','z, y, -x','-z, -y, x',
                                              'z, -y, x','-z, y, -x','-z, y, x','z, -y, -x','-z, -y, -x','z, y, x','x, y+1/2, z+1/2',
                                              '-x, -y+1/2, -z+1/2','-x, -y+1/2, z+1/2','x, y+1/2, -z+1/2','-x, y+1/2, -z+1/2',
                                              'x, -y+1/2, z+1/2','x, -y+1/2, -z+1/2','-x, y+1/2, z+1/2','z, x+1/2, y+1/2','-z, -x+1/2, -y+1/2',
                                              'z, -x+1/2, -y+1/2','-z, x+1/2, y+1/2','-z, -x+1/2, y+1/2','z, x+1/2, -y+1/2','-z, x+1/2, -y+1/2',
                                              'z, -x+1/2, y+1/2','y, z+1/2, x+1/2','-y, -z+1/2, -x+1/2','-y, z+1/2, -x+1/2','y, -z+1/2, x+1/2',
                                              'y, -z+1/2, -x+1/2','-y, z+1/2, x+1/2','-y, -z+1/2, x+1/2','y, z+1/2, -x+1/2','y, x+1/2, -z+1/2',
                                              '-y, -x+1/2, z+1/2','-y, -x+1/2, -z+1/2','y, x+1/2, z+1/2','y, -x+1/2, z+1/2','-y, x+1/2, -z+1/2',
                                              '-y, x+1/2, z+1/2','y, -x+1/2, -z+1/2','x, z+1/2, -y+1/2','-x, -z+1/2, y+1/2','-x, z+1/2, y+1/2',
                                              'x, -z+1/2, -y+1/2','-x, -z+1/2, -y+1/2','x, z+1/2, y+1/2','x, -z+1/2, y+1/2','-x, z+1/2, -y+1/2',
                                              'z, y+1/2, -x+1/2','-z, -y+1/2, x+1/2','z, -y+1/2, x+1/2','-z, y+1/2, -x+1/2','-z, y+1/2, x+1/2',
                                              'z, -y+1/2, -x+1/2','-z, -y+1/2, -x+1/2','z, y+1/2, x+1/2','x+1/2, y, z+1/2','-x+1/2, -y, -z+1/2',
                                              '-x+1/2, -y, z+1/2','x+1/2, y, -z+1/2','-x+1/2, y, -z+1/2','x+1/2, -y, z+1/2','x+1/2, -y, -z+1/2',
                                              '-x+1/2, y, z+1/2','z+1/2, x, y+1/2','-z+1/2, -x, -y+1/2','z+1/2, -x, -y+1/2','-z+1/2, x, y+1/2',
                                              '-z+1/2, -x, y+1/2','z+1/2, x, -y+1/2','-z+1/2, x, -y+1/2','z+1/2, -x, y+1/2','y+1/2, z, x+1/2',
                                              '-y+1/2, -z, -x+1/2','-y+1/2, z, -x+1/2','y+1/2, -z, x+1/2','y+1/2, -z, -x+1/2','-y+1/2, z, x+1/2',
                                              '-y+1/2, -z, x+1/2','y+1/2, z, -x+1/2','y+1/2, x, -z+1/2','-y+1/2, -x, z+1/2','-y+1/2, -x, -z+1/2',
                                              'y+1/2, x, z+1/2','y+1/2, -x, z+1/2','-y+1/2, x, -z+1/2','-y+1/2, x, z+1/2','y+1/2, -x, -z+1/2',
                                              'x+1/2, z, -y+1/2','-x+1/2, -z, y+1/2','-x+1/2, z, y+1/2','x+1/2, -z, -y+1/2','-x+1/2, -z, -y+1/2',
                                              'x+1/2, z, y+1/2','x+1/2, -z, y+1/2','-x+1/2, z, -y+1/2','z+1/2, y, -x+1/2','-z+1/2, -y, x+1/2',
                                              'z+1/2, -y, x+1/2','-z+1/2, y, -x+1/2','-z+1/2, y, x+1/2','z+1/2, -y, -x+1/2','-z+1/2, -y, -x+1/2',
                                              'z+1/2, y, x+1/2','x+1/2, y+1/2, z','-x+1/2, -y+1/2, -z','-x+1/2, -y+1/2, z','x+1/2, y+1/2, -z',
                                              '-x+1/2, y+1/2, -z','x+1/2, -y+1/2, z','x+1/2, -y+1/2, -z','-x+1/2, y+1/2, z','z+1/2, x+1/2, y',
                                              '-z+1/2, -x+1/2, -y','z+1/2, -x+1/2, -y','-z+1/2, x+1/2, y','-z+1/2, -x+1/2, y','z+1/2, x+1/2, -y',
                                              '-z+1/2, x+1/2, -y','z+1/2, -x+1/2, y','y+1/2, z+1/2, x','-y+1/2, -z+1/2, -x','-y+1/2, z+1/2, -x',
                                              'y+1/2, -z+1/2, x','y+1/2, -z+1/2, -x','-y+1/2, z+1/2, x','-y+1/2, -z+1/2, x','y+1/2, z+1/2, -x',
                                              'y+1/2, x+1/2, -z','-y+1/2, -x+1/2, z','-y+1/2, -x+1/2, -z','y+1/2, x+1/2, z','y+1/2, -x+1/2, z',
                                              '-y+1/2, x+1/2, -z','-y+1/2, x+1/2, z','y+1/2, -x+1/2, -z','x+1/2, z+1/2, -y','-x+1/2, -z+1/2, y',
                                              '-x+1/2, z+1/2, y','x+1/2, -z+1/2, -y','-x+1/2, -z+1/2, -y','x+1/2, z+1/2, y','x+1/2, -z+1/2, y',
                                              '-x+1/2, z+1/2, -y','z+1/2, y+1/2, -x','-z+1/2, -y+1/2, x','z+1/2, -y+1/2, x','-z+1/2, y+1/2, -x',
                                              '-z+1/2, y+1/2, x','z+1/2, -y+1/2, -x','-z+1/2, -y+1/2, -x','z+1/2, y+1/2, x'],
    
                                #'P4(2)/nmc #137':['x,y,z', '-x,-y,z', '-y+1/2,x+1/2,z+1/2', 'y+1/2,-x+1/2,z+1/2', '-x+1/2,y+1/2,-z+1/2', 
                                #                  'x+1/2,-y+1/2,-z+1/2', 'y,x,-z', '-y,-x,-z', '-x+1/2,-y+1/2,-z+1/2', 'x+1/2,y+1/2,-z+1/2', 
                                #                  'y,-x,-z', '-y,x,-z', 'x,-y,z', '-x,y,z', '-y+1/2,-x+1/2,z+1/2', 'y+1/2,x+1/2,z+1/2'],
        
        
                                #'Im-3 (#204)':['x,y,z', '-x,-y,z', '-x,y,-z', 'x,-y,-z', 'z,x,y', 'z,-x,-y', '-z,-x,y', '-z,x,-y',
                                #                  'y,z,x', '-y,z,-x', 'y,-z,-x', '-y,-z,x','-x,-y,-z', 'x,y,-z', 'x,-y,z' ,'-x,y,z',
                                #                  '-z,-x,-y', '-z,x,y', 'z,x,-y', 'z,-x,y','-y,-z,-x', 'y,-z,x', '-y,z,x', 'y,z,-x'],
        
                        
                                'Pnma (#62)':['x,y,z', '-x+1/2,-y,z+1/2', '-x,y+1/2,-z',	'x+1/2,-y+1/2,-z+1/2', '-x,-y,-z', 
                                              'x+1/2,y,-z+1/2',	'x,-y+1/2,z','-x+1/2,y+1/2,z+1/2'],
                                'Pbnm (#62)':['x,y,z', 'x+1/2,-y+1/2,-z', '-x,-y,z+1/2',	'-x+1/2,y+1/2,-z+1/2', '-x,-y,-z', '-x+1/2,y+1/2,z',	
                                              'x,y,-z+1/2','x+1/2,-y+1/2,z+1/2'],
        
                                #'Cmcm (#63)':['x, y, z', '-x, -y, -z','-x, -y, z+1/2','x, y, -z+1/2','-x, y, -z+1/2','x, -y, z+1/2','x, -y, -z','-x, y, z',
                                #              'x+1/2, y+1/2, z','-x+1/2, -y+1/2, -z','-x+1/2, -y+1/2, z+1/2','x+1/2, y+1/2, -z+1/2','-x+1/2, y+1/2, -z+1/2',
                                #              'x+1/2, -y+1/2, z+1/2','x+1/2, -y+1/2, -z','-x+1/2, y+1/2, z'],
                                'Cmcm (#63)':['x,y,z','-x,-y,1/2+z','x,-y,-z','-x,y,1/2-z','-x,-y,-z',
                                            'x,y,1/2-z','-x,y,z','x,-y,1/2+z','1/2+x,1/2+y,z','1/2-x,1/2-y,1/2+z',
                                            '1/2+x,1/2-y,-z','1/2-x,1/2+y,1/2-z','1/2-x,1/2-y,-z','1/2+x,1/2+y,1/2-z',
                                            '1/2-x,1/2+y,z','1/2+x,1/2-y,1/2+z'],                                
    
                                'Imma (#74)':['x, y, z','-x, -y, -z','-x, -y+1/2, z','x, y+1/2, -z','-x, y+1/2, -z','x, -y+1/2, z','x, -y, -z',
                                              '-x, y, z','x+1/2, y+1/2, z+1/2','-x+1/2, -y+1/2, -z+1/2','-x+1/2, -y, z+1/2','x+1/2, y, -z+1/2',
                                              '-x+1/2, y, -z+1/2','x+1/2, -y, z+1/2','x+1/2, -y+1/2, -z+1/2','-x+1/2, y+1/2, z+1/2'],
                                'P4/mbm (#127)':['x,y,z', '-x,-y,z', '-y,x,z', 'y,-x,z', '-x+1/2,y+1/2,-z', 'x+1/2,-y+1/2,-z', 'y+1/2,x+1/2,-z', 
                                                 '-y+1/2,-x+1/2,-z', '-x,-y,-z', 'x,y,-z', 'y,-x,-z', '-y,x,-z','x+1/2,-y+1/2,z', '-x+1/2,y+1/2,z', 
                                                 '-y+1/2,-x+1/2,z', 'y+1/2,x+1/2,z'],
                                'I4/mcm (#140)':['x, y, z','-x, -y, -z','-x, -y, z','x, y, -z','-y, x, z','y, -x, -z','y, -x, z','-y, x, -z','-x, y, -z+1/2',
                                                 'x, -y, z+1/2','x, -y, -z+1/2','-x, y, z+1/2','y, x, -z+1/2','-y, -x, z+1/2','-y, -x, -z+1/2','y, x, z+1/2',
                                                 'x+1/2, y+1/2, z+1/2','-x+1/2, -y+1/2, -z+1/2','-x+1/2, -y+1/2, z+1/2','x+1/2, y+1/2, -z+1/2','-y+1/2, x+1/2, z+1/2',
                                                 'y+1/2, -x+1/2, -z+1/2','y+1/2, -x+1/2, z+1/2','-y+1/2, x+1/2, -z+1/2','-x+1/2, y+1/2, -z','x+1/2, -y+1/2, z',
                                                 'x+1/2, -y+1/2, -z','-x+1/2, y+1/2, z','y+1/2, x+1/2, -z','-y+1/2, -x+1/2, z','-y+1/2, -x+1/2, -z','y+1/2, x+1/2, z'],
                                'Pm-3m (#221)':['x,y,z', '-x,-y,z', '-x,y,-z', 'x,-y,-z','z,x,y', 'z,-x,-y', '-z,-x,y', '-z,x,-y', 'y,z,x', 
                                                '-y,z,-x', 'y,-z,-x', '-y,-z,x', 'y,x,-z', '-y,-x,-z', 'y,-x,z', '-y,x,z', 'x,z,-y', '-x,z,y', '-x,-z,-y', 
                                                'x,-z,y', 'z,y,-x',	'z,-y,x', '-z,y,x', '-z,-y,-x', '-x,-y,-z', 'x,y,-z', 'x,-y,z',	'-x,y,z',
                                                '-z,-x,-y',	'-z,x,y', 'z,x,-y', 'z,-x,y', '-y,-z,-x', 'y,-z,x', '-y,z,x', 'y,z,-x', '-y,-x,z', 
                                                'y,x,z', '-y,x,-z', 'y,-x,-z', '-x,-z,y', 'x,-z,-y', 'x,z,y', '-x,z,-y', '-z,-y,x', '-z,y,-x', 'z,-y,-x', 'z,y,x'],
                                #'Im-3 (#204)':['x, y, z','-x, -y, -z','-x, -y, z','x, y, -z','-x, y, -z',
                                #                'x, -y, z','x, -y, -z','-x, y, z','z, x, y','-z, -x, -y','z, -x, -y',
                                #                '-z, x, y','-z, -x, y','z, x, -y','-z, x, -y','z, -x, y','y, z, x',
                                #                '-y, -z, -x','-y, z, -x','y, -z, x','y, -z, -x','-y, z, x','-y, -z, x',
                                #                'y, z, -x','x+1/2, y+1/2, z+1/2','-x+1/2, -y+1/2, -z+1/2','-x+1/2, -y+1/2, z+1/2',
                                #                'x+1/2, y+1/2, -z+1/2','-x+1/2, y+1/2, -z+1/2','x+1/2, -y+1/2, z+1/2','x+1/2, -y+1/2, -z+1/2',
                                #                '-x+1/2, y+1/2, z+1/2','z+1/2, x+1/2, y+1/2','-z+1/2, -x+1/2, -y+1/2','z+1/2, -x+1/2, -y+1/2',
                                #                '-z+1/2, x+1/2, y+1/2','-z+1/2, -x+1/2, y+1/2','z+1/2, x+1/2, -y+1/2','-z+1/2, x+1/2, -y+1/2',
                                #                'z+1/2, -x+1/2, y+1/2','y+1/2, z+1/2, x+1/2','-y+1/2, -z+1/2, -x+1/2','-y+1/2, z+1/2, -x+1/2',
                                #                'y+1/2, -z+1/2, x+1/2','y+1/2, -z+1/2, -x+1/2','-y+1/2, z+1/2, x+1/2','-y+1/2, -z+1/2, x+1/2',
                                #                'y+1/2, z+1/2, -x+1/2'],
                                'Im-3 (#204)':['x,y,z','z,x,y','y,z,x','-y,-z,x','z,-x,-y',
                                                '-y,z,-x','-z,-x,y','-z,x,-y','y,-z,-x','-x,-y,z',
                                                'x,-y,-z','-x,y,-z','-x,-y,-z','-z,-x,-y','-y,-z,-x',
                                                'y,z,-x','-z,x,y','y,-z,x','z,x,-y','z,-x,y','-y,z,x',
                                                'x,y,-z','-x,y,z','x,-y,z','1/2+x,1/2+y,1/2+z','1/2+z,1/2+x,1/2+y',
                                                '1/2+y,1/2+z,1/2+x','1/2-y,1/2-z,1/2+x','1/2+z,1/2-x,1/2-y',
                                                '1/2-y,1/2+z,1/2-x','1/2-z,1/2-x,1/2+y','1/2-z,1/2+x,1/2-y',
                                                '1/2+y,1/2-z,1/2-x','1/2-x,1/2-y,1/2+z','1/2+x,1/2-y,1/2-z',
                                                '1/2-x,1/2+y,1/2-z','1/2-x,1/2-y,1/2-z','1/2-z,1/2-x,1/2-y',
                                                '1/2-y,1/2-z,1/2-x','1/2+y,1/2+z,1/2-x','1/2-z,1/2+x,1/2+y',
                                                '1/2+y,1/2-z,1/2+x','1/2+z,1/2+x,1/2-y','1/2+z,1/2-x,1/2+y',
                                                '1/2-y,1/2+z,1/2+x','1/2+x,1/2+y,1/2-z','1/2-x,1/2+y,1/2+z',
                                                '1/2+x,1/2-y,1/2+z'],   
                                'P4(2)/nmc #137':['x,y,z','1/2-y,x,1/2+z','1/2-x,1/2-y,z','y,1/2-x,1/2+z',
                                                    '1/2+x,-y,-z','-x,1/2+y,-z','1/2+y,1/2+x,1/2-z','-y,-x,1/2-z',
                                                    '-x,-y,-z','1/2+y,-x,1/2-z','1/2+x,1/2+y,-z','-y,1/2+x,1/2-z',
                                                    '1/2-x,y,z','x,1/2-y,z','1/2-y,1/2-x,1/2+z','y,x,1/2+z'],                             
                                                
                                #'P4(2)/nmc (#137)':['x, y, z','-x, -y, -z','-x+1/2, -y+1/2, z','x+1/2, y+1/2, -z',
                                #                    '-y+1/2, x, z+1/2','y+1/2, -x, -z+1/2','y, -x+1/2, z+1/2',
                                #                    '-y, x+1/2, -z+1/2','-x, y+1/2, -z','x, -y+1/2, z',
                                #                    'x+1/2, -y, -z','-x+1/2, y, z','y+1/2, x+1/2, -z+1/2',
                                #                    '-y+1/2, -x+1/2, z+1/2','-y, -x, -z+1/2','y, x, z+1/2'],
                            
                                #'Pnma (#62)':['x, y, z','-x, -y, -z','-x+1/2, -y, z+1/2','x+1/2, y, -z+1/2',
                                #'-x, y+1/2, -z','x, -y+1/2, z','x+1/2, -y+1/2, -z+1/2','-x+1/2, y+1/2, z+1/2'],
                            
                                'C2/c (#15)':['x, y, z','-x, -y, -z','-x, y, -z+1/2','x, -y, z+1/2','x+1/2, y+1/2, z',
                                '-x+1/2, -y+1/2, -z','-x+1/2, y+1/2, -z+1/2','x+1/2, -y+1/2, z+1/2'],
                            
                                'R-3c (#167)':['x, y, z','-x, -y, -z','-y, x-y, z','y, -x+y, -z','-x+y, -x, z','x-y, x, -z',
                                'y, x, -z+1/2','-y, -x, z+1/2','x-y, -y, -z+1/2','-x+y, y, z+1/2','-x, -x+y, -z+1/2',
                                'x, x-y, z+1/2','x+2/3, y+1/3, z+1/3','-x+2/3, -y+1/3, -z+1/3','-y+2/3, x-y+1/3, z+1/3',
                                'y+2/3, -x+y+1/3, -z+1/3','-x+y+2/3, -x+1/3, z+1/3','x-y+2/3, x+1/3, -z+1/3','y+2/3, x+1/3, -z+5/6',
                                '-y+2/3, -x+1/3, z+5/6','x-y+2/3, -y+1/3, -z+5/6','-x+y+2/3, y+1/3, z+5/6','-x+2/3, -x+y+1/3, -z+5/6',
                                'x+2/3, x-y+1/3, z+5/6','x+1/3, y+2/3, z+2/3','-x+1/3, -y+2/3, -z+2/3','-y+1/3, x-y+2/3, z+2/3',
                                'y+1/3, -x+y+2/3, -z+2/3','-x+y+1/3, -x+2/3, z+2/3','x-y+1/3, x+2/3, -z+2/3','y+1/3, x+2/3, -z+1/6',
                                '-y+1/3, -x+2/3, z+1/6','x-y+1/3, -y+2/3, -z+1/6','-x+y+1/3, y+2/3, z+1/6','-x+1/3, -x+y+2/3, -z+1/6',
                                'x+1/3, x-y+2/3, z+1/6'],
                            
                                #'I4/mmm (#139)':['x, y, z','-x, -y, -z','-x, -y, z','x, y, -z','-y, x, z','y, -x, -z',
                                #'y, -x, z','-y, x, -z','-x, y, -z','x, -y, z','x, -y, -z','-x, y, z','y, x, -z',
                                #'-y, -x, z','-y, -x, -z','y, x, z','x+1/2, y+1/2, z+1/2','-x+1/2, -y+1/2, -z+1/2','-x+1/2, -y+1/2, z+1/2',
                                #'x+1/2, y+1/2, -z+1/2','-y+1/2, x+1/2, z+1/2','y+1/2, -x+1/2, -z+1/2','y+1/2, -x+1/2, z+1/2','-y+1/2, x+1/2, -z+1/2',
                                #'-x+1/2, y+1/2, -z+1/2','x+1/2, -y+1/2, z+1/2','x+1/2, -y+1/2, -z+1/2','-x+1/2, y+1/2, z+1/2',
                                #'y+1/2, x+1/2, -z+1/2','-y+1/2, -x+1/2, z+1/2','-y+1/2, -x+1/2, -z+1/2','y+1/2, x+1/2, z+1/2']
                                'I4/mmm (#139)':['x,y,z','-y,x,z','-x,-y,z','y,-x,z','x,-y,-z','-x,y,-z',
                                                'y,x,-z','-y,-x,-z','-x,-y,-z','y,-x,-z','x,y,-z','-y,x,-z',
                                                '-x,y,z','x,-y,z','-y,-x,z','y,x,z','1/2+x,1/2+y,1/2+z',
                                                '1/2-y,1/2+x,1/2+z','1/2-x,1/2-y,1/2+z','1/2+y,1/2-x,1/2+z',
                                                '1/2+x,1/2-y,1/2-z','1/2-x,1/2+y,1/2-z','1/2+y,1/2+x,1/2-z',
                                                '1/2-y,1/2-x,1/2-z','1/2-x,1/2-y,1/2-z','1/2+y,1/2-x,1/2-z',
                                                '1/2+x,1/2+y,1/2-z','1/2-y,1/2+x,1/2-z','1/2-x,1/2+y,1/2+z',
                                                '1/2+x,1/2-y,1/2+z','1/2-y,1/2-x,1/2+z','1/2+y,1/2+x,1/2+z']}
                            
                                #'Cmcm (#63)':['x, y, z','-x, -y, -z','-x, -y, z+1/2','x, y, -z+1/2','-x, y, -z+1/2',
                                #'x, -y, z+1/2','x, -y, -z','-x, y, z','x+1/2, y+1/2, z','-x+1/2, -y+1/2, -z',
                                #'-x+1/2, -y+1/2, z+1/2','x+1/2, y+1/2, -z+1/2','-x+1/2, y+1/2, -z+1/2','x+1/2, -y+1/2, z+1/2',
                                #'x+1/2, -y+1/2, -z','-x+1/2, y+1/2, z'],
                            
                                #'Imma (#74)':['x, y, z','-x, -y, -z','-x, -y+1/2, z','x, y+1/2, -z','-x, y+1/2, -z',
                                #'x, -y+1/2, z','x, -y, -z','-x, y, z','x+1/2, y+1/2, z+1/2','-x+1/2, -y+1/2, -z+1/2',
                                #'-x+1/2, -y, z+1/2','x+1/2, y, -z+1/2','-x+1/2, y, -z+1/2','x+1/2, -y, z+1/2',
                                #'x+1/2, -y+1/2, -z+1/2','-x+1/2, y+1/2, z+1/2'],
                            
                                #'P4/mbm (#127)':['x, y, z','-x, -y, -z','-x, -y, z','x, y, -z','-y, x, z',
                                #'y, -x, -z','y, -x, z','-y, x, -z','-x+1/2, y+1/2, -z','x+1/2, -y+1/2, z','x+1/2, -y+1/2, -z',
                                #'-x+1/2, y+1/2, z','y+1/2, x+1/2, -z','-y+1/2, -x+1/2, z','-y+1/2, -x+1/2, -z','y+1/2, x+1/2, z'],
                            
                                #'I4/mcm (#140)':['x, y, z','-x, -y, -z','-x, -y, z','x, y, -z','-y, x, z','y, -x, -z',
                                #'y, -x, z','-y, x, -z','-x, y, -z+1/2','x, -y, z+1/2','x, -y, -z+1/2','-x, y, z+1/2',
                                #'y, x, -z+1/2','-y, -x, z+1/2','-y, -x, -z+1/2','y, x, z+1/2','x+1/2, y+1/2, z+1/2',
                                #'-x+1/2, -y+1/2, -z+1/2','-x+1/2, -y+1/2, z+1/2','x+1/2, y+1/2, -z+1/2','-y+1/2, x+1/2, z+1/2',
                                #'y+1/2, -x+1/2, -z+1/2','y+1/2, -x+1/2, z+1/2','-y+1/2, x+1/2, -z+1/2','-x+1/2, y+1/2, -z',
                                #'x+1/2, -y+1/2, z','x+1/2, -y+1/2, -z','-x+1/2, y+1/2, z','y+1/2, x+1/2, -z','-y+1/2, -x+1/2, z',
                                #'-y+1/2, -x+1/2, -z','y+1/2, x+1/2, z'],
                            
                                #'Pm-3m (#221)':['x, y, z','-x, -y, -z','-x, -y, z','x, y, -z','-x, y, -z',
                                #'x, -y, z','x, -y, -z','-x, y, z','z, x, y','-z, -x, -y','z, -x, -y',
                                #'-z, x, y','-z, -x, y','z, x, -y','-z, x, -y','z, -x, y','y, z, x',
                                #'-y, -z, -x','-y, z, -x','y, -z, x','y, -z, -x','-y, z, x','-y, -z, x',
                                #'y, z, -x','y, x, -z','-y, -x, z','-y, -x, -z','y, x, z','y, -x, z',
                                #'-y, x, -z','-y, x, z','y, -x, -z','x, z, -y','-x, -z, y','-x, z, y',
                                #'x, -z, -y','-x, -z, -y','x, z, y','x, -z, y','-x, z, -y','z, y, -x',
                                #'-z, -y, x','z, -y, x','-z, y, -x','-z, y, x','z, -y, -x','-z, -y, -x',
                                #'z, y, x']}                                   
        
        elcounts = {}
        for d in struct_data:
            for el in self.ellist:
                if d[0] == el:
                    if el not in elcounts.keys():
                        elcounts[el] = 1
                    else:
                        elcounts[el]+=1
        for k,v in elcounts.items():
            count = 1
            ii=0
            for d in struct_data:
                if d[0] == k:
                    if v >1:
                        label = k+str(count)
                        struct_data[ii].insert(0,label)
                        count+=1
                    else:
                        label = k
                        struct_data[ii].insert(0,label)
                ii+=1
                        
        
        data = {}
        data['_cell_length_a'] = a
        data['_cell_length_b'] = b
        data['_cell_length_c'] = c
        data['_cell_angle_alpha'] = alpha
        data['_cell_angle_beta'] = beta
        data['_cell_angle_gamma'] = gamma
        data['_space_group_name_H-M_alt'] = space_group_name_HM_and_number[key]['name']
        data['_symmetry_Int_tables_number'] = space_group_name_HM_and_number[key]['number']
        data['_symmetry_cell_setting'] = space_group_sym_cell[key]
        
        data['_space_group_symop_operation_xyz'] = space_group_symops[key]
        
        data['_atom_type_symbol'] = self.ellist
        data['_atom_type_oxidation_number'] = self.oxilist
        
        
        data['_atom_site_label'] = []
        data['_atom_site_type_symbol'] = []
        data['_atom_site_symmetry_multiplicity'] = []
        data['_atom_site_Wycoff_symbol'] = []
        data['_atom_site_fract_x'] = []
        data['_atom_site_fract_y'] = []
        data['_atom_site_fract_z'] = []
        data['_atom_site_occupancy'] = []
        
        for d in struct_data:
            data['_atom_site_label'].append(d[0]) 
            data['_atom_site_type_symbol'].append(d[1])
            data['_atom_site_symmetry_multiplicity'].append(d[2])
            data['_atom_site_Wycoff_symbol'].append(d[3])
            data['_atom_site_fract_x'].append(d[4])
            data['_atom_site_fract_y'].append(d[5])
            data['_atom_site_fract_z'].append(d[6])
            #data['_atom_site_occupancy'].append(d[7])
            data['_atom_site_occupancy'].append(1)
    
        d = OrderedDict()
        d[self.formula] = CifBlock(data,loops,header)
        cf = CifFile(d)
        with zopen(self.cif_file, "wt") as f:
            f.write(cf.__str__())     
                
    def store_output(self):
        subdir = str(self.formula)
        path = opj(self.store_dir,subdir)
        if os.path.exists(path) == False:
            os.mkdir(path)
        #else:
            #warnings.warn(subdir+' file exists. Existing files will be overwritten.')
            
        move('gii.txt',opj(path,'gii.txt'))
        move('output.txt',opj(path,str(self.tilt)+'_output.txt'))
        move(self.cif_file ,opj(path,self.cif_file))
        
        
    def store_best_poscar(self):
        
        def get_cif_file_name(self):
            tilt = self.tilt
            if self.menu == 1:
                if tilt == 3:
                    sg = 'Im-3'
                if tilt == 5:
                    sg = 'P42_nmc'
                if tilt == 10:
                    sg = 'Pnma'
                if tilt == 13:
                    sg = 'C2_c'
                if tilt == 14:
                    sg = 'R-3c'
                if tilt == 16:
                    sg = 'I4_mmm'
                if tilt == 17:
                    sg = 'Cmcm'
                if tilt == 20:
                    sg = 'Imma'
                if tilt == 21:
                    sg = 'P4_mbm'
                if tilt == 22:
                    sg = 'I4_mcm'
                if tilt == 23:
                    sg = 'Pm-3m'
                self.cif_file = 'output_'+sg+'.cif'
                
        
        def get_min_O_bond(structure):
            dist_matrix = structure.distance_matrix
            O_inds = [i for i,specie in enumerate(structure.species) if specie.symbol == 'O']
            minO_dist = 100
            for i in O_inds:
                for j in O_inds:
                    if j>i:
                        if dist_matrix[i][j] < minO_dist:
                            minO_dist = dist_matrix[i][j]
            return minO_dist
        
        top = opj(self.pwd,self.store_dir)
        first = opj(top,self.formula)
        second = opj(first,str(self.tilt)+'_tilt')
        pos_store = opj(second,'POSCAR')
        if os.path.exists(top) == False:
            os.mkdir(top)
        if os.path.exists(first) == False:
            os.mkdir(first)
        if os.path.exists(second) == False:
            os.mkdir(second)
        #get_cif_file_name(self)
        #print(self.cif_file)
        #atoms = read(self.cif_file)
        structure = CifParser(self.cif_file).get_structures()[0]
        """
        #print('cif',structure)
        try:
            
            dict_structs = OD.apply_transformation(structure,return_ranked_list=1)
            #print(list_structs)
            structure = dict_structs[0]['structure']
        #print('od',structure)
        except:
            try:
                structure.make_supercell([2,2,2])
                dict_structs = OD.apply_transformation(structure,return_ranked_list=1)
            #print(list_structs)
                structure = dict_structs[0]['structure']
            except:               
                pass
        prim = SpacegroupAnalyzer(structure).find_primitive()
        print(prim)
        """
        prim = structure
        minO = get_min_O_bond(prim)
        if minO < 1.5:
            print(self.formula,str(self.tilt)+'_tilt','Warning: O-O bond length = '+str(round(minO,3))+' ang')
            with open(opj(second,'Structure_Warning'),'w') as f:
                f.write('Warning: O-O bond length = '+str(round(minO,3))+' ang')
                f.close()
            
        Poscar(prim).write_file(pos_store)
        move(self.cif_file,opj(second,self.cif_file))
        move('output.txt',opj(second,str(self.tilt)+'_output.txt'))
        move('gii.txt',opj(first,'gii.txt'))
        

       
               
if __name__ == '__main__':
    
    As = {'Sn':2}
    #Bs = {'Mn':4,'Ce':4}
    Bs = {'Co':4}
    X = {'O':2}
    
    Model = SPuDS(As,Bs,X,store_dir = 'A_B_switch')
    """
    for tilt in Model.allowed_tilts:   
        #Model.gen_default_input(tilt)
        Model.gen_JT_input(tilt)
        Model.run()
        #print(tilt)
        Model.make_cif()
        Model.store_best_poscar()
    """
        
    tilt = 10
    #Model.gen_default_input(tilt)
    Model.gen_JT_input(tilt)
    Model.run()
    Model.make_cif()
    Model.store_best_poscar()
   
    




































    
    
   



































    
    
   



































    
    
   

































