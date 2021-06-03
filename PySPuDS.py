# -*- coding: utf-8 -*-
"""
Created on Fri Apr 16 18:56:04 2021

@author: zaba1157
"""

import os
from shutil import move
import subprocess
import re
from pymatgen.io.cif import CifBlock,CifFile
from collections import OrderedDict
from monty.io import zopen
import json
import pandas as pd

module_dir = os.getcwd()
SPuDS_install_directory = os.getcwd()

#"""Edit path to your SPuDS installation directory."""
#SPuDS_install_directory = '/path/to/your/SPuDS/install_directory'


class SPuDS_BVparams():
    """
    #***************************************************************
    # COPYRIGHT NOTICE
    # This table may be used and distributed without fee for
    # non-profit purposes providing 
    # 1) that this copyright notice is included and 
    # 2) no fee is charged for the table and 
    # 3) details of any changes made in this list by anyone other than
    # the copyright owner are suitably noted in the _audit_update record
    # Please consult the copyright owner regarding any other uses.
    #
    # The copyright is owned by I. David Brown, Brockhouse Institute for
    # Materials Research, McMaster University, Hamilton, Ontario Canada.
    # idbrown@mcmaster.ca
    #
    #*****************************DISCLAIMER************************
    #
    # The values reported here are taken from the literature and
    # other sources and the author does not warrant their correctness
    # nor accept any responsibility for errors.  Users are advised to
    # consult the primary sources. 
    #
    #***************************************************************
    """    
    def __init__(self):
        """
        Reads table of SPuDS provided bond valence parameters.
        """
        bvfile = os.path.join(SPuDS_install_directory, "bvparm20.cif")
        params = pd.read_csv(bvfile, sep='\s+',
                                  header=None,
                                  names=['Atom1', 'Atom1_valence',
                                         'Atom2', 'Atom2_valence',
                                         'Ro', 'B', 'ref_id', 'details'],
                                  skiprows=179,
                                  skipfooter=1,
                                  index_col=False,
                                  engine="python")
        self.params = params        
        
    def get_bv_params(self, cation, anion, cat_val, an_val):
        """
        Retrieves bond valence parameters from SPuDS table.       
        Args:
            cation str(el): cation
            anion str(el): anion
            cat_val int: cation formal oxidation state
            an_val int: anion formal oxidation state
        Returns:
            bvlist: dataframe of bond valence parameters            
        """
        bv_data = self.params
        bvlist = self.params.loc[(bv_data['Atom1'] == str(cation)) \
                                & (bv_data['Atom1_valence'] == cat_val) \
                                & (bv_data['Atom2'] == str(anion)) \
                                & (bv_data['Atom2_valence'] == an_val)]
        return bvlist.iloc[0] # First entry if multiple values exist.


class SPuDS():
    """
    Args:
        A: dict of A cation and oxidation state, {Ca:2}
        B: dict of B cation and oxidation state, {Ti:4}
        X: dict of X anion and oxidation state, {O:-2}
        
    Implemented SPuDS Menus:
        1. ABX3
    
    #TODO Add Aditional SPuDS Menus:   
        3. AA'B2X6
        6. A2BB'X6
        7. AA'BB'X6     
    """  
    
    def __init__(self, As, Bs, Xs, store_dir='SPuDS_output'):
        
        # Read in ABX3 space group symmetry operations
        with open(os.path.join(module_dir, 'ABX3_SPuDS_symops.json'),'r') as f:
            self.symops_dict = json.load(f)
            f.close()
            
        # Load dataframe of tabulated bond valence parameters    
        self.bv = SPuDS_BVparams()
        
        # Assign __init__ args to model
        self.As = As
        self.Bs = Bs
        self.Xs = Xs
        self.store_dir = store_dir
        self.SPuDS_dir = SPuDS_install_directory        
        
        # SPuDS compatible nomenclature for single letter elements (i.e. K_) 
        spuds_single_cats = ['K','P','U','V','W','Y']
        for el in As.keys():
            if el in spuds_single_cats:
                As[el+'_'] = As.pop(el)                    
        for el in Bs.keys():
            if el in spuds_single_cats:
                Bs[el+'_'] = Bs.pop(el)        
                
        # Assign elements and oxidation states to model
        self.Ael = [el for el in As.keys()][0]
        self.Bel = [el for el in Bs.keys()][0]
        self.Xel = [el for el in Xs.keys()][0]
        self.n = len(As.keys()) + len(Bs.keys())
        
        els = [el for el in As.keys()]
        els.extend([el for el in Bs.keys()])
        els.extend([el for el in Xs.keys()])
        self.ellist = els
        
        oxis = [oxi for oxi in As.values()]
        oxis.extend([oxi for oxi in Bs.values()])
        oxis.extend([oxi for oxi in Xs.values()])
        self.oxilist = oxis
        
        # Assign SPuDS interface capabilities to model        
        self.allowed_menus = [1]
        self.menu = 1
        self.allowed_tilts = [3,5,10,13,14,16,17,20,21,22,23]  
        # ABX3 formula and storage nameing scheme          
        self.formula = "".join(self.ellist)+"3"
        self.store_form = self.formula+'__'+"".join([str(x)+'('+str(oxis[i])+')'
                                                     for i,x in enumerate(els)])
        # Some minimal error checking
        if len(X.keys()) > 1:
            raise Exception('Only single X anion supported.')            
        if len(As.keys()) > 1 or len(Bs.keys()) > 1:
            raise Exception('Only ABX3 compositions currently supported.')           
        # Precheck A-site in tabulated bond valence parameters
        try:
            self.bv.get_bv_params(self.Ael,self.Xel,As[self.Ael],Xs[self.Xel])
        except:
            raise Exception('A-site '+self.Ael+'('+str(As[self.Ael])+')-'+
                            self.Xel+'('+str(Xs[self.Xel])+')'+
                            ' not in tabulated bv parameters')
        # Precheck B-site in tabulated bond valence parameters
        try:
            self.bv.get_bv_params(self.Bel,self.Xel,Bs[self.Bel],Xs[self.Xel])
        except:
            raise Exception('B-site '+self.Bel+'('+str(Bs[self.Bel])+')-'+
                            self.Xel+'('+str(Xs[self.Xel])+')'+
                            ' not in tabulated bv parameters')            
        

    def write_default_input(self,tilt):
        """    
        Generates and writes the SPuDS input.txt file from default settings
        
        Args: 
            tilt: Glazer tilt number (int)
        
        Defaults: 
            No Jahn-Teller distortions and T = 298 K (i.e. dR/dT = 0.0)      
        """
        # Check if specified tilt in self.allowed_tilts
        if tilt not in self.allowed_tilts:
            raise Exception('Specified tilt not in allowed_tilts')
        # Assign defaults to model
        self.tilt = tilt        
        self.JT = 0
        self.JT_params = ['-0.1','-0.04','+0.18']
        self.T = 298
        self.dRdT = [0] * self.n
        # Get ordered list of write lines for input.txt
        writelines = [self.menu, tilt, self.JT, 0]
        for j in self.JT_params:
            writelines.append(j)            
        writelines.append(self.formula)        
        for oxi in self.oxilist:
            writelines.append(abs(oxi))            
        writelines.append(self.T)        
        for drdt in self.dRdT:
            writelines.append(drdt)
        # Write the input.txt file
        with open(os.path.join(self.SPuDS_dir,'input.txt'),'w') as f:
            for line in writelines:
                f.write(str(line)+"\n")
            f.close()

    def run(self):
        """
        Runs SPuDS using input parameters specified in the input.txt file
        """
        inputs = "\n".join(["0","0"])            
        result = subprocess.run(os.path.join(self.SPuDS_dir,'spuds'), 
                                input=inputs, capture_output=True, text=True)
        # Print SPuDS run stdout
        print(result.stdout)

    def write_cif(self):
        """
        Generates and writes a structure file (.cif) from SPuDS output.txt 
        """
        
        def parse_spuds_out(self): 
            """
            Read SPuDS output.txt and get predicted structure info
            Returns:
                site_list: list of lists of site info 
                           list([element, multiplicity, Wycoff label,
                           x-coord, y-coord, z-coord, site occupancy])
                a_len: SPuDS a lattice length (Ang)
                b_len: SPuDS b lattice length (Ang)
                c_len: SPuDS c lattice length (Ang) 
                alp: SPuDS alpha lattice angle (deg)
                bet: SPuDS beta lattice angle (deg)
                gam: SPuDS gamma lattice angle (deg)                
                
            """
            
            a_lat,b_lat,c_lat,alp,bet,gam = None,None,None,None,None,None
            dvdr = ['***********************************'+
                    '***********************************']
            site_list = []
            with open(os.path.join(self.SPuDS_dir,'output.txt'),'r') as f:
                linecount = 0
                count = 0
                for line in f:
                    linelist = line.split()
                    if linelist == dvdr:
                        count += 1                            
                    if count == 2:
                        if len(linelist) == 4:
                            if linelist[0] == 'Space' and linelist[1] == 'group':
                                # Get dict key for retrieving symmetry operations
                                self.symops_key = " ".join([linelist[2],linelist[3]])
                        if len(linelist) == 3:
                            # Get lattice lengths (Ang)
                            if linelist[0] == 'a':
                                a_lat = linelist[2]
                            if linelist[0] == 'b':
                                b_lat = linelist[2]
                            if linelist[0] == 'c':
                                c_lat = linelist[2] 
                            # Get lattice angles (deg)
                            if linelist[0] == 'alpha':
                                alp = linelist[2]
                            if linelist[0] == 'beta':
                                bet = linelist[2]
                            if linelist[0] == 'gamma':
                                gam = linelist[2]                    
                        if len(linelist) > 3:
                            # Get site_list
                            for el in self.ellist:
                                if linelist[0] == el:                                    
                                    if el == self.Xel:
                                        # Anion should be fully occupied
                                        occup = 1.00
                                    else:
                                        # Cations can be partially occupied
                                        occup = linelist[5]                                        
                                    multsite = re.split('(\d+)',linelist[1])
                                    mult = multsite[1]
                                    wycoff = multsite[2]
                                    x = linelist[2]
                                    y = linelist[3]
                                    z = linelist[4]
                                    # Append site info to site_list
                                    site_list.append([el,mult,wycoff,
                                                        x,y,z,occup])
                    linecount += 1
                f.close()
            
            # Add count labels to site_list
            elcounts = {}
            for d in site_list:
                for el in self.ellist:
                    if d[0] == el:
                        if el not in elcounts.keys():
                            elcounts[el] = 0
                        elcounts[el] += 1  
                          
            for k,v in elcounts.items():
                count = 1
                ii = 0
                for d in site_list:
                    if d[0] == k:
                        if v > 1:
                            label = k+str(count)
                            site_list[ii].insert(0,label)
                            count += 1
                        else:
                            label = k
                            site_list[ii].insert(0,label)
                    ii += 1


            return site_list, a_lat, b_lat, c_lat, alp, bet, gam
        
        def map_lattice_menu_1(self,a_len,b_len,c_len,alp,bet,gam):
            """
            Map to .cif compatible lattice parameters from SPuDS output.txt
            Args: 
                a_len: SPuDS a lattice length (Ang)
                b_len: SPuDS b lattice length (Ang)
                c_len: SPuDS c lattice length (Ang) 
                alp: SPuDS alpha lattice angle (deg)
                bet: SPuDS beta lattice angle (deg)
                gam: SPuDS gamma lattice angle (deg)
            Returns: 
                .cif compatible lattice parameters:
                a: a lattice length (Ang)
                b: b lattice length (Ang)
                c: c lattice length (Ang)
                alpha: alpha lattice angle (deg)
                beta: beta lattice angle (deg)
                gamma: gamma lattice angle (deg)                
            """
            # Lattice angles (deg)
            if self.tilt in [3,5,10,16,17,20,21,22,23]:
                alpha = 90
                beta = 90
                gamma = 90                 
            elif self.tilt == 13:
                alpha = 90
                beta = bet
                gamma = 90
            elif self.tilt == 14:
                alpha = 90
                beta = 90
                gamma = 120 
            # Lattice lengths (Ang)    
            if self.tilt in [3,23]:
                a = a_len
                b = a_len
                c = a_len
            elif self.tilt in [5,14,16,21,22]:
                a = a_len
                b = a_len
                c = c_len
            elif self.tilt in [10,13,17,20]:
                a = a_len
                b = b_len
                c = c_len                
 
            return a,b,c,alpha,beta,gamma    
        
        def make_cif(self):
            """
            Generates a pymatgen CifFile object using structure info parsed
            from SPuDS output.txt.
            Returns:
                cf: pymatgen CifFile object
            """
            # SPuDS ouput structure info
            site_list,a_lat,b_lat,c_lat,alp,bet,gam = parse_spuds_out(self)
            # Mapped lattice parameters to .cif compatibility
            a,b,c,alpha,beta,gamma = map_lattice_menu_1(self,
                                                        a_lat,b_lat,c_lat,
                                                        alp,bet,gam) 
            symd = self.symops_dict[self.symops_key] # symops dict data         
            # Create dict of .cif parameters                   
            data = {}
            data['_cell_length_a'] = a
            data['_cell_length_b'] = b
            data['_cell_length_c'] = c
            data['_cell_angle_alpha'] = alpha
            data['_cell_angle_beta'] = beta
            data['_cell_angle_gamma'] = gamma
            data['_space_group_name_H-M_alt'] = symd['name']
            data['_symmetry_Int_tables_number'] = symd['number']
            data['_symmetry_cell_setting'] = symd['latsym']        
            data['_space_group_symop_operation_xyz'] = symd['symops']        
            data['_atom_type_symbol'] = self.ellist
            data['_atom_type_oxidation_number'] = self.oxilist
                    
            data['_atom_site_label'] = [d[0] for d in site_list]
            data['_atom_site_type_symbol'] = [d[1] for d in site_list]
            data['_atom_site_symmetry_multiplicity'] = [d[2] for d in site_list]
            data['_atom_site_Wycoff_symbol'] = [d[3] for d in site_list]
            data['_atom_site_fract_x'] = [d[4] for d in site_list]
            data['_atom_site_fract_y'] = [d[5] for d in site_list]
            data['_atom_site_fract_z'] = [d[6] for d in site_list]
            data['_atom_site_occupancy'] = [d[7] for d in site_list]
            # .cif file header
            cif_header = 'SPuDS'
            # .cif file loops
            cif_loops = [['_space_group_symop_operation_xyz'],                        
                        ['_atom_type_symbol','_atom_type_oxidation_number'],                      
                        ['_atom_site_label','_atom_site_type_symbol',
                          '_atom_site_symmetry_multiplicity',
                          '_atom_site_Wycoff_symbol','_atom_site_fract_x',
                          '_atom_site_fract_y','_atom_site_fract_z',
                          '_atom_site_occupancy']]
            # Create CifFile object
            d = OrderedDict()
            d[self.formula] = CifBlock(data,cif_loops,cif_header)
            cf = CifFile(d)
            
            return cf

        # .cif file nameing scheme
        self.cif_file = str(self.tilt)+'_'+self.formula+'.cif'
        # Generate pymatgen CifFile object
        cf = make_cif(self)
        # Write created .cif file to SPuDS_dir
        with zopen(os.path.join(self.SPuDS_dir,self.cif_file), "wt") as f:
            f.write(cf.__str__())
            f.close()
            
    def store_results(self):
        """
        Move SPuDS output files and generated .cif structure file to 
        specified storage directory (store_dir). Storage nameing scheme moves
        files to /store_dir/store_form/files
        """
        store_path = os.path.join(self.store_dir,self.store_form)

        # Create storage directory if does not already exist
        if os.path.exists(self.store_dir) == False:
            os.mkdir(self.store_dir)        
        # Create storage subdirectory if does not already exist
        if os.path.exists(store_path) == False:
            os.mkdir(store_path)
            
        # Move files to /store_dir/store_form/    
        move(os.path.join(self.SPuDS_dir,'gii.txt'), 
             os.path.join(store_path,'gii.txt'))
        move(os.path.join(self.SPuDS_dir,'output.txt'),
             os.path.join(store_path,str(self.tilt)+'_output.txt'))
        move(os.path.join(self.SPuDS_dir,self.cif_file),
             os.path.join(store_path,self.cif_file))

    
if __name__ == '__main__':
    
    A = {'Ca':2}
    B = {'Ti':4}
    X = {'O':-2}
    
       
    Model = SPuDS(A,B,X,store_dir = 'PySPuDS_results')
    
    for tilt in Model.allowed_tilts:   
        Model.write_default_input(tilt)
        Model.run()
        Model.write_cif()
        Model.store_results()
