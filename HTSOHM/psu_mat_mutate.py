#ifrom numpy.random import choice
import numpy as np
#from np.random import *
import pylab as py
import os
from numpy.random import sample
from random import randrange, random, choice
import shutil


# function for finding "closest" distance over periodic boundaries
def closestDist(x_o, x_r):  
    a = 1 - x_r + x_o
    b = abs(x_r - x_o)
    c = 1 - x_o + x_r
    dx = min(a, b, c)

    return dx

# given intiral and "random" x-fraction, returns new x fraction
def deltax(x_o, x_r, strength):
    dx = closestDist(x_o, x_r)
    
    if (x_o > x_r 
        and (x_o - x_r) > 0.5
        and x_o != 1.0 
        and x_r != 1.0):
        xfrac = round((x_o + strength * random() * dx) % 1., 4)
        
    if (x_o < x_r 
        and (x_r - x_o) > 0.5
        and x_o != 1.0 
        and x_r != 1.0):
        xfrac = round((x_o - strength * random() * dx) % 1., 4)
        
    if x_o == 1.0 or x_r == 1.0:
        xfrac = round(x_o + strength * random() * dx, 4)
    
    if (x_o > x_r
        and (x_o - x_r) < 0.5
        and not (x_o == 1.0 or x_r == 1.0)):
        xfrac = round(x_o - strength * random() * dx, 4)
        
    if (x_o < x_r
        and (x_r - x_o) < 0.5
        and not (x_o == 1.0 or x_r == 1.0)):
        xfrac = round(x_o + strength * random() * dx, 4)
        
    if x_o == x_r:
        xfrac = x_o
        
    return xfrac

def mutate(data_path, atom_types, number_of_children, strength, 
           bins, bin_dimensions=['HV', 'SA', 'CH4']):

    #boundaries
    xmin = 25.6
    xmax = 51.2
    ymin = xmin
    ymax = xmax
    zmin = xmin
    zmax = xmax
    ndenmin = 0.0000013905
    ndenmax = 0.04302
    epmin = 1.2580
    epmax = 513.264
    sigmin = 1.052342
    sigmax = 6.549291

    if 'HV' in bin_dimensions:
    
        HV_labels = np.genfromtxt(data_path + '/HVdata2col.txt', usecols=0,
                                  dtype=str)
        HV_values = np.genfromtxt(data_path + '/HVdata2col.txt', usecols=1,
                                  dtype=float)
   
        HV_binmin = 0.
        HV_binmax = 1.
        HV_step = HV_binmax/bins
        edges_HV = np.arange(HV_binmin, HV_binmax + HV_step, HV_step)

        HV_IDs = np.empty(bins).tolist()
        for i in range(bins):
            HV_IDs[i] = []
    
        for i in range(bins):
            for j in range(len(HV_labels)):
                if HV_values[j] >= edges_HV[i] and HV_values[j] <= edges_HV[i + 1]:
                    HVID = HV_labels[j]
                    HV_IDs[i].append(HVID)

    if 'SA' in bin_dimensions:
    
        SA_labels = np.genfromtxt(data_path + '/SAdata_m2_cc.txt', usecols=0,
                                  dtype=str)
        SA_values = np.genfromtxt(data_path + '/SAdata_m2_cc.txt', usecols=1,
                                  dtype=float)

        SA_binmin = 0.
        SA_binmax = 4500.
        SA_step = SA_binmax/bins
        edges_SA = np.arange(SA_binmin, SA_binmax + SA_step, SA_step)

        SA_IDs = np.empty(bins).tolist()
        for i in range(bins):
            SA_IDs[i] = []
    
        for i in range(bins):
            for j in range(len(SA_labels)):
                if SA_values[j] >= edges_SA[i] and SA_values[j] <= edges_SA[i + 1]:
                    SAID = SA_labels[j]
                    SA_IDs[i].append(SAID)

    if 'CH4' in bin_dimensions:
    
        CH4_labels = np.genfromtxt(data_path + '/ch4_abs_cc_cc.txt', usecols=0,
                                   dtype=str)
        CH4_values = np.genfromtxt(data_path + '/ch4_abs_cc_cc.txt', usecols=1,
                                   dtype=float)

        CH4_binmin = 0.
        CH4_binmax = 400.
        CH4_step = CH4_binmax/bins
        edges_CH4 = np.arange(CH4_binmin, CH4_binmax + CH4_step, CH4_step)

        CH4_IDs = np.empty(bins).tolist()
        for i in range(bins):
            CH4_IDs[i] = []
    
        for i in range(bins):
            for j in range(len(CH4_labels)):
                if CH4_values[j] >= edges_CH4[i] and CH4_values[j] <= edges_CH4[i + 1]:
                    CH4ID = CH4_labels[j]
                    CH4_IDs[i].append(CH4ID)   

#################################
#################################
###                           ###
###  ONE DIMENSIONAL BINNING  ###
###                           ###
#################################
#################################
    
    if len(bin_dimensions) == 1:
        freq = np.empty(bins)
    
        if 'CH4' in bin_dimensions:
            IDs = CH4_IDs
        if 'HV' in bin_dimensions:
            IDs = HV_IDs
        if 'SA' in bin_dimensions:
            IDs = SA_IDs
    
        for i in range(bins):
            freq[i] = len(IDs[i])
        
        mat_per_bin = np.zeros(bins)
        for i in range(bins):
            mat_per_bin[i] = freq.sum() / freq[i]
    
        mat_per_bin = mat_per_bin * freq.sum() / mat_per_bin.sum()
    
        for i in range(bins):
            mat_per_bin[i] = round(mat_per_bin[i])
        
        weights = mat_per_bin / mat_per_bin.sum()

        w_list = weights
        ID_list = IDs

#################################
#################################
###                           ###
###  TWO DIMENSIONAL BINNING  ###
###                           ###
#################################
#################################

    if len(bin_dimensions) == 2:
        IDs = np.empty([bins,bins]).tolist()
        freq = np.empty([bins,bins])
        for i in range(bins):
            for j in range(bins):
                IDs[i][j] = []
    
        if 'CH4' in bin_dimensions and 'HV' in bin_dimensions:
            for i in range(bins):
                for j in range(bins):
                    CH4_list = CH4_IDs[i]
                    HV_list = HV_IDs[j]
                    for k in CH4_list:
                        for l in HV_list:
                            if k == l:
                                IDs[i][j].append(k)
                            
        if 'CH4' in bin_dimensions and 'SA' in bin_dimensions:
            for i in range(bins):
                for j in range(bins):
                    CH4_list = CH4_IDs[i]
                    SA_list = SA_IDs[j]
                    for k in CH4_list:
                        for l in SA_list:
                            if k == l:
                                IDs[i][j].append(k)

        if 'SA' in bin_dimensions and 'HV' in bin_dimensions:
            for i in range(bins):
                for j in range(bins):
                    SA_list = SA_IDs[i]
                    HV_list = HV_IDs[j]
                    for k in SA_list:
                        for l in HV_list:
                            if k == l:
                                IDs[i][j].append(k)
        
        for i in range(bins):
            for j in range(bins):
                freq[i,j] = len(IDs[i][j])
            
        mat_per_bin = np.zeros([bins,bins])
        for i in range(bins):
            for j in range(bins):
                if freq[i][j] != 0.:
                    mat_per_bin[i][j] = freq.sum() / freq[i][j]

        mat_per_bin = mat_per_bin * freq.sum() / mat_per_bin.sum() 

        for i in range(bins):
            for j in range(bins):
                mat_per_bin[i][j] = round(mat_per_bin[i][j])
    
        weights = mat_per_bin / mat_per_bin.sum()
    
        w_list = []
        ID_list = []
        for i in range(bins):
            w_list = np.concatenate([w_list, weights[i,:]])
            for j in range(bins):
                ID_list = ID_list + [IDs[i][j]]
            
#################################
#################################
###                           ###
### THREE DIMENSIONAL BINNING ###
###                           ###
#################################
#################################

    if len(bin_dimensions) == 3:
    
        IDs = np.empty([bins,bins,bins]).tolist()
        freq = np.empty([bins,bins,bins])
        for i in range(bins):
            for j in range(bins):
                for k in range(bins):
                    IDs[i][j][k] = []
                
        if 'CH4' in bin_dimensions and 'HV' in bin_dimensions and 'SA' in bin_dimensions:
            for i in range(bins):
                for j in range(bins):
                    for k in range(bins):
                        CH4_list = CH4_IDs[i]    
                        HV_list = HV_IDs[j]
                        SA_list = SA_IDs[k]
                        for l in CH4_list:
                            for m in HV_list:
                                if l == m:
                                    for n in SA_list:
                                        if m == n:
                                            IDs[i][j][k].append(n)
        
            for i in range(bins):
                for j in range(bins):
                    for k in range(bins):
                        freq[i,j,k] = len(IDs[i][j][k])
                    
        mat_per_bin = np.zeros([bins,bins,bins])
        for i in range(bins):
            for j in range(bins):
                for k in range(bins):
                    if freq[i][j][k] != 0.:
                        mat_per_bin[i][j][k] = freq.sum() / freq[i][j][k]

        mat_per_bin = mat_per_bin * freq.sum() / mat_per_bin.sum() 

        for i in range(bins):
            for j in range(bins):
                for k in range(bins):
                    mat_per_bin[i][j][k] = round(mat_per_bin[i][j][k])
    
        weights = mat_per_bin / mat_per_bin.sum()
    
        w_list = []
        ID_list = []
        for i in range(bins):
            for j in range(bins):
                w_list = np.concatenate([w_list, weights[i,j,:]])
                for k in range(bins):
                    ID_list = ID_list + [IDs[i][j][k]]

    print 'original distribution :  \n',freq
    print '\ntotal materials       :   ', freq.sum()
    print '\n\nbinned distribution   :   \n', mat_per_bin
    print '\ntotal materials       :   ', mat_per_bin.sum()


    print 'weights :   \n', weights

    top_path = str(atom_types) + 'atmtyp_' + str(number_of_children) + 'mutants_strength' + str(int(strength*100))

    if not os.path.exists(top_path):
        os.mkdir(top_path)

    for i in range(number_of_children):
        print 'making MAT-' + str(i) + '...'        

        ID_bin = np.random.sample(ID_list, p=w_list)
        parent_ID = np.random.sample(ID_bin)

        # creating child directory...
        mat_path = top_path + '/MAT-' + str(i)
        os.mkdir(mat_path)
        
        # copying parent data files...
        shutil.copy(data_path + '/' + parent_ID + 
                    '/force_field_mixing_rules.def', mat_path + 
                    '/old_mixing_rules.txt')
        
        shutil.copy(data_path + '/' + parent_ID + 
                    '/pseudo_atoms.def', mat_path) 
        shutil.copy(data_path + '/' + parent_ID + 
                    '/force_field.def', mat_path) 
        shutil.copy(data_path + '/' + parent_ID + '/' + parent_ID
                    + '.cif', mat_path + '/old_cif.txt')                        
        
        # importing values from parent data files...
        n1, n2, ep_o, sig_o = np.genfromtxt(mat_path + 
                                            '/old_mixing_rules.txt',
                                            unpack=True,
                                            skip_header=7,
                                            skip_footer=9)
        
        cif_atyp = np.genfromtxt(mat_path + '/old_cif.txt', usecols=0,
                                 dtype=str, skip_header=19)
        
        n1, n2, x_o, y_o, z_o, q_o = np.genfromtxt(mat_path + 
                                                   '/old_cif.txt',
                                                   unpack=True,
                                                   skip_header=19)
        # opening child data files ...            
        cif_file = open(os.path.abspath(mat_path) + '/MAT-' + str(i) + 
                        '.cif', 'w')
        
        mixing_file = open(os.path.abspath(mat_path) + 
                           '/force_field_mixing_rules.def', 'w')

        # perturbing crystal lattice parameters (a,b,c)...        
        a_footer = len(x_o)+12
        b_footer = len(x_o)+11
        c_footer = len(x_o)+10
        
        n1, a_o = np.genfromtxt(mat_path + '/old_cif.txt', unpack=True,
                                skip_header=5, skip_footer=a_footer)
        
        n1, b_o = np.genfromtxt(mat_path + '/old_cif.txt', unpack=True,
                                skip_header=6, skip_footer=b_footer)
        
        n1, c_o = np.genfromtxt(mat_path + '/old_cif.txt', unpack=True,
                                skip_header=7, skip_footer=c_footer)
        
        a_r = randrange(xmin * 10.**4, xmax * 10.**4) / 10.**4
        b_r = randrange(xmin * 10.**4, xmax * 10.**4) / 10.**4
        c_r = randrange(xmin * 10.**4, xmax * 10.**4) / 10.**4

        aval = round(a_o + strength * random() * (a_r - a_o), 4)
        bval = round(b_o + strength * random() * (b_r - b_o), 4)
        cval = round(c_o + strength * random() * (c_r - c_o), 4)

        # writing new crystal lattice parameters to file...                    
        cif_head = ('material' + str(i) + '\n' +
                    '\n' +
                    'loop_\n' +
                    '_symmetry_equiv_pos_as_xyz\n' +
                    '  x,y,z\n' +
                    '_cell_length_a          ' + str(aval) + '\n' +
                    '_cell_length_b          ' + str(bval) + '\n' +
                    '_cell_length_c          ' + str(cval) + '\n' +
                    '_cell_angle_alpha       90.0000\n' +
                    '_cell_angle_beta        90.0000\n' +
                    '_cell_angle_gamma       90.0000\n' +
                    'loop_\n' +
                    '_atom_site_label\n' +
                    '_atom_site_type_symbol\n' +
                    '_atom_site_fract_x\n' +
                    '_atom_site_fract_y\n' +
                    '_atom_site_fract_z\n' +
                    '_atom_site_charge\n' +
                    '#NET CHARGE: ???\n')
        cif_file.write(cif_head)
        

        # perturbing number density...
        n_o = len(x_o)
        nden_o = n_o / (a_o * b_o * c_o)
                    
        nden_r = randrange(ndenmin * 10.**16, ndenmax * 10.**16) / 10.**16
                    
        nden = nden_o + strength * random() * (nden_r - nden_o)
        n = int(nden * aval * bval * cval)

        # removing pseudo-atoms from unit cell (if necessary)...
        omit_n = 0
        if n < n_o:
            omit_n = n_o - n

        # perturbing atomic positions (xfrac, yfrac, zfrac)...
        for l in range(n_o - omit_n):
                        
            x_r = random()
            y_r = random()
            z_r = random()

            xfrac = deltax(x_o[l], x_r, strength)
            yfrac = deltax(y_o[l], y_r, strength)
            zfrac = deltax(z_o[l], z_r, strength)

            charge = q_o[l]
         
            cif_line = (cif_atyp[l] + '     C     ' + str(xfrac) +
                        '     ' + str(yfrac) + '     ' + 
                        str(zfrac) + '    ' + str(charge) + '\n')
            cif_file.write(cif_line)
        
        # adding pseudo-atoms to unit cell (if necessary)...
        if n > n_o:
            add_n = n - n_o
            for m in range(add_n):

                atyp = choice(cif_atyp)
                charge = choice(q_o)                            
 
                xfrac = round(random(), 4)
                yfrac = round(random(), 4)
                zfrac = round(random(), 4)

                new_line = (cif_atyp[l] + '     C     ' + str(xfrac) +
                            '     ' + str(yfrac) + '     ' + 
                            str(zfrac) + '    ' + str(charge) + '\n')
                cif_file.write(new_line)
        
        mixing_head = ('# general rule for shifted vs truncated\n' +
                       'shifted\n' +
                       '# general rule tailcorrections\n' +
                       'no\n' +
                       '# number of defined interactions\n' +
                       str(atom_types + 8) + '\n' +
                       '# type interaction, parameters.    IMPORTANT:' + 
                       ' define shortest matches first, so that more ' +
                       'specific ones overwrites these\n')
        mixing_file.write(mixing_head)
        
        # perturbing LJ parameters (sigma, epsilon)...
        for o in range(atom_types):
            
            ep_r = randrange(epmin * 10.**16, epmax * 10.**16) / 10.**16
            sig_r = randrange(sigmin * 10.**16, sigmax * 10.**16) / 10.**16

            epsilon = round(ep_o[o] + strength * random() * (ep_r - ep_o[o]), 4) 
            sigma = round(sig_o[o] + strength * random() * (sig_r - sig_o[o]), 4)
        
            mixing_line = ('A_' + str(o) + '   lennard-jones   ' + 
                                       str(epsilon) + '   ' + str(sigma) + '\n')
            mixing_file.write(mixing_line)
                    
        mixing_foot = ('N_n2        lennard-jones   36.0     3.31\n' +
                       'N_com       none\n' +
                       'C_co2       lennard-jones   27.0     2.80\n' +
                       'O_co2       lennard-jones   79.0     3.05\n' +
                       'CH4_sp3     lennard-jones   158.5    3.72\n' +
                       'He          lennard-jones   10.9     2.64\n' +
                       'H_h2        none\n' +
                       'H_com       lennard-jones   36.7     2.958\n' +
                       '# general mixing rule for Lennard-Jones\n' +
                       'Lorentz-Berthelot')
        mixing_file.write(mixing_foot)

        cif_file.close() 
        mixing_file.close()

