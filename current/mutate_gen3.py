import numpy as np
import pylab as py
import os
from random import choice, randrange, random
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

def mutate(atom_types, strength):

    #boundaries
    xmin = 14.8193
    xmax = 52.3940
    ymin = xmin
    ymax = xmax
    zmin = xmin
    zmax = xmax
    ndenmin = 0.000013905
    ndenmax = 0.04243
    epmin = 2.0128
    epmax = 410.6115
    sigmin = 1.6838
    sigmax = 5.2394

    #atom_types = 4
    data_path = 'gen3'

    SA_labels = np.genfromtxt(data_path + '/SAdata_m2_cc.txt', usecols=0,
                              dtype=str)
    SA_values = np.genfromtxt(data_path + '/SAdata_m2_cc.txt', usecols=1,
                              dtype=float)

    HV_labels = np.genfromtxt(data_path + '/HVdata2col.txt', usecols=0,
                              dtype=str)
    HV_values = np.genfromtxt(data_path + '/HVdata2col.txt', usecols=1,
                              dtype=float)

    freq,SAedges,HVedges = np.histogram2d(SA_values, HV_values, bins=10)

    IDs = np.empty(np.shape(freq)).tolist()
    SA_IDs = np.empty(len(freq)).tolist()
    HV_IDs = np.empty(len(freq)).tolist()

    for i in range(len(freq)):
        SA_IDs[i] = []
        HV_IDs[i] = []
        for j in range(len(freq)):
            IDs[i][j] = []

    for i in range(len(freq)):
    
        for j in range(len(SA_labels)):
            if SA_values[j] >= SAedges[i] and SA_values[j] <= SAedges[i + 1]:
            #SA_IDs[i] = SA_IDs[i] + [SA_labels[j]]
                SAID = SA_labels[j]
                SA_IDs[i].append(SAID)

        for j in range(len(HV_labels)):
            if HV_values[j] >= HVedges[i] and HV_values[j] <= HVedges[i + 1]:
                HVID = HV_labels[j]
                HV_IDs[i].append(HVID)

    for i in range(len(freq)):
        for j in range(len(freq)):
            print i, j
            SAlist = SA_IDs[i]
            HVlist = HV_IDs[j]
            for k in SAlist:
                for l in HVlist:
                    if k == l:
                        IDs[i][j].append(k)


    top_path = 'gen3_strength' + str(int(strength*100))

    if not os.path.exists(top_path):
        os.mkdir(top_path)

    print freq

    k = 0


#    mat_per_bin = (1 - freq / np.max(freq)) * np.max(freq)
#    print mat_per_bin, type(mat_per_bin), mat_per_bin.sum()
#    scale_factor = freq.sum() / mat_per_bin.sum()
#    mat_per_bin = mat_per_bin * scale_factor    
#    print mat_per_bin, type(mat_per_bin), mat_per_bin.sum()

    freq_percent = freq / freq.sum() 
    bin_percent = freq_percent

#* np.max(freq)
    for i in range(len(freq)):
        for j in range(len(freq)):
            if freq[i][j] != 0.:
                bin_percent[i][j] = 1 - freq_percent[i][j]
            if freq[i][j] == 0.:
                bin_percent[i][j] == 0.
    bin_percent = bin_percent / bin_percent.sum()

    print bin_percent, bin_percent.sum()

    mat_per_bin = bin_percent
    for i in range(len(freq)):
        for j in range(len(freq)):
            mat_per_bin[i][j] = round(bin_percent[i][j] * freq.sum())

    print mat_per_bin, mat_per_bin.sum()

    for i in range(len(IDs)):
        for j in range(len(IDs)):
            if len(IDs[i][j]) != 0:    
                for n in range(int(mat_per_bin[i,j])):
    #                parent_ID = choice(IDs[i][j])
    
                    parent_ID = choice(IDs[i][j])
                    #parent_ID = IDs[0][0][index]
    
                    mat_path = top_path + '/MAT-' + str(k)
                    os.mkdir(mat_path)
                    shutil.copy(data_path + '/' + parent_ID + 
                                '/force_field_mixing_rules.def', mat_path + 
                                '/old_mixing_rules.txt')
        
                    shutil.copy(data_path + '/' + parent_ID + 
                                '/pseudo_atoms.def', mat_path) 
                    shutil.copy(data_path + '/' + parent_ID + 
                                '/force_field.def', mat_path) 
                    shutil.copy(data_path + '/' + parent_ID + '/' + parent_ID
                                + '.cif', mat_path + '/old_cif.txt')                        
        #+ '.cif', mat_path + '/MAT-' + str(k) + '.cif') 
        
                    print 'creating MAT-'+ str(k) + '...'
        
                    n1, n2, ep_o, sig_o = np.genfromtxt(mat_path + 
                                                            '/old_mixing_rules.txt',
                                                            unpack=True,
        #                                            usecols=(2, 3),
                                                            skip_header=7,
                                                            skip_footer=9)
        
                    cif_atyp = np.genfromtxt(mat_path + '/old_cif.txt', usecols=0,
                                             dtype=str, skip_header=19)
        
                    n1, n2, x_o, y_o, z_o, q_o = np.genfromtxt(mat_path + 
                                                               '/old_cif.txt',
                                                               unpack=True,
        #                                          usecols=(2, 3, 4),
                                                               skip_header=19)
                    
                    cif_file = open(os.path.abspath(mat_path) + '/MAT-' + str(k) + 
                                    '.cif', 'w')
        
                    mixing_file = open(os.path.abspath(mat_path) + 
                                    '/force_field_mixing_rules.def', 'w')
        
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
                    
                    cif_head = ('material' + str(k) + '\n' +
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
                    
                    n_o = len(x_o)
                    nden_o = n_o / (a_o * b_o * c_o)
                    
                    nden_r = randrange(ndenmin * 10.**9, ndenmax * 10.**9) / 10.**9
                    
                    nden = nden_o + strength * random() * (nden_r - nden_o)
                    n = int(nden * aval * bval * cval)

                    omit_n = 0
                    if n < n_o:
                        omit_n = n_o - n

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
        
                    for o in range(atom_types):
        
                        ep_r = randrange(epmin * 10.**4, epmax * 10.**4) / 10.**4
                        sig_r = randrange(sigmin * 10.**4, sigmax * 10.**4) / 10.**4

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
    
                    k = k + 1

