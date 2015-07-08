import numpy as np
import pylab as py
import os
from random import randrange, random
from numpy.random import choice
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

def mutate(atom_types, number_of_parents, number_of_children, strength, 
           bins, bin_dimensions=['HV', 'SA', 'CH4']):

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
    data_path = str(number_of_parents) + 'mat_' + str(atom_types) + 'atmtyp'

    if 'HV' in bin_dimensions:
    
        HV_labels = np.genfromtxt(data_path + '/HVdata2col.txt', usecols=0,
                                  dtype=str)
        HV_values = np.genfromtxt(data_path + '/HVdata2col.txt', usecols=1,
                                  dtype=float)
    
        freq_HV, edges_HV = np.histogram(HV_values, bins=bins)

        HV_IDs = np.empty(bins).tolist()
        for i in range(bins):
            HV_IDs[i] = []
    
        for i in range(bins):
            for j in range(len(HV_labels)):
                if HV_values[j] >= edges_HV[i] and HV_values[j] <= edges_HV[i + 1]:
                    HVID = HV_labels[j]
                    HV_IDs[i].append(HVID)
                             
    #print freq_HV, edges_HV, np.shape(HV_IDs)
    #print HV_IDs[0][0]

    if 'SA' in bin_dimensions:
    
        SA_labels = np.genfromtxt(data_path + '/SAdata_m2_cc.txt', usecols=0,
                                  dtype=str)
        SA_values = np.genfromtxt(data_path + '/SAdata_m2_cc.txt', usecols=1,
                                  dtype=float)

        freq_SA, edges_SA = np.histogram(SA_values, bins=bins)

        SA_IDs = np.empty(bins).tolist()
        for i in range(bins):
            SA_IDs[i] = []
    
        for i in range(bins):
            for j in range(len(SA_labels)):
                if SA_values[j] >= edges_SA[i] and SA_values[j] <= edges_SA[i + 1]:
                    SAID = SA_labels[j]
                    SA_IDs[i].append(SAID)
                             
    #print freq_SA, edges_SA, np.shape(SA_IDs)
    #print SA_IDs[0][0]

    if 'CH4' in bin_dimensions:
    
        CH4_labels = np.genfromtxt(data_path + '/ch4_abs_cc_cc.txt', usecols=0,
                                   dtype=str)
        CH4_values = np.genfromtxt(data_path + '/ch4_abs_cc_cc.txt', usecols=1,
                                   dtype=float)
    
        freq_CH4, edges_CH4 = np.histogram(CH4_values, bins=bins)

        CH4_IDs = np.empty(bins).tolist()
        for i in range(bins):
            CH4_IDs[i] = []
    
        for i in range(bins):
            for j in range(len(CH4_labels)):
                if CH4_values[j] >= edges_CH4[i] and CH4_values[j] <= edges_CH4[i + 1]:
                    CH4ID = CH4_labels[j]
                    CH4_IDs[i].append(CH4ID)
                             
    #print freq_CH4, edges_CH4, np.shape(CH4_IDs)
    #print CH4_IDs[0][0]    

#################################
#################################
###                           ###
###  ONE DIMENSIONAL BINNING  ###
###                           ###
#################################
#################################
    
    if len(bin_dimensions) == 1:
    #IDs = np.empty(bins).tolist()
    #for i in range(bins):
    #    for j in range(bins):
    #        IDs[i] = []
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
            
#        print np.shape(IDs)
    
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

#print np.shape(IDs)

    print 'original distribution :  \n',freq
    print '\ntotal materials       :   ', freq.sum()
    print '\n\nbinned distribution   :   \n', mat_per_bin
    print '\ntotal materials       :   ', mat_per_bin.sum()


    print 'weights :   \n', weights

#    if SA in dims:
#        SA_labels = np.genfromtxt(data_path + '/SAdata_m2_cc.txt', usecols=0,
#                                  dtype=str)
#        SA_values = np.genfromtxt(data_path + '/SAdata_m2_cc.txt', usecols=1,
#                                  dtype=float)

#    if HV in dims:
#        HV_labels = np.genfromtxt(data_path + '/HVdata2col.txt', usecols=0,
#                                  dtype=str)
#        HV_values = np.genfromtxt(data_path + '/HVdata2col.txt', usecols=1,
#                                  dtype=float)
#    if CH4 in dims:
#        CH4_labels = np.genfromtxt(data_path + '/CH4_abs_cc_cc.txt', usecols=0,
#                                  dtype=str)
#        CH4_values = np.genfromtxt(data_path + '/CH4_abs_cc_cc.txt', usecols=1,
#                                  dtype=float)
        
#    freq,SAedges,HVedges = np.histogram2d(SA_values, HV_values, bins=5)

#    IDs = np.empty(np.shape(freq)).tolist()
#    SA_IDs = np.empty(len(freq)).tolist()
#    HV_IDs = np.empty(len(freq)).tolist()

#    for i in range(len(freq)):
#        SA_IDs[i] = []
#        HV_IDs[i] = []
#        for j in range(len(freq)):
#            IDs[i][j] = []

#    for i in range(len(freq)):
    
#        for j in range(len(SA_labels)):
#            if SA_values[j] >= SAedges[i] and SA_values[j] <= SAedges[i + 1]:
            #SA_IDs[i] = SA_IDs[i] + [SA_labels[j]]
#                SAID = SA_labels[j]
#                SA_IDs[i].append(SAID)

#        for j in range(len(HV_labels)):
#            if HV_values[j] >= HVedges[i] and HV_values[j] <= HVedges[i + 1]:
#                HVID = HV_labels[j]
#                HV_IDs[i].append(HVID)

#    for i in range(len(freq)):
#        for j in range(len(freq)):
#            print i, j
#            SAlist = SA_IDs[i]
#            HVlist = HV_IDs[j]
#            for k in SAlist:
#                for l in HVlist:
#                    if k == l:
#                        IDs[i][j].append(k)


    top_path = str(atom_types) + 'atmtyp_' + str(number_of_children) + 'mutants_strength' + str(int(strength*100))

    if not os.path.exists(top_path):
        os.mkdir(top_path)

#    print freq, freq.sum()

#    k = 0


#    mat_per_bin = (1 - freq / np.max(freq)) * np.max(freq)
#    print mat_per_bin, type(mat_per_bin), mat_per_bin.sum()
#    scale_factor = freq.sum() / mat_per_bin.sum()
#    mat_per_bin = mat_per_bin * scale_factor    
#    print mat_per_bin, type(mat_per_bin), mat_per_bin.sum()

#    freq_percent = freq / freq.sum() 
#    bin_percent = freq_percent

#* np.max(freq)
#    for i in range(len(freq)):
#        for j in range(len(freq)):
#            if freq[i][j] != 0.:
#                bin_percent[i][j] = 1 - freq_percent[i][j]
#            if freq[i][j] == 0.:
#                bin_percent[i][j] == 0.
#    bin_percent = bin_percent / bin_percent.sum()

#    print bin_percent, bin_percent.sum()

#    mat_per_bin = bin_percent
#    for i in range(len(freq)):
#        for j in range(len(freq)):
#            mat_per_bin[i][j] = round(bin_percent[i][j] * freq.sum())

#    mat_per_bin = np.zeros([len(freq),len(freq)])
#    for i in range(len(freq)):
#        for j in range(len(freq)):
#            if freq[i][j] != 0.:
#                mat_per_bin[i][j] = freq.sum() / freq[i][j]

#    mat_per_bin = mat_per_bin * freq.sum() / mat_per_bin.sum() 

#    for i in range(len(freq)):
#        for j in range(len(freq)):
#            mat_per_bin[i][j] = round(mat_per_bin[i][j])

#    print 'original distribution :  \n',freq
#    print '\ntotal materials       :   ', freq.sum()
#    print '\n\nbinned distribution   :   \n', mat_per_bin
#    print '\ntotal materials       :   ', mat_per_bin.sum()

 #   weights = mat_per_bin / mat_per_bin.sum()

#    print 'weights :   \n', weights
#    k = 0
#    w_list = []
#ID_list = np.zeros([25,])
#    ID_list = []
#    for i in range(len(weights)):
#        w_list = np.concatenate([w_list, weights[i,:]])
#        for j in range(len(weights)):
        #ID_list[k] = IDs[i][j]
#            ID_list = ID_list + [IDs[i][j]]
#            k = k + 1
    #ID_list = np.concatenate([ID_list, IDs[i,:]])

    #number_of_children = 100

    for i in range(number_of_children):
    #print '\n\n cycle :   ', i
        print 'making MAT-' + str(i) + '...'        

        ID_bin = choice(ID_list, p=w_list)
        parent_ID = choice(ID_bin)

        #print ID_bin
        #print parent_ID


#    print mat_per_bin, mat_per_bin.sum()

    #for i in range(len(IDs)):
     #   for j in range(len(IDs)):
      #      if len(IDs[i][j]) != 0:    
       #         for n in range(int(mat_per_bin[i,j])):
    #                parent_ID = choice(IDs[i][j])
    
        #            parent_ID = choice(IDs[i][j])
                    #parent_ID = IDs[0][0][index]
    
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
        #+ '.cif', mat_path + '/MAT-' + str(k) + '.cif') 
        
           #         print 'creating MAT-'+ str(k) + '...'
        
        # importing values from parent data files...
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
        

        # perturbing number density...
        n_o = len(x_o)
        nden_o = n_o / (a_o * b_o * c_o)
                    
        nden_r = randrange(ndenmin * 10.**9, ndenmax * 10.**9) / 10.**9
                    
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
    
        #k = k + 1
        cif_file.close() 
        mixing_file.close()

