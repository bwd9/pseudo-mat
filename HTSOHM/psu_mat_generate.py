from random import choice, random, randrange
from math import fsum
import os
import numpy as np
# import functools need in Python 3 or greater

#define the greatest common divisor for use in charge calculation
def GCD(a,b):
	while b:
		a, b = b, a % b
	return a
	
#define the lowest common multiple (two args) for use in charge calculation
def LCM(a,b):
	return a * b // GCD(a,b)

#define lowest common multiple (multiple args) for use in charge calculation
def LCMM(*args):
	return reduce(LCM, args)

#main function (two inputs: N (number of materials) and ATOM_TYPES (number of atom types)	
def generate(
N, ATOM_TYPES, 
ndenmax=0.004302, ndenmin=0.000013905, 
xmax=52.3940, xmin=14.8193, 
ymax=52.3940, ymin=14.8193, 
zmax=52.3940, zmin=14.8193,    #unit cell parameters
epmax=410.6115, epmin=2.0128,  #LJ potentials (epsilon and sigma)
sigmax=5.2394, sigmin=1.6834,  #maximum charge (change to 3.0 for Blake's simulations
qmax=6.0,
elem_charge=0.0001):           #define an elementary charge to control sig figs

#epmax DEFINED WRT TO X-Y-Z LIMITS?
#max number density based on that of pure Iron
#max unit cell dimensions based on PCN-777 cages size
#max LJ parameters (for now using 1.5x highest values in GenericMOFs)
#max charge... UFF?

#    ATOM_TYPES = 10

    if type(N) != int:
        print 'N must be an integer.'
   
    Ntag = str(N)
    ntag = str(ndenmax)
    xtag = str(xmax)
    ytag = str(xmax)
    ztag = str(xmax)
    eptag = str(xmax)
    sigtag = str(xmax)
    qtag = str(xmax)
    
    top_path = (Ntag + 'mat_' + str(ATOM_TYPES) + 'atmtyp')
   
    if not os.path.exists(top_path):
        os.mkdir(top_path) 

    #open mat_stats.txt, to track material data    
    mat_stats = open(os.path.abspath(top_path)+ '/mat_stats.txt', 'w')
    mat_stat_heading = ('\nBOUNDARIES\nNumber of particles: ' + Ntag +
                    	'\nnumber density:   ' + ntag + '\nx-coordinate: ' +
			xtag + '\ny-coordinate: ' + ytag + '\nz-coordinate: ' +
			 ztag + '\nEpsilon: ' + eptag + '\nSigma: ' + sigtag 
			+ '\nCharge: ' + qtag + '\n\n' +
			'#name     number density     xdim     ydim     '+
			'zdim     total particles     net charge\n')
    mat_stats.write(mat_stat_heading)
    
    #MAT-XXX loop...
    for i in range(N + 1):
       
        mat_name = 'MAT-' + str(i)

 
	#make MAT-XXX directory
	os.mkdir(top_path+'/'+mat_name)
	
	#open .cif file
        cif_file = open(os.path.abspath(top_path) + '/'+mat_name + '/' + 
			mat_name+'.cif', 'w')
	
	#open force_field_mixing_rules.def
	mixing_rules = open(os.path.abspath(top_path) + '/'+mat_name +
			'/force_field_mixing_rules.def', 'w')
        
	#open pseudo_atoms.def
        pseudo_atoms = open(os.path.abspath(top_path) + '/'+mat_name + 
			'/pseudo_atoms.def', 'w')
	
	#open force_field.def
        force_field = open(os.path.abspath(top_path) + '/'+mat_name +
			'/force_field.def', 'w')

	xdim_ = round(random() * (xmax - xmin) + xmin, 4)
        ydim_ = round(random() * (ymax - ymin) + ymin, 4)
        zdim_ = round(random() * (zmax - zmin) + zmin, 4)
        
        Nmax = int(ndenmax * xdim_ * ydim_ * zdim_)
        n_ = randrange(2, Nmax, 1)
        nden_ = round(n_ / (xdim_ * ydim_ * zdim_))
        #N_ = xdim_ * ydim_ * zdim_ * nden_
        #n_ = int(N_)

        cif_heading = ('material' + str(i) + 
			'\n\nloop_\n' +
			'_symmetry_equiv_pos_as_xyz\n' +
			'  x,y,z\n' +
			'_cell_length_a          ' + str(xdim_) +
			'\n_cell_length_b          ' + str(ydim_) +
			'\n_cell_length_c          ' + str(zdim_) + 
			'\n_cell_angle_alpha       90.0000\n' +
			'_cell_angle_beta        90.0000\n' +
			'_cell_angle_gamma       90.0000\n' +
			'loop_\n' +
			'_atom_site_label\n' +
			'_atom_site_type_symbol\n' +
			'_atom_site_fract_x\n' +
			'_atom_site_fract_y\n' +
			'_atom_site_fract_z\n')
	cif_file.write(cif_heading)

        mixing_heading = ('# general rule for shifted vs truncated\n' +
                          'shifted\n' +
                          '# general rule tailcorrections\n' +
                          'no\n' +
                          '# number of defined interactions\n' +
                          str(ATOM_TYPES + 10) + '\n' +
                          '# type interaction, parameters.    IMPORTANT: define shortest matches first, so that more specific ones overwrites these\n')
        mixing_rules.write(mixing_heading)
        
        pseudo_heading = ('#number of pseudo atoms\n' + str(ATOM_TYPES + 10) + 
			'\n#type          print    as     chem     oxidation' +
			'     mass       charge     polarization     ' +
			'B-factor     radii    connectivity     anisotropic' +
			'   anisotrop-type  tinker-type\n')
        pseudo_atoms.write(pseudo_heading)
        
        #LJ parameters

	ep = []
	sig = []
        q = []
        for i in range(ATOM_TYPES):
            epsilon = round(random() * (epmax - epmin) + epmin, 4)
            ep.append(epsilon)
            sigma = round(random() * (sigmax -sigmin) + sigmin, 4)
            sig.append(sigma)
            charge = 0
            q.append(charge)
       
        ep_ = np.asarray(ep)
        sig_ = np.asarray(sig)
        q_ = np.asarray(q)
        ID_ = np.asarray(range(0,ATOM_TYPES))

        ep = ep_.reshape(-1,1)
        sig = sig_.reshape(-1,1)
        q = q_.reshape(-1,1)
        ID = ID_.reshape(-1,1)

        atoms = np.hstack((ID, ep, sig, q))


	#Begin charge calculation
        n_atoms = np.empty([0, 4])
        for i in range(n_):
            atomtype = choice(range(ATOM_TYPES))
            n_atoms = np.vstack([n_atoms, atoms[atomtype, :]])
    
        IDs = n_atoms[:,0]
		count = []
		for i in range(ATOM_TYPES):
            if i in IDs:
                count_i = list(IDs).count(i)
				count.append(count_i)
		
		temp_LCM = LCMM(*count)
		
		cm_max = floor(qmax/(temp_LCM*elem_charge/min(count)))
		
		cm_list = []
		for i in range(ATOM_TYPES - 1):
			cm_i = np.randrange(-1 * cm_max, cm_max, 1)
			cm_list.append(cm_i)
			
			atoms[i, 3] = cm_i * temp_LCM * elem_charge / count[i]
		
		cm_f = -1 * sum(cm_list)
		atoms[ATOM_TYPES, 3] = cm_f * temp_LCM * elem_charge / count[ATOM_TYPES]
	
	#old charging method	
#        for i in range(ATOM_TYPES):
#@            if i in IDs:
#               charge = round(random() * (qmax - qmin) + qmin, 4)
                # weight_i = list(IDs).count(i)
                # k = choice(IDs)
                # weight_k = list(IDs).count(k)
                # for j in range(n_):
                    # if n_atoms[j,0] == i:
                        # n_atoms[j,3] = n_atoms[j,3] + charge * int(weight_k)
			# atoms[i,3] = n_atoms[j,3] + charge * int(weight_k)
                    # if n_atoms[j,0] == k:
                        # n_atoms[j,3] = n_atoms[j,3] - charge * int(weight_i)
                        # atoms[k,3] = n_atoms[j,3] - charge * int(weight_i)

        mat_charge = str(sum(n_atoms[:,3]))
	cif_file.write('#NET CHARGE: ' + mat_charge + '\n')
	mat_X_stats = (mat_name + '     ' + str(nden_) + '     ' + str(xdim_) + '     ' + str(ydim_) +
			'     ' + str(zdim_) + '     ' + str(n_) + '     ' + 
			mat_charge + '\n')
	mat_stats.write(mat_X_stats)
	
        eps = n_atoms[:,1]
        sigs = n_atoms[:,2]
        qs = n_atoms[:,3]
	
	#writing mixing_rules, pseudo_atoms...
        for i in range(ATOM_TYPES):

            atom_X_pseudo = ('A_' + str(int(atoms[i,0])) + '   yes   C   C   0   ' +
                             '12.0   ' + str(atoms[i,3]) + '   0.0   0.0   ' +
                             '1.0  1.00   0   0  absolute   0\n')
            pseudo_atoms.write(atom_X_pseudo)

            atom_X_mixing = ('A_' + str(int(atoms[i,0])) + ' ' +
                             'lennard-jones ' + str(atoms[i,1]) + ' '
                             + str(atoms[i,2]) + '\n')
            mixing_rules.write(atom_X_mixing)                    

	#writing cif...

        for i in range(n_):
#FIX THIS TO ALLOW FOR NON-INT VALUES?
            x = round(random(), 4)
            y = round(random(), 4)
            z = round(random(), 4)
            
            atom_X_cif = ('A_' + str(int(n_atoms[i,0])) + '     ' + 'C     ' + 
                          str(x) + '     ' + str(y) + '     ' + str(z) + '\n')    
            cif_file.write(atom_X_cif)
 
#SUPPORTED ADSORBATES
# name         pseudo-atoms
# N2       :   N_n2; N_com
# CO2      :   C_co2; O_co2
# methane  :   CH4_sp3
# helium   :   He
# hydrogen :   H_h2; H_com
# H2       :   H_h2; H_com
# O2        :   O_o2; O_com


        adsorbate_mixing = ('N_n2        lennard-jones   36.0     3.31\n' +
                            'N_com       none\n' +
                            'C_co2       lennard-jones   27.0     2.80\n' +
                            'O_co2       lennard-jones   79.0     3.05\n' +
                            'CH4_sp3     lennard-jones   158.5    3.72\n' +
                            'He          lennard-jones   10.9     2.64\n' +
                            'H_h2        none\n' +
                            'H_com       lennard-jones   36.7     2.958\n' +
                            'O_o2        lennard-jones   49.0     3.02\n' +
                            'O_com       none\n' +
                            '# general mixing rule for Lennard-Jones\n' +
                            'Lorentz-Berthelot')
        mixing_rules.write(adsorbate_mixing)

        adsorbate_pseudo = ('N_n2     yes   N   N   0   14.00674   -0.4048' +
			'   0.0   1.0   0.7   0   0   relative   0\n' +
			'N_com    no    N   -   0   0.0         0.8096' +
			'   0.0   1.0   0.7   0   0   relative   0\n' +
			'C_co2    yes   C   C   0   12.0        0.70' +
			'     0.0   1.0   0.720 0   0   relative   0\n' +
			'O_co2    yes   O   O   0   15.9994    -0.35' +
			'     0.0   1.0   0.68  0   0   relative   0\n' +
			'CH4_sp3  yes   C   C   0   16.04246    0.0' +
			'      0.0   1.0   1.00  0   0   relative   0\n' +
			'He       yes   He  He  0   4.002602    0.0' +
			'      0.0   1.0   1.0   0   0   relative   0\n' +
			'H_h2     yes   H   H   0   1.00794     0.468' +
			'    0.0   1.0   0.7   0   0   relative   0\n' +
			'H_com    no    H   H   0   0.0        - 0.936' +
			'   0.0   1.0   0.7   0   0   relative   0\n'
			'O_o2     yes   O   O   0   15.9994   -0.112' +
			   0.0  1.0   0.7   0   0   relative   0\n' +
			'O_com   no   O   O   0   0.0   0.224' +
			'   0.0   1.0   0.7   0   0   relative   0\n')
        pseudo_atoms.write(adsorbate_pseudo)
        
        force_field_rules = ('# rules to overwrite\n0\n' +
				'# number of defined interactions\n0\n' +
				'# mixing rules to overwrite\n0')
	force_field.write(force_field_rules)
        
    cif_file.close()
    mixing_rules.close()
    pseudo_atoms.close()
    force_field.close()
    mat_stats.close()
