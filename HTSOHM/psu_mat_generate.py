from random import choice, random, randrange
from math import fsum
import os
import numpy as np

def generate(N, ATOM_TYPES, ndenmax=0.04302, ndenmin=0.0000013905, xmax=51.2, xmin=25.6, ymax=51.2, ymin=25.6,
zmax=51.2, zmin=25.6, epmax=513.264, epmin=1.2580, sigmax=6.549291, sigmin=1.052342, qmax=0.0, qmin=0.0):
#epmax DEFINED WRT TO X-Y-Z LIMITS?
#max number density based on that of pure Iron
#max unit cell dimensions based on PCN-777 cages size
#max LJ parameters (for now using 1.5x highest values in GenericMOFs)
#max charge... UFF?

    #ATOM_TYPES = 4

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
    
    top_path = ('materials' + '_' + Ntag + '.' + ntag + '_' + xtag + '.' + ytag
		 + '.' + ztag + '_' + eptag + '.' + sigtag + '_' + qtag)
    
    if not os.path.exists(top_path):
        os.mkdir(top_path) 

#    def drange(start, stop, step):
#        r = start
#        while r < stop:
#            yield r
#            r+= step
    
#    nden0 = drange(1, ndenmax*10000, ndenp*10000)
#    ndendim = [nden for nden in nden0]
    
#    x0 = drange(0, xmax + xp, xp)
#    xdim = [x for x in x0]
    
#    y0 = drange(0, ymax + yp, yp)
#    ydim = [y for y in y0]
    
#    z0 = drange(0, zmax + zp, zp)
#    zdim = [z for z in z0]
    
#    ep0 = drange(0, epmax + epp, epp)
#    epdim = [ep for ep in ep0]
#    sig0 = drange(0, sigmax + sigp, sigp)
#    sigdim = [sig for sig in sig0]   


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

 	#nden_ = choice(ndendim)/10000.
        #xdim_ = choice(xdim)
        #ydim_ = choice(ydim)
        #zdim_ = choice(zdim)
        #nden_ = randrange(0.0001, ndenmax, 1)
        #xdim_ = randrange(15., xmax, 0.1)
	#ydim_ = randrange(15., ymax, 0.1)
        #zdim_ = randrange(15., zmax, 0.1)
        #N_ = xdim_ * ydim_ * zdim_ * nden_
        #n_ = int(N_)        
        nden_ = round(random() * (ndenmax - ndenmin) + ndenmin, 6)
	xdim_ = round(random() * (xmax - xmin) + xmin, 4)
        ydim_ = round(random() * (ymax - ymin) + ymin, 4)
        zdim_ = round(random() * (zmax - zmin) + zmin, 4)
        N_ = xdim_ * ydim_ * zdim_ * nden_
        n_ = int(N_)

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
			'_atom_site_fract_z\n' +
			'_atom_site_charge\n')
	cif_file.write(cif_heading)

#        mixing_heading = ('# general rule for shifted vs truncated\nshifted\n' +
#			'# general rule for tailcorrections\nno\n' +
#			'# number of defined interactions\n' + str(108) +  #check these + XXX values
#			'\n# type interaction\n')

        mixing_heading = ('# general rule for shifted vs truncated\n' +
                          'shifted\n' +
                          '# general rule tailcorrections\n' +
                          'no\n' +
                          '# number of defined interactions\n' +
                          str(ATOM_TYPES + 8) + '\n' +
                          '# type interaction, parameters.    IMPORTANT: define shortest matches first, so that more specific ones overwrites these\n')
        mixing_rules.write(mixing_heading)
        
        pseudo_heading = ('#number of pseudo atoms\n' + str(ATOM_TYPES + 8) + 
			'\n#type          print    as     chem     oxidation' +
			'     mass       charge     polarization     ' +
			'B-factor     radii    connectivity     anisotropic' +
			'   anisotrop-type  tinker-type\n')
        pseudo_atoms.write(pseudo_heading)
        
        ##make charges
        #q = []
       	#for k in range(n_ + 1):
        #    q.append(0)
        #for l in range(5*(n_ + 1)):
        #    m = choice(range(n_ + 1))
        #    n = choice(range(n_ + 1))
        #    if m == n:
        #        n = choice(range(n_ + 1))
        #    dq = random() * qmax
        #    if q[m] + dq <= qmax and q[n] - dq >= -1 * qmax:
        #        q[m] = float(float(q[m]) + dq)
        #        q[n] = float(float(q[n]) - dq)
        #    if q[m] > qmax or q[n] < -1 * qmax:
        #        q[m] = q[m] - dq
        #        q[n] = q[n] + dq
        #for o in range(5*(n_ + 1)):
        #    m = choice(range(n_ + 1))
        #    n = choice(range(n_ + 1))
        #    if m == n:
        #        n = choice(range(n_ + 1))
        #    dq = random() * qmax
        #    if q[m] + dq <= qmax and q[n] - dq >= -1 * qmax:
        #        q[m] = float(float(q[m]) + dq)
        #        q[n] = float(float(q[n]) - dq)
        #    if q[m] > qmax or q[n] < -1 * qmax:
        #        q[m] = q[m] - dq
        #        q[n] = q[n] + dq
        #p = choice(range(n_ + 1))
        #q[p] = q[p] - sum(q)
        #if sum(q) != 0.000000000000000000000:
        #    for l in range(5*(n_ + 1)):
        #        m = choice(range(n_ + 1))
        #        n = choice(range(n_ + 1))
        #        if m == n:
        #            n = choice(range(n_ + 1))
        #        dq = random() * qmax
        #        if q[m] + dq <= qmax and q[n] - dq >= -1 * qmax:
        #            q[m] = float(float(q[m]) + dq)
        #            q[n] = float(float(q[n]) - dq)
        #        if q[m] > qmax or q[n] < -1 * qmax:
        #            q[m] = q[m] - dq
        #            q[n] = q[n] + dq
        #    for o in range(5*(n_ + 1)):
        #        m = choice(range(n_ + 1))
        #        n = choice(range(n_ + 1))
        #        if m == n:
        #            n = choice(range(n_ + 1))
        #        dq = random() * qmax
        #        if q[m] + dq <= qmax and q[n] - dq >= -1 * qmax:
        #            q[m] = float(float(q[m]) + dq)
        #            q[n] = float(float(q[n]) - dq)
        #        if q[m] > qmax or q[n] < -1 * qmax:
        #            q[m] = q[m] - dq
        #            q[n] = q[n] + dq
        #        p = choice(range(n_ + 1))
        #        q[p] = q[p] - sum(q)
	
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

        n_atoms = np.empty([0, 4])
        for i in range(n_):
            atomtype = choice(range(ATOM_TYPES))
            n_atoms = np.vstack([n_atoms, atoms[atomtype, :]])
    
        IDs = n_atoms[:,0]

        for i in range(ATOM_TYPES):
            if i in IDs:
                charge = round(random() * (qmax - qmin) + qmin, 4)
                weight_i = list(IDs).count(i)
                k = choice(IDs)
                weight_k = list(IDs).count(k)
                for j in range(n_):
                    if n_atoms[j,0] == i:
                        n_atoms[j,3] = n_atoms[j,3] + charge * int(weight_k)
			atoms[i,3] = n_atoms[j,3] + charge * int(weight_k)
                    if n_atoms[j,0] == k:
                        n_atoms[j,3] = n_atoms[j,3] - charge * int(weight_i)
                        atoms[k,3] = n_atoms[j,3] - charge * int(weight_i)

#        for i in range(100):
#            atoms[i,3] = round(atoms[i,3], 4)

#        for i in range(n_):
#            n_atoms[i,3] = round(n_atoms[i,3], 4)



#        net_charge = sum(n_atoms[:,3])
#        if net_charge != 0:
#            atomID = choice(range(100))
#            weight = list(IDs).count(atomID)
#            atoms[atomID,3] = atoms[atomID,3] - net_charge/weight
#            for i in range(n_):
#                if n_atoms[i,0] == atomID:
#                    n_atoms[atomID,3] = n_atoms[atomID,3] - net_charge/weight


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
            x = choice(range(int(xdim_ + 1)))
            y = choice(range(int(ydim_ + 1)))
            z = choice(range(int(zdim_ + 1)))
            
            atom_X_cif = ('A_' + str(int(n_atoms[i,0])) + '     ' + 'C     ' + 
                          str(round(x/xdim_, 4)) + '     ' + str(round(y/ydim_, 4)) + 
                          '     ' + str(round(z/zdim_, 4)) + '    ' +
                          str(n_atoms[i,3]) + '\n')    
            cif_file.write(atom_X_cif)



        #    #ep = choice(epdim)
        #    #sig = choice(sigdim)
        #    epval = ep[atomtype]
        #    sigval = sig[atomtype]
        #    charge = q[n_]
        #    #if charge < 0:
        #    atom_X_cif = ('A' + str(atomtype) + '     ' + 'C     ' + 
	#			str(x/xdim_) + '     ' + str(y/ydim_) + 
	#			'     ' + str(z/zdim_) + '    ' +
	#			str(charge) + '\n')    
        #    cif_file.write(atom_X_cif)
        #    for k in range(100):
        #        if k != atomtype:
        #            atom_X_pseudo = ('A' + str(k) + '   yes   C   C   0   12.0   0' +
	#			     '   0.0   0.0   1.0  1.00   0   ' +
	#			     '0  absolute   0\n')
        #        if k == atomtype:
        #            atom_X_pseudo = ('A' + str(k) + '   yes   C   C   0   12.0   ' +
	#			     str(q[n_]) + '   0.0   0.0   1.0  1.00   0   ' +
	#			     '0  absolute   0\n')
        #            
        #            pseudo_atoms.write(atom_X_pseudo)
        #             
        #            atom_X_mixing = ('A' + str(k) + '          LENNARD_JONES     ' +
	#			     str(ep[k]) + '     ' + str(sig[k]) + '\n')
        #            mixing_rules.write(atom_X_mixing)

 

            #if charge >= 0:
            #    atom_X_cif = ('A' + str(atomtype) + '     ' + str(x) + '     ' +
		#		str(y) + '     ' + str(z) + '     ' +
		#		str(charge) + '\n')
		#cif_file.write(atom_X_cif)
             	
        #for i in range(100):

         #   atom_X_mixing = ('A' + str(i) + '          LENNARD_JONES     ' +
	#			str(ep[i]) + '     ' + str(sig[i]) + '\n')
         #   mixing_rules.write(atom_X_mixing)
#
 #           atom_X_pseudo = ('A' + str(i) + '   yes   C   C   0   12.0   ' +
#				str(q[i]) + '   0.0   0.0   1.0  1.00   0   ' +
#				'0  absolute   0\n')
 ##           pseudo_atoms.write(atom_X_pseudo)
 
#SUPPORTED ADSORBATES
# name         pseudo-atoms
# N2       :   N_n2; N_com
# CO2      :   C_co2; O_co2
# methane  :   CH4_sp3
# helium   :   He
# hydrogen :   H_h2; H_com
# H2       :   H_h2; H_com

        #adsorbate_mixing = ('N_n2        LENNARD_JONES   36.0     3.31\n' +
	#		'N_com       none\n' +
	#		'C_co2       LENNARD_JONES   27.0     2.80\n' +
	#		'O_co2       LENNARD_JONES   79.0     3.05\n' +
	#		'CH4_sp3     LENNARD_JONES   158.5    3.72\n' +
	#		'He          LENNARD_JONES   10.9     2.64\n' +
	#		'H_h2        none\n' +
	#		'H_com       LENNARD_JONES   36.7     2.958\n' +
	#		'# general mixing rule for Lennard-Jones\n' +
	#		'Lorentz-Berthlot')
        adsorbate_mixing = ('N_n2        lennard-jones   36.0     3.31\n' +
                            'N_com       none\n' +
                            'C_co2       lennard-jones   27.0     2.80\n' +
                            'O_co2       lennard-jones   79.0     3.05\n' +
                            'CH4_sp3     lennard-jones   158.5    3.72\n' +
                            'He          lennard-jones   10.9     2.64\n' +
                            'H_h2        none\n' +
                            'H_com       lennard-jones   36.7     2.958\n' +
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
