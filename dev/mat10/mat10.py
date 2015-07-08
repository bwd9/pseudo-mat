from random import choice, random
import os
def mat10(N, ndenmax, gas, ndenp=0.01, xmax=20, xp=0.01, ymax=20, yp=0.01, zmax=20, zp=0.01, epmax=1, epp=0.01, sigmax=1, sigp=0.01, qmax=1, cycles=500, initcycles=100, printevery=100, unitcells=(2,2,2), temperature=298, pressure='1e5 2.5e5 5e5 1e6 2.5e6 5e6'):
    #
    #_max refers to the maximum upper booundary for a given paramter, minimum boundaries are zero, except  
    #in the case of charges where the minimum boundary is qmax * -1
    # 
    #_p refers to the precision for a given parameter. acceptable parameter values are determined by the upper boundary
    #and precison; that is acceptable values are all within the upper and lower boundaries subdivided by the precision.
    #
    #N [=] number of hypothetical materials
    #n [=] number of particles in material
    #x [=] x coordinate
    #y [=] y coordinate
    #z [=] z coordinate
    #ep [=] epsilon, LJ parameter=0.l
    #sig [=] sigma, LJ parameter
    #q [=] charge
    #
    #UNITS???
    #
    if type(N) != int:
        print 'N must be an integer.'
    if type(gas) != str:
        print 'please input \'gas\' as \'string\''
    if gas == 'ch4':
        gas = 'methane'
    if gas == 'n2':
        gas = 'N2'
    #if gas != 'co2' or 'n2' or 'methane':
    #    print 'undefined adsorbant; please choose \'co2\', \'n2\', or \'ch4\'/\'methane\'' 
    #^^fix this error message
    Ntag = str(N)
    ntag = str(ndenmax)
    xtag = str(xmax)
    ytag = str(xmax)
    ztag = str(xmax)
    eptag = str(xmax)
    sigtag = str(xmax)
    qtag = str(xmax)
    run_path = 'materials' + '_' + Ntag + '.' + ntag + '_' + xtag + '.' + ytag + '.' + ztag + '_' + eptag + '.' + sigtag + '_' + qtag
    if not os.path.exists(run_path):
        os.mkdir(run_path) 
    def drange(start, stop, step):
        r = start
        while r < stop:
            yield r
            r+= step
    
    nden0 = drange(1, ndenmax, ndenp)
    ndendim = [nden for nden in nden0]
    x0 = drange(0, xmax + xp, xp)
    xdim = [x for x in x0]
    y0 = drange(0, ymax + yp, yp)
    ydim = [y for y in y0]
    z0 = drange(0, zmax + zp, zp)
    zdim = [z for z in z0]
    ep0 = drange(0, epmax + epp, epp)
    epdim = [ep for ep in ep0]
    sig0 = drange(0, sigmax + sigp, sigp)
    sigdim = [sig for sig in sig0]
    runshell_file = open(os.path.join(os.path.abspath(run_path), 'run_shell'), 'wb')
    runshell_file.write('echo \"...\"\nwhile [ -f \"$(find ./ -name shell -type f)\" ];\ndo\n    cd -- \"$(find ./ -name shell -type f -printf \'%h\' -quit)\"\n    source shell\ndone\necho \"...done!\"')
    
    
  #  numden_file = open(os.path.join(os.path.abspath(run_path), 'numden'), 'wb')
    runshell_file.close()

    for i in range(N + 1):
	
        mat_path = 'MAT-' + str(i)
        mat_path_name = 'Adsorption_of_' + gas + '_in_' + mat_path
        #clean this^ up later
        os.mkdir(os.path.join(run_path, mat_path_name))
	
        cif_filename = mat_path + '.cif'
        mixing_filename = 'force_field_mixing_rules.def'
        pseudo_filename = 'pseudo_atoms.def'
        run_filename = 'run'
        sim_filename = 'simulation.input'
        shell_filename = 'shell'
        force_filename = 'force_field.def'
	
        sim_output = open(os.path.join(os.path.abspath(run_path)+'/Adsorption_of_' + gas + '_in_' +'MAT-'+str(i), sim_filename), 'wb')
        cif_output = open(os.path.join(os.path.abspath(run_path)+'/Adsorption_of_' + gas + '_in_' +'MAT-'+str(i), cif_filename), 'wb')
        mixing_output = open(os.path.join(os.path.abspath(run_path)+'/Adsorption_of_' + gas + '_in_' +'MAT-'+str(i), mixing_filename), 'wb')
        pseudo_output = open(os.path.join(os.path.abspath(run_path)+'/Adsorption_of_' + gas + '_in_' +'MAT-'+str(i), pseudo_filename), 'wb')
        run_output = open(os.path.join(os.path.abspath(run_path)+'/Adsorption_of_' + gas + '_in_' +'MAT-'+str(i), run_filename), 'wb')
        force_output = open(os.path.join(os.path.abspath(run_path)+'/Adsorption_of_' + gas + '_in_' +'MAT-'+str(i), force_filename), 'wb')
        shell_output = open(os.path.join(os.path.abspath(run_path)+'/Adsorption_of_' + gas + '_in_' +'MAT-'+str(i), shell_filename), 'wb')
	
        cif_header = 'material' + str(i) + '\n\nBOUNDARIES\nNumber of particles: ' + ntag + '\nx-coordinate: ' + xtag + '\ny-coordinate: ' + ytag + '\nz-coordinate: ' + ztag + '\nEpsilon: ' + eptag + '\nSigma: ' + sigtag + '\nCharge: ' + qtag + '\n\nloop_\n_symmetry_equiv_pos_as_xyz\n  x,y,z\n_cell_length_a          ' + xtag + '\n_cell_length_b          ' + ytag + '\n_cell_length_c          ' + ztag + '\n_cell_angle_alpha       90.0000\n_cell_angle_beta        90.0000\n_cell_angle_gamma       90.0000\nloop_\n_atom_site_label\n_atom_site_fract_x\n_atom_site_fract_y\n_atom_site_fract_z\n_atom_site_charge\n'
	
        ndendim_ = choice(ndendim)
        xdim_ = choice(xdim)
        ydim_ = choice(ydim)
        zdim_ = choice(zdim)
 #       n_mat = int(xdim_ * ydim_ * zdim_ / ndendim_)
	densityinfo = 'mat '+str(i)+'  den  :   '+str(ndendim_)+'   x    :   '+str(xdim_)+'   y    :   '+str(ydim_)+'   z    :   '+str(zdim_)+'   vol  :   '+str(xdim_*ydim_*zdim_)+'   n    :   '+str(n_mat)+'\n'
#	numden_file.write(densityinfo)	

        mixing_header = '# general rule for shifted vs truncated\nshifted\n# general rule for tailcorrections\nno\n# number of defined interactions\n' + str(n_mat+4) + '\n# type interaction\n'
        cif_output.write(cif_header)
        mixing_output.write(mixing_header)
        
        pseudo_header = '#number of pseudo atoms\n' + str(n_mat+9) + '\n#type          print    as     chem     oxidation     mass       charge     polarization     B-factor     radii    connectivity     anisotropic   anisotrop-type  tinker-type\n'
        pseudo_output.write(pseudo_header)
        	
        sim_output.write('SimulationType                     MonteCarlo\nNumberofCycles                     ' + str(cycles) + '\nNumberofInitializationCycles       ' + str(initcycles) + '\nPrintEvery                         ' + str(printevery) + '\n\nForcefield                         GenericMOFs\nChargeMethod                       EWald\nCutOff                             15.0\n\nFramework  0\nFrameworkName  ' + str(mat_path) + '\nUnitCells  ' + str(unitcells[0]) + ' ' + str(unitcells[1]) + ' ' + str(unitcells[2]) + '\nExternalTemperature  ' + str(temperature) + '\nExternalPressure  ' + str(pressure) + '\n\nComponent 0  MoleculeName                      ' + str(gas) + '\n             MoleculeDefinition                TraPPE\n             TranslationProbability            1.0\n             ReinsertionProbability            1.0\n             SwapProbability                   1.0\n             CreateNumberOfMolecules           0')
        
        #make charges
        q = []
        for k in range(n_mat + 1):
            q.append(0)
        for l in range(5*(n_mat + 1)):
            m = choice(range(n_mat + 1))
            n = choice(range(n_mat + 1))
            if m == n:
                n = choice(range(n_mat + 1))
            dq = random() * qmax
            if q[m] + dq <= qmax and q[n] - dq >= -1 * qmax:
                q[m] = float(float(q[m]) + dq)
                q[n] = float(float(q[n]) - dq)
            if q[m] > qmax or q[n] < -1 * qmax:
                q[m] = q[m] - dq
                q[n] = q[n] + dq
        for o in range(5*(n_mat + 1)):
            m = choice(range(n_mat + 1))
            n = choice(range(n_mat + 1))
            if m == n:
                n = choice(range(n_mat + 1))
            dq = random() * qmax
            if q[m] + dq <= qmax and q[n] - dq >= -1 * qmax:
                q[m] = float(float(q[m]) + dq)
                q[n] = float(float(q[n]) - dq)
            if q[m] > qmax or q[n] < -1 * qmax:
                q[m] = q[m] - dq
                q[n] = q[n] + dq
        p = choice(range(n_mat + 1))
        q[p] = q[p] - sum(q)
        if sum(q) != 0.000000000000000000000:
            for l in range(5*(n_mat + 1)):
                m = choice(range(n_mat + 1))
                n = choice(range(n_mat + 1))
                if m == n:
                    n = choice(range(n_mat + 1))
                dq = random() * qmax
                if q[m] + dq <= qmax and q[n] - dq >= -1 * qmax:
                    q[m] = float(float(q[m]) + dq)
                    q[n] = float(float(q[n]) - dq)
                if q[m] > qmax or q[n] < -1 * qmax:
                    q[m] = q[m] - dq
                    q[n] = q[n] + dq
            for o in range(5*(n_mat + 1)):
                m = choice(range(n_mat + 1))
                n = choice(range(n_mat + 1))
                if m == n:
                    n = choice(range(n_mat + 1))
                dq = random() * qmax
                if q[m] + dq <= qmax and q[n] - dq >= -1 * qmax:
                    q[m] = float(float(q[m]) + dq)
                    q[n] = float(float(q[n]) - dq)
                if q[m] > qmax or q[n] < -1 * qmax:
                    q[m] = q[m] - dq
                    q[n] = q[n] + dq
                p = choice(range(n_mat + 1))
                q[p] = q[p] - sum(q)
                if sum(q) != 0.000000000000000000000:
                    for l in range(5*(n_mat + 1)):
                        m = choice(range(n_mat + 1))
                        n = choice(range(n_mat + 1))
                        if m == n:
                            n = choice(range(n_mat + 1))
                        dq = random() * qmax
                        if q[m] + dq <= qmax and q[n] - dq >= -1 * qmax:
                            q[m] = float(float(q[m]) + dq)
                            q[n] = float(float(q[n]) - dq)
                        if q[m] > qmax or q[n] < -1 * qmax:
                            q[m] = q[m] - dq
                            q[n] = q[n] + dq
                    for o in range(5*(n_mat + 1)):
                        m = choice(range(n_mat + 1))
                        n = choice(range(n_mat + 1))
                        if m == n:
                            n = choice(range(n_mat + 1))
                        dq = random() * qmax
                        if q[m] + dq <= qmax and q[n] - dq >= -1 * qmax:
                            q[m] = float(float(q[m]) + dq)
                            q[n] = float(float(q[n]) - dq)
                        if q[m] > qmax or q[n] < -1 * qmax:
                            q[m] = q[m] - dq
                            q[n] = q[n] + dq
                    p = choice(range(n_mat + 1))
                    q[p] = q[p] - sum(q)
        cif_output.write('#NET CHARGE: ' + str(sum(q)) + '\n')
	#WRITE PARAMETERS TO FILE
        for j in range(n_mat+1):
#FIX THIS TO ALLOW FOR NON-INT VALUES
            x = choice(range(int(xdim_+1)))
            y = choice(range(int(ydim_+1)))
            z = choice(range(int(zdim_+1)))
            ep = choice(epdim)
            sig = choice(sigdim)
            charge = q[j]
            if charge < 0:
                cif_output.write('A' + str(j) + '     ' + str(x) + '     ' + str(y) + '     ' + str(z) + '    ' + str(charge) + '\n')    
            if charge >= 0:
                cif_output.write('A' + str(j) + '     ' + str(x) + '     ' + str(y) + '     ' + str(z) + '     ' + str(charge) + '\n')
            #also account for values with fewer sigfigs???
            mixing_output.write('A' + str(j) + '          LENNARD_JONES     ' + str(sig) + '     ' + str(ep) + '\n')
            pseudo_output.write('A' + str(j) + '   yes   C   C   0   12.0   ' + str(q[j]) + '   0.0   0.0   1.0  1.00   0   0  absolute   0')
        mixing_output.write('CH4_sp3      LENNARD_JONES          148.0          3.73\nO_co2       LENNARD_JONES          79.0          3.05\nC_co2       LENNARD_JONES          27.0          2.80\nN_n2        LENNARD_JONES          36.0          3.31\n# general mixing rule for Lennard-Jones\nLorentz-Berthelot')
        pseudo_output.write('CH4-sp3       yes      C      C        0             16.04246   0.0        0.0              1.0          1.00     0                0             absolute        0\nCH3_sp3       yes      C      C        0             15.03452   0.0        0.0              1.0          1.00     0                0             absolute        0\nCH2_sp3       yes      C      C        0             14.02658   0.0        0.0              1.0          1.00     0                0             absolute        0\nCH_sp3        yes      C      C        0             13.01864   0.0        0.0              1.0          1.00     0                0             absolute        0\nC_sp3         yes      C      C        0             12.0       0.0        0.0              1.0          1.00     0                0             absolute        0\nC_co2         yes      C      C        0             12.0107    0.7        1.508            1.0          0.720    0                0             absolute        0\nO_co2         yes      O      O        0             15.9994   -0.35       0.9475           1.0          0.68     0                0             relative        0\nN_n2          yes      N      N        0             14.00674  -0.482      0.0              1.0          0.7      0                0             relative        0\nN_com         no       N      N        0             0.0        0.964      0.0              1.0          0.7      0                0             relative        0\n')
        #if framework pseudo atoms are not created, consider changing line2 of pseudo_atoms.def 
        run_output.write('#! /bin/sh -f\nexport RASPA_DIR=${HOME}/RASPA-2.0-CF/\n$RASPA_DIR/src/simulate $1')
        force_output.write('# rules to overwrite\n0\n# number of defined interactions\n0\n# mixing rules to overwrite\n0')
        shell_output.write('source run\ncd Output/System_0\ngrep \'Average loading absolute\' $PATTERN * >> ~/RASPA-2.0-CF/mat6_output.txt\ngrep \'WARNING\' $PATTERN * >> ~/RASPA-2.0-CF/mat6_output.txt\ncd ..\ncd ..\nrm -rf Output\nrm shell\ncd ~/RASPA-2.0-CF')
        sim_output.close()
        cif_output.close()
        mixing_output.close()
        pseudo_output.close()
        run_output.close()
        force_output.close()
        shell_output.close()
#	numden_file.close()
# run file may vary depending on installation

#check what atoms are present in ch4 structure file
