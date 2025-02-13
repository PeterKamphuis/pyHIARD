#!/usr/local/bin/ python3
# -*- coding: utf-8 -*-
#This program is a pyhton script to create a data base of artificial galaxies.
#It first creates a model with tirific at a high resolution and then runs it through casa to get obtain a realistic observation.
# once a numerical list is set in length we can convert it to a numpy array in order to do operations faster.
# first we import numpy
from pyHIARD.constants import G_agc, H_0, c_kms, HI_rest_freq, c
from multiprocessing import get_context
from scipy.ndimage import gaussian_filter
from scipy import integrate
from scipy import interpolate
from pyHIARD.AGC.base_galaxies import Base_Galaxy
from astropy.io import fits

import copy
import numpy as np
import os
import shutil
import pyHIARD.common_functions as cf
import subprocess
import sys
import time
import warnings

with warnings.catch_warnings():
     warnings.simplefilter("ignore")
     import matplotlib
     matplotlib.use('pdf')
     import matplotlib.pyplot as plt
     import matplotlib.font_manager as mpl_fm
     from matplotlib.ticker import AutoMinorLocator


#Some errors
class TirificRunError(Exception):
    pass
#Some errors

class InputError(Exception):
    pass

class RunningError(Exception):
    pass
#------------------------------!!!!!!!!!!!!!!!!!!!!!!!!!!!!!Main for creating the AGC!!!!!!!!!!!!!!!!!!!!!----------------------


def AGC(cfg):
    '''Create realistic artificial galaxies.'''
    # Let's give an overview of the database that will be created
    print(
        f"We will create a database in the directory {cfg.general.main_directory}.\n")
    if cfg.agc.delete_existing:
        print("All previous models will be removed prior to the build. \n")
    else:
        print("We will retain previously build models. \n")
    if cfg.agc.inhomogenous:
        print(
            f"We will use {cfg.agc.corruption_method} as corruption method and inhomogeneities will be added.\n")
    else:
        print(
            f"We will use {cfg.agc.corruption_method} as corruption method.\n")

    if 'Base' in cfg.agc.variables_to_vary:
        print(
            f"We will produce the base galaxies: {','.join([str(e) for e in cfg.agc.base_galaxies])}.\n")
    else:
        print(
            f"We will use the following galaxies as base: {','.join([str(e) for e in cfg.agc.base_galaxies])}.\n")
    print("We will vary the following parameters")
    if 'Inclination' in cfg.agc.variables_to_vary:
        print(
            f"We vary the inclination with the following values: {','.join([str(e) for e in cfg.agc.inclination])}.\n")
    if 'PA' in cfg.agc.variables_to_vary:
        print(
            f"We vary the PA with the following values: {','.join([str(e) for e in cfg.agc.pa])}.\n")
    if 'Beams' in cfg.agc.variables_to_vary:
        print(
            f"We create the model with {','.join([str(e) for e in cfg.agc.beams])} beams across the major axis.\n")
    if 'Radial_Motions' in cfg.agc.variables_to_vary:
        print(
            f"Inject radial motions with speeds of {','.join([str(e) for e in cfg.agc.radial_motions])} km/s.\n")
    if 'Flare' in cfg.agc.variables_to_vary:
        print(f"We will swap the flaring of the base galaxies")
    if 'Arms' in cfg.agc.variables_to_vary:
        print(f"We will swap the inclusion of the arms in the base galaxies")
    if 'Bar' in cfg.agc.variables_to_vary:
        print(f"We will swap the inclusion of the bar in the base galaxies")
    if 'Mass' in cfg.agc.variables_to_vary:
        print(
            f"We add the following masses to each base set: {','.join([str(e) for e in cfg.agc.masses])}.\n")
    if 'Channelwidth' in cfg.agc.variables_to_vary:
        print(
            f"Varying the channel width with: {','.join([str(e) for e in cfg.agc.channelwidth])} km/s.\n")
    if 'SNR' in cfg.agc.variables_to_vary:
        print(
            f"Varying the signal to noise ratio with: {','.join([str(e) for e in cfg.agc.snr])}.\n")
    if 'Warp' in cfg.agc.variables_to_vary: 
        print(
            f"Varying the theta angle of the angular momentum vector with: {','.join([str(e[0]) for e in cfg.agc.warp])}.\n")
        print(
            f"Varying the phi angle of the angular momentum vector with: {','.join([str(e[1]) for e in cfg.agc.warp])}.\n")
    if 'Beam_Size' in cfg.agc.variables_to_vary:
        print(
            f"Varying the beam size with: {','.join([str(e) for e in cfg.agc.beam_size])}.\n")
   
    # Let's just make 1 catalogue with everything we need and adjust
    # the fitting program to have this one catalogue as input

    Catalogue = f"{cfg.general.main_directory}Output_AGC_Summary.txt"
    # If we are making new models we want to ensure this is a new file
    number_models = cf.get_created_models(Catalogue, cfg.agc.delete_existing)

    #Copy a fits file from the WHISP data base to use as template if it is not ther3 yet
    # Check for the existence of a template fits file
    templatethere = os.path.isfile(cfg.general.main_directory+'/Input.fits')
    #templatethere =False
    # If it doesn't exist copy it into the directory
    if not templatethere:
        create_template_fits(cfg.general.main_directory)


    # All this went well
    templatethere = os.path.isfile(cfg.general.main_directory+'/Input.fits')
    if templatethere:
        print(" The Input Template is found and all is well")
    else:
        print(" The Input Template is NOT found !!! ABORTING")
        sys.exit()

    # start a loop over the various base galaxies
    set_done = [1024]
    plot_ax = []

    # If we make new models delete everything in the directory
    if cfg.agc.delete_existing:
        masses_to_delete = []
        if 'Mass' in cfg.agc.variables_to_vary:
            #Because python is the dumbest language ever
            masses_to_delete = copy.deepcopy(cfg.agc.masses)
        for galaxy in cfg.agc.base_galaxies:
            if galaxy == 7:
                continue
            else:
                masses_to_delete.append(Base_Galaxy(galaxy).Mass)
        to_delete = f"rm -R {' '.join([f'{cfg.general.main_directory}Mass{x:.1e}-*rm*' for x in masses_to_delete])} {cfg.general.main_directory}Fractions_and_Masses.txt"
        print("All previous models of the requested base galaxies will be removed prior to the build. \n")
        print(f"The command {to_delete} will be run.")
        cfg.agc.delete_existing = cf.get_bool(
            "Are you sure you want to do this? (Yes/No, default=No): ", default=False)
        if cfg.agc.delete_existing:
            os.system(to_delete)

    colors = iter(plt.cm.rainbow(np.linspace(
        0, 1, len(cfg.agc.base_galaxies)+len(cfg.agc.masses))))
    max_rad = 0.
    All_Galaxies = []
    Gauss_Galaxies = []
    Casa_Galaxies = []
    created = []
    rc_counter=0
    next_ucmodel = 6
    next_casamodel = 5
    for base in range(len(cfg.agc.base_galaxies)):
        base_defined = False
        # We want to keep the center constant per base galaxy, for easy comparison as well as to be able to investigate how center determination is affected
        if cfg.agc.corruption_method == 'Gaussian':
            RAdeg = np.random.uniform()*360
            DECdeg = (np.arccos(2*np.random.uniform()-1)*(360./(2.*np.pi)))-90
        else:
            RAdeg = np.random.uniform()*360
            DECdeg = -60
            while DECdeg < -20.:
                DECdeg = (np.arccos(2*np.random.uniform()-1)
                          * (360./(2.*np.pi)))-90
        # From here we go into a loop to adjust variables over the bases
        for ix in range(len(cfg.agc.variables_to_vary)):
            if cfg.agc.variables_to_vary[ix] == 'Inclination': numloops = len(
                cfg.agc.inclination)
            elif cfg.agc.variables_to_vary[ix] == 'PA': numloops = len(cfg.agc.pa)
            elif cfg.agc.variables_to_vary[ix] in ['Flare', 'Arms', 'Bar', 'Base']: numloops = 1
            elif cfg.agc.variables_to_vary[ix] == 'Warp': numloops = len(cfg.agc.warp)
            elif cfg.agc.variables_to_vary[ix] == 'Beams': numloops = len(cfg.agc.beams)
            elif cfg.agc.variables_to_vary[ix] == 'SNR': numloops = len(cfg.agc.snr)
            elif cfg.agc.variables_to_vary[ix] == 'Channelwidth': numloops = len(cfg.agc.channelwidth)
            elif cfg.agc.variables_to_vary[ix] == 'Beam_Size': numloops = len(cfg.agc.beam_size)
            elif cfg.agc.variables_to_vary[ix] == 'Radial_Motions': numloops = len(cfg.agc.radial_motions)
            elif cfg.agc.variables_to_vary[ix] == 'Dispersion': numloops = len(cfg.agc.dispersion)
            elif cfg.agc.variables_to_vary[ix] == 'Mass': numloops = len(cfg.agc.masses)
            else:
                print("This is not a supported parameter")
                exit()
            for jx in range(numloops):
                if cfg.agc.base_galaxies[base] > 6:
                    if not base_defined:
                        Current_Galaxy = Base_Galaxy(
                            cfg.agc.base_galaxies[base])
                        Current_Galaxy_Base = copy.deepcopy(Current_Galaxy)
                        base_defined = True
                        to_delete = f"rm -R {' '.join([f'{cfg.general.main_directory}Mass{Current_Galaxy.Mass:.1e}-*rm*'])}"
                        print(
                            "All previous models of the requested base galaxy will be removed prior to the build. \n")
                        print(f"The command {to_delete} will be run.")
                        cfg.agc.delete_existing = cf.get_bool(
                            "Are you sure you want to do this? (Yes/No, default=No): ", default=False)
                        if cfg.agc.delete_existing:
                            os.system(to_delete)
                    else:
                        Current_Galaxy = copy.deepcopy(Current_Galaxy_Base)
                else:
                    Current_Galaxy = Base_Galaxy(cfg.agc.base_galaxies[base])
                if cfg.agc.variables_to_vary[ix] == 'Inclination':
                    Current_Galaxy.Inclination = cfg.agc.inclination[jx]
                elif cfg.agc.variables_to_vary[ix] == 'PA': 
                    Current_Galaxy.PA = cfg.agc.pa[jx]
                elif cfg.agc.variables_to_vary[ix] == 'Flare':
                    if Current_Galaxy.Flare == 'Flare':
                        Current_Galaxy.Flare = "No_Flare"
                    else:
                        Current_Galaxy.Flare = "Flare"
                elif cfg.agc.variables_to_vary[ix] == 'Warp': 
                    Current_Galaxy.Warp = [cfg.agc.warp[jx][0], cfg.agc.warp[jx][1]]
                elif cfg.agc.variables_to_vary[ix] == 'Beams': 
                    Current_Galaxy.Beams = cfg.agc.beams[jx]
                elif cfg.agc.variables_to_vary[ix] == 'SNR': 
                    Current_Galaxy.SNR = cfg.agc.snr[jx]
                    
                elif cfg.agc.variables_to_vary[ix] == 'Channelwidth': 
                    Current_Galaxy.Channelwidth = cfg.agc.channelwidth[jx]
                elif cfg.agc.variables_to_vary[ix] == 'Beam_Size': 
                    Current_Galaxy.Res_Beam = [cfg.agc.beam_size[jx][0], cfg.agc.beam_size[jx][1], cfg.agc.beam_size[jx][2]]
                elif cfg.agc.variables_to_vary[ix] == 'Arms':
                    if Current_Galaxy.Arms == 'Arms':
                        Current_Galaxy.Arms = "No_Arms"
                    else:
                        Current_Galaxy.Arms = "Arms"
                elif cfg.agc.variables_to_vary[ix] == 'Bar':
                    if Current_Galaxy.Bar == 'Bar':
                        Current_Galaxy.Bar = "No_Bar"
                    else:
                        Current_Galaxy.Bar = "Bar"
                elif cfg.agc.variables_to_vary[ix] == 'Radial_Motions':
                    Current_Galaxy.Radial_Motions = cfg.agc.radial_motions[jx]
                elif cfg.agc.variables_to_vary[ix] == 'Dispersion': 
                    Current_Galaxy.Dispersion = cfg.agc.dispersion[jx]
                elif cfg.agc.variables_to_vary[ix] == 'Mass': 
                    Current_Galaxy.Mass = cfg.agc.masses[jx]
                elif cfg.agc.variables_to_vary[ix] == 'Base': 
                    pass
                else: print("This is not a supported parameter")
                Current_Galaxy.Res_Beam[0:1] = np.sort(
                    Current_Galaxy.Res_Beam[0:1])
                setattr(Current_Galaxy, "Coord", [RAdeg, DECdeg])
                setattr(Current_Galaxy, "Model_Number", number_models)

                number_models += 1

                if (cfg.agc.corruption_method == 'Casa_5' and \
                    int(number_models/5.) == number_models/5.) or\
                    (cfg.agc.corruption_method == 'Tres' and \
                    number_models == next_casamodel)\
                    or (cfg.agc.corruption_method == 'Casa_Sim'):
                        setattr(Current_Galaxy, "Corruption", "Casa_Sim")
                        if cfg.agc.corruption_method == 'Tres':
                            if next_ucmodel == number_models:
                                next_ucmodel += 1
                            next_casamodel = int(number_models+5+np.random.randint(5))
                elif cfg.agc.corruption_method == 'No_Corrupt' or \
                    (cfg.agc.corruption_method == 'Tres' \
                    and number_models == next_ucmodel):
                    if not (cfg.agc.variables_to_vary[ix] == 'SNR'):
                        setattr(Current_Galaxy, "Corruption", "No_Corrupt")
                        next_ucmodel = int(number_models+5+np.random.randint(5))
                    else:
                        setattr(Current_Galaxy, "Corruption", "Gaussian")
                        next_ucmodel += 1                     
                else:
                    setattr(Current_Galaxy, "Corruption", "Gaussian")
                
                Achieved = copy.deepcopy(Current_Galaxy)
                #Check if galaxy already exists
                name = set_name(Current_Galaxy)
                print(f"{name} is the name we attach to the current galaxy")

                # Check for the existence of the directory
                constructstring = f"mkdir {cfg.general.main_directory}{name}"
                checkdir = False
                galaxy_dir = f"{cfg.general.main_directory}{name}/"
                galaxy_exists = os.path.isdir(galaxy_dir)
                if not galaxy_exists:
                    os.system(constructstring)
                    created.append(name)
                else:
                    time.sleep(0.1)
                    # Do we have a cube
                    galaxy_cube_exist = os.path.isfile(
                        f"{galaxy_dir}Convolved_Cube.fits")
                    if not galaxy_cube_exist:
                        galaxy_cube_exist = os.path.isfile(
                            f"{galaxy_dir}Convolved_Cube_CS.fits")
                    if not galaxy_cube_exist:
                        galaxy_cube_exist = os.path.isfile(
                            f"{galaxy_dir}Convolved_Cube_Gauss.fits")
                    if galaxy_cube_exist:
                        print("This galaxy appears fully produced")
                        checkdir = True
                        continue
                    else:
                        if name in created:
                            number_models -= 1
                            continue
                        else:
                            print(
                                "The directory was made but there is no full cube available")
                            #print("Reproducing the galaxy. Be aware of Double Table entries")
                            print("This is too dangerous. Breaking the code.")
                            exit()
                if Current_Galaxy.Corruption in ['No_Corrupt', 'Gaussian']:
                    Gauss_Galaxies.append((cfg, Current_Galaxy, Achieved))
                else:
                    Casa_Galaxies.append((cfg, Current_Galaxy, Achieved))
                All_Galaxies.append(name)
                if Current_Galaxy.Mass not in set_done:
                    SBRprof, Rad, sclength, MHI, Rad_HI, Vrot, sub_ring, molecular_profile, DynMass = build_sbr_prof(
                        Current_Galaxy, symmetric=cfg.agc.symmetric, no_template=True)  # Column densities,Raii in kpc, Opt_scalelength in kpc, HI mass in M_solar
                    set_done, max_rad, colors, plot_ax = plot_RC(
                        set_done, Current_Galaxy.Mass, Rad, Vrot, colors, max_rad, sub_ring, plot_ax,font_file=cfg.general.font_file)
                    rc_counter += 1
    if rc_counter > 0:
        plt.figure(59)
        plot_ax.set_ylim(ymin=0)
        plot_ax.set_xlim(xmin=0, xmax=max_rad)
        plt.legend(loc='lower right', fontsize=12)
        plt.savefig('Rotation_Curves.pdf', bbox_inches='tight')
        # this is not closing properly
        plot_ax.cla()
        plt.clf()
        plt.close()
                #print(f"This is the parameter to vary {cfg.agc.variables_to_vary[ix]}.")
    if len(Casa_Galaxies) > 0:
        #with open(f"{cfg.general.main_directory}/Casa_Noise_Statistics.txt", 'w') as file:
        #    file.write(f"{'Req SNR':10s} {'Achieved SNR':10s} {'Inc. Fact':10s} {'Factor Check':10s} {'Mean Flux':10s} {'Input Noise':10s} {'Output Noise':10s} {'Input Noise':10s} \n")
        #    file.write(
        #        f"{' ':10s} {' ':10s} {' ':10s} {' ':10s} {'Jy/Beam':10s} {'Jy/Beam':10s} {'Jy/Beam':10s} {'Jy':10s} \n")
#
        #tclean is parallel inmplemented so we run all casa corruption singular
        results_casa = []
        for the_galaxy in Casa_Galaxies:
            single_result = one_galaxy(*the_galaxy)
            results_casa.append(single_result)
    #Create All Uncoorupted and Gaussian Galaxies
    if len(Gauss_Galaxies) > 0:
        if cfg.general.multiprocessing:
            no_process = cfg.general.ncpu
            if no_process > len(Gauss_Galaxies):
                no_process = len(Gauss_Galaxies)
            with get_context("spawn").Pool(processes=no_process) as pool:
                results_gauss = pool.starmap(one_galaxy, Gauss_Galaxies)
        else:
            results_gauss = []
            for the_galaxy in Gauss_Galaxies:
                single_result = one_galaxy(*the_galaxy)
                results_gauss.append(single_result)

    results = ['empty']*len(All_Galaxies)
    if len(Gauss_Galaxies) > 0:
        results = sort_results(All_Galaxies, results_gauss, results)
    if len(Casa_Galaxies) > 0:
        results = sort_results(All_Galaxies, results_casa, results)

    os.system(f'rm -f {cfg.general.main_directory}/Input.fits')
    print("We created {} models".format(number_models))

    #Add the results to the catalogue
    with open(Catalogue, 'a') as cat:
        for line in results:
            cat.write(line)


AGC.__doc__ = f'''
NAME:
   AGC

PURPOSE:
    Create realistic artificial galaxies.

CATEGORY:
    agc

INPUTS:
    cfg = OmegaConf Configuration file

OPTIONAL INPUTS:

OUTPUTS:
    A set of artificial galaxies in a setup ready for FAT fitting

OPTIONAL OUTPUTS:

PROCEDURES CALLED:
   Unspecified

NOTE:
'''


def calc_vc_NFW(DM_mass, MHI, m_star, rad):
    r200 = (DM_mass*G_agc/(100. * H_0**2)*1e12)**(1./3.)
    v200square = DM_mass*G_agc/r200
    xr = rad/(r200/10**3)
    # From A. Dutton 2014
    c200 = 10**(0.905-0.101*np.log10(DM_mass/(10**12*100./H_0)))
    NFWvflat = np.sqrt(v200square*(1./xr)*((np.log(1.+c200*xr)
                       - (c200*xr)/(1+c200*xr))/(np.log(1+c200)-c200/(1+c200))))
    v_star = np.sqrt(m_star*G_agc/(rad*10**3))
    v_HI = np.sqrt(MHI*1.4*G_agc/(rad*10**3))
    return np.sqrt(NFWvflat**2+v_star**2+v_HI**2)


calc_vc_NFW.__doc__ = f'''
NAME:
   calc_vc_NFW

PURPOSE:
    Calculate the total circular velocity bases on the DM mass with a NFW profile, the stellar mass and the HI mass.

CATEGORY:
   agc

INPUTS:
   DM_mass = the DM mass
   MHI = the HI mass
   m_star = the stellar mass
   rad = the radius at which to calculate the circular velocity

OPTIONAL INPUTS:

OUTPUTS:
    The total circular velocity at the provided radius

OPTIONAL OUTPUTS:

PROCEDURES CALLED:
   Unspecified

NOTE:
'''

# First we define some useful routines for converting and copying
# A routine to copy disks


def copy_disk(olddisk, newdisk):
    '''Routine to copy a disk in the tirific template file'''
    start = 0
    startlast = 0.
    if int(Template["NDISKS"].split('=')[1]) < newdisk:
        Template["NDISKS"] = "NDISKS = {:d}".format(newdisk)
    copkeys = "Empty"
    last = 'RADI'
    for key in Template.keys():
        if start == 2:
            if copkeys == "Empty":
                copkeys = [key]
            else:
                copkeys.append(key)
        if key == 'RADI':
            start += 1
            if olddisk == 1:
                start += 1
        if '_' in key:
            key_ext = key.split('_')
            if key_ext[1] == str(olddisk) and start == 1:
                start += 1
                if olddisk > 1:
                    copkeys = [key]
            elif (key_ext[1] != str(olddisk)) and start == 2:
                del copkeys[-1]
                start += 1
            if key_ext[1] == str(newdisk-1) and startlast == 0:
                startlast = 1
            if key_ext[1] == str(newdisk-1) and startlast == 1:
                last = key
            if key_ext[1] != str(newdisk-1) and startlast == 1:
                startlast += 1
        if key == 'CONDISP' and (start == 2 or startlast == 1):
            if start == 2:
                del copkeys[-1]
            startlast += 1
            start += 1
    for key in reversed(copkeys):
        var_name = key.split('_')[0]
        Template.insert(last, var_name+"_{:d}".format(newdisk), var_name
                        + "_{:d} =".format(newdisk)+Template[key].split('=')[1])
    Template.insert("CFLUX_{:d}".format(newdisk-1), "CFLUX_{:d}".format(newdisk),
                    "CFLUX_{:d} =".format(newdisk)+Template["CFLUX_{:d}".format(newdisk-1)].split('=')[1])
    if '_' in copkeys[-1]:
        return copkeys[-1].split('_')[0]+"_{:d}".format(newdisk)
    else:
        return copkeys[-1]+"_{:d}".format(newdisk)


copy_disk.__doc__ = f'''
NAME:
   copy_disk

PURPOSE:
    Routine to copy a disk in the tirific template file

CATEGORY:
   agc

INPUTS:
   Olddisk = the # of the disk to be copied in def file notation, i.e 1 = paramaters without extension, 2 with extension _2 and so on
   newdisk = the # of the disk to be copied in def file notation, i.e 1 = paramaters without extension, 2 with extension _2 and so on

OPTIONAL INPUTS:

OUTPUTS:
    the variable names of the new disk

OPTIONAL OUTPUTS:

PROCEDURES CALLED:
   Unspecified

NOTE:
'''


# Obtaining the derivative at any give point of a function
def derivative(x, func):
    '''Obtaining the derivative at any give point of a function'''
    h = x/1000.
    der = (func(x+h)-func(x-h))/(2*h)
    return der


derivative.__doc__ = f'''
NAME:
    derivative

PURPOSE:
    Obtaining the derivative at any give point of a function

CATEGORY:
   agc

INPUTS:
    x = point where to obtain the derivative
    func = the function to evaluate

OPTIONAL INPUTS:

OUTPUTS:
    the derivative

OPTIONAL OUTPUTS:

PROCEDURES CALLED:
   Unspecified

NOTE:
'''


# Then the functions that create the different components of the galaxy

# function to get the surface brightness profile based on the last element in the rotation curve
def build_sbr_prof(Current_Galaxy, symmetric=False, no_template=False):
    '''function to get the surface brightness profile based on the last element in the rotation curve'''
    # First we use the mass to build a rotation curve and sbr profile
    # Avanti used 	V_Ropt = (200*(l**0.41))/(0.80 +0.49*np.log10(l) +(0.75*np.exp(-0.4*l))/(0.47 + 2.25*(l**0.4)))**0.5 (Persic & Sallucci)
    # We will not use the URC as it get's silly at high mass and remains showing flat parts at low mass
    # We will use the parameterisation of Courteau 1997
    v_circ, HIrad, MHI, MK, DynMass = get_HI_disk(Current_Galaxy.Mass)

        #First we set the radii at with a 10 elements per rings out to 1.5 times HIrad. However we always need at least 25 rings for transition purposes
    if Current_Galaxy.Beams < 5.:
        sub_ring = int(25./Current_Galaxy.Beams)
    else:
        sub_ring = 5
    Rad = np.linspace(0., 1.5*HIrad, int(sub_ring*Current_Galaxy.Beams))
    # Finaaly the courteau presciption
    v_flat, R_0, beta, gamma = get_sparcs_fit(v_circ, HIrad)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        x = R_0/Rad
    x[0] = 1e-8

    vrot = v_flat*(1+x)**beta*(1+x**gamma)**(-1./gamma)
    vrot[0] = 0.
    # from these we use the the prescription of Martinsson to get an HI profile.
    # As done by A. Gogate in initialparam.py
    # First we set the radii at with a 5 elements per rings out to 1.5 times HIrad
    a = (Rad) * 10 ** 3  # in parsecs
    fracHIrad = (HIrad - np.floor(HIrad))
    # hole in the center should be amplified because of the conversion to H_2
    # with a turn over at ~ 120 V_max where it should start adding back
    # the H_2 follows an exponential distribution similar to the stellar disk (ref?)
    # The scale length is h=Vmax^2/(0.88**2*np.pi*G*sig0) (Freeman 1970)

    # The scale length for the molecular disk is then from Graham 2014 rederived relation in cal_scl.py
    h_r = (-4.13422991 - 0.31576291 * MK) * 1000.  # parsec

    if not symmetric:
        if fracHIrad < 0.5:
            Rhi_p = np.array([(HIrad - HIrad / 15.) * 10 ** 3,
                             (HIrad + HIrad / 15.) * 10 ** 3])
            # Rhi_p2 = HIrad + HIrad / 20. * 10 ** 3
            h_r = np.array([h_r + h_r / 5, h_r - h_r / 5])
            # h_r2 = h_r - h_r / 20.
        else:
            Rhi_p = np.array([(HIrad + HIrad / 15.) * 10 ** 3,
                             (HIrad - HIrad / 15.) * 10 ** 3])
            # Rhi_p2 = HIrad - HIrad / 20. * 10 ** 3
            h_r = np.array([h_r - h_r / 5., h_r + h_r / 5.])
            # h_r2 = h_r + h_r / 20.
    else:
        Rhi_p = np.array([HIrad, HIrad]*10**3, dtype=float)
        h_r = np.array([h_r, h_r], dtype=float)
    Hiradindex = [0, 0]
    Hiradindex[0] = np.where((Rad > (Rhi_p[0] - Rhi_p[0] / (sub_ring * Current_Galaxy.Beams)) / 10 ** 3) & (
            Rad < (Rhi_p[0] + Rhi_p[0] / (sub_ring * Current_Galaxy.Beams)) / 10 ** 3))[0][0]
    Hiradindex[1] = np.where((Rad > (Rhi_p[1] - Rhi_p[1] / (sub_ring * Current_Galaxy.Beams)) / 10 ** 3) & (
            Rad < (Rhi_p[1] + Rhi_p[1] / (sub_ring * Current_Galaxy.Beams)) / 10 ** 3))[0][0]

    # print("This the HI Radius in kpc {}".format(HIrad))
    # std deviation or dispersion of gaussian
    s1 = 0.36 * Rhi_p
    # s2 = 0.36 * Rhi_p2
    # We assume that at 120 km/s v_max the H_2 disk imprint on the HI disk is adequately described by the prescription of Martinsson.
    # Lower some of the disk is added back to the HI higher the central hole is more pronounced
    I_cen = ((v_circ / 120.) ** 0.5 - 1)
    #print("This the scale length {} and central brightness {}".format(h_r, I_cen))
    # So our molecular profile is
    Exp = np.zeros([len(a), 2])
    Sig2 = np.zeros([len(a), 2])
    Sigma = np.zeros([len(a), 2])
    new = np.zeros(2)
    for i in [0, 1]:
        Exp[:, i] = I_cen * np.exp(-a / h_r[i])
        # gaussian2 From Martinsson 2015
        Sig2[:, i] = np.exp(-(a - 0.4 * Rhi_p[i]) ** 2 / (2 * (s1[i]) ** 2))
        # total
        Sigma[:, i] = Sig2[:, i] - Exp[:, i]
        Sigma[Sigma[:, i] < 0.,
            i] = 0.  # for negative sigma max, does not include negative values
        # scale Sigma such that it is one at HI rad
        new[i] = 1. / Sigma[Hiradindex[i], i]
        Sigma[:, i] = new[i] * Sigma[:, i]
    Sigma[0:4, 0] = (Sigma[0:4, 0]+Sigma[0:4, 1])/2.
    Sigma[0:4, 1] = Sigma[0:4, 0]
    # get the HI Mass in the profile
    OutHIMass = integrate.simpson(
        (np.pi * a) * Sigma[:, 0], x=a) + integrate.simpson((np.pi * a) * Sigma[:, 1], x=a)
    # And check that it matches the required HI mas within 5%
    counter = 1
    while np.absolute(MHI - OutHIMass) > MHI * 0.02:
        # if not rescale sigma1
        if MHI - OutHIMass > 0:
            s1 = (0.36 - (0.0025 * counter)) * Rhi_p
        else:
            s1 = (0.36 + (0.0025 * counter)) * Rhi_p
        # and recalculate the profile
        for i in [0, 1]:
            Sig2[:, i] = np.exp(-(a - 0.4 * Rhi_p[i])
                                ** 2 / (2 * (s1[i]) ** 2))
            Sigma[:, i] = Sig2[:, i] - Exp[:, i]
            Sigma[Sigma[:, i] < 0., i] = 0.
            new[i] = 1. / Sigma[Hiradindex[i], i]
            Sigma[:, i] = new[i] * Sigma[:, i]
        Sigma[0:4, 0] = (Sigma[0:4, 0]+Sigma[0:4, 1])/2.
        Sigma[0:4, 1] = Sigma[0:4, 0]
        OutHIMass = integrate.simpson(
            (np.pi * a) * Sigma[:, 0], x=a) + integrate.simpson((np.pi * a) * Sigma[:, 1], x=a)
        counter += 1
    #print("The requested HI mass = {:.2e} and the retrieved HI mass = {:.2e}".format(MHI, OutHIMass))
    # final HI radial distribution by renormalisation
    S = Sigma * (1.24756e+20)
    # Where A. Gogate contribution stops
    # S is column densities but tirific takes Jy * km/s/arcsec^2 so
    conv_column_arsec = 605.7383 * 1.823E18 * \
        (2. * np.pi / (np.log(256.)))  # S(mJy/beam)*conv_column_arcsec=N_HI
    sbr_prof = S / (conv_column_arsec * 1000.)
    # Let's write these to the Template immediately
    # The surface brightness profile, which is still symmetric
    if not no_template:
        Template["SBR"] = "SBR = " + " ".join(str(e) for e in sbr_prof[:, 0])
        Template["SBR_2"] = "SBR_2 = " + \
            " ".join(str(e) for e in sbr_prof[:, 1])

    # as the rest on of the code is based on a single SBR profile let's average the value
    HIrad = [np.mean(Rhi_p) / 10 ** 3, Rhi_p[0] / 10 ** 3,
                     Rhi_p[1] / 10 ** 3]  # HI radii in kpc mean,appr, rec
    molecular_profile = []
    for i in [0, 1]:
        molecular_profile.append(
            Exp[:, i]*new[i]*1.24756e+20/(conv_column_arsec*1000.))  # Jy/arcsec^2 km/s
    molecular_profile = np.array(molecular_profile)
    h_r = np.mean(h_r)
    average_sbr_profile = np.zeros(len(sbr_prof[:, 0]))
    for i in range(len(sbr_prof)):
        average_sbr_profile[i] = np.mean(sbr_prof[i, :])
    return average_sbr_profile, Rad, h_r/1000., OutHIMass, HIrad, vrot, sub_ring, molecular_profile, DynMass


build_sbr_prof.__doc__ = f'''
NAME:
   build_sbr_prof

PURPOSE:
    Take the input galaxy and from the mass build the SBR profile, the RC,
    the scale length and calculated the number of subring the model requires

CATEGORY:
   agc

INPUTS:
   Current_Galaxy = Galaxy class input
   symmetric = Boolean do we want symmetric profiles?

OPTIONAL INPUTS:

OUTPUTS:
    average_sbr_profile = an average profile for the model in Jy Km/s arcsec^-2.
                    The assymetric profile is written to Template whioch is a global
    Rad = the Rad in kpc
    h_r/1000. = the scale length in kpc
    OutHIMass = the final retrieved HI mass
    HIrad= the symmetric HI radius at 1 M_solar/pc^2
    vrot = The rotation curve
    sub_ring = the number sub rings  per beam in the model.


OPTIONAL OUTPUTS:

PROCEDURES CALLED:
   Unspecified

NOTE:
'''


def create_arms(velocity, Radii, disk_brightness, disk=1, WarpStart=-1,\
                    Bar="No_Bar"):
    '''Create the spiral Arms'''
    if WarpStart == -1: WarpStart = Radii[-1]
    max_vrot = np.max(velocity)
    max_rad = Radii[-1]
    #The pattern speed at a given radius is vrot/radii
    V_Rot = interpolate.interpolate.interp1d(
        Radii, velocity, fill_value="extrapolate")
    # The radius of co-ration can be approximated by the extend of the visible disk ~ Warpstart (Roberts et al. 1975)
    Omega_CR = V_Rot(WarpStart)/WarpStart
    # From this we can estimate the inner and outer Lindblad resonances (Eq 38 Dobbs & Baba 2014)
     #The epicyclic frequency k^2=R*dOmega^2/dR+4*Omega^2
    # f'(x) = (f(x+h)-f(x-h))/2h
    h = WarpStart/1000.
    derive = (V_Rot(float(WarpStart+h))**2/(WarpStart+h)**2
              - V_Rot(float(WarpStart-h))**2/(WarpStart-h)**2)/(2*h)

    k_CR = (WarpStart * derive+4*Omega_CR**2)**0.5
    # So the ILR =
    if Bar == "Barred":
        num_arms = 2
    else:
        #other wise 2 arms when a bulge is present and 4 arms when not
        if Radii[np.where(max_vrot == velocity)[0]] < np.mean(Radii):
            num_arms = 2
        else:
            num_arms = 4

    LLR = Omega_CR-k_CR/num_arms
    ULR = Omega_CR+k_CR/num_arms
    Radii[0] = 0.1
    om = interpolate.interpolate.interp1d(
        Radii, velocity/Radii, fill_value="extrapolate")
    Radii[0] = 0.
    r_cur = Radii[1]
    while om(r_cur) > ULR and r_cur < max_rad:
        r_cur += 0.1
    ILR = 0.75*r_cur
    r_cur = Radii[1]
    while om(r_cur) > LLR and r_cur < max_rad:
        r_cur += 0.1
    OLR = 0.75*r_cur

    # From Seigar et al. 2006 we get the relation between shear (S) and pitch angle
    S = 0.5*(1-Radii[1:]/velocity[1:]*derivative(Radii[1:], V_Rot))
    pitch2 = 64.25-73.24*S
    #As we assume a constant pitch angle we will take the average between ILR and OLR as the pitch angle
    tmp = np.where((Radii[1:] > ILR) & (Radii[1:] < OLR))[0]
    pitch = np.sum(pitch2[tmp])/len(tmp)

    #print("This is the average pitch angle {}".format(pitch))
    #Using Kennicut's prescription.This prescription incorrectly states cot(P) instead of tan(P) see Davis et. al 2012
    # The arms start at the inner Lindblad Resonance and hence the phase is 0 there
    phase = np.log(Radii[1:]/ILR)/np.tan(pitch*np.pi/180.)*180/np.pi
    phase = np.hstack((phase[0], phase))

    #How many arms do we make
    # with bar it is always a grand design

    # we take a brightness in the arms 1/no_arms the brightness of the disk
    brightness = 1./np.sqrt(num_arms)*disk_brightness
    # and only between the resonances
    index = np.where((Radii < ILR) | (Radii > OLR))[0]
    brightness[index] = 0.
    # with a ten ring transition
    brightness[tmp[0]:tmp[0]+10] = brightness[tmp[0]:tmp[0]+10]*(1-1/np.arange(1, 11))
    brightness[tmp[-1]-10:tmp[-1]] = brightness[tmp[-1]
        - 10:tmp[-1]]*(1/np.arange(1, 11))
    # For the width we take 15% of the full circle but scaled to the total galaxy size a good description would be swell.
    width = 0.15*2.*np.pi*Radii
    if WarpStart < 10.:
        width = width*10./WarpStart
    ndisk = int(Template["NDISKS"].split('=', 1)[1])
    ndisk += 1
    last_add = copy_disk(disk, ndisk)
    Template[f"SBR_{ndisk:d}"] = f"SBR_{ndisk:d} = 0."
    # To simulate strems towrds the arms we spin this up by 10 km/s
    Template[f"VROT_{ndisk:d}"] = f"VROT_{ndisk:d} = {' '.join(str(e+20.) for e in velocity)}"
    phaseshift = 360./num_arms
    #we offset the phase by 37 degrees
    for i in range(num_arms):
        Template.insert(last_add, "GA{:d}A_{:d}".format(
            i+1, ndisk), "GA{:d}A_{:d} =".format(i+1, ndisk)+" ".join(str(e) for e in brightness))
        Template.insert("GA{:d}A_{:d}".format(i+1, ndisk), "GA{:d}P_{:d}".format(i+1, ndisk),
                        "GA{:d}P_{:d} =".format(i+1, ndisk)+" ".join(str(e+i*phaseshift+37) for e in phase))
        Template.insert("GA{:d}P_{:d}".format(i+1, ndisk), "GA{:d}D_{:d}".format(
            i+1, ndisk), "GA{:d}D_{:d} =".format(i+1, ndisk)+" ".join(str(e) for e in width))
    # and we add radial motions of 40 km/s
    Template.insert("VROT_{:d}".format(ndisk), "VRAD_{:d}".format(
        ndisk), "VRAD_{:d} = -40.".format(ndisk))

    return phase, brightness, width


create_arms.__doc__ = f'''
NAME:
   create_arms

PURPOSE:
    Create the spiral Arms if requested

CATEGORY:
   agc

INPUTS:
    velocity = Rotation curve
    Radii= the radii
    disk_brightness
    disk=1,WarpStart=-1,Bar="No_Bar"

OPTIONAL INPUTS:
    disk = 1
    Are we creating in the approaching (1) or receding side (2)

    WarpStart=-1
    radius of the start of the warp

    Bar="No_Bar"
    Boolean whether bar is included or not

OUTPUTS:
    phase = the phases of the arms
    brightness = the brightness amplitude of the arms
    width =the used width in deg

OPTIONAL OUTPUTS:

PROCEDURES CALLED:
   Unspecified

NOTE:
'''

# Routine to create a bar


def create_bar(velocity, Radii, disk_brightness, Template, disk=1, WarpStart=-1):
    '''Routine to create a Bar'''
    if WarpStart == -1: WarpStart = Radii[-1]
    max_rad = Radii[-1]
    bar_width = 1+0.05*WarpStart  # kpc + 5% of the total optical disk
    max_vrot = np.max(velocity)
    #The pattern speed at a given radius is vrot/radii
    V_Rot = interpolate.interpolate.interp1d(
        Radii, velocity, fill_value="extrapolate")
    # The radius of co-ration can be approximated by the extend of the visible disk ~ Warpstart (Roberts et al. 1975)
    Omega_CR = V_Rot(WarpStart)/WarpStart
    # From this we can estimate the inner and outer Lindblad resonances (Eq 38 Dobbs & Baba 2012)
     #The epicyclic frequency k^2=R*d/drOmega^2+4*Omega^2
    # f'(x) = (f(x+h)-f(x-h))/2h
    h = WarpStart/1000.
    derive = (V_Rot(float(WarpStart+h))**2/(WarpStart+h)**2
              - V_Rot(float(WarpStart-h))**2/(WarpStart-h)**2)/(2*h)
    k_CR = (WarpStart * derive+4*Omega_CR**2)**0.5
    # So the ILR =
    LLR = Omega_CR-k_CR/2.
    ULR = Omega_CR+k_CR/2.
    Radii[0] = 0.1
    om = interpolate.interpolate.interp1d(
        Radii, velocity/Radii, fill_value="extrapolate")
    Radii[0] = 0.
    r_cur = Radii[1]
    while om(r_cur) > ULR and r_cur < max_rad:
        r_cur += 0.1
    ILR = 0.75*r_cur
    r_cur = Radii[1]
    while om(r_cur) > LLR and r_cur < max_rad:
        r_cur += 0.1
    # We set the full brightness to the maximum of the disk
    brightness = np.zeros(len(disk_brightness))
    brightness[np.where(Radii < ILR)[0]] = np.max(disk_brightness)
    # The width has to be 180 when R < width and 180*width/(pi*r)
    width = np.zeros(len(disk_brightness))
    width[Radii <= bar_width] = 180.
    # the angle made up of the radius and width *2.
    width[Radii > bar_width] = 360./np.pi * \
        np.arcsin(bar_width/Radii[Radii > bar_width])
    # Get the number of disks present
    ndisk = int(Template["NDISKS"].split('=', 1)[1])
    ndisk += 1
    #print("We are adding disk no {:d}".format(ndisk))
    # we also need streaming motions
    vrad_bar = np.zeros(len(disk_brightness))
    vrad_bar[:] = -50.
    #copy the input disk
    last_add = copy_disk(disk, ndisk)
    #We offset by 37 deg.
    Template.insert("VSYS_{:d}".format(ndisk), "AZ1P_{:d}".format(
        ndisk), "AZ1P_{:d} = 37.".format(ndisk))
    Template.insert("AZ1P_{:d}".format(ndisk), "AZ1W_{:d}".format(
        ndisk), "AZ1W_{:d} = ".format(ndisk)+" ".join(str(e) for e in width))
    Template.insert("AZ1W_{:d}".format(ndisk), "AZ2P_{:d}".format(
        ndisk), "AZ2P_{:d} = 217.".format(ndisk))
    Template.insert("AZ2P_{:d}".format(ndisk), "AZ2W_{:d}".format(
        ndisk), "AZ2W_{:d} = ".format(ndisk)+" ".join(str(e) for e in width))
    # And we add streaming motions to the bar km/s
    Template.insert("VROT_{:d}".format(ndisk), "VRAD_{:d}".format(
        ndisk), "VRAD_{:d} = ".format(ndisk)+" ".join(str(e) for e in vrad_bar))
    Template["SBR_{:d}".format(ndisk)] = "SBR_{:d} =".format(
        ndisk)+" ".join(str(e) for e in brightness)
    return ILR


create_bar.__doc__ = f'''
NAME:
   create_bar

PURPOSE:
    Create a bar

CATEGORY:
   agc

INPUTS:
    velocity,Radii,disk_brightness,Template, disk=1,WarpStart=-1
    velocity = Rotation curve
    Radii= the radii
    disk_brightness = the SBR profile of the disk
    disk=1,WarpStart=-1,Bar="No_Bar"

OPTIONAL INPUTS:
    disk = 1
    Are we creating in the approaching (1) or receding side (2)

    WarpStart=-1
    radius of the start of the warp

    Bar="No_Bar"
    Boolean whether bar is included or not

OUTPUTS:
    phase = the phases of the arms
    brightness = the brightness amplitude of the arms
    width =the used width in deg

OPTIONAL OUTPUTS:

PROCEDURES CALLED:
   Unspecified

NOTE:
'''

# Thi function creates the flarin g which is based on the rotation curve and dispersion according Puche et al. 1992


def create_flare(Radii, velocity, dispersion, flare, Max_Rad, sub_ring, \
                    distance=1.):
    '''This function creates the scale height and  dispersion profiles according to Puche et al. 1992'''
    #make sure we have arrays so we can do numerical operations
    Radii = np.array(Radii)
    velocity = np.array(velocity)
    # first we need the dispersion
    # if there is no change from inside to outside it is easy
    if dispersion[0] == dispersion[1]:
       disp = np.full(len(Radii), dispersion[0])
    else:
        # else we use a tangent change with the center at halfway
       disp = -1*np.arctan((Radii-np.mean(Max_Rad/2.))/(Radii[-1]/10.))/np.pi*np.absolute(
           dispersion[0]-dispersion[1])+np.mean(dispersion)
    # We recalculate the Dynamical Mass and densities
    # Dynamical Mass

    Dynmass = (Radii*10**3)*velocity**2/G_agc
    # This is in a Volume of
    Volume = (4./3.)*np.pi*Radii**3
    Dynmass[0] = 1.
    Volume[0] = 1.
    Density = Dynmass/Volume * \
        2.  # In M_solar/kpc^3 the two comes from Puche 1992 but seems random
    # Then we check wether we want a flare or not
    G2 = G_agc/(3.086e+13**2)  # pc^3 M_sol^-1 s^-2
    halfint = int((len(Radii[Radii < Max_Rad])+10)/2.)
    if flare.lower() == 'flare':
        flare = disp/((4.*np.pi*G2*Density/1000.**3)**0.5*3.086e+16)  # in kpc
        flare[:halfint-10] = flare[halfint-10]
        fact = np.arange(1/21, 1, 1./21)
        flare[halfint-10:halfint
            + 10] = (1-fact)*flare[halfint-10]+fact*flare[halfint-10:halfint+10]
    elif flare.lower() == 'no_flare':
        flare = np.full(len(
            Radii), disp[halfint]/((4.*np.pi*G2*Density[halfint]/1000.**3)**0.5*3.086e+16))
    else:
        print(f"{flare} is not an option for the flare. Choose Flare or No_Flare")
        sys.exit()

    flare[0] = flare[1]

    # convert the scale heights to arcsec
    h_z_arcsec = cf.convertskyangle(flare, distance=distance, physical=True)
    # and write both to the Template
    Template["SDIS"] = "SDIS = "+" ".join(str(e) for e in disp)
    Template["SDIS_2"] = "SDIS_2 = "+" ".join(str(e) for e in disp)
    # The scaleheight
    Template["Z0"] = "Z0 = "+" ".join(str(e) for e in h_z_arcsec)
    Template["Z0_2"] = "Z0_2 = "+" ".join(str(e) for e in h_z_arcsec)
    return flare, disp


create_flare.__doc__ = f'''
NAME:
   create_flare

PURPOSE:
    This function creates the scale height and  dispersion profiles according to Puche et al. 1992

CATEGORY:
   agc

INPUTS:
    Radii = The radii at which the RC is defined
    velocity = the RC
    dispersion = the requested inner and outer dispersion
    flare = Flare or No_Flare
    Max_Rad = the HI radius
    sub_ring = the sub rings used

OPTIONAL INPUTS:
    distance= distance to the galaxy


OUTPUTS:
    flare = the scale height profile
    disp = the set dispersion profile

OPTIONAL OUTPUTS:

PROCEDURES CALLED:
   Unspecified

NOTE:
'''

#A function to create the template fits file
def create_template_fits(directory):
    from pyHIARD.Templates.Template_Header import sample_header_dictionary
    new_cube = np.full([120,500,500],1.)
    new_cube[60, 250, 250] = 5.
    new_list = fits.PrimaryHDU(new_cube)
    new_cube = []
    header_keys = [x for x in new_list.header]
    for key, value in sample_header_dictionary.items() :
        if key not in header_keys:
            new_list.header[key] = value

    dummy = fits.HDUList([new_list])

    header_sets = {'BMAJ': 0., 'BMIN': 0., 'BPA': 0.,
                   'CRPIX1': 200., 'CRPIX2': 200., 'CRPIX3': 60.,
                   'CDELT1': -4./3600., 'CDELT2': 4./3600., 'CDELT3': 4.,
                   'CUNIT2': 'M/S', 'CRVAL3': 500.,
                   'CTYPE1': 'RA---SIN',  'CTYPE2': 'DEC--SIN', 'CTYPE3': 'VELO-HEL',
                   'BITPIX': -32, 'NAXIS': 3,
                   }
    for key, value in header_sets.items():
        dummy[0].header[key] = value

    dummy.writeto(f'{directory}/Input.fits',\
                    output_verify='silentfix+ignore', overwrite=True)
    dummy.close()

# A function for varying the PA and inclination as function of radius
def create_warp(Radii,PA,inclination,warp_change,
                warp_radii,disk=1,debug=False):
    if debug:
        print(f'''We are starting create_warp.
warp_change = {warp_change}
warp_radii = {warp_radii}''')
    Radii = np.array(Radii)
    if ((np.sum(warp_change) != 0) and (warp_radii[0] < warp_radii[1])).all():
        # First we need to obtain the vector that constitutes the inner area
        #it runs exactly counter to inclination
        inclination = 90-inclination
        # For this the PA has to be between 0-90
        mult = np.floor(PA/90.)
        inPA = PA-mult*90.
        #avoid singularities
        if inPA == 0.:
            inPA = 0.001
        if inclination < 0.001:
            inclination = 0.001
        # define the angular momentum vector of the plane and the outer most ring
        theta = np.arctan(np.tan(inclination*(np.pi/180.))
                          * np.tan(inPA*(np.pi/180.)))
        phi = np.arctan(np.tan(inPA*(np.pi/180.))/np.sin(theta))
        # and the maximum values at Rad_HI
        thetamax = theta+warp_change[0]
        phimax = phi+warp_change[1]
                # indices of the warp start and the Rad_HI
        start_index = int(np.sum(Radii < warp_radii[0]))
        end_index = int(np.sum(Radii < warp_radii[1]))
        # step size of theta and phi
        # As we will increase the step size triangular we need the total number of point in the sequence
        warprange = end_index-start_index
        if warprange < 2:
            thetamax = theta
            phimax = phi
            warprange = 1
        increasetheta = (thetamax-theta)/(0.5*warprange*(warprange+1))
        increasephi = (phimax-phi)/(0.5*warprange*(warprange+1))
        #print(warprange,thetamax,phimax,Radii[1]-Radii[0],warp_radii[1]-warp_radii[0])
        # calculate theta
        thetarings = np.array(np.full(len(Radii), theta))
        index_array = np.arange(start_index, len(Radii))-start_index
        thetarings[start_index:] = theta+0.5 * \
            index_array*(index_array+1)*increasetheta
        #calculate phi
        phirings = np.array(np.full(len(Radii), phi))
        phirings[start_index:] = phi+0.5 * \
            index_array*(index_array+1)*increasephi
        # return to PA

        if (phirings[0] < np.pi/2.) and (phirings[-1] > np.pi/2) and (inclination < 5.):
            PA = np.arctan(np.sin(thetarings)*np.tan(phirings-np.pi/2.))*(360./(2*np.pi)) + \
                           mult*90+np.arctan(np.sin(theta)
                                             * np.tan(phi))*(360./(2*np.pi))
        else:
            PA = np.arctan(np.sin(thetarings)*np.tan(phirings)) * \
                           (360./(2*np.pi))+mult*90
            PA[phirings > 0.5*np.pi] = PA[phirings > 0.5*np.pi]+180.

        # return inclination
        inc = 90-np.arctan(1./(np.cos(thetarings)
                           * np.tan(phirings)))*(360./(2*np.pi))
        # return inclination boundary adjustements
        inc[np.where(inc > 90.)] = 180 - inc[np.where(inc > 90.)]
        inc[np.where(inc < 0.)] = -1 * inc[np.where(inc < 0.)]
         # return a correct quadrant phirings
        phirings = phirings+mult/2.*np.pi
    else:
        #if there is no warp then all values are the same
        PA = np.full(len(Radii), PA)
        inc = np.full(len(Radii), inclination)
        phirings = np.full(len(Radii), 0)
        thetarings = np.full(len(Radii), 0)

    #write to our template file
    # let's see if we can retrace intrinsic phi with the formula's from Peters
    #theta_test = np.arctan((np.sin(PA*np.pi/180.)*np.sin(inc*np.pi/180.)-np.cos(PA*np.pi/180.)*np.sin(inc*np.pi/180.))/(np.cos(PA*np.pi/180.)*np.cos(inc*np.pi/180.)-np.sin(PA*np.pi/180.)*np.cos(inc*np.pi/180.)))*180./np.pi
    phirings = phirings*180./np.pi
    # thetarings=thetarings*180./np.pi
    # This seems to work mostly but not at some extremes exactly for some reason
    # According to josh tan is required, and yes that makes it work at large angles as well.
    #angle_adjust = np.tan(
    #    (PA[0]-PA)*np.cos(inc*np.pi/180.)*np.pi/180)*180/np.pi
    angle_adjust = np.array(np.degrees(
        np.tan(np.radians((PA[0]-PA)*np.cos(np.radians(inc))))), dtype=float)

    if disk == 1:
        Template["INCL"] = "INCL = "+" ".join(str(e) for e in inc)
        Template["PA"] = "PA = "+" ".join(str(e) for e in PA)
        try:
            phase = float(Template["AZ1P"].split('=')[1])
        except KeyError:
            phase = 0.
        Template.insert(
            "AZ1W", "AZ1P", f"AZ1P = {' '.join([f'{phase+e:10.5f}' for e in angle_adjust])}")
    else:
        Template["INCL_{:d}".format(disk)] = "INCL_{:d} =".format(
            disk)+" ".join(str(e) for e in inc)
        Template["PA_{:d}".format(disk)] = "PA_{:d} =".format(
            disk)+" ".join(str(e) for e in PA)
        try:
            phase = float(Template["AZ1P_{:d}".format(disk)].split('=')[1])
        except KeyError:
            phase = 0.
        Template.insert("AZ1W_{:d}".format(disk), "AZ1P_{:d}".format(
            disk), "AZ1P_{:d} = 180.".format(disk))
        Template.insert("AZ1W_{:d}".format(disk), "AZ1P_{:d}".format(
            disk), "AZ1P_{:d} =".format(disk)+" ".join(str(phase+(e)) for e in angle_adjust))
    return PA, inc, phirings


create_warp.__doc__ = f'''
NAME:
   create_warp

PURPOSE:
    Bend the PA and the INCL in the way a linear change in the angular momentum vector would cause.

CATEGORY:
   agc

INPUTS:
    Radii = radii of the rings
    PA = base PA
    inclination = base inclination
    warp_change = change in the angular momentum vector
    warp_radii = start and end radius of the warp vecotor input.
   
OPTIONAL INPUTS:
    disk = 1
    Are we creating in the approaching (1) or receding side (2)


OUTPUTS:

OPTIONAL OUTPUTS:

PROCEDURES CALLED:
   Unspecified

NOTE:
'''

# afunction for creating inhomogeneities acros the disk intrinsicially


def create_inhomogeneity(mass, SNR, disks=1.):
    if len(disks) == 0:
        disks = [disks]

    for i in range(len(disks)):
        ndisk = int(Template["NDISKS"].split('=', 1)[1])
        ndisk += 1
        last_add = copy_disk(disks[i], ndisk)
        sbr = np.array(Template["SBR_{:d}".format(ndisk)].split(
            '=', 1)[1].strip().split(" "), dtype=np.float64)
        rad = np.array(Template["RADI"].split('=', 1)[
                       1].strip().split(" "), dtype=np.float64)
        rad = rad*4.*np.pi/max(rad)
        # This part is tricky as we do not want negative emission but cutting out the negatives would result in added flux
        #Hence these profiles should always be smaller than the sbr
        newprof = np.sin(rad)*sbr*np.log10(mass)/11.*8./SNR*0.9
        #print("We are modifying by factor {} for mass {} and SNR {}".format(np.log10(mass)/11.*8./SNR, np.log10(mass),SNR))
        nextprof = -1*newprof

        req_flux = abs(np.sum(newprof*np.pi*rad*2)/200.)
        if req_flux == 0:
            req_flux = 1e-5
        Template["SBR_{:d}".format(ndisk)] = "SBR_{:d} = ".format(
            ndisk)+" ".join(str(e) for e in newprof)
        Template["CFLUX_{:d}".format(ndisk)] = "CFLUX_{:d} = {}" .format(
            ndisk, req_flux)
        ndisk = ndisk+1
        last_add = copy_disk(disks[i], ndisk)
        Template["SBR_{:d}".format(ndisk)] = "SBR_{:d} = ".format(
            ndisk)+" ".join(str(e) for e in nextprof)
        Template["CFLUX_{:d}".format(ndisk)] = "CFLUX_{:d} = {}" .format(
            ndisk, req_flux)
    #print("We will create inhomogeneities on top of disk(s) ="+" ".join(str(e) for e in [disks]))

def clean_up_casa(work_dir,sim_observe_graphics = 'none'):
    print('We are starting to remove the unnecessary products of the casa corrupt sequence')
    os.mkdir(f'{work_dir}Casa_Log')
    #check that it is created
    succes = os.path.isdir(f'{work_dir}Casa_Log')
    if not succes:
        time.sleep(1)
        succes = os.path.isdir(f'{work_dir}Casa_Log')
        if not succes:
            raise RunningError(
                f'We failed to create the directory {work_dir}Casa_Log')
    if sim_observe_graphics == 'file':
        os.system(f'mv {work_dir}sim_data/*.png {work_dir}Casa_Log/')
    os.system(f'mv {work_dir}casa*.log {work_dir}Casa_Log/')
    shutil.move(f'{work_dir}Observation_Overview.txt',f'{work_dir}Casa_Log/Observation_Overview.txt')
    shutil.move(f'{work_dir}casa_corrupt.py',f'{work_dir}Casa_Log/casa_corrupt.py')

    files_to_remove = ['testsmooth_cube.fits']

    for files in files_to_remove:
        os.remove(f'{work_dir}{files}')

    directories_to_remove = ['sim_data','in_cube.image','Final_Cube.image','mask.image','mask.image.tmp','mask.image.tmp2', 'Uni_Cube.image',\
        'Uni_Cube.psf','Uni_Cube.sumwt','Uni_Cube.mask','Uni_Cube.model','Uni_Cube.pb','Uni_Cube.residual']
    for dir in directories_to_remove:
        succes = os.path.isdir(f'{work_dir}{dir}')
        if succes:
            cf.delete_directory(f'{work_dir}{dir}')
        succes = os.path.isdir(f'{work_dir}{dir}')
        if succes:
            time.sleep(1)
            succes = os.path.isdir(f'{work_dir}{dir}')
            if succes:
                raise RunningError(
                    f'We failed to delete the directory {work_dir}{dir}')
clean_up_casa.__doc__ = f'''
NAME:
   clean_up_casa

PURPOSE:
   delete all the casa stuff

CATEGORY:
   agc

INPUTS:
   work_dir = directory to clean


OPTIONAL INPUTS:

OUTPUTS:

OPTIONAL OUTPUTS:

PROCEDURES CALLED:
   Unspecified

NOTE: T
'''


def create_mask(work_dir, beam, casa=False):
    # First open the model cube

    dummy = fits.open(work_dir+'/Unconvolved_Cube.fits', uint=False,
                      do_not_scale_image_data=True, ignore_blank=True)
        #fits.writeto(work_dir+'/Unconvolved_Cube.fits',dummy[0].data,original_hdr, overwrite = True)
    # downgrade the velocity resolution to the one we want
    tmp_cube = np.zeros([int(len(dummy[0].data[:, 0, 0])/3),
                        len(dummy[0].data[0, :, 0]), len(dummy[0].data[0, 0, :])])
    for n in range(1, len(dummy[0].data[:, 0, 0]), 3):
        tmp_cube[int((n-1)/3), :,
                     :] = np.mean(dummy[0].data[n-1:n+2, :, :], axis=0)
    #This line updates the NAXiS3 keyword
    dummy[0].data = tmp_cube
    # reset the header
    dummy[0].header['CDELT3'] = dummy[0].header['CDELT3']*3.
    dummy[0].header['CRPIX3'] = dummy[0].header['CRPIX3']/3.
    dummy[0].header['ALTRPIX'] = dummy[0].header['CRPIX3']
    #We do not want to do the casa corruption at high resolution as it takes too long
    if casa:

        fits.writeto(f'{work_dir}Unconvolved_Cube.fits',
                     dummy[0].data, dummy[0].header, overwrite=True)
        pixperbeam = cf.get_beam_area_in_pixels(
            dummy[0].header, beam=[x/3600. for x in beam])
    # In order to create a cleaning mask we smooth to the beam size and cut at 1e-5 Jy/beam
    #   Calculate the sigma's from the required beam size
    sigma = [(beam[0]/abs(dummy[0].header['CDELT1']*3600.))/(2*np.sqrt(2*np.log(2))),
              (beam[1]/abs(dummy[0].header['CDELT2']*3600.))/(2*np.sqrt(2*np.log(2)))]
    #
    #our BPA is set to 0 which means the major axis smoothing should be in DEC axis and the minor on the RA
    smooth = gaussian_filter(dummy[0].data, sigma=(
        0, sigma[0], sigma[1]), order=0)
    if casa:
        #If we are doing the casa corruption we want to do the proper conservation of brightness
        #In casa of gaussian corruption this is done after combining the model and the noise
        smooth = smooth*pixperbeam

        # We want the mean signal in the smoothed cube
        fits.writeto(f'{work_dir}testsmooth_cube.fits',
                     smooth, dummy[0].header, overwrite=True)
    mean_signal = cf.get_mean_flux(smooth)
    mask = cf.get_mask(smooth)
    mask_header = copy.deepcopy(dummy[0].header)
    mask = np.array(mask, dtype='float32')
    mask_header['BITPIX'] = -32

    fits.writeto(work_dir+'/mask.fits', mask, mask_header, overwrite=True)
    hdr = copy.deepcopy(dummy[0].header)
    data = copy.deepcopy(dummy[0].data)
    dummy.close()

    return mean_signal, hdr, data


create_mask.__doc = f'''
NAME:
   create_mask

PURPOSE:
   Create the mask for the final cube from the unsmoothed data cube and calculate the average signal.

CATEGORY:
   agc

INPUTS:
   work_dir = directory where to put the mask
   beam = size of the required beam

OPTIONAL INPUTS:
   casa = When the casa corruption method is used the unconvolved cube is require at low velocity resolution.
          Thus when True it is written to work_dir

OUTPUTS:
    the mask is written to disk and the final uncorrupted signal is returned in Jy/beam
    mean_signal = uncorrupted signal is returned in Jy/beam with the correct beam
    hdr = hdr is the header of the mask such that it can be used for further calculations
    data = is the unconvolved signal cube but the velocity axis is cut down to the final axis.

OPTIONAL OUTPUTS:

PROCEDURES CALLED:
   Unspecified

NOTE:
'''


def vel_to_freq(hdr):
    '''obatin the frequency of the m/s cube'''
    central_freq = HI_rest_freq*(1-float(hdr['CRVAL3'])/c)
    cdelt_freq = -HI_rest_freq*float(hdr['CDELT3']/1000.)/c_kms
    top_freq = central_freq+(hdr['NAXIS3']-1-hdr['CRPIX3']-1)*cdelt_freq
    low_freq = central_freq-(hdr['CRPIX3']-1)*cdelt_freq
    return central_freq, cdelt_freq, top_freq, low_freq

def corrupt_casa(work_dir, beam, SNR, maindir,sim_observe_graphics = 'none'):

    # the casa tasks ooze memory so they should be  run in a subprocess
    '''Corrupt our artifical galaxy with a casa's sim observe routines'''

    mean_signal, hdr, data = create_mask(work_dir, beam, casa=True)
    #This next line should be commented as it is merely for testing
    # In order to corrupt we need to know the average signal.
    # we do this taking the mean in each chaneel above a tenth of the max and then take the mean of that profile
    # This is the noise in the final cube
    #!!!! Is this correct or should we conserve for brightness??
    noise = mean_signal/SNR
    print(f'We use this noise estimate {noise}')
    # We need to know the location in the sky dec > 30 we use WSRT 30 > dec > -30 --> VLA -30 > dec  --> Atca
    RA, DEC = cf.convertRADEC(
        hdr['CRVAL1']+hdr['CDELT1'], hdr['CRVAL2']+hdr['CDELT2'])
    #mask sure previously used products are gone
    os.system('rm -Rf in_cube')
    vel_frame = 'LSRK'
    # flat noise = weighted by sqrt(weight), flatsky = weighted by weight, pbsquare=no normalisation
    skynormalistation = 'flatnoise'
    #weight is  sum ( square ( pb ) ), pb =  Primary beam calculated as sqrt ( xxx.weight )
    inframe = fits.getval(f'{work_dir}Unconvolved_Cube.fits', 'SPECSYS')
    fits.setval(f'{work_dir}Unconvolved_Cube.fits', 'SPECSYS', value=vel_frame)
    # make some names to use
    msname = f'{work_dir}sim_data.ms'
    casa_image = f'{work_dir}in_cube.image'
    casa_mask = f'{work_dir}mask.image'
    #get the header frequency informatio
    freq, cdelt, up_band, low_band = vel_to_freq(hdr)
    print(freq, cdelt, up_band, low_band)

    #!!!!!!!!!!!!!!!!!!!
    with open(f'{work_dir}casa_corrupt.py','w') as casa_program:
        casa_program.write(f'''
# Read the model into a casa format
def main():
    from casatools import simulator,  ctsys, measures, table, msmetadata, synthesisutils, ms
    from casatasks import tclean,  imhead, exportfits, flagdata, simanalyze, \\
                          importfits, listobs, vishead, imsmooth, simobserve
    from casatasks.private import simutil
    import os
    import numpy as np
    importfits(fitsimage='{work_dir}Unconvolved_Cube.fits', imagename='{casa_image}',            # Name of input image FITS file # Name of output CASA image # If fits image has multiple coordinate ,# If its file contains multiple images,Set blanked pixels to zero (not NaN)
                whichrep=0, whichhdu=-1, zeroblanks=True,
                overwrite=True, defaultaxes=True, defaultaxesvalues=['{RA}', '{DEC}','{up_band}Hz', 'I'])


    # We have options to use different telescopes but lets's stick with the VLA
    if {hdr['CRVAL2']} > 90:
        ant_list = 'WSRT.cfg'
    elif {hdr['CRVAL2']} > -30:
        #If the requested beam is small we need to use a configuration
        if {beam[1]} < 8.:
            #nat beam ~ 2."
            ant_list = 'vla.a.cfg'
            ant_ext = 'vla.a'
        elif {beam[1]} < 22:
            #nat beam ~ 6"
            ant_list = 'vla.b.cfg'
            ant_ext = 'vla.b'
        elif {beam[1]} < 50:
            #nat beam ~ 20"
            ant_list = 'vla.c.cfg'
            ant_ext = 'vla.c'
        else:
            #nat beam ~75"
            ant_list = 'vla.d.cfg'
            ant_ext = 'vla.d'
    else:
        #if at low declination we need to use atca
        ant_list = 'atca_6c.cfg'

    print(f"We selected {{ant_list}} as the antenna layout")
    # Instantiate all the required tools for the simulation
    sm = simulator()
    msmd = msmetadata()
    me = measures()
    ms = ms()
    #ms =measurementset()
    mysu = simutil.simutil()
    su = synthesisutils()
    tb = table()
    print(f"starting the simulation")
    # First clean up
    os.system(f'rm -rf {msname}')
    # Open the simulator
    # We will use simobserve as it works better
    #simobserve need to work in the directory
    sm.setvp(dovp=False)


    simobserve(
        #  simobserve :: visibility simulation task
        project='sim_data',  # root prefix for output file names
        skymodel='{casa_image}',  # model image to observe
        complist='',  # componentlist to observe
        setpointings=True,
        integration='20s',  # integration (sampling) time
        inbright='{np.max(data)}Jy/pixel',
        indirection=f'J2000 {RA} {DEC}',
        incenter='{freq}Hz',
        inwidth='{cdelt}Hz',
        incell='{hdr["CDELT2"]*3600.}arcsec',
        # observation mode to simulate [int(interferometer)|sd(singledish)|""(none)]
        obsmode='int',
        outframe='{vel_frame}',
        antennalist=f'{{ant_list}}',  # interferometer antenna position file
        # hour angle of observation center e.g. "-3:00:00", "5h", "-4.5" (a number without units will be
        hourangle='transit',
        #   interpreted as hours), or "transit"
        totaltime='12h',  # total time of observation or number of repetitions
        thermalnoise='',  # No noise to be able to add the correct noise""]
        # display graphics at each stage to [screen|file|both|none]. Have to be off to be able to run in screen.
        graphics= '{sim_observe_graphics}',
        verbose=False,
        overwrite=True)  # overwrite files starting with $project
    #There is a memeory leak in simobserve in  1109 sm.predict(imagename=newmodel,complist=complist)

    antennalist = os.path.join(ctsys.resolve("alma/simmos"), ant_list)
    ## Fictitious telescopes can be simulated by specifying x, y, z, d, an, telname, antpos.
    ##     x,y,z are locations in meters in ITRF (Earth centered) coordinates.
    ##     d, an are lists of antenna diameter and name.
    ##     telname and obspos are the name and coordinates of the observatory.
    (x, y, z, d, an, an2, telname, obspos) = mysu.readantenna(antennalist)
    #calculate the number of baselines
    no_integrations = 12.*3600./20.
    source = '{work_dir.split('/')[-2]}'

    sm.setvp(dovp=False)

    image_size = [su.getOptimumSize({int(hdr['NAXIS1'])}),
                    su.getOptimumSize({int(hdr['NAXIS2'])})]
    tclean(vis=f'{work_dir}sim_data/sim_data.{{ant_ext}}.ms',
        usemask='user',
        restart=False,
        imagename='{work_dir}Uni_Cube',
        niter=0,
        threshold='1Jy/beam',
        mask='{casa_mask}',
        normtype='{skynormalistation}',
        datacolumn='observed',
        imsize=image_size,
        smallscalebias=0.6,
        cell=["{abs(hdr['CDELT1']*3600.)}arcsec",
                     "{abs(hdr['CDELT2']*3600.)}arcsec"],
        pblimit=-1.0,
        pbmask=0.0,
        pbcor=False,
        restoringbeam='common',
        specmode='cube',
        start=0,
        cfcache='sim_predict.cfcache',
        outframe='{vel_frame}',
        width=1,
        calcres=False,
        calcpsf=True,
        restfreq='{HI_rest_freq}Hz',
        veltype='radio',
        field='0',
        weighting='natural',
        gridder='standard',
        deconvolver='clark',
        interactive=False,
        parallel=True
    )
    #Tclean also has a memory leak.
    #Get the uniform weighted beam size

    summary = imhead(imagename='{work_dir}Uni_Cube.psf', mode='summary')

    uniform_beam = [summary['perplanebeams']['beams']['*4']['*0']['major']['value'],
                    summary['perplanebeams']['beams']['*4']['*0']['minor']['value'],
                    summary['perplanebeams']['beams']['*4']['*0']['positionangle']['value']]
    print(
        f'We find the following beam from the uniform weighted image {{uniform_beam}}')
    print(f"Adding noise to the observations")

    print(f' For the beam {beam} We are increasing the original noise({noise}) with {{ np.mean(uniform_beam[:2])/{np.mean(beam[:2])} }} to obtain untapered noise.')
    no_baselines = len(an)*(len(an)-1.)/2.
    no_pol = 2
    noise_at_uni = {noise}*np.mean(uniform_beam[:2])/np.mean({beam[:2]})
    visnoise = noise_at_uni * \
        np.sqrt(no_pol)*np.sqrt(no_baselines)*np.sqrt(no_integrations)*1.2

    print(f'Corrupting the visbilities with noise of {{visnoise}} Jy')
    sm.openfromms(f'{work_dir}sim_data/sim_data.{{ant_ext}}.ms');
    sm.setvp(dovp=False)
    sm.setseed(50)
    sm.setnoise(mode='simplenoise', simplenoise=f'{{visnoise}}Jy');
    #sm.setgain(mode='fbm',amplitude=0.025) #Minimal amplitude errors
    sm.corrupt();
    sm.close();

    # read mask
    importfits(fitsimage='{work_dir}mask.fits', imagename='{casa_mask}',            # Name of input image FITS file # Name of output CASA image
                # If fits image has multiple coordinate ,# If its file contains multiple images,Set blanked pixels to zero (not NaN)
                whichrep=0, whichhdu=-1, zeroblanks=True,
                overwrite=True, defaultaxes=True,
                defaultaxesvalues=['{RA}', '{DEC}', '{up_band}Hz', 'I'])
    #Create the final Cube
    imhead(imagename='{casa_mask}', mode='put', hdkey='telescope', hdvalue=telname)
    imhead(imagename='{casa_mask}', mode='put',
           hdkey='date-obs', hdvalue='2019/10/4/00:00:00')
    # use simanalyze to produce an image
    #cleaning our visbilities
    listobs(vis=f'{work_dir}sim_data/sim_data.{{ant_ext}}.ms',
            listfile=f'{work_dir}Observation_Overview.txt', verbose=True, overwrite=True)

    tclean(vis=f'{work_dir}sim_data/sim_data.{{ant_ext}}.ms',
        usemask='user',
        restart=True,
        imagename=f'{work_dir}Uni_Cube',
        #imagename=f'{work_dir}Final_Cube',
        niter=1000,
        threshold=f'{{noise_at_uni/4.}}Jy/beam',
        mask='{casa_mask}',
        normtype='{skynormalistation}',
        datacolumn='observed',
        imsize=image_size,
        smallscalebias=0.6,
        cell=[f"{abs(hdr['CDELT1']*3600.)}arcsec",
                     f"{abs(hdr['CDELT2']*3600.)}arcsec"],
        pblimit=-1.0,
        pbmask=0.0,
        pbcor=False,
        restoringbeam='common',
        specmode='cube',
        start=0,
        cfcache='sim_predict.cfcache',
        outframe='{vel_frame}',
        width=1,
        calcres=True,
        calcpsf=False,
        restfreq=f'{HI_rest_freq}Hz',
        veltype='radio',
        field='0',
        weighting='natural',
        gridder='standard',
        deconvolver='clark',
        interactive=False,
        parallel=True
    )
    if {beam[0]} < uniform_beam[0] or {beam[1]} < uniform_beam[1]:
        print(f'!!!!!!!!!!!!!!!!!!!!!!We can not make the beam as small as you want it. We are simply copying the naturally weighted cube.')
        print(f'cp -r {work_dir}Uni_Cube.image {work_dir}Final_Cube.image')
        shutil.move(f'{work_dir}Uni_Cube.image', f'{work_dir}Uni_Cube.image')
    else:
        #imsmooth(imagename=f'{work_dir}sim_data/sim_data.{{ant_ext}}.image', outfile=f'{work_dir}Final_Cube.image',
        imsmooth(imagename=f'{work_dir}Uni_Cube.image', outfile=f'{work_dir}Final_Cube.image',
                    overwrite=True, major=f'{beam[0]}arcsec', minor=f'{beam[1]}arcsec',
                    pa=f'{beam[2]}deg', targetres=True)

    imhead(imagename=f'{work_dir}Final_Cube.image',
           mode='put', hdkey='object', hdvalue='AGC Galaxy')

    exportfits(imagename=f'{work_dir}Final_Cube.image',  # Name of input CASA image
               # Name of output image FITS file
               fitsimage=f'{work_dir}Convolved_Cube.fits',
               # Use velocity (rather than frequency) as spectral axis
               velocity=True,
               # Use the optical (rather than radio) velocity convention
               optical=False,
               bitpix=-32,  # Bits per pixel
               # Minimum pixel value (if minpix > maxpix, value is automatically determined)
               minpix=0,
               # Maximum pixel value (if minpix > maxpix, value is automatically determined)
               maxpix=-1,
               overwrite=True,  # Overwrite pre-existing imagename
               dropstokes=True,  # Drop the Stokes axis?
               stokeslast=False,  # Put Stokes axis last in header?
               history=True,  # Write history to the FITS image?
               dropdeg=True)
if __name__ == '__main__':
    main()
        ''')
    del data
    current_run = subprocess.Popen(['python', "casa_corrupt.py"],
                           stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                           cwd=f"{work_dir}", universal_newlines=True)
    casa_run, casa_warnings_are_annoying = current_run.communicate()
    print(casa_run)
    print(casa_warnings_are_annoying)

    fits.setval(f'{work_dir}Unconvolved_Cube.fits', 'SPECSYS', value=inframe)
    dummy = fits.open(work_dir+'/Convolved_Cube.fits', uint=False,
                      do_not_scale_image_data=True, ignore_blank=True)

    #We need to regrid to our desired resolution
    for key in hdr:
        try:
            print(
                f'For {key} we have in the header {hdr[key]} and {dummy[0].header[key]} in the original')
        except:
            pass
    newdummy = cf.regrid_array(dummy[0].data, Out_Shape=((int(hdr['NAXIS3']),
        int(hdr['NAXIS2']), int(hdr['NAXIS1']))))
    for ext in ['1', '2']:
        achieved_regrid = dummy[0].data.shape[int(
            ext)] / newdummy.shape[int(ext)]
        dummy[0].header[f'CDELT{ext}'] = dummy[0].header[f'CDELT{ext}'] * \
            achieved_regrid
        dummy[0].header[f'CRPIX{ext}'] = dummy[0].header[f'CRPIX{ext}'] / \
            achieved_regrid
        dummy[0].header[f'NAXIS{ext}'] = newdummy.shape[int(ext)]

    # We cut 20 pixels around the edges
    #fits.writeto(work_dir+'/Convolved_Cube.fits',newdummy,hdr, overwrite = True)
    #Cut empty channels at the start and beginning
    dummymask = fits.open(work_dir+'/mask.fits', uint=False,
                          do_not_scale_image_data=True, ignore_blank=True)
    clean_smooth = fits.open(f'{work_dir}testsmooth_cube.fits')
    while dummymask[0].data.shape[0] > newdummy.shape[0]:
        print(f'We are removing a back channel from the mask')
        dummymask[0].data = dummymask[0].data[:-1]
        clean_smooth[0].data = clean_smooth[0].data[:-1]

    while np.sum(newdummy[-1]) == 0:
        newdummy = newdummy[:-1]
        clean_smooth[0].data = clean_smooth[0].data[:-1]
        dummymask[0].data = dummymask[0].data[:-1]
        dummy[0].header['NAXIS3'] = dummy[0].header['NAXIS3']-1

    while np.sum(newdummy[0]) == 0:
        newdummy = newdummy[1:]
        dummymask[0].data = dummymask[0].data[1:]
        clean_smooth[0].data = clean_smooth[0].data[1:]
        dummy[0].header['NAXIS3'] = dummy[0].header['NAXIS3']-1
        dummy[0].header['CRPIX3'] = dummy[0].header['CRPIX3']-1
    dummy[0].header['SPECSYS'] = inframe
    outnoise = (np.std(newdummy[0])+np.std(newdummy[-1]))/2.
    getscaling = copy.deepcopy(clean_smooth[0].data)
    mean_signal = cf.get_mean_flux(
        clean_smooth[0].data, Mask=dummymask[0].data)
    getscaling[dummymask[0].data < 0.5] = -0.1
    getscaling[getscaling < 0.5*mean_signal] = -0.1
    #And scale the final cube to make sure that we maintained the total flux
    #fits.writeto(work_dir+'/Scaling_Mask.fits', getscaling,
    #             dummy[0].header, overwrite=True)
    #fits.writeto(f'{work_dir}testsmooth_cube.fits', clean_smooth[0].data,
    #             dummy[0].header, overwrite=True)
    mean_scale = np.mean(
        clean_smooth[0].data[getscaling > 0.]/newdummy[getscaling > 0.])

    print(f'We scale the final cube by {mean_scale} to retain the M_HI')
    newdummy *= mean_scale
    fits.writeto(work_dir+'/Convolved_Cube.fits', newdummy,
                 dummy[0].header, overwrite=True)
    dummy[0].header['BITPIX'] = 16
    dummymask[0].data = np.array(dummymask[0].data, dtype=int)
    fits.writeto(work_dir+'/mask.fits',
                 dummymask[0].data, dummy[0].header, overwrite=True)
    dummy.close()
    dummymask.close()
    del hdr
    del newdummy
    del getscaling

    #clean up the mess
    clean_up_casa(work_dir,sim_observe_graphics=sim_observe_graphics)

    # As these are a lot of system operations on big files let's give the system time to catch up
    finished = False
    counter = 0
    while not finished:
        ms_exists = os.path.isdir(f'{work_dir}Final_Cube.image')
        final_cube_exists = os.path.isfile(f'{work_dir}Convolved_Cube.fits')
        if not ms_exists and final_cube_exists:
            finished = True
        else:
            counter += 1
            time.sleep(1)
        if counter in [30,60,90]:
            #try again
            clean_up_casa(work_dir,sim_observe_graphics=sim_observe_graphics)
        if counter > 120:
            raise RunningError(
                f'Something failed in CASA corruption. please check {work_dir} and the log there.')

corrupt_casa.__doc__ = f'''
NAME:
   corrupt_casa

PURPOSE:
   Corrupt our artifical galaxy with a casa's sim observe routines

CATEGORY:
   agc

INPUTS:
   work_dir = directory where to put the mask
   beam = size of the required beam
   SNR = is the requested SNR

OPTIONAL INPUTS:

OUTPUTS:
   the Convolved and corrupted cube is writen to disk


OPTIONAL OUTPUTS:

PROCEDURES CALLED:
   Unspecified

NOTE: The final SNR is a guesstimate.
'''

def corrupt_gauss(work_dir, beam, SNR, no_corrupt=False):
    '''Corrupt our artifical galaxy with the correct noise cube such that the average SNR over the galaxy is the required one'''
    mean_signal, hdr, data = create_mask(work_dir, beam)
    # Calculate the area of the beam in arcsec
    #
    #beamarea=(np.pi*abs(beam[0]*beam[1]))/(4.*np.log(2.))
    # Convert arcsec to pixels
    #pixperbeam=beamarea/(abs(hdr['CDELT1']*3600.)*abs(hdr['CDELT2']*3600.))
    pixperbeam = cf.get_beam_area_in_pixels(hdr, beam=[x/3600. for x in beam])
    # Make an the size of the model with random values distributed as the noise
    # The noise we want is the mean signal divided by the signal to noise ratio
    # that value needs to be deconvolved so from https://en.wikipedia.org/wiki/Gaussian_blur
    # The formula is uncited but works
    sigma = [(beam[0]/abs(hdr['CDELT1']*3600.))/(2*np.sqrt(2*np.log(2))),
              (beam[1]/abs(hdr['CDELT2']*3600.))/(2*np.sqrt(2*np.log(2)))]
    noisescl = (mean_signal/SNR*np.mean(sigma)*2*np.sqrt(np.pi))
    #print("This our mean signal {} and noise {} in Jy/pixel. ".format(mean_signal,noisescl))

    if beam[2] != 0.:
        #As we are going to rotatet the cube we should first extend it
        shift = int(abs(np.sin(np.radians(beam[2])))*hdr['NAXIS1']/2.+5)
        Pix_Extend = [shift, shift]
        data = np.pad(data, [[0, 0], Pix_Extend, Pix_Extend], 'constant')
        data = cf.rotateCube(
            data, beam[2], [hdr['CRPIX1']+shift, hdr['CRPIX2']+shift], order=1)

    # combine the two cubes
    if no_corrupt:
        noisedcube = data
    else:
        cuberms = np.random.normal(scale=noisescl, size=np.shape(data))
        noisedcube = data+cuberms
    # Smooth to the requred resolution

    #our BPA is set to 0 which means the major axis smoothing should be in DEC axis and the minor on the RA
    final = gaussian_filter(noisedcube, sigma=(0, sigma[0], sigma[1]), order=0)

    if beam[2] != 0.:
        #rotate back
        final_tmp = cf.rotateCube(
            final, -1*(beam[2]), [hdr['CRPIX1']+shift, hdr['CRPIX2']+shift], order=1)
        final = copy.deepcopy(final_tmp[:,
            shift:hdr['NAXIS2']
            + shift, shift:hdr['NAXIS1']
             + shift])
        final_tmp = []

    # to preserve brightness temperature this should be multiplied with the increase in beam area
    # which is the same as the amount of pixels in the beam as we go from 1 pixel area to an area the size of the beam which is assuming two gaussians so we need to correct with a difference factor of the area of the. Last factor is to correct for the fact that the unconvolved cube hassquare pixels not a circular beam.
    final = final*pixperbeam
    # And write this final cube to the directory
    hdr['BMAJ'] = beam[0]/3600.
    hdr['BMIN'] = beam[1]/3600.
    hdr['BPA'] = beam[2]
    hdr['BUNIT'] = 'Jy/beam'
    fits.writeto(work_dir+'/Convolved_Cube.fits', final, hdr, overwrite=True)


corrupt_gauss.__doc__ = f'''
NAME:
   corrupt_gauss

PURPOSE:
   Corrupt our artifical galaxy with the correct noise cube such that the average SNR over the galaxy is the required one

CATEGORY:
   agc

INPUTS:
   work_dir = directory where to put the mask
   beam = size of the required beam
   SNR = is the requested SNR

OPTIONAL INPUTS:

OUTPUTS:
   the Convolved and corrupted cube is writen to disk

OPTIONAL OUTPUTS:

PROCEDURES CALLED:
   Unspecified

NOTE:
'''


def get_HI_disk(Mass, output_directory=None):
    # About 18.5% is cosmic baryon fraction from planck (Omb=0.048, OmDM = 0.2589). In a galactic halo this is roughly 70%-90% of the cosmic mean (See Crain et al. 2006, Pezzulli 2019)
    bary_frac = 0.1854
    diff = 0.9
    counter = 0
    # if the baryon mass is higher the 2* a typical gass mass disk (i.e. M_HI = 1e9 Msol)
    # then we subtract this disk from the baryonic fraction and count the remainder as the initial stellar mass guess
    if bary_frac*Mass > 2.8*10**9:
        m_star = bary_frac*Mass-1.4*10**9
    else:
        # else half baryonic in stars
        m_star = bary_frac*Mass/2.

    # and convert stellar mass to an HI Mass following van Driel (2016)
    m_star_prev = 1
    while abs(m_star_prev-m_star)/m_star > 0.1:
        M_HI = m_star*10**(-0.43*np.log10(m_star)+3.75)
        # If the total gas mass is higher than the calculated baryonic mass then modify
        if 1.4*M_HI > bary_frac*Mass:
            if (1.4*M_HI+m_star) > Mass:
                #MHI is the total mass  - the expected cosmic mean stellar mass
                M_HI = (Mass-m_star)/1.4
                bary_frac = 1.
            else:
                #The baryonic fraction is raised to match the current one
                bary_frac = (1.4*M_HI+m_star)/(Mass)
        m_star_prev = copy.deepcopy(m_star)
        m_star = bary_frac * Mass-M_HI*1.4

    #print("This is MHI {:.2e} and mstar {:.2e}".format(MHI,m_star,m_star_prev))
    # the mass leads to the radius at which we need to hit 1 M/pc^2 from Wang (2016) in kpc
    R_HI = 10**(0.506*np.log10(M_HI)-3.293)/2.

    # McGaugh 2014  M_*/Lk = 0.6
    M_K = np.log10(m_star/0.6)*-2.5+3.29
    # Ponomareva 2017 has
    # log(LT,b,i) = (3.7 ± 0.11) × log(2Vflat) + 1.3 ± 0.3 for 3.6 mu  and then in Table 3 the values for Ks band
    v_circ_TF = 10**((np.log10(m_star/0.6)-1.22)/3.81)/2.

    #  Let's check that this roughly gives our DM mass at the virial radius in pc
    v_circ_NFW = calc_vc_NFW(Mass, M_HI, m_star, R_HI)
    #our final velocity is an average between TF and NFW
    V_HI = (v_circ_TF+v_circ_NFW)/2.

    DynMass = R_HI*10**3*V_HI**2/G_agc
    if output_directory:
        DynMassNFW = R_HI*10**3*v_circ_NFW**2/G_agc
        with open(f'{output_directory}Fractions_and_Masses.txt', 'a') as file:
            file.write(f'''The input Mass = {Mass: .2e} and the retrieved NFW Dynamical mass = {DynMassNFW: .2e} and Dynamical Mass based on v_circ = {DynMass: .2e}.
The current the baryon fraction = {bary_frac: .5f}\n''')
    return V_HI, R_HI, M_HI, M_K, DynMass


get_HI_disk.__doc__ = f'''
NAME:
    get_HI_disk

PURPOSE:
    Obtaining the V_HI, R_HI and M_HI to build the HI disk from.

CATEGORY:
   agc

INPUTS:
    Mass = the requested DM mass

OPTIONAL INPUTS:

OUTPUTS:
     V_HI, R_HI and M_HI and the absolute magnitude of the stellar disk.

OPTIONAL OUTPUTS:

PROCEDURES CALLED:
   Unspecified

NOTE:
'''


def get_sparcs_fit(v_end, radius):
    '''Calculate the R_0, gamma and beta parameters'''

    v_flat = 11.+0.76*v_end
    R_0 = 10*np.exp(-1*radius**2/10.**2)+2.25+0.08*radius
    beta = 0.5*(1.14-0.0006*v_end)+0.5*(0.7*np.arctan((radius - 32.)/25.))
    if beta < 0.25:
        beta = 0.25
    gamma = 0.5*(gaussian_function(v_end, 1.5, 170., 60.)) + \
                 0.5*(3.7 - 0.003*radius)
    return v_flat, R_0, beta, gamma


get_sparcs_fit.__doc__ = f'''
NAME:
   get_sparcs_fit

PURPOSE:
    Calculate the R_0, gamma and beta parameters for Courteau's parameterization of the RC
    based on a the fitting of all SPARCS RCs

CATEGORY:
   agc

INPUTS:
   v_flat = the flat velocity of the mass distribution
   radius = the HI radius, used in R_0
OPTIONAL INPUTS:

OUTPUTS:
    R_0, gamma and beta

OPTIONAL OUTPUTS:

PROCEDURES CALLED:
   Unspecified

NOTE:
'''


def gaussian_function(x, amp, center, sigma):
    return amp*np.exp(-(x-center)**2/(2.*sigma**2))


gaussian_function.__doc__ = f'''
NAME:
   gaussian_function

PURPOSE:
    return a value on a gaussian function

CATEGORY:
   agc

INPUTS:
   x input value
   amp = amplitude of the gaussian function
   center = center of the gaussian
   sigma = sigma of the gaussian function

OPTIONAL INPUTS:

OUTPUTS:
    y value of the gaussian

OPTIONAL OUTPUTS:

PROCEDURES CALLED:
   Unspecified

NOTE:
'''

def one_galaxy(cfg, Current_Galaxy, Achieved):
    # Build a name and a directory where to stor the specific output
    print(f'''This is the galaxy we are creating.''')
    cf.print_base_galaxy(Current_Galaxy)
    name = set_name(Current_Galaxy)
    #print("This is the Current Mass = {:.1e}".format(Current_Galaxy.Mass))
    # Then we copy the original Template def
    global Template
    #dictionaries are mutible and hence there should not have to be declared global
    # But since python is the most alogical programming language ever invented
    # it does have to be declared because we  are in a function not in Global
    #Template=copy.deepcopy(Template_in)
    shutil.copy(f'{cfg.general.main_directory}/Input.fits',f'{cfg.general.main_directory}/{name}/Input.fits')
    #make sure it arrives
    counter = 0
    found = False
    while not found:
        input_exists = os.path.isfile(
            f'{cfg.general.main_directory}/{name}/Input.fits')
        if input_exists:
            found = True
        else:
            time.sleep(0.1)
            counter += 1
        if counter in [30,60,90]:
            shutil.copy(f'{cfg.general.main_directory}/Input.fits',f'{cfg.general.main_directory}/{name}/Input.fits')
        if counter > 120:
            raise RunningError(
                f'Something went wrong copying the input file to {name}. Please file an issue on Github')


    Template=cf.read_template_file('Template.def')

     # First set the beam to 0.
    Template["BMAJ"]="BMAJ = 0."
    Template["BMIN"]="BMIN = 0."
    Template["BPA"]="BPA = 0."

    #Then we need to build the Surface Brightnes profile
    SBRprof, Rad, sclength, MHI, Rad_HI, Vrot, sub_ring, molecular_profile, DynMass=build_sbr_prof(
        Current_Galaxy, symmetric = cfg.agc.symmetric)  # Column densities,Raii in kpc, Opt_scalelength in kpc, HI mass in M_solar
    setattr(Achieved, "HI_Radius", Rad_HI)
    setattr(Current_Galaxy, "HI_Mass", MHI)
    Achieved.Mass=DynMass
    setattr(Achieved, "Scalelength", sclength)
    #We need central coordinates the vsys will come from the required distance and hubble flow. The RA and dec should not matter hance it will be the only random component in the code as we do want to test variations of them
    Sky_Size=np.radians(
        Current_Galaxy.Res_Beam[0]*Current_Galaxy.Beams/3600.)
    Distance=(Rad_HI[0]/(np.tan(Sky_Size/2.)))/1000.

  
    vsys=Distance*H_0
    # RAhr, DEChr=cf.convertRADEC(*Current_Galaxy.Coord)
    # print("It's central coordinates are RA={} DEC={} vsys={} km/s".format(RAhr,DEChr,vsys))
    # Write them to the template
    for ext in ['', '_2']:
        Template[f"XPOS{ext}"]=f"XPOS{ext}= {Current_Galaxy.Coord[0]}"
        Template[f"YPOS{ext}"]=f"YPOS{ext}= {Current_Galaxy.Coord[1]}"
        Template[f"VSYS{ext}"]=f"VSYS{ext}= {vsys}"

    # The Distance, this is not required but can be usefull
    Template["DISTANCE"]="DISTANCE= {}".format(Distance)
    # With the distance we can also convert our radii
     # we need our radii in the corresponding arcsec
    Rad_arcsec=cf.convertskyangle(Rad, distance = Distance, physical = True)
     # then the number of rings to the total number of rings
    Template["NUR"]="NUR= {}".format(len(Rad))
    # The radii in arcsec
    Template["RADI"]="RADI = "+" ".join(str(e) for e in Rad_arcsec)
    # Then we need to get a starting radius for the warp.
    # the warp should start at the edge of the optical radius which is the HI scale length/0.6
    # which are about ~ 4 * h_r
    WarpStart=4.*sclength
    WarpEnd=Rad[np.where(SBRprof >= 4.98534620064e-05)[0][-1]]
    #If we demand a Warp but it is too small to take effect we lower the start radius
    if np.sum(Current_Galaxy.Warp) != 0. and (WarpEnd-WarpStart) < 0.1*WarpEnd:
        WarpStart = 0.95*WarpEnd
        WarpEnd = 1.1*WarpEnd
    setattr(Achieved, "Warp_Radius", WarpStart)
    # Write it to the Template
    Template["VROT"]="VROT = "+" ".join(str(e) for e in Vrot)
    Template["VROT_2"]="VROT_2 = "+" ".join(str(e) for e in Vrot)
    # We need a scale height and dispersion for each ring. They are coupled and hence they are both created in create_flare
    h_z, dispersion=create_flare(Rad, Vrot, Current_Galaxy.Dispersion,
                                   Current_Galaxy.Flare, Rad_HI[0], sub_ring, distance = Distance)

    # Finally we need to set the warping
    PA, inc, phirings=create_warp(Rad, Current_Galaxy.PA, Current_Galaxy.Inclination,\
                                   Current_Galaxy.Warp, [WarpStart, WarpEnd],
                                   debug=cfg.general.debug)
    if cfg.agc.symmetric:
        PA_2, inc_2, phirings_2=create_warp(Rad, Current_Galaxy.PA,\
            Current_Galaxy.Inclination, Current_Galaxy.Warp,\
            [WarpStart, WarpEnd], disk = 2)
    else:
        # we want an assymetric warp so we redo the PA and inc but swap the variation
        PA_2, inc_2, phirings_2=create_warp(Rad, Current_Galaxy.PA, \
            Current_Galaxy.Inclination, [Current_Galaxy.Warp[0]-Current_Galaxy.Warp[1],\
            Current_Galaxy.Warp[1]/2.+Current_Galaxy.Warp[0]/2.],\
            [WarpStart, WarpEnd], disk = 2)

    #If we want Radial Motions then they need to be inserted
    if Current_Galaxy.Radial_Motions != 0.:
          Template.insert("VROT", "VRAD", "VRAD = {}".format(
              Current_Galaxy.Radial_Motions))
          Template.insert("VROT_2", "VRAD_2", "VRAD_2 = {}".format(
              Current_Galaxy.Radial_Motions))

    # This comes from FAT. If I remember correctly this is the sine response to the channels *1.2/(2*SQRT(2*ALOG(2.))))
    # However in our input we want independent channels which means we should set this to 0.
    # Template["CONDISP"]="CONDISP = 0."
    if cfg.agc.channel_dependency == 'sinusoidal':
        Template["CONDISP"]="CONDISP = {}".format(
            Current_Galaxy.Channelwidth*1.2/(2*np.sqrt(2*np.log(2.))))
    elif cfg.agc.channel_dependency == 'hanning':
        Template["CONDISP"] = "CONDISP = {}".format(
            Current_Galaxy.Channelwidth*2./(2*np.sqrt(2*np.log(2.))))
    elif cfg.agc.channel_dependency == 'independent':
        Template["CONDISP"] = "CONDISP = 0."
    else:
        raise InputError(
            f"{cfg.agc.channel_dependency} is not an option for the channel dependency")
    setattr(Achieved, "Channel_Dep", cfg.agc.channel_dependency)
    # We need to set the input and output cube
    Template["INSET"] = "INSET = Input.fits"
    Template["OUTSET"] = "OUTSET = Unconvolved_Cube.fits"
    #Some tirific varaiables
    Template["LOOPS"] = "LOOPS = 0 "
    # We need models with about 3 million particles but not more as it takes too long
    TotFlux = MHI/(2.36e5*(Distance)**2)
    Template["CFLUX"] = "CFLUX = {}".format(TotFlux/5e6)
    Template["CFLUX_2"] = "CFLUX_2 = {}".format(TotFlux/5e6)
    Template.insert("RMS", "NDISKS", "NDISKS = 2")
    Template["TIRDEF"] = "TIRDEF = ModelInput.def"
    Template["GR_DEVICE"] = "GR_DEVICE = "
    Template["GR_CONT"] = "GR_CONT = "
    Template.insert("INSET", "ACTION", "ACTION = 1")
    Template["PROGRESSLOG"] = "PROGRESSLOG = "
    Template["LOGNAME"] = "LOGNAME= "

    #-------------------------------This finishes the basic disk the following are optional components----------------------

    # The possible arms
    if Current_Galaxy.Arms == 'Arms':
        phase, arm_brightness, arm_width = create_arms(
            Vrot, Rad, SBRprof, WarpStart=WarpStart, Bar=Current_Galaxy.Bar)
        phase, arm_brightness, arm_width = create_arms(
            Vrot, Rad, SBRprof, disk=2, WarpStart=WarpStart, Bar=Current_Galaxy.Bar)
        Achieved.Arms = 'Arms'
    else:
        Achieved.Arms = 'No_Arms'
    # A possible Bar
    if Current_Galaxy.Bar == 'Bar':
        bar_length = create_bar(
            Vrot, Rad, SBRprof, Template, WarpStart=WarpStart)
        Achieved.Bar = 'Bar'
    else:
        Achieved.Bar = 'No_Bar'
    # and possible inhomogeneities
    if cfg.agc.inhomogenous:
        inhomogeneity_amp = create_inhomogeneity(
            MHI, Current_Galaxy.SNR, disks=[1, 2])
    # we write the def files
    with open(f"{cfg.general.main_directory}{name}/tirific.def", 'w') as file:
        file.writelines([Template[key]+"\n" for key in Template])

    # So we need to  modify the input file to the correct coordinates else we'll get an empty cube

    dummy = fits.open(f"{cfg.general.main_directory}{name}/Input.fits",
                      uint=False, do_not_scale_image_data=True, ignore_blank=True)
    #first we do a check wether the BMAJ
    #is written correctly for the python
    #program. We set a bunch of generic header values.
    #First we make sure that it fits

    size = 2.*Rad_arcsec[-1] + (Current_Galaxy.Res_Beam[0])
    if Current_Galaxy.Beams < 6.:
        size = 2.*Rad_arcsec[-1] + (6.*Current_Galaxy.Res_Beam[0])
    #pix_size = (Rad_arcsec[1]-Rad_arcsec[0])
    pix_size = Current_Galaxy.Res_Beam[1]/5.
    if Current_Galaxy.Corruption == 'Casa_Sim':
        size += 20*pix_size
        #we want to ensure 4 pixels in the natural beam
        if Current_Galaxy.Res_Beam[1] < 8.:
            #nat beam ~2"
            if pix_size > 0.5:
                pix_size = 0.5
        elif Current_Galaxy.Res_Beam[1] < 22:
            if pix_size > 1.5:
                pix_size = 1.5
            #nat beam ~ 6"
        elif Current_Galaxy.Res_Beam[1] < 80:
            if pix_size > 5:
                pix_size = 5
            #nat beam ~ 20"
        else:
            #nat beam ~75"
            if pix_size > 18:
                pix_size = 18

    required_pixels = int(np.ceil(size/pix_size))
    vel_max = 2.5*np.max([Vrot*np.sin(inc*np.pi/180.)+4
                         * dispersion, Vrot*np.sin(inc_2*np.pi/180.)+4*dispersion])
    velpix = int(np.ceil(vel_max/Current_Galaxy.Channelwidth)*3)
    dummy[0].header['CRPIX1'] = np.floor(required_pixels/2.)
    dummy[0].header['CRPIX2'] = np.floor(required_pixels/2.)
    dummy[0].header['CRPIX3'] = np.floor(velpix/2.)
    # Stupid astropy doesn't account for minus in header of cdelt1 then cdelt has different precision
    tmp = int(pix_size/3600.*1e15)/1e15
    dummy[0].header['CDELT1'] = -1*tmp
    dummy[0].header['CDELT2'] = tmp

    dummy[0].header['CDELT3'] = Current_Galaxy.Channelwidth*1000./3.
    dummy[0].header['CRVAL1'] = Current_Galaxy.Coord[0]
    dummy[0].header['CRVAL2'] = Current_Galaxy.Coord[1]
    dummy[0].header['CRVAL3'] = vsys*1000.
    dummy[0].header['NAXIS1'] = required_pixels
    dummy[0].header['NAXIS2'] = required_pixels
    dummy[0].header['NAXIS3'] = velpix
    dummy[0].header['BMAJ'] = 0.
    dummy[0].header['BMIN'] = 0.

    try:
        del dummy[0].header['BLANK']
    except:
        pass
    # make the cube a typical size
    dummy2 = np.zeros(
        (velpix, required_pixels, required_pixels), dtype=np.float32)
    dummy2[int(np.floor(velpix/2.)), int(np.floor(required_pixels/2.)),
               int(np.floor(required_pixels/2.))] = 5
    fits.writeto(f"{cfg.general.main_directory}{name}/Input.fits", dummy2,
                 dummy[0].header, output_verify='silentfix+ignore', overwrite=True)
    dummy.close()
    dummy2 = []

    current_run = subprocess.Popen([cfg.general.tirific, "DEFFILE=tirific.def", "ACTION=1"],
                           stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                           cwd=f"{cfg.general.main_directory}{name}", universal_newlines=True)
    tirific_run, tirific_warnings_are_annoying = current_run.communicate()
    #print(tirific_run)
    if current_run.returncode == 1:
        pass
    else:
        print(tirific_warnings_are_annoying)
        raise TirificRunError(
            "AGC:Tirific did not execute properly. See screen for details")
    tirific_run=[]
    tirific_warnings_are_annoying = []

    # Make sure that our tirific output confirms with the FITs standard header
    fits_to_modify = f"{cfg.general.main_directory}{name}/Unconvolved_Cube.fits"
    unconvolved_cube = fits.open(fits_to_modify,
                     uint=False, do_not_scale_image_data=True, ignore_blank=True)
    hdr = unconvolved_cube[0].header
    data = unconvolved_cube[0].data

    freq, cdelt, up_band, low_band = vel_to_freq(hdr)
    hdr['BITPIX'] = -32
    hdr['SPECSYS'] = 'BARYCENT'
    hdr['OBJECT'] = 'AGC_GALAXY'
    hdr['INSTRUME'] = 'AGC'
    hdr['CTYPE3']='VELO-HEL'
    hdr['BMAJ']=0.
    hdr['BMIN']=0.
    hdr['BPA']=0.
    hdr['RESTFRQ']=HI_rest_freq
    hdr['BUNIT']='Jy/pixel'
    hdr['ALTRVAL']=freq
    hdr['ALTRPIX']=hdr['CRPIX3']


    fits.writeto(fits_to_modify,data, hdr, overwrite = True)


    unconvolved_cube.close()
    data= []
    hdr = []
    # Now we want to corrupt this cube with some realistic noise
    # For this we first want to get the noise we want in terms of Jansky per beam
    # we will define the SNR as the mean(Intensity)/noiselevel hence noise =mean(In)/SNR

    if Current_Galaxy.Corruption == 'Casa_Sim':
        fits.setval(fits_to_modify, 'CTYPE3', value = 'VRAD')
        from contextlib import redirect_stdout
        if cfg.agc.sim_observe_graphics:
            graphs= 'file'
        else:
            graphs= 'none'
        with open(f'{cfg.general.main_directory}{name}/Casa_Log.txt', 'w') as f:
            with redirect_stdout(f):
                corrupt_casa(f"{cfg.general.main_directory}{name}/",
                             Current_Galaxy.Res_Beam, Current_Galaxy.SNR, cfg.general.main_directory,
                             sim_observe_graphics=graphs)
        os.system(
            f"mv {cfg.general.main_directory}{name}/Casa_Log.txt {cfg.general.main_directory}{name}/Casa_Log/")
        fits.setval(fits_to_modify, 'CTYPE3', value = 'VELO-HEL')
        fits.setval(
            f"{cfg.general.main_directory}{name}/mask.fits", 'CTYPE3', value = 'VELO-HEL')
        fits.setval(
            f"{cfg.general.main_directory}{name}/Convolved_Cube.fits", 'CTYPE3', value = 'VELO-HEL')
        Achieved.Corruption='Casa_Sim'

    elif Current_Galaxy.Corruption == 'Gaussian' or Current_Galaxy.Corruption == 'No_Corrupt':
        no_corrupt= False
        if Current_Galaxy.Corruption == 'No_Corrupt':
            no_corrupt= True
        corrupt_gauss(f"{cfg.general.main_directory}{name}/",
                      Current_Galaxy.Res_Beam, Current_Galaxy.SNR,no_corrupt = no_corrupt)
        Achieved.Corruption=Current_Galaxy.Corruption
    else:
        raise RunningError(f'These boots (Corruption = {Current_Galaxy.Corruption}) are not made for running pyHIARD')


    mask=fits.open(f"{cfg.general.main_directory}{name}/mask.fits")
    Cube=fits.open(f"{cfg.general.main_directory}{name}/Convolved_Cube.fits")
    #Here we have no control over the BPA it is what it is.

    Current_Galaxy.Res_Beam[2]=Cube[0].header['BPA']
    maskr=mask[0].data
    sigma=(np.std(Cube[0].data[0])+np.std(Cube[0].data[-1]))/2.
    Cube_Clean=Cube[0].data
    Cube_Clean[maskr < 0.5]=0.
    maskr = []
    mask.close()
    pixperbeam=cf.get_beam_area_in_pixels(Cube[0].header)
    totalsignal=np.sum(Cube_Clean)/pixperbeam
    mass=2.36E5*Distance**2*totalsignal*Cube[0].header['CDELT3']/1000.
    mean_signal=cf.get_mean_flux(Cube_Clean)
    Cube_Clean = []
    Cube.close()
    SNRachieved=float(mean_signal)/(float(sigma))
    Achieved.Res_Beam=[Cube[0].header['BMAJ']*3600.,
                         Cube[0].header['BMIN']*3600.,
                         Cube[0].header['BPA']]

    if Current_Galaxy.Corruption == 'Casa_Sim':
        catalog_cube_name= 'Convolved_Cube_CS'
    elif Current_Galaxy.Corruption == 'Gaussian':
        catalog_cube_name='Convolved_Cube_Gauss'
    elif Current_Galaxy.Corruption == 'No_Corrupt':
        catalog_cube_name='Convolved_Cube_UC'
    else:
        catalog_cube_name='Convolved_Cube'

    if  catalog_cube_name != 'Convolved_Cube':
        os.system(
            f"mv {cfg.general.main_directory}{name}/Convolved_Cube.fits {cfg.general.main_directory}{name}/{catalog_cube_name}.fits")


    if not cfg.agc.retain_unconvolved_model:
        os.remove(f'{cfg.general.main_directory}{name}/Unconvolved_Cube.fits')
        #if we have

        # We'll create a little text file with an Overview of all the parameters
    Template['BMAJ'] = f'BMAJ= {Cube[0].header["BMAJ"]*3600.}'
    Template['BMIN'] = f'BMIN= {Cube[0].header["BMIN"]*3600.}'
    Template['BPA'] = f'BPA= {Cube[0].header["BPA"]}'
    Achieved.HI_Mass= mass
    Achieved.SNR= SNRachieved

    setattr(Achieved, "Mean_Signal", mean_signal )
    setattr(Achieved, "Noise", sigma )
    setattr(Achieved, "Pixel_Beam", pixperbeam )
    setattr(Achieved, "Channelwidth", Cube[0].header['CDELT3']/1000. )

    write_overview_file(
        f"{cfg.general.main_directory}{name}/{name}-Info.txt", Current_Galaxy,Template,Achieved)


    # We also want a file that contains initial estimates for all the parameters. We scramble them with gaussian variations
    cf.scrambled_initial(f"{cfg.general.main_directory}{name}/",Template)
    cf.plot_input(f"{cfg.general.main_directory}{name}/",Template, \
                Title=f'DM Mass in = {Current_Galaxy.Mass:.2e}',RHI=Rad_HI
                ,add_sbr=molecular_profile, WarpR=[WarpStart,WarpEnd],Distance=Distance,\
                font_file= cfg.general.font_file )

    Template.clear()


    # and cleanup
    os.system(f"rm -f {cfg.general.main_directory}{name}/ModelInput.def")
    os.system(f"rm -f {cfg.general.main_directory}{name}/Input.fits")
    os.system(
        f"mv {cfg.general.main_directory}{name}/tirific.def {cfg.general.main_directory}{name}/ModelInput.def")
    return f'{int(Current_Galaxy.Model_Number):d}|{Distance:.2f}|{name}|{catalog_cube_name}\n'

one_galaxy.__doc__= f'''
NAME:
   one_galaxy(cfg, Current_Galaxy)

PURPOSE:
    Routine to build a singlle galaxy on a single processor

CATEGORY:
   agc

INPUTS:
    cfg = Configuration file
    Current_Galaxy = Base Galaxy class to build the galaxy from

OPTIONAL INPUTS:

OUTPUTS:
    summary line to at to catalog

OPTIONAL OUTPUTS:

PROCEDURES CALLED:
   Unspecified

NOTE:
'''

def plot_RC(set_done,Mass,Rad,Vrot,colors,max_rad,sub_ring,ax,\
                font_file='empty.ttf'):
    '''dd the RC to the overview plot and return updated tracker'''
    if set_done[0] == 1024:
        set_done= [Mass]
        try:
            mpl_fm.fontManager.addfont(font_file)
            font_name = mpl_fm.FontProperties(fname=font_file).get_name()
        except FileNotFoundError:
            font_name = 'DejaVu Sans'
        labelfont= {'family': font_name,
                 'weight':'normal',
                 'size':22}
        plt.rc('font',**labelfont)
        c = next(colors)
        plt.figure(59, figsize=(8, 8), dpi=300, facecolor='w', edgecolor='k')
        ax = plt.subplot(1, 1, 1)
        plt.plot(Rad,Vrot,c=c)
        #plt.plot(Rad,Vrot,'ko',label='M$_{\odot}$ = {:.1e}'.format(Current_Galaxy.Mass))
        plt.plot(Rad,Vrot,'o',c=c,
                 label='M$_{\odot}$ = '+' {:.1e}'.format(Mass))
        plt.ylabel('V$_{rot}$ (km s$^{-1}$)',**labelfont)
        plt.xlabel('Radius (kpc)',**labelfont)
        ax.yaxis.set_minor_locator(AutoMinorLocator(4))
        ax.xaxis.set_minor_locator(AutoMinorLocator(4))
        for axis in ['top','bottom','left','right']:
            ax.spines[axis].set_linewidth(1.5)
        plt.tick_params(axis='both', which='minor',
                        bottom=True,left=True,length=3)
        plt.tick_params(axis='both', which='major', labelsize=17, length=6)
        plt.tick_params(axis='both', which='both', direction = 'in',
                        width=1.5 , bottom=True,left=True ,right =True, top=True)
        max_rad = np.max(Rad)+1

    else:
        set_done.append(Mass)
        plt.figure(59)
        c = next(colors)
        if np.max(Rad)+1 > max_rad:
            max_rad =  np.max(Rad)+1
        plt.plot(Rad,Vrot,c=c)
        #plt.plot(Rad,Vrot,'o',label=r'M$_{\odot}$ = {:.1e}'.format(Current_Galaxy.Mass),c=c)
        plt.plot(Rad,Vrot,'o',
                 label='M$_{\odot}$ = '+' {:.1e}'.format(Mass),c=c)

    return set_done,max_rad,colors,ax
plot_RC.__doc__ = f'''
NAME:
   plot_RC

PURPOSE:
    Add the RC to the overview plot and return updated tracker

CATEGORY:
   agc

INPUTS:
   set_done = tracker on which masses are plotted
   Mass = mass of the current galaxy
   Rad = radii of the RC
   Vrot = RC
   colors = tracker which colors are used
   max_rad = the maximum radius in the plot

OPTIONAL INPUTS:

OUTPUTS:
    set_done = updated tracker on which masses are plotted
    max_rad = updated the maximum radius in the plot
    colors = updated tracker which colors are used

OPTIONAL OUTPUTS:

PROCEDURES CALLED:
   Unspecified

NOTE:
'''

def set_name(Current_Galaxy):
    '''set the name for the galaxy'''
    #name=""
    #for key,value in Current_Galaxy.items:
    name=f"Mass{float(Current_Galaxy.Mass):.1e}-i{Current_Galaxy.Inclination:.1f}d{Current_Galaxy.Dispersion[0]:.1f}-{Current_Galaxy.Dispersion[1]:.1f}"
    name=f"{name}pa{Current_Galaxy.PA:.1f}w{Current_Galaxy.Warp[0]:.2f}-{Current_Galaxy.Warp[1]:.2f}-"
    name=f"{name}{Current_Galaxy.Flare}-ba{Current_Galaxy.Beams:.1f}SNR{Current_Galaxy.SNR:.1f}"
    name=f"{name}bm{Current_Galaxy.Res_Beam[0]:.1f}-{Current_Galaxy.Res_Beam[1]:.1f}ch{Current_Galaxy.Channelwidth:.1f}"
    name=f"{name}-{Current_Galaxy.Arms}-{Current_Galaxy.Bar}-rm{Current_Galaxy.Radial_Motions:.1f}"

    return name

set_name.__doc__ = f'''
NAME:
   set_name

PURPOSE:
    set the name for the galaxy

CATEGORY:
   agc

INPUTS:
   Current_Galaxy = the Current Gakaxy that is processed

OPTIONAL INPUTS:

OUTPUTS:
    name = string with the name
OPTIONAL OUTPUTS:

PROCEDURES CALLED:
   Unspecified

NOTE: I'm sure this can be done in a more pythonic fashion
'''
def sort_results(All_Galaxies,Result_Galaxies,results):
    for line in Result_Galaxies:
        split_line=line.split('|')
        index_no = np.where(split_line[2] == np.array(All_Galaxies))[0][0]
        results[index_no]=line
    return results
sort_results.__doc__ = f'''
NAME:
   sort_results(All_Galaxies, Result_Galaxies, results)

PURPOSE:
    sort the results into the final output

CATEGORY:
    agc

INPUTS:
    All_Galaxies = a list with the name of all galaxies
    Result_Galaxies = Part of the results to arrange
    results = final sorted array where to plant Result_Galaxies

OPTIONAL INPUTS:

OUTPUTS:
    An organized results list

OPTIONAL OUTPUTS:

PROCEDURES CALLED:
   Unspecified

NOTE:
'''

def write_overview_file(filename,Current_Galaxy,Template,Achieved):
    RAhr,DEChr= cf.convertRADEC(*Current_Galaxy.Coord)
    RAhra,DEChra= cf.convertRADEC(*Achieved.Coord)
    with open(filename, 'w') as overview:
            overview.write(f'''  # This file contains the basic parameters of this galaxy. For the radial dependencies look at Overview.png or ModelInput.def.
#{'Variable':<14s} {'Requested':<15s} {'Achieved':<15s} {'Units':<15s}
{'Inclination':<15s} {Current_Galaxy.Inclination:<15.2f} {float(Template['INCL'].split('=')[1].split()[0]):<15.2f} {'degree':<15s}
{'PA':<15s} {Current_Galaxy.PA:<15.2f} {float(Template['PA'].split('=')[1].split()[0]):<15.2f} {'degree':<15s}
{'Sys. Velocity':<15s} {'':<15s} {float(Template['VSYS'].split('=')[1].split()[0]):<15.2f} {'km/s':<15s}
{'RA':<15s} {RAhr.strip():<15s} {RAhra.strip():<15s} {'':<15s}
{'Declination':<15s} {DEChr.strip():<15s} {DEChra.strip():<15s} {'':<15s}
{'Dispersion':<15s} {f"{Current_Galaxy.Dispersion[0]:.2f}-{Current_Galaxy.Dispersion[1]:.2f}":<15s} {f"{float(Template['SDIS'].split('=')[1].split()[0]):.2f}-{float(Template['SDIS'].split('=')[1].split()[-1]):.2f}":<15s} {'km/s':<15s}
{'Scale height':<15s} {'':<15s} {f"{float(Template['Z0'].split('=')[1].split()[0]):.2f}-{float(Template['Z0'].split('=')[1].split()[-1]):.2f}":<15s} {'arcsec':<15s}
{'Warp':<15s} {f"{Current_Galaxy.Warp[0]:.2f}-{Current_Galaxy.Warp[1]:.2f}":<15s} {f"{Achieved.Warp[0]:.2f}-{Achieved.Warp[1]:.2f}":<15s} {'radian':<15s}
{'Warp Radius':<15s} {'':<15s} {Achieved.Warp_Radius:<15.3f} {'kpc':<15s}
{'HI Radius':<15s} {'':<15s} {Achieved.HI_Radius[0]:<15.3f} {'kpc':<15s}
{'Scalelength':<15s} {'':<15s} {Achieved.Scalelength:<15.3f} {'kpc':<15s}
{'Maj Axis':<15s} {Current_Galaxy.Beams:<15.3f} {Achieved.Beams:<15.3f} {'Beams'}
{'Total Mass':<15s} {Current_Galaxy.Mass:<15.1e} {Achieved.Mass:<15.1e} {'M_solar':<15s}
{'HI Mass':<15s} {Current_Galaxy.HI_Mass:<15.1e} {Achieved.HI_Mass:<15.1e} {'M_solar':<15s}
{'Chan. Width':<15s} {Current_Galaxy.Channelwidth:<15.3f} {Achieved.Channelwidth:<15.3f} {'km/s':<15.2s}
{'Chan. Dep.':<15s} {'':<15s} {Achieved.Channel_Dep:<15s}
{'SNR':<15s} {Current_Galaxy.SNR:<15.3f} {Achieved.SNR:<15.3f}
{'Mean Signal':<15s} {'':<15s} {Achieved.Mean_Signal*1000.:<15.3f} {'mJy/beam':<15s}
{'Noise':<15s} {'':<15s} {Achieved.Noise*1000.:<15.3f} {'mJy/beam':<15s}
{'Distance':<15s} {'':<15s} {float(Template['DISTANCE'].split('=')[1]):<15.3f} {'Mpc':<15s}
{'Maj. FWHM':<15s} {Current_Galaxy.Res_Beam[0]:<15.2f} {Achieved.Res_Beam[0]:<15.2f} {'arcsec':<15s}
{'Min. FWHM':<15s} {Current_Galaxy.Res_Beam[1]:<15.2f} {Achieved.Res_Beam[1]:<15.2f} {'arcsec':<15s}
{'Beam BPA':<15s} {Current_Galaxy.Res_Beam[2]:<15.2f} {Achieved.Res_Beam[2]:<15.2f} {'degree':<15s}
{'Pixel per Beam':<15s} {'':<15s} {Achieved.Pixel_Beam:<15.2f} {'pixel':<15s}
{'Flare':<15s} {Current_Galaxy.Flare:<15s} {Achieved.Flare:<15s}
{'Arms':<15s} {Current_Galaxy.Arms:<15s} {Achieved.Arms:<15s}
{'Bar':<15s} {Current_Galaxy.Bar:<15s} {Achieved.Bar:<15s}
{'Corruption':<15s} {Current_Galaxy.Corruption:<15s} {Achieved.Corruption:<15s}''')

write_overview_file.__doc__ = f'''
NAME:
   write_overview_file

PURPOSE:
    write the overview file with requested and achieved values

CATEGORY:
   agc

INPUTS:
   Current_Galaxy = the Current Gakaxy that is processed(The requested values)
   Template = The tirific template that is used to create the galaxy
   Achieved = Copy of Current_Galaxy that contains the values measured from the cube

OPTIONAL INPUTS:

OUTPUTS:
    file with a table with the requested and achieved values
OPTIONAL OUTPUTS:

PROCEDURES CALLED:
   Unspecified

NOTE: For now the Achieved.Flare and Achieved.Warp are merely copies of the input simply assumed to be implemented correctly
'''
