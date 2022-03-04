#!/usr/local/bin/ python3
# -*- coding: utf-8 -*-
#This program is a pyhton script to create a data base of artificial galaxies.
#It first creates a model with tirific at a high resolution and then runs it through casa to get obtain a realistic observation.
# once a numerical list is set in length we can convert it to a numpy array in order to do operations faster.
# first we import numpy
import copy
import numpy as np
import os
import psutil
import pyHIARD.common_functions as cf
import resource
import scipy.ndimage as ndimage
import subprocess
import sys
import time
import warnings

with warnings.catch_warnings():
     warnings.simplefilter("ignore")
     import matplotlib
     matplotlib.use('pdf')
     import matplotlib.pyplot as plt
     from matplotlib.ticker import AutoMinorLocator
try:
    import importlib.resources as import_res
except ImportError:
    # Try backported to PY<37 `importlib_resources`.
    import importlib_resources as import_res


from astropy.io import fits
from astropy.io.fits.verify import VerifyWarning
from multiprocessing import Pool,get_context
from pyHIARD.AGC.base_galaxies import Base_Galaxy
from pyHIARD.constants import G_agc,H_0
from scipy import interpolate
from scipy import integrate



#Some errors
class TirificRunError(Exception):
    pass

#------------------------------!!!!!!!!!!!!!!!!!!!!!!!!!!!!!Main for creating the AGC!!!!!!!!!!!!!!!!!!!!!----------------------
def AGC(cfg):
    '''Create realistic artificial galaxies.'''
    # Let's give an overview of the database that will be created
    print(f"We will create a database with {len(cfg.agc.base_galaxies)} basic sets in the directory {cfg.general.main_directory}.\n")
    if cfg.agc.delete_existing:
        print("All previous models will be removed prior to the build. \n")
    else:
        print("We will retain previously build models. \n")
    if cfg.agc.inhomogenous:
        print(f"We will use {cfg.agc.corruption_method} as corruption method and inhomogeneities will be added.\n")
    else:
        print(f"We will use {cfg.agc.corruption_method} as corruption method.\n")
    print("We will vary the following parameters")
    if 'Inclination' in cfg.agc.variables_to_vary:
        print(f"We vary the inclination with the following values: {','.join([str(e) for e in cfg.agc.inclination])}.\n")
    if 'PA' in cfg.agc.variables_to_vary:
        print(f"We vary the PA with the following values: {','.join([str(e) for e in cfg.agc.pa])}.\n")
    if 'Beams' in cfg.agc.variables_to_vary:
        print(f"We create the model with {','.join([str(e) for e in cfg.agc.beams])} beams across the major axis.\n")
    if 'Radial_Motions' in cfg.agc.variables_to_vary:
        print(f"Inject radial motions with speeds of {','.join([str(e) for e in cfg.agc.radial_motions])} km/s.\n")
    if 'Flare' in cfg.agc.variables_to_vary:
        print(f"We will swap the flaring of the base galaxies")
    if 'Arms' in cfg.agc.variables_to_vary:
        print(f"We will swap the inclusion of the arms in the base galaxies")
    if 'Bar' in cfg.agc.variables_to_vary:
        print(f"We will swap the inclusion of the bar in the base galaxies")
    if 'Mass' in cfg.agc.variables_to_vary:
        print(f"We add the following masses to each base set: {','.join([str(e) for e in cfg.agc.masses])}.\n")
    if 'Channelwidth' in cfg.agc.variables_to_vary:
        print(f"Varying the channel width with: {','.join([str(e) for e in cfg.agc.channelwidth])} km/s.\n")
    if 'SNR' in cfg.agc.variables_to_vary:
        print(f"Varying the signal to noise ratio with: {','.join([str(e) for e in cfg.agc.snr])}.\n")
    if 'Warp' in cfg.agc.variables_to_vary:
        print(f"Varying the theta angle of the angular momentum vector with: {','.join([str(e) for e in cfg.agc.warp[:][0]])}.\n")
        print(f"Varying the phi angle of the angular momentum vector with: {','.join([str(e) for e in cfg.agc.warp[:][1]])}.\n")
    if 'Beam_Resolution' in cfg.agc.variables_to_vary:
        print(f"Varying the beam size with: {','.join([str(e) for e in cfg.agc.beam_size])}.\n")



    # Let's just make 1 catalogue with everything we need and adjust
    # the fitting program to have this one catalogue as input

    Catalogue=f"{cfg.general.main_directory}Output_AGC_Summary.txt"
    # If we are making new models we want to ensure this is a new file
    number_models= cf.get_created_models(Catalogue,cfg.agc.delete_existing)

    #Copy a fits file from the WHISP data base to use as template if it is not ther3 yet
    # Check for the existence of a template fits file
    templatethere= os.path.isfile(cfg.general.main_directory+'/Input.fits')
    #templatethere =False
    # If it doesn't exist copy it into the directory
    if not templatethere:
        my_resources = import_res.files('pyHIARD.Templates')
        data = (my_resources / 'Input.fits').read_bytes()
        with open(f"{cfg.general.main_directory}/Input.fits",'w+b') as tmp:
            tmp.write(data)
        # let's make sure it has BMAJ, BMIN and BPA
        with warnings.catch_warnings():
            warnings.simplefilter("ignore",category=VerifyWarning)
            dummy = fits.open(cfg.general.main_directory+'/Input.fits',uint = False, do_not_scale_image_data=True,ignore_blank = True)
        sizex = 500
        sizey = 500
        sizez= 120
        header_sets = {'BMAJ': 0., 'BMIN': 0., 'BPA': 0.,\
                       'CRPIX1': 200., 'CRPIX2': 200.,'CRPIX3': 60.,\
                       'CDELT1': -4./3600.,'CDELT2': 4./3600.,'CDELT3': 4.,\
                       'CUNIT2': 'M/S', 'CRVAL3': 500.,\
                       'CTYPE1': 'RA---SIN',  'CTYPE2': 'DEC--SIN', 'CTYPE3': 'VELO-HEL',\
                       'BITPIX': -32,'NAXIS': 3,\
                       'NAXIS1':sizex,'NAXIS2':sizey,'NAXIS3':sizez}
        for key,value in header_sets.items():
            dummy[0].header[key] = value
        header_deletes = ['CROTA1','CROTA2','DRVAL3','DTYPE3','DUNIT3','HISTORY',\
                            'BMMAJ','BMMIN','BMPA','MAPLAB','BLANK']
        for key in header_deletes:
            del dummy[0].header[key]
        # make the cube a typical size
        dummy2 = np.zeros((sizez,sizey,sizex))
        dummy2[60,250,250] = 5
        fits.writeto(cfg.general.main_directory+'/Input.fits',dummy2,dummy[0].header, output_verify='silentfix+ignore', overwrite = True)
        dummy.close()


    # All this went well
    templatethere= os.path.isfile(cfg.general.main_directory+'/Input.fits')
    if templatethere:
        print(" The Input Template is found and all is well")
    else:
        print(" The Input Template is NOT found !!! ABORTING")
        sys.exit()

    # start a loop over the various base galaxies
    set_done= [1024]
    plot_ax =[]


    # If we make new models delete everything in the directory
    if cfg.agc.delete_existing:
        masses_to_delete= []
        if 'Mass' in cfg.agc.variables_to_vary:
            #Because python is the dumbest language ever
            masses_to_delete= copy.deepcopy(cfg.agc.masses)
        for galaxy in cfg.agc.base_galaxies:
            if galaxy == 7:
                continue
            else:
                masses_to_delete.append(Base_Galaxy(galaxy).Mass)
        to_delete= f"rm -R {' '.join([f'{cfg.general.main_directory}Mass{x:.1e}-*rm*' for x in masses_to_delete])} {cfg.general.main_directory}Fractions_and_Masses.txt"
        print("All previous models of the requested base galaxies will be removed prior to the build. \n")
        print(f"The command {to_delete} will be run.")
        cfg.agc.delete_existing = cf.get_bool("Are you sure you want to do this? (Yes/No, default=No): ",default=False)
        if cfg.agc.delete_existing:
            os.system(to_delete)





    colors=iter(plt.cm.rainbow(np.linspace(0,1,len(cfg.agc.base_galaxies)+len(cfg.agc.masses))))
    max_rad = 0.
    All_Galaxies = []
    Gauss_Galaxies = []
    Casa_Galaxies = []
    created = []

    for base in range(len(cfg.agc.base_galaxies)):
        base_defined=False
        # We want to keep the center constant per base galaxy, for easy comparison as well as to be able to investigate how center determination is affected
        if cfg.agc.corruption_method == 'Gaussian':
            RAdeg=np.random.uniform()*360
            DECdeg=(np.arccos(2*np.random.uniform()-1)*(360./(2.*np.pi)))-90
        else:
            RAdeg = np.random.uniform()*360
            DECdeg=-60
            while DECdeg < -20.:
                DECdeg = (np.arccos(2*np.random.uniform()-1)*(360./(2.*np.pi)))-90
        # From here we go into a loop to adjust variables over the bases
        for ix in range(len(cfg.agc.variables_to_vary)):
            if cfg.agc.variables_to_vary[ix] == 'Inclination':numloops=len(cfg.agc.inclination)
            elif cfg.agc.variables_to_vary[ix] == 'PA': numloops=len(cfg.agc.pa)
            elif cfg.agc.variables_to_vary[ix] in ['Flare','Arms','Bar','Base']: numloops=1
            elif cfg.agc.variables_to_vary[ix] == 'Warp': numloops=len(cfg.agc.warp)
            elif cfg.agc.variables_to_vary[ix] == 'Beams': numloops=len(cfg.agc.beams)
            elif cfg.agc.variables_to_vary[ix] == 'SNR': numloops=len(cfg.agc.snr)
            elif cfg.agc.variables_to_vary[ix] == 'Channelwidth': numloops=len(cfg.agc.channelwidth)
            elif cfg.agc.variables_to_vary[ix] == 'Beam_Resolution': numloops=len(cfg.agc.beam_size)
            elif cfg.agc.variables_to_vary[ix] == 'Radial_Motions': numloops=len(cfg.agc.radial_motions)
            elif cfg.agc.variables_to_vary[ix] == 'Dispersion': numloops=len(cfg.agc.dispersion)
            elif cfg.agc.variables_to_vary[ix] == 'Mass': numloops=len(cfg.agc.masses)
            else:
                print("This is not a supported parameter")
                exit()
            for jx in range (numloops):
                if cfg.agc.base_galaxies[base] > 6:
                    if not base_defined:
                        Current_Galaxy = Base_Galaxy(cfg.agc.base_galaxies[base])
                        Current_Galaxy_Base = copy.deepcopy(Current_Galaxy)
                        base_defined= True
                        to_delete= f"rm -R {' '.join([f'{cfg.general.main_directory}Mass{Current_Galaxy.Mass:.1e}-*rm*'])}"
                        print("All previous models of the requested base galaxy will be removed prior to the build. \n")
                        print(f"The command {to_delete} will be run.")
                        cfg.agc.delete_existing = cf.get_bool("Are you sure you want to do this? (Yes/No, default=No): ",default=False)
                        if cfg.agc.delete_existing:
                            os.system(to_delete)
                    else:
                        Current_Galaxy = copy.deepcopy(Current_Galaxy_Base)
                else:
                    Current_Galaxy = Base_Galaxy(cfg.agc.base_galaxies[base])
                if cfg.agc.variables_to_vary[ix] == 'Inclination':Current_Galaxy.Inclination = cfg.agc.inclination[jx]
                elif cfg.agc.variables_to_vary[ix] == 'PA': Current_Galaxy.PA = cfg.agc.pa[jx]
                elif cfg.agc.variables_to_vary[ix] == 'Flare':
                    if Current_Galaxy.Flare == 'Flare':
                        Current_Galaxy.Flare = "No_Flare"
                    else:
                        Current_Galaxy.Flare = "Flare"
                elif cfg.agc.variables_to_vary[ix] == 'Warp': Current_Galaxy.Warp = [cfg.agc.warp[jx][0],cfg.agc.warp[jx][1]]
                elif cfg.agc.variables_to_vary[ix] == 'Beams':Current_Galaxy.Beams = cfg.agc.beams[jx]
                elif cfg.agc.variables_to_vary[ix] == 'SNR': Current_Galaxy.SNR = cfg.agc.snr[jx]
                elif cfg.agc.variables_to_vary[ix] == 'Channelwidth': Current_Galaxy.Channelwidth = cfg.agc.channelwidth[jx]
                elif cfg.agc.variables_to_vary[ix] == 'Beam_Resolution': Current_Galaxy.Res_Beam = [cfg.agc.beam_size[jx][0],cfg.agc.beam_size[jx][1],cfg.agc.beam_size[jx][2]]
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
                elif cfg.agc.variables_to_vary[ix] == 'Radial_Motions': Current_Galaxy.Radial_Motions = cfg.agc.radial_motions[jx]
                elif cfg.agc.variables_to_vary[ix] == 'Dispersion': Current_Galaxy.Dispersion = cfg.agc.dispersion[jx]
                elif cfg.agc.variables_to_vary[ix] == 'Mass': Current_Galaxy.Mass = cfg.agc.masses[jx]
                elif cfg.agc.variables_to_vary[ix] == 'Base': pass
                else:print("This is not a supported parameter")
                Current_Galaxy.Res_Beam[0:1] = np.sort(Current_Galaxy.Res_Beam[0:1])
                setattr(Current_Galaxy,"Coord",[RAdeg,DECdeg])
                setattr(Current_Galaxy,"Model_Number",number_models)
                number_models += 1
                if cfg.agc.corrupt_models:
                    if (cfg.agc.corruption_method == 'Casa_5' and (int(number_models/5.) == number_models/5.)) or (cfg.agc.corruption_method == 'Casa_Sim'):
                        setattr(Current_Galaxy,"Corruption","Casa_Sim")
                    elif (cfg.agc.corruption_method == 'Gaussian' or cfg.agc.corruption_method == 'Casa_5'):
                        setattr(Current_Galaxy,"Corruption","Gaussian")
                    else:
                      print("!!!!!!!This corruption method is unknown, leaving the cube uncorrupted and unconvolved!!!!!!!!")
                else:
                    setattr(Current_Galaxy,"Corruption","Uncorrupted")
                #Check if galaxy already exists
                name = set_name(Current_Galaxy)
                print(f"{name} is the name we attach to the current galaxy")

                # Check for the existence of the directory
                constructstring=f"mkdir {cfg.general.main_directory}{name}"
                checkdir=False
                galaxy_dir =f"{cfg.general.main_directory}{name}/"
                galaxy_exists= os.path.isdir(galaxy_dir)
                if not galaxy_exists:
                    os.system(constructstring)
                    created.append(name)
                else:
                    time.sleep(0.1)
                    # Do we have a cube
                    galaxy_cube_exist = os.path.isfile(f"{galaxy_dir}Convolved_Cube.fits")
                    if not galaxy_cube_exist:
                        galaxy_cube_exist = os.path.isfile(f"{galaxy_dir}Convolved_Cube_CS.fits")
                    if galaxy_cube_exist:
                        print("This galaxy appears fully produced")
                        checkdir = True
                        continue
                    else:
                        if name in created:
                            continue
                        else:
                            print("The directory was made but there is no full cube avalaible")
                            #print("Reproducing the galaxy. Be aware of Double Table entries")
                            print("This is too dangerous. Breaking the code.")
                            exit()
                if Current_Galaxy.Corruption in ['Uncorrupted','Gaussian']:
                    Gauss_Galaxies.append((cfg,Current_Galaxy))
                else:
                    Casa_Galaxies.append((cfg,Current_Galaxy))
                All_Galaxies.append(name)
                if Current_Galaxy.Mass not in set_done:
                    SBRprof,Rad,sclength,MHI,Rad_HI,Vrot,sub_ring,molecular_profile = build_sbr_prof(Current_Galaxy,symmetric=cfg.agc.symmetric,no_template=True) #Column densities,Raii in kpc, Opt_scalelength in kpc, HI mass in M_solar
                    set_done,max_rad,colors,plot_ax = plot_RC(set_done,Current_Galaxy.Mass,Rad,Vrot,colors,max_rad,sub_ring,plot_ax)

                #print(f"This is the parameter to vary {cfg.agc.variables_to_vary[ix]}.")

    if len(Casa_Galaxies) > 0:
        if cfg.general.multiprocessing:
            available_memory = psutil.virtual_memory().total/2**30
            no_process = int(np.floor(available_memory/8.))
            if no_process > len(Casa_Galaxies):
                no_process =  len(Casa_Galaxies)
            #tclean is parallel inmplemented but simobserve is not, need betterhandling in casa for multprocessing
            #The problem is memory limits combined with CPU limits for simobserve. No easy solution
            # Casa recomend 8Gb per CPU so we limit the no processes per 8Gb memory for now.
            with get_context("spawn").Pool(processes=no_process) as pool:
                results_casa = pool.starmap(one_galaxy, Casa_Galaxies)
        else:
            results_casa = []
            for the_galaxy in Casa_Galaxies:
                single_result = one_galaxy(*the_galaxy)
                results_casa.append(single_result)
    #Create All Uncoorupted and Gaussian Galaxies
    if len(Gauss_Galaxies) > 0:
        if cfg.general.multiprocessing:
            no_process = cfg.general.ncpu
            if no_process > len(Gauss_Galaxies):
                no_process =  len(Gauss_Galaxies)
            with get_context("spawn").Pool(processes=no_process) as pool:
                results_gauss = pool.starmap(one_galaxy, Gauss_Galaxies)
        else:
            results_gauss = []
            for the_galaxy in Gauss_Galaxies:
                single_result = one_galaxy(*the_galaxy)
                results_gauss.append(single_result)

    results = ['empty']*len(All_Galaxies)
    if len(Gauss_Galaxies) > 0:
        results= sort_results(All_Galaxies,results_gauss,results)
    if len(Casa_Galaxies) > 0:
        results= sort_results(All_Galaxies,results_casa,results)


    os.system(f'rm -f {cfg.general.main_directory}/Input.fits')
    print("We created {} models".format(number_models))

    #Add the results to the catalogue
    with open(Catalogue, 'a') as cat:
        for line in results:
            cat.write(line)

    plt.figure(59)

    plot_ax.set_ylim(ymin=0)
    plot_ax.set_xlim(xmin=0, xmax=max_rad)
    plt.legend(loc='lower right',fontsize=12)
    plt.savefig('Rotation_Curves.pdf', bbox_inches='tight')
    plt.close()
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

def calc_vc_NFW(DM_mass,MHI,m_star, rad):
    r200= (DM_mass*G_agc/(100.* H_0**2)*1e12)**(1./3.)
    v200square=DM_mass*G_agc/r200
    xr=rad/(r200/10**3)
    c200=10**(0.905-0.101*np.log10(DM_mass/(10**12*100./H_0))) #From A. Dutton 2014
    NFWvflat=np.sqrt(v200square*(1./xr)*((np.log(1.+c200*xr)-(c200*xr)/(1+c200*xr))/(np.log(1+c200)-c200/(1+c200))))
    v_star=np.sqrt(m_star*G_agc/(rad*10**3))
    v_HI=np.sqrt(MHI*1.4*G_agc/(rad*10**3))
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
def copy_disk(olddisk,newdisk):
    '''Routine to copy a disk in the tirific template file'''
    start = 0
    startlast = 0.
    if int(Template["NDISKS"].split('=')[1]) < newdisk:
        Template["NDISKS"] = "NDISKS = {:d}".format(newdisk)
    copkeys ="Empty"
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
            elif (key_ext[1] != str(olddisk))  and start == 2:
                del copkeys[-1]
                start += 1
            if key_ext[1] == str(newdisk-1) and startlast == 0:
                startlast = 1
            if key_ext[1] == str(newdisk-1) and startlast == 1:
                last = key
            if key_ext[1] != str(newdisk-1) and startlast == 1:
                startlast +=1
        if key == 'CONDISP' and (start == 2 or startlast ==1):
            if start == 2:
                del copkeys[-1]
            startlast += 1
            start += 1
    for key in reversed(copkeys):
        var_name = key.split('_')[0]
        Template.insert(last,var_name+"_{:d}".format(newdisk),var_name+"_{:d} =".format(newdisk)+Template[key].split('=')[1])
    Template.insert("CFLUX_{:d}".format(newdisk-1),"CFLUX_{:d}".format(newdisk),"CFLUX_{:d} =".format(newdisk)+Template["CFLUX_{:d}".format(newdisk-1)].split('=')[1])
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
def derivative(x,func):
    '''Obtaining the derivative at any give point of a function'''
    h=x/1000.
    der=(func(x+h)-func(x-h))/(2*h)
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
def build_sbr_prof(Current_Galaxy,symmetric = False,no_template=False):
    '''function to get the surface brightness profile based on the last element in the rotation curve'''
    # First we use the mass to build a rotation curve and sbr profile
    # Avanti used 	V_Ropt = (200*(l**0.41))/(0.80 +0.49*np.log10(l) +(0.75*np.exp(-0.4*l))/(0.47 + 2.25*(l**0.4)))**0.5 (Persic & Sallucci)
    # We will not use the URC as it get's silly at high mass and remains showing flat parts at low mass
    # We will use the parameterisation of Courteau 1997
    v_circ, HIrad, MHI,MK = get_HI_disk(Current_Galaxy.Mass)

        #First we set the radii at with a 10 elements per rings out to 1.5 times HIrad. However we always need at least 25 rings for transition purposes
    if Current_Galaxy.Beams < 5.:
        sub_ring = int(25./Current_Galaxy.Beams)
    else:
        sub_ring = 5
    Rad= np.linspace(0.,1.5*HIrad, int(sub_ring*Current_Galaxy.Beams))
    # Finaaly the courteau presciption
    v_flat,R_0,beta,gamma = get_sparcs_fit(v_circ,HIrad)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        x = R_0/Rad
    x[0]=1e-8

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
    h_r = (-4.13422991 - 0.31576291 * MK) * 1000. #parsec

    if not symmetric:
        if fracHIrad < 0.5:
            Rhi_p = np.array([(HIrad - HIrad / 15.) * 10 ** 3, (HIrad + HIrad / 15.) * 10 ** 3])
            # Rhi_p2 = HIrad + HIrad / 20. * 10 ** 3
            h_r = np.array([h_r + h_r / 5, h_r - h_r / 5])
            # h_r2 = h_r - h_r / 20.
        else:
            Rhi_p = np.array([(HIrad + HIrad / 15.) * 10 ** 3, (HIrad - HIrad / 15.) * 10 ** 3])
            # Rhi_p2 = HIrad - HIrad / 20. * 10 ** 3
            h_r = np.array([h_r - h_r / 5., h_r + h_r / 5.])
            # h_r2 = h_r + h_r / 20.
    else:
        Rhi_p = np.array([HIrad,HIrad]*10**3,dtype=float)
        h_r= np.array([h_r,h_r],dtype=float)
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
    for i in [0,1]:
        Exp[:, i] = I_cen * np.exp(-a / h_r[i])
        # gaussian2 From Martinsson 2015
        Sig2[:, i] = np.exp(-(a - 0.4 * Rhi_p[i]) ** 2 / (2 * (s1[i]) ** 2))
        # total
        Sigma[:, i] = Sig2[:, i] - Exp[:, i]
        Sigma[Sigma[:,i] < 0.,i] = 0.  # for negative sigma max, does not include negative values
        # scale Sigma such that it is one at HI rad
        new[i] = 1. / Sigma[Hiradindex[i], i]
        Sigma[:, i] = new[i] * Sigma[:, i]
    Sigma[0:4,0]=(Sigma[0:4,0]+Sigma[0:4,1])/2.
    Sigma[0:4, 1] = Sigma[0:4, 0]
    # get the HI Mass in the profile
    OutHIMass = integrate.simps((np.pi * a) * Sigma[:, 0], a) + integrate.simps((np.pi * a) * Sigma[:, 1], a)
    # And check that it matches the required HI mas within 5%
    counter = 1
    while np.absolute(MHI - OutHIMass) > MHI *0.02:
        # if not rescale sigma1
        if MHI - OutHIMass > 0:
            s1 = (0.36 - (0.0025 * counter)) * Rhi_p
        else:
            s1 = (0.36 + (0.0025 * counter)) * Rhi_p
        # and recalculate the profile
        for i in [0,1]:
            Sig2[:, i] = np.exp(-(a - 0.4 * Rhi_p[i]) ** 2 / (2 * (s1[i]) ** 2))
            Sigma[:,i] = Sig2[:, i] - Exp[:, i]
            Sigma[Sigma[:,i] < 0.,i] = 0.
            new[i] = 1. / Sigma[Hiradindex[i], i]
            Sigma[:, i] = new[i] * Sigma[:, i]
        Sigma[0:4,0]=(Sigma[0:4,0]+Sigma[0:4,1])/2.
        Sigma[0:4, 1]=Sigma[0:4,0]
        OutHIMass = integrate.simps((np.pi * a) * Sigma[:, 0], a) + integrate.simps((np.pi * a) * Sigma[:, 1], a)
        counter += 1
    #print("The requested HI mass = {:.2e} and the retrieved HI mass = {:.2e}".format(MHI, OutHIMass))
    # final HI radial distribution by renormalisation
    S = Sigma * (1.24756e+20)
    # Where A. Gogate contribution stops
    # S is column densities but tirific takes Jy * km/s/arcsec^2 so
    conv_column_arsec = 605.7383 * 1.823E18 * (2. * np.pi / (np.log(256.)))  # S(mJy/beam)*conv_column_arcsec=N_HI
    sbr_prof = S / (conv_column_arsec * 1000.)
    # Let's write these to the Template immediately
    # The surface brightness profile, which is still symmetric
    if not no_template:
        Template["SBR"] = "SBR = " + " ".join(str(e) for e in sbr_prof[:, 0])
        Template["SBR_2"] = "SBR_2 = " + " ".join(str(e) for e in sbr_prof[:, 1])


    # as the rest on of the code is based on a single SBR profile let's average the value
    HIrad = [np.mean(Rhi_p) / 10 ** 3,Rhi_p[0] / 10 ** 3,Rhi_p[1] / 10 ** 3] #HI radii in kpc mean,appr, rec
    molecular_profile = []
    for i in [0,1]:
        molecular_profile.append(Exp[:,i]*new[i]*1.24756e+20/(conv_column_arsec*1000.)) # Jy/arcsec^2 km/s
    molecular_profile=np.array(molecular_profile)
    h_r = np.mean(h_r)
    average_sbr_profile = np.zeros(len(sbr_prof[:, 0]))
    for i in range(len(sbr_prof)):
        average_sbr_profile[i] = np.mean(sbr_prof[i, :])
    return average_sbr_profile,Rad,h_r/1000.,OutHIMass, HIrad,vrot,sub_ring,molecular_profile
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


def create_arms(velocity,Radii,disk_brightness, disk=1,WarpStart=-1,Bar="No_Bar"):
    '''Create the spiral Arms'''
    if WarpStart == -1: WarpStart = Radii[-1]
    max_vrot=np.max(velocity)
    max_rad=Radii[-1]
    #The pattern speed at a given radius is vrot/radii
    V_Rot = interpolate.interpolate.interp1d(Radii,velocity,fill_value="extrapolate")
    # The radius of co-ration can be approximated by the extend of the visible disk ~ Warpstart (Roberts et al. 1975)
    Omega_CR=V_Rot(WarpStart)/WarpStart
    # From this we can estimate the inner and outer Lindblad resonances (Eq 38 Dobbs & Baba 2014)
     #The epicyclic frequency k^2=R*dOmega^2/dR+4*Omega^2
    # f'(x) = (f(x+h)-f(x-h))/2h
    h=WarpStart/1000.
    derive=(V_Rot(float(WarpStart+h))**2/(WarpStart+h)**2-V_Rot(float(WarpStart-h))**2/(WarpStart-h)**2)/(2*h)

    k_CR = (WarpStart * derive+4*Omega_CR**2)**0.5
    # So the ILR =
    if Bar == "Barred":
        num_arms=2
    else:
        #other wise 2 arms when a bulge is present and 4 arms when not
        if Radii[np.where(max_vrot == velocity)[0]] < np.mean(Radii):
            num_arms=2
        else:
            num_arms=4


    LLR = Omega_CR-k_CR/num_arms
    ULR = Omega_CR+k_CR/num_arms
    Radii[0]=0.1
    om= interpolate.interpolate.interp1d(Radii,velocity/Radii,fill_value="extrapolate")
    Radii[0]=0.
    r_cur=Radii[1]
    while om(r_cur) > ULR and r_cur < max_rad:
        r_cur += 0.1
    ILR = 0.75*r_cur
    r_cur= Radii[1]
    while om(r_cur) > LLR and r_cur < max_rad:
        r_cur += 0.1
    OLR = 0.75*r_cur



    # From Seigar et al. 2006 we get the relation between shear (S) and pitch angle
    S = 0.5*(1-Radii[1:]/velocity[1:]*derivative(Radii[1:],V_Rot))
    pitch2= 64.25-73.24*S
    #As we assume a constant pitch angle we will take the average between ILR and OLR as the pitch angle
    tmp = np.where((Radii[1:] > ILR) & (Radii[1:] < OLR))[0]
    pitch = np.sum(pitch2[tmp])/len(tmp)

    #print("This is the average pitch angle {}".format(pitch))
    #Using Kennicut's prescription.This prescription incorrectly states cot(P) instead of tan(P) see Davis et. al 2012
    # The arms start at the inner Lindblad Resonance and hence the phase is 0 there
    phase=np.log(Radii[1:]/ILR)/np.tan(pitch*np.pi/180.)*180/np.pi
    phase=np.hstack((phase[0],phase))

    #How many arms do we make
    # with bar it is always a grand design

    # we take a brightness in the arms 1/no_arms the brightness of the disk
    brightness=1./np.sqrt(num_arms)*disk_brightness
    # and only between the resonances
    index= np.where((Radii < ILR) | (Radii > OLR))[0]
    brightness[index]=0.
    # with a ten ring transition
    brightness[tmp[0]:tmp[0]+10]=brightness[tmp[0]:tmp[0]+10]*(1-1/np.arange(1,11))
    brightness[tmp[-1]-10:tmp[-1]]=brightness[tmp[-1]-10:tmp[-1]]*(1/np.arange(1,11))
    # For the width we take 15% of the full circle but scaled to the total galaxy size a good description would be swell.
    width=0.15*2.*np.pi*Radii
    if WarpStart < 10.:
        width= width*10./WarpStart
    ndisk=int(Template["NDISKS"].split('=',1)[1])
    ndisk+=1
    last_add = copy_disk(disk,ndisk)
    Template[f"SBR_{ndisk:d}"] = f"SBR_{ndisk:d} = 0."
    # To simulate strems towrds the arms we spin this up by 10 km/s
    Template[f"VROT_{ndisk:d}"] = f"VROT_{ndisk:d} = {' '.join(str(e+20.) for e in velocity)}"
    phaseshift=360./num_arms
    #we offset the phase by 37 degrees
    for i in range(num_arms):
        Template.insert(last_add,"GA{:d}A_{:d}".format(i+1,ndisk),"GA{:d}A_{:d} =".format(i+1,ndisk)+" ".join(str(e) for e in brightness))
        Template.insert("GA{:d}A_{:d}".format(i+1,ndisk),"GA{:d}P_{:d}".format(i+1,ndisk),"GA{:d}P_{:d} =".format(i+1,ndisk)+" ".join(str(e+i*phaseshift+37) for e in phase))
        Template.insert("GA{:d}P_{:d}".format(i+1,ndisk),"GA{:d}D_{:d}".format(i+1,ndisk),"GA{:d}D_{:d} =".format(i+1,ndisk)+" ".join(str(e) for e in width))
    # and we add radial motions of 40 km/s
    Template.insert("VROT_{:d}".format(ndisk),"VRAD_{:d}".format(ndisk),"VRAD_{:d} = -40.".format(ndisk))

    return phase,brightness,width
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
def create_bar(velocity,Radii,disk_brightness,Template, disk=1,WarpStart=-1):
    '''Routine to create a Bar'''
    if WarpStart == -1: WarpStart = Radii[-1]
    max_rad = Radii[-1]
    bar_width = 1+0.05*WarpStart # kpc + 5% of the total optical disk
    max_vrot=np.max(velocity)
    #The pattern speed at a given radius is vrot/radii
    V_Rot = interpolate.interpolate.interp1d(Radii,velocity,fill_value="extrapolate")
    # The radius of co-ration can be approximated by the extend of the visible disk ~ Warpstart (Roberts et al. 1975)
    Omega_CR=V_Rot(WarpStart)/WarpStart
    # From this we can estimate the inner and outer Lindblad resonances (Eq 38 Dobbs & Baba 2012)
     #The epicyclic frequency k^2=R*d/drOmega^2+4*Omega^2
    # f'(x) = (f(x+h)-f(x-h))/2h
    h=WarpStart/1000.
    derive=(V_Rot(float(WarpStart+h))**2/(WarpStart+h)**2-V_Rot(float(WarpStart-h))**2/(WarpStart-h)**2)/(2*h)
    k_CR = (WarpStart * derive+4*Omega_CR**2)**0.5
    # So the ILR =
    LLR = Omega_CR-k_CR/2.
    ULR = Omega_CR+k_CR/2.
    Radii[0]=0.1
    om= interpolate.interpolate.interp1d(Radii,velocity/Radii,fill_value="extrapolate")
    Radii[0]=0.
    r_cur=Radii[1]
    while om(r_cur) > ULR and r_cur < max_rad:
        r_cur += 0.1
    ILR = 0.75*r_cur
    r_cur= Radii[1]
    while om(r_cur) > LLR and r_cur < max_rad:
        r_cur += 0.1
    # We set the full brightness to the maximum of the disk
    brightness= np.zeros(len(disk_brightness))
    brightness[np.where(Radii < ILR)[0]]=np.max(disk_brightness)
    # The width has to be 180 when R < width and 180*width/(pi*r)
    width =  np.zeros(len(disk_brightness))
    width[Radii <= bar_width] =180.
    width[Radii > bar_width]=360./np.pi*np.arcsin(bar_width/Radii[Radii > bar_width]) #the angle made up of the radius and width *2.
    # Get the number of disks present
    ndisk=int(Template["NDISKS"].split('=',1)[1])
    ndisk+=1
    #print("We are adding disk no {:d}".format(ndisk))
    # we also need streaming motions
    vrad_bar=np.zeros(len(disk_brightness))
    vrad_bar[:]=-50.
    #copy the input disk
    last_add = copy_disk(disk,ndisk)
    #We offset by 37 deg.
    Template.insert("VSYS_{:d}".format(ndisk), "AZ1P_{:d}".format(ndisk), "AZ1P_{:d} = 37.".format(ndisk))
    Template.insert("AZ1P_{:d}".format(ndisk), "AZ1W_{:d}".format(ndisk), "AZ1W_{:d} = ".format(ndisk)+" ".join(str(e) for e in width))
    Template.insert("AZ1W_{:d}".format(ndisk), "AZ2P_{:d}".format(ndisk), "AZ2P_{:d} = 217.".format(ndisk))
    Template.insert("AZ2P_{:d}".format(ndisk), "AZ2W_{:d}".format(ndisk), "AZ2W_{:d} = ".format(ndisk)+" ".join(str(e) for e in width))
    # And we add streaming motions to the bar km/s
    Template.insert("VROT_{:d}".format(ndisk), "VRAD_{:d}".format(ndisk), "VRAD_{:d} = ".format(ndisk)+" ".join(str(e) for e in vrad_bar))
    Template["SBR_{:d}".format(ndisk)]="SBR_{:d} =".format(ndisk)+" ".join(str(e) for e in brightness)
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
def create_flare(Radii,velocity,dispersion,flare,Max_Rad,sub_ring,distance=1.):
    '''This function creates the scale height and  dispersion profiles according to Puche et al. 1992'''
    #make sure we have arrays so we can do numerical operations
    Radii = np.array(Radii)
    velocity = np.array(velocity)
    # first we need the dispersion
    # if there is no change from inside to outside it is easy
    if dispersion[0] == dispersion[1]:
       disp=np.full(len(Radii),dispersion[0])
    else:
        # else we use a tangent change with the center at halfway
       disp=-1*np.arctan((Radii-np.mean(Max_Rad/2.))/(Radii[-1]/10.))/np.pi*np.absolute(dispersion[0]-dispersion[1])+np.mean(dispersion)
    # We recalculate the Dynamical Mass and densities
    # Dynamical Mass

    Dynmass=(Radii*10**3)*velocity**2/G_agc
    # This is in a Volume of
    Volume=(4./3.)*np.pi*Radii**3
    Dynmass[0]=1.
    Volume[0]=1.
    Density=Dynmass/Volume*2.  #  In M_solar/kpc^3 the two comes from Puche 1992 but seems random
    # Then we check wether we want a flare or not
    G2=G_agc/(3.086e+13**2) #pc^3 M_sol^-1 s^-2
    halfint=int((len(Radii[Radii < Max_Rad])+10)/2.)
    if flare.lower() == 'flare':
        flare =disp/((4.*np.pi*G2*Density/1000.**3)**0.5*3.086e+16) # in kpc
        flare[:halfint-10] = flare[halfint-10]
        fact=np.arange(1/21,1,1./21)
        flare[halfint-10:halfint+10] = (1-fact)*flare[halfint-10]+fact*flare[halfint-10:halfint+10]
    elif flare.lower() == 'no_flare':
        flare = np.full(len(Radii),disp[halfint]/((4.*np.pi*G2*Density[halfint]/1000.**3)**0.5*3.086e+16))
    else:
        print(f"{flare} is not an option for the flare. Choose Flare or No_Flare")
        sys.exit()

    flare[0]=flare[1]

    # convert the scale heights to arcsec
    h_z_arcsec = cf.convertskyangle(flare,distance=distance,physical = True)
    # and write both to the Template
    Template["SDIS"]="SDIS = "+" ".join(str(e) for e in disp)
    Template["SDIS_2"]="SDIS_2 = "+" ".join(str(e) for e in disp)
    # The scaleheight
    Template["Z0"]="Z0 = "+" ".join(str(e) for e in h_z_arcsec)
    Template["Z0_2"]="Z0_2 = "+" ".join(str(e) for e in h_z_arcsec)
    return flare,disp

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

# A function for varying the PA and inclination as function of radius
def create_warp(Radii,
                PA,
                inclination,
                warp_change,
                warp_radii,
                sub_ring,
                disk=1):
    '''A function for varying the PA and inclination as function of radius'''
    Radii =np.array(Radii)
    if ((np.sum(warp_change) != 0) and (warp_radii[0] < warp_radii[1])).all():
        # First we need to obtain the vector that constitutes the inner area
        #it runs exactly counter to inclination
        inclination=90-inclination
        # For this the PA has to be between 0-90
        mult=np.floor(PA/90.)
        inPA= PA-mult*90.
        #avoid singularities
        if inPA == 0.:
            inPA = 0.001
        if inclination < 0.001 :
            inclination = 0.001
        # define the angular momentum vector of the plane and the outer most ring
        theta=np.arctan(np.tan(inclination*(np.pi/180.))*np.tan(inPA*(np.pi/180.)))
        phi=np.arctan(np.tan(inPA*(np.pi/180.))/np.sin(theta))
        # and the maximum values at Rad_HI
        thetamax=theta+warp_change[0]
        phimax=phi+warp_change[1]
                # indices of the warp start and the Rad_HI
        start_index = int(np.sum(Radii < warp_radii[0]))
        end_index = int(np.sum(Radii < warp_radii[1]))
        # step size of theta and phi
        # As we will increase the step size triangular we need the total number of point in the sequence
        warprange = end_index-start_index
        if warprange < 2:
            thetamax=theta
            phimax=phi
            warprange=1
        increasetheta=(thetamax-theta)/(0.5*warprange*(warprange+1))
        increasephi=(phimax-phi)/(0.5*warprange*(warprange+1))
        #print(warprange,thetamax,phimax,Radii[1]-Radii[0],warp_radii[1]-warp_radii[0])
        # calculate theta
        thetarings = np.array(np.full(len(Radii),theta))
        index_array=np.arange(start_index,len(Radii))-start_index
        thetarings[start_index:] = theta+0.5*index_array*(index_array+1)*increasetheta
        #calculate phi
        phirings = np.array(np.full(len(Radii),phi))
        phirings[start_index:] = phi+0.5*index_array*(index_array+1)*increasephi
        # return to PA

        if (phirings[0] < np.pi/2.) and (phirings[-1] > np.pi/2) and (inclination < 5.):
            PA= np.arctan(np.sin(thetarings)*np.tan(phirings-np.pi/2.))*(360./(2*np.pi))+mult*90+np.arctan(np.sin(theta)*np.tan(phi))*(360./(2*np.pi))
        else:
            PA= np.arctan(np.sin(thetarings)*np.tan(phirings))*(360./(2*np.pi))+mult*90
            PA[phirings > 0.5*np.pi]=PA[phirings > 0.5*np.pi]+180.

        # return inclination
        inc=90-np.arctan(1./(np.cos(thetarings)*np.tan(phirings)))*(360./(2*np.pi))
        # return inclination boundary adjustements
        inc[np.where(inc > 90.)] = 180 - inc[np.where(inc > 90.)]
        inc[np.where(inc < 0.)] = -1 * inc[np.where(inc < 0.)]
         # return a correct quadrant phirings
        phirings=phirings+mult/2.*np.pi
    else:
        #if there is no warp then all values are the same
        PA =np.full(len(Radii), PA)
        inc =np.full(len(Radii), inclination)
        phirings=np.full(len(Radii),0)
        thetarings=np.full(len(Radii),0)



    #write to our template file
    # let's see if we can retrace intrinsic phi with the formula's from Peters
    #theta_test = np.arctan((np.sin(PA*np.pi/180.)*np.sin(inc*np.pi/180.)-np.cos(PA*np.pi/180.)*np.sin(inc*np.pi/180.))/(np.cos(PA*np.pi/180.)*np.cos(inc*np.pi/180.)-np.sin(PA*np.pi/180.)*np.cos(inc*np.pi/180.)))*180./np.pi
    phirings=phirings*180./np.pi
    # thetarings=thetarings*180./np.pi
    # This seems to work mostly but not at some extremes exactly for some reason
    # According to josh tan is required, and yes that makes it work at large angles as well.
    angle_adjust=np.tan((PA[0]-PA)*np.cos(inc*np.pi/180.)*np.pi/180)*180/np.pi
    if disk == 1:
        Template["INCL"]="INCL = "+" ".join(str(e) for e in inc)
        Template["PA"]="PA = "+" ".join(str(e) for e in PA)
        try:
            phase =  Template["AZ1P"].split('=')[1]
        except KeyError:
            phase =  0.
            Template.insert("AZ1W","AZ1P","AZ1P = "+" ".join(str(phase+(e)) for e in angle_adjust))
    else:
        Template["INCL_{:d}".format(disk)]="INCL_{:d} =".format(disk)+" ".join(str(e) for e in inc)
        Template["PA_{:d}".format(disk)]="PA_{:d} =".format(disk)+" ".join(str(e) for e in PA)
        try:
            phase =  float(Template["AZ1P_{:d}".format(disk)].split('=')[1])
        except KeyError:
            phase =  0.
        Template.insert("AZ1W_{:d}".format(disk),"AZ1P_{:d}".format(disk),"AZ1P_{:d} = 180.".format(disk))
        Template.insert("AZ1W_{:d}".format(disk),"AZ1P_{:d}".format(disk),"AZ1P_{:d} =".format(disk)+" ".join(str(phase+(e)) for e in angle_adjust))
    return PA,inc,phirings
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
    warp_radii = ??
    sub_ring = the sub rings

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
def create_inhomogeneity(mass,SNR,disks=1.):
    if len(disks) == 0:
        disks=[disks]

    for i in range(len(disks)):
        ndisk=int(Template["NDISKS"].split('=',1)[1])
        ndisk+=1
        last_add = copy_disk(disks[i],ndisk)
        sbr=np.array(Template["SBR_{:d}".format(ndisk)].split('=',1)[1].strip().split(" "),dtype=np.float64)
        rad=np.array(Template["RADI"].split('=',1)[1].strip().split(" "),dtype=np.float64)
        rad=rad*4.*np.pi/max(rad)
        # This part is tricky as we do not want negative emission but cutting out the negatives would result in added flux
        #Hence these profiles should always be smaller than the sbr
        newprof=np.sin(rad)*sbr*np.log10(mass)/11.*8./SNR*0.9
        #print("We are modifying by factor {} for mass {} and SNR {}".format(np.log10(mass)/11.*8./SNR, np.log10(mass),SNR))
        nextprof=-1*newprof

        req_flux=abs(np.sum(newprof*np.pi*rad*2)/200.)
        if req_flux == 0:
            req_flux = 1e-5
        Template["SBR_{:d}".format(ndisk)]="SBR_{:d} = ".format(ndisk)+" ".join(str(e) for e in newprof)
        Template["CFLUX_{:d}".format(ndisk)]="CFLUX_{:d} = {}" .format(ndisk,req_flux)
        ndisk=ndisk+1
        last_add = copy_disk(disks[i],ndisk)
        Template["SBR_{:d}".format(ndisk)]="SBR_{:d} = ".format(ndisk)+" ".join(str(e) for e in nextprof)
        Template["CFLUX_{:d}".format(ndisk)]="CFLUX_{:d} = {}" .format(ndisk,req_flux)
    #print("We will create inhomogeneities on top of disk(s) ="+" ".join(str(e) for e in [disks]))
def create_mask(work_dir,beam,casa=False):
    # First open the model cube

    dummy = fits.open(work_dir+'/unconvolved_cube.fits',uint = False, do_not_scale_image_data=True,ignore_blank = True)
    # Add to the hedear some info
    dummy[0].header.append('RESTFRQ')
    dummy[0].header['RESTFRQ'] = 1.420405752E+09
    # downgrade the velocity resolution to the one we want
    tmp_cube=np.zeros([int(len(dummy[0].data[:,0,0])/3),len(dummy[0].data[0,:,0]),len(dummy[0].data[0,0,:])])
    for n in range(1,len(dummy[0].data[:,0,0]),3):
        tmp_cube[int((n-1)/3),:,:] = np.mean(dummy[0].data[n-1:n+2,:,:],axis=0)
    dummy[0].data = tmp_cube
    #print(dummy[0].data.shape)
    # reset the header
    dummy[0].header['CDELT3'] = dummy[0].header['CDELT3']*3.
    dummy[0].header['CRPIX3'] = dummy[0].header['CRPIX3']/3.
    dummy[0].header['NAXIS3'] = int(dummy[0].header['NAXIS3']/3.)

    if casa:
        fits.writeto(work_dir+'/unconvolved_cube.fits',dummy[0].data,dummy[0].header, overwrite = True)
    # In order to create a cleaning mask we smooth to the beam size and cut at 1e-5 Jy/beam
    #   Calculate the sigma's from the required beam size
    sigma=[(beam[0]/abs(dummy[0].header['CDELT1']*3600.))/(2*np.sqrt(2*np.log(2))),(beam[1]/abs(dummy[0].header['CDELT2']*3600.))/(2*np.sqrt(2*np.log(2)))]
    #
    reals = dummy[0].data[dummy[0].data > 1e-5]
    # Remember that python is a piece of shiit so z,y,x
    cutoff = (np.mean(reals)/2.)
    # correct the cutoff for the smoothing
    cutoff=cutoff/(0.5*(2*np.sqrt(np.pi)*sigma[0]+sigma[1]*np.sqrt(np.pi)*2.))
    # smooth the image
    #our BPA is set to 0 which means the major axis smoothing should be in DEC axis and the minor on the RA
    smooth = ndimage.gaussian_filter(dummy[0].data, sigma=(0,sigma[0], sigma[1]), order=0)

    # We want the mean signal in the smoothed cube
    # !!!!! This one is not correct for conserved surface brightness temperature but we want it to estimate the noise per pixel!!!!!!!!
    #But tcorrection is applied to the final cube so we need this one
    smooth[smooth < cutoff]=0
    mean_signal = cf.get_mean_flux(smooth)
    #mean_signal= np.mean(smooth[smooth > cutoff])
    # let's check that the minimum is less than the noise.
    #if np.min(smooth) < mean_signal/SNR*-3:
    #    print(np.min(smooth),mean_signal/SNR*-3)
    #    exit()
    #print("This is the minimum signal {} and comparable noise {}".format( np.min(smooth),mean_signal/SNR*-3))
    # Then create the mask
    smooth[smooth > cutoff]=1
    # Write mask
    fits.writeto(work_dir+'/mask.fits',smooth,dummy[0].header, overwrite = True)
    hdr = copy.deepcopy(dummy[0].header)
    data = copy.deepcopy(dummy[0].data)
    dummy.close()
    return mean_signal,hdr,data

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



# this is a routine to use the casa task simobserve to corrupt the observations.
def corrupt_casa(work_dir,beam,SNR,casa_call='casa'):
    '''Corrupt our artifical galaxy with a casa's sim observe routines'''

    Template_Casa = cf.read_casa_template('Template_Casa.py')
    mean_signal,hdr,data= create_mask(work_dir,beam,casa=True)
    # In order to corrupt we need to know the average signal.
    # we do this taking the mean in each chaneel above a tenth of the max and then take the mean of that profile
    # This is the noise in the final cube
    #!!!! Is this correct or should we conserve for brightness??
    noise=mean_signal/SNR
    #*(abs(dummy[0].header['CDELT1']*3600.)*abs(dummy[0].header['CDELT2']*3600.))/1.1330900354567984


    # We need to know the location in the sky dec > 30 we use WSRT 30 > dec > -30 --> VLA -30 > dec  --> Atca

    # Full synthesis 12 hrs leads to noise levels 12: 0.000701111 ,24:  0.00065377   ,48:0.000499886 ,96:0.000436001 ,300:0.000213708 ,600:0.000178363   ,900:0.000149931 1200:0.000119707 2400:9.85791e-05 , 4800: 8.18873e-05 # Based on a simulation with no taper
    noisin = [0.000701111,0.00065377 ,0.000499886,0.000436001,0.000213708,0.000178363,0.000149931,0.000119707 ,9.85791e-05 , 8.18873e-05]
    timein = [12.,24,48,96,300,600,900,1200,2400,4800]
    a1,a2 =np.polyfit(noisin,timein,1)
    totaltime = a1*noise+a2
    if totaltime < 12:
        totaltime = 12.
    RA,DEC =cf.convertRADEC(hdr['CRVAL1'],hdr['CRVAL2'])
    tri = open(work_dir+'/pntings.txt', 'w')
    tri.writelines("#Epoch     RA          DEC      TIME(optional) \n ")
    tri.writelines("J2000     {}          {}      {} ".format(RA,DEC,str(int(12*3600.))))

    tri.close()

    if  hdr['CRVAL2'] > 90:
        Template_Casa['simobserve_antennalist']="antennalist = 'WSRT.cfg'  #  interferometer antenna position file"
        Template_Casa['imhead_hdvalue'] = "hdvalue = 'WSRT'"
        Template_Casa['tclean_vis'] = "vis = 'simulated/simulated.WSRT.noisy.ms'"
    elif hdr['CRVAL2'] > -30:
        Template_Casa['simobserve_antennalist']="antennalist = 'vla.b.cfg'  #  interferometer antenna position file"
        Template_Casa['imhead_hdvalue'] = "hdvalue = 'VLA'"
        Template_Casa['tclean_vis'] = "vis = 'simulated/simulated.vla.b.noisy.ms'"
    else:
        Template_Casa['simobserve_antennalist']="antennalist = 'atca_6c.cfg'  #  interferometer antenna position file"
        Template_Casa['imhead_hdvalue'] = "hdvalue = 'ATCA'"
        Template_Casa['tclean_vis'] = "vis = 'simulated/simulated.atca_6c.noisy.ms'"
    #let's assure the same cube size
    # We want 2000 integrations
    Template_Casa['simobserve_integration'] = "integration    = '{:d}s'".format(int(totaltime*3600/2000.))
    Template_Casa['simobserve_totaltime'] = "totaltime    = '{:d}'".format(int(totaltime/12.))
    Template_Casa['tclean_cell']= "cell = ['{}arcsec','{}arcsec']".format(abs(hdr['CDELT1']*3600.),abs(hdr['CDELT2']*3600.))
    Template_Casa['tclean_imsize'] = "imsize=[{:d},{:d}]".format(int(abs(hdr['NAXIS1'])),int(abs(hdr['NAXIS2'])))
    Template_Casa['tclean_scales'] = "scales=[0,{:d},{:d}]".format(int(2.*beam[0]/abs(hdr['CDELT1']*3600.)),int(5.*beam[1]/abs(hdr['CDELT1']*3600.)))
    Template_Casa['tclean_threshold'] = "threshold = '{}Jy/beam'".format(noise/2.)
    Template_Casa['tclean_uvtaper'] = "uvtaper = ['{}arcsec','{}arcsec']".format(beam[0],beam[1])

    # In order to create a cleaning mask we smooth to twice the beam size and cut at 1e-5 Jy/beam
    with open(work_dir+'/run_casa.py', 'w') as tri:
        tri.writelines([Template_Casa[key]+"\n" for key in Template_Casa])

    os.chdir(work_dir)
    bla = subprocess.call([casa_call,'--nologger','--nogui','-c','run_casa.py'])

    dummy = fits.open(work_dir+'/Convolved_Cube.fits',uint = False, do_not_scale_image_data=True,ignore_blank = True)
    # We cut 20 pixels around the edges
    newsize = int(np.shape(dummy[0].data)[2]-60.)
    newdummy=np.zeros((np.shape(dummy[0].data)[0],newsize,newsize))
    newdummy[:,:,:]=dummy[0].data[:,int(np.floor(dummy[0].header['CRPIX2']-newsize/2.)):int(np.floor(dummy[0].header['CRPIX2']+newsize/2.)),int(np.floor(dummy[0].header['CRPIX1']-newsize/2.)):int(np.floor(dummy[0].header['CRPIX1']+newsize/2.))]
    dummy[0].header['NAXIS1']=newsize+1
    dummy[0].header['NAXIS2']=newsize+1

    dummy[0].header['CRPIX1']=dummy[0].header['CRPIX1']-np.floor(dummy[0].header['CRPIX1']-newsize/2.)
    dummy[0].header['CRPIX2']=dummy[0].header['CRPIX2']-np.floor(dummy[0].header['CRPIX2']-newsize/2.)
    fits.writeto(work_dir+'/Convolved_Cube.fits',newdummy,dummy[0].header, overwrite = True)
    # Also the mask then

    dummy = fits.open(work_dir+'/mask.fits',uint = False, do_not_scale_image_data=True,ignore_blank = True)
    # We cut 20 pixels around the edges
    newsize = int(np.shape(dummy[0].data)[2]-60.)
    newdummy=np.zeros((np.shape(dummy[0].data)[0],newsize,newsize))
    newdummy[:,:,:]=dummy[0].data[:,int(np.floor(dummy[0].header['CRPIX2']-newsize/2.)):int(np.floor(dummy[0].header['CRPIX2']+newsize/2.)),int(np.floor(dummy[0].header['CRPIX1']-newsize/2.)):int(np.floor(dummy[0].header['CRPIX1']+newsize/2.))]
    dummy[0].header['NAXIS1']=newsize+1
    dummy[0].header['NAXIS2']=newsize+1

    dummy[0].header['CRPIX1']=dummy[0].header['CRPIX1']-np.floor(dummy[0].header['CRPIX1']-newsize/2.)
    dummy[0].header['CRPIX2']=dummy[0].header['CRPIX2']-np.floor(dummy[0].header['CRPIX2']-newsize/2.)
    fits.writeto(work_dir+'/mask.fits',newdummy,dummy[0].header, overwrite = True)
corrupt_casa.__doc__=f'''
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
   Template_Casa = The template for the sim_observe python file for casa

OPTIONAL INPUTS:
   casa_call = casa command line command

OUTPUTS:
   the Convolved and corrupted cube is writen to disk


OPTIONAL OUTPUTS:

PROCEDURES CALLED:
   Unspecified

NOTE: The final SNR is a guesstimate as the actuall observation is based on the time. This part could be improved.
'''

def corrupt_gauss(work_dir,beam,SNR):
    '''Corrupt our artifical galaxy with the correct noise cube such that the average SNR over the galaxy is the required one'''
    mean_signal,hdr,data= create_mask(work_dir,beam)
    # Calculate the area of the beam in arcsec
    beamarea=(np.pi*abs(beam[0]*beam[1]))/(4.*np.log(2.))
    # Convert arcsec to pixels
    pixperbeam=beamarea/(abs(hdr['CDELT1']*3600.)*abs(hdr['CDELT2']*3600.))
    # Make an the size of the model with random values distributed as the noise
    # The noise we want is the mean signal divided by the signal to noise ratio
    # that value needs to be deconvolved so from https://en.wikipedia.org/wiki/Gaussian_blur
    # The formula is uncited but works
    sigma=[(beam[0]/abs(hdr['CDELT1']*3600.))/(2*np.sqrt(2*np.log(2))),(beam[1]/abs(hdr['CDELT2']*3600.))/(2*np.sqrt(2*np.log(2)))]

    noisescl=(mean_signal/SNR*sigma[0]*2*np.sqrt(np.pi))
    #print("This our mean signal {} and noise {} in Jy/pixel. ".format(mean_signal,noisescl))

    if beam[2] !=0. :
        #As we are going to rotatet the cube we should first extend it
        shift = int(abs(np.sin(np.radians(beam[2])))*hdr['NAXIS1']/2.+5)
        Pix_Extend= [shift,shift]
        data = np.pad(data, [[0, 0], Pix_Extend, Pix_Extend], 'constant')
        data = cf.rotateCube(data,beam[2],[hdr['CRPIX1']+shift,hdr['CRPIX2']+shift],order=1)

    cuberms = np.random.normal(scale=noisescl,size=np.shape(data))
    # combine the two cubes

    noisedcube=data+cuberms
    # Smooth to the requred resolution

    #our BPA is set to 0 which means the major axis smoothing should be in DEC axis and the minor on the RA
    final = ndimage.gaussian_filter(noisedcube, sigma=(0,sigma[0], sigma[1]), order=0)

    if beam[2] != 0.:
        #rotate back
        final_tmp = cf.rotateCube(final,-1*(beam[2]),[hdr['CRPIX1']+shift,hdr['CRPIX2']+shift],order=1)
        final = copy.deepcopy(final_tmp[:, \
            shift:hdr['NAXIS2'] + \
            shift,shift:hdr['NAXIS1']\
             + shift])
        final_tmp =[]

    # to preserve brightness temperature this should be multiplied with the increase in beam area
    # which is the same as the amount of pixels in the beam as we go from 1 pixel area to an area the size of the beam which is assuming two gaussians so we need to correct with a difference factor of the area of the. Last factor is to correct for the fact that the unconvolved cube hassquare pixels not a circular beam.
    final=final*pixperbeam
    # And write this final cube to the directory
    hdr['BMAJ']=beam[0]/3600.
    hdr['BMIN']=beam[1]/3600.
    hdr['BPA'] = beam[2]
    fits.writeto(work_dir+'/Convolved_Cube.fits',final,hdr, overwrite = True)
corrupt_gauss.__doc__=f'''
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

def get_HI_disk(Mass,output_directory=None):
    # About 18.5% is cosmic baryon fraction from planck (Omb=0.048, OmDM = 0.2589). In a galactic halo this is roughly 70%-90% of the cosmic mean (See Crain et al. 2006, Pezzulli 2019)
    bary_frac = 0.1854
    diff=0.9
    counter = 0
    # if the baryon mass is higher the 2* a typical gass mass disk (i.e. M_HI = 1e9 Msol)
    # then we subtract this disk from the baryonic fraction and count the remainder as the initial stellar mass guess
    if bary_frac*Mass > 2.8*10**9:
        m_star=bary_frac*Mass-1.4*10**9
    else:
        # else half baryonic in stars
        m_star=bary_frac*Mass/2.

    # and convert stellar mass to an HI Mass following van Driel (2016)
    m_star_prev=1
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
        m_star_prev=copy.deepcopy(m_star)
        m_star=bary_frac* Mass-M_HI*1.4

    #print("This is MHI {:.2e} and mstar {:.2e}".format(MHI,m_star,m_star_prev))
    # the mass leads to the radius at which we need to hit 1 M/pc^2 from Wang (2016) in kpc
    R_HI=10**(0.506*np.log10(M_HI)-3.293)/2.

    # McGaugh 2014  M_*/Lk = 0.6
    M_K = np.log10(m_star/0.6)*-2.5+3.29
    # Ponomareva 2017 has
    # log(LT,b,i) = (3.7 ± 0.11) × log(2Vflat) + 1.3 ± 0.3 for 3.6 mu  and then in Table 3 the values for Ks band
    v_circ_TF = 10**((np.log10(m_star/0.6)-1.22)/3.81)/2.

    #  Let's check that this roughly gives our DM mass at the virial radius in pc
    v_circ_NFW = calc_vc_NFW(Mass,M_HI,m_star,R_HI)
    #our final velocity is an average between TF and NFW
    V_HI=(v_circ_TF+v_circ_NFW )/2.

    if output_directory:
        DynMass=R_HI*10**3*V_HI**2/G_agc
        DynMassNFW= R_HI*10**3*v_circ_NFW**2/G_agc
        with open(f'{output_directory}Fractions_and_Masses.txt','a') as file:
            file.write(f'''The input Mass = {Mass:.2e}  and the retrieved NFW Dynamical mass = {DynMassNFW:.2e} and Dynamical Mass based on v_circ = {DynMass:.2e}.
The current the baryon fraction = {bary_frac:.5f}\n''')
    return V_HI, R_HI, M_HI,M_K
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

def get_sparcs_fit(v_end,radius):
    '''Calculate the R_0, gamma and beta parameters'''

    v_flat=11.+0.76*v_end
    R_0 = 10*np.exp(-1*radius**2/10.**2)+2.25+0.08*radius
    beta = 0.5*(1.14-0.0006*v_end)+0.5*(0.7*np.arctan((radius - 32.)/25.))
    if beta < 0.25:
        beta = 0.25
    gamma = 0.5*(gaussian_function(v_end,1.5,170.,60.))+0.5*(3.7- 0.003*radius)
    return v_flat,R_0,beta,gamma

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

def gaussian_function(x,amp,center,sigma):
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

def one_galaxy(cfg,Current_Galaxy):
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
    os.system(f"cp {cfg.general.main_directory}/Input.fits {cfg.general.main_directory}/{name}/Input.fits  ")
    Template = cf.read_template_file('Template.def')

     # First set the beam to 0.
    Template["BMAJ"]= "BMAJ = 0."
    Template["BMIN"]= "BMIN = 0."
    Template["BPA"]= "BPA = 0."


    #Then we need to build the Surface Brightnes profile
    SBRprof,Rad,sclength,MHI,Rad_HI,Vrot,sub_ring,molecular_profile = build_sbr_prof(Current_Galaxy,symmetric=cfg.agc.symmetric) #Column densities,Raii in kpc, Opt_scalelength in kpc, HI mass in M_solar
    #We need central coordinates the vsys will come from the required distance and hubble flow. The RA and dec should not matter hance it will be the only random component in the code as we do want to test variations of them
    Sky_Size = np.radians(Current_Galaxy.Res_Beam[0]*Current_Galaxy.Beams/3600.)
    Distance = (Rad_HI[0]/(np.tan(Sky_Size/2.)))/1000.
    #print("The Distance is {:5.2f} Mpc".format(Distance))
    vsys = Distance*H_0
    #if cfg.agc.corruption_method == 'Gaussian' or (cfg.agc.corruption_method == 'Casa_5' and (int(number_models/5.) != number_models/5.)):
    #    RAdeg=np.random.uniform()*360
    #    DECdeg=(np.arccos(2*np.random.uniform()-1)*(360./(2.*np.pi)))-90
    #else:
    #    RAdeg = np.random.uniform()*360
    #    DECdeg=-60
    #    while DECdeg < -20.:
    #        DECdeg = (np.arccos(2*np.random.uniform()-1)*(360./(2.*np.pi)))-90
    RAhr,DEChr= cf.convertRADEC(*Current_Galaxy.Coord)
    #print("It's central coordinates are RA={} DEC={} vsys={} km/s".format(RAhr,DEChr,vsys))
    # Write them to the template
    for ext in ['','_2']:
        Template[f"XPOS{ext}"]=f"XPOS{ext}= {Current_Galaxy.Coord[0]}"
        Template[f"YPOS{ext}"]=f"YPOS{ext}= {Current_Galaxy.Coord[1]}"
        Template[f"VSYS{ext}"]=f"VSYS{ext}= {vsys}"

    # The Distance, this is not required but can be usefull
    Template["DISTANCE"]="DISTANCE= {}".format(Distance)
    # With the distance we can also convert our radii
     # we need our radii in the corresponding arcsec
    Rad_arcsec = cf.convertskyangle(Rad,distance=Distance,physical = True)
     # then the number of rings to the total number of rings
    Template["NUR"]="NUR= {}".format(len(Rad))
    # The radii in arcsec
    Template["RADI"]="RADI = "+" ".join(str(e) for e in Rad_arcsec)
    # Then we need to get a starting radius for the warp.
    # the warp should start at the edge of the optical radius which is the HI scale length/0.6
    # which are about ~ 4 * h_r
    WarpStart = 4.*sclength
    WarpEnd=Rad[np.where(SBRprof >= 4.98534620064e-05)[0][-1]]
    # Write it to the Template
    Template["VROT"]="VROT = "+" ".join(str(e) for e in Vrot)
    Template["VROT_2"]="VROT_2 = "+" ".join(str(e) for e in Vrot)
    # We need a scale height and dispersion for each ring. They are coupled and hence they are both created in create_flare
    h_z,dispersion=create_flare(Rad,Vrot,Current_Galaxy.Dispersion,Current_Galaxy.Flare,Rad_HI[0],sub_ring,distance=Distance)

    # Finally we need to set the warping
    PA,inc,phirings = create_warp(Rad,Current_Galaxy.PA,Current_Galaxy.Inclination,Current_Galaxy.Warp,[WarpStart,WarpEnd],sub_ring)
    if cfg.agc.symmetric:
        PA_2,inc_2,phirings_2 = create_warp(Rad,Current_Galaxy.PA,Current_Galaxy.Inclination,Current_Galaxy.Warp,[WarpStart,WarpEnd],sub_ring,disk=2)
    else:
        # we want an assymetric warp so we redo the PA and inc but swap the variation
        PA_2,inc_2,phirings_2 = create_warp(Rad,Current_Galaxy.PA,Current_Galaxy.Inclination,[Current_Galaxy.Warp[0]-Current_Galaxy.Warp[1],Current_Galaxy.Warp[1]/2.+Current_Galaxy.Warp[0]/2.],[WarpStart,WarpEnd],sub_ring,disk=2)

    #If we want Radial Motions then they need to be inserted
    if Current_Galaxy.Radial_Motions != 0.:
          Template.insert("VROT","VRAD","VRAD = {}".format(Current_Galaxy.Radial_Motions))
          Template.insert("VROT_2","VRAD_2","VRAD_2 = {}".format(Current_Galaxy.Radial_Motions))

    # This comes from FAT. If I remember correctly this is the sine response to the channels *1.2/(2*SQRT(2*ALOG(2.))))
    # However in our input we want independent channels which means we should set this to 0.
    # Template["CONDISP"]="CONDISP = 0."
    if cfg.agc.channel_dependency == 'sinusoidal':
        Template["CONDISP"]="CONDISP = {}".format(Current_Galaxy.Channelwidth*1.2/(2*np.sqrt(2*np.log(2.))))
    elif cfg.agc.channel_dependency == 'hanning':
        Template["CONDISP"]="CONDISP = {}".format(Current_Galaxy.Channelwidth*2./(2*np.sqrt(2*np.log(2.))))
    elif cfg.agc.channel_dependency == 'independent':
        Template["CONDISP"]="CONDISP = 0."
    else:
        raise InputError(f"{cfg.agc.channel_dependency} is not an option for the channel dependency")
    # We need to set the input and output cube
    Template["INSET"]="INSET = Input.fits"
    Template["OUTSET"]="OUTSET = unconvolved_cube.fits"
    #Some tirific varaiables
    Template["LOOPS"]="LOOPS = 0 "
    # We need models with about 3 million particles but not more as it takes too long
    TotFlux=MHI/(2.36e5*(Distance)**2)
    Template["CFLUX"]="CFLUX = {}".format(TotFlux/5e6)
    Template["CFLUX_2"]="CFLUX_2 = {}".format(TotFlux/5e6)
    Template.insert("RMS","NDISKS","NDISKS = 2")
    Template["TIRDEF"]="TIRDEF = ModelInput.def"
    Template["GR_DEVICE"]="GR_DEVICE = "
    Template["GR_CONT"]="GR_CONT = "
    Template.insert("INSET","ACTION","ACTION = 1")
    Template["PROGRESSLOG"]="PROGRESSLOG = "
    Template["LOGNAME"]= "LOGNAME= "

    #-------------------------------This finishes the basic disk the following are optional components----------------------

    # The possible arms
    if Current_Galaxy.Arms == 'Arms':
        phase,arm_brightness,arm_width = create_arms(Vrot,Rad,SBRprof,WarpStart=WarpStart, Bar=Current_Galaxy.Bar)
        phase,arm_brightness,arm_width = create_arms(Vrot,Rad,SBRprof,disk=2,WarpStart=WarpStart, Bar=Current_Galaxy.Bar)
    # A possible Bar
    if Current_Galaxy.Bar == 'Barred':
        bar_length = create_bar(Vrot,Rad,SBRprof,Template,WarpStart=WarpStart)
    # and possible inhomogeneities
    if cfg.agc.inhomogenous:
        inhomogeneity_amp = create_inhomogeneity(MHI,Current_Galaxy.SNR,disks=[1,2])
    # we write the def files
    with open(f"{cfg.general.main_directory}{name}/tirific.def", 'w') as file:
        file.writelines([Template[key]+"\n" for key in Template])

    # So we need to  modify the input file to the correct coordinates else we'll get an empty cube

    dummy = fits.open(f"{cfg.general.main_directory}{name}/Input.fits",uint = False, do_not_scale_image_data=True,ignore_blank = True)
    #first we do a check wether the BMAJ
    #is written correctly for the python
    #program. We set a bunch of generic header values.
    #First we make sure that it fits

    size = 2.*Rad_arcsec[-1] +(Current_Galaxy.Res_Beam[0])
    if  Current_Galaxy.Beams < 6.:
        size = 2.*Rad_arcsec[-1] +(6.*Current_Galaxy.Res_Beam[0])
    #pix_size = (Rad_arcsec[1]-Rad_arcsec[0])
    pix_size = Current_Galaxy.Res_Beam[1]/5.
    if Current_Galaxy.Corruption =='Casa_Sim':
        size += 60*pix_size
    required_pixels=int(np.ceil(size/pix_size))
    vel_max=2.5*np.max([Vrot*np.sin(inc*np.pi/180.)+4*dispersion,Vrot*np.sin(inc_2*np.pi/180.)+4*dispersion])
    velpix=int(np.ceil(vel_max/Current_Galaxy.Channelwidth)*3 )
    dummy[0].header['CRPIX1'] = np.floor(required_pixels/2.)
    dummy[0].header['CRPIX2'] = np.floor(required_pixels/2.)
    dummy[0].header['CRPIX3'] = np.floor(velpix/2.)
    # Stupid astropy doesn't account for minus in header of cdelt1 then cdelt has different precision
    tmp=int(pix_size/3600.*1e15)/1e15
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
    dummy[0].header['OBJECT'] = f'AGC_GALAXY'
    dummy[0].header['INSTRUME'] =f'AGC'
    try:
        del  dummy[0].header['BLANK']
    except:
        pass
    # make the cube a typical size
    dummy2 = np.zeros((velpix,required_pixels,required_pixels),dtype= np.float32)
    dummy2[int(np.floor(velpix/2.)),int(np.floor(required_pixels/2.)),int(np.floor(required_pixels/2.))] = 5
    fits.writeto(f"{cfg.general.main_directory}{name}/Input.fits",dummy2,dummy[0].header, output_verify='silentfix+ignore', overwrite = True)
    dummy.close()

    current_run = subprocess.Popen([cfg.general.tirific,"DEFFILE=tirific.def","ACTION=1"],\
                           stdout = subprocess.PIPE, stderr = subprocess.PIPE,\
                           cwd=f"{cfg.general.main_directory}{name}",universal_newlines = True)
    tirific_run, tirific_warnings_are_annoying = current_run.communicate()
    #print(tirific_run)
    if current_run.returncode == 1:
        pass
    else:
        print(tirific_warnings_are_annoying)
        raise TirificRunError("AGC:Tirific did not execute properly. See screen for details")
    #os.chdir(f"{cfg.general.main_directory}{name}")
    #os.system(f"{cfg.general.tirific} deffile=tirific.def")

    #os.chdir(cfg.general.main_directory)


    # Now we want to corrupt this cube with some realistic noise
    # For this we first want to get the noise we want in terms of Jansky per beam
    # we will define the SNR as the mean(Intensity)/noiselevel hence noise =mean(In)/SNR
    if Current_Galaxy.Corruption == 'Casa_Sim':
        corrupt_casa(f"{cfg.general.main_directory}{name}/",Current_Galaxy.Res_Beam,Current_Galaxy.SNR,casa_call=cfg.general.casa)
        os.chdir(cfg.general.main_directory)
        mask = fits.open(f"{cfg.general.main_directory}{name}/mask.fits",uint = False, do_not_scale_image_data=True,ignore_blank = True)
        Cube = fits.open(f"{cfg.general.main_directory}{name}/Convolved_Cube.fits",uint = False, do_not_scale_image_data=True,ignore_blank = True)
        #Here we have no control over the BPA it is what it is.
        Current_Galaxy.Res_Beam[2]=Cube[0].header['BPA']
        maskr = mask[0].data[1:]
        sigma = (np.std(Cube[0].data[0])+np.std(Cube[0].data[-1]))/2.
        Cube_Clean = Cube[0].data
        Cube_Clean[maskr < 0.5] = 0.
        beamarea=(np.pi*abs(Cube[0].header['BMAJ']*3600.*Cube[0].header['BMIN']*3600.))/(4.*np.log(2.))
        pixperbeam=beamarea/(abs(Cube[0].header['CDELT1']*3600.)*abs(Cube[0].header['CDELT2']*3600.))
        totalsignal = np.sum(Cube_Clean)/pixperbeam
        mass = 2.36E5*Distance**2*totalsignal*Cube[0].header['CDELT3']/1000.
        #totsig=np.zeros(len(Cube_Clean[:]))
        #for j in range(len(totsig)):
        #    if  len(Cube_Clean[j][Cube_Clean[j] > 0.]) > 0:
        #        totsig[j]=np.mean(Cube_Clean[j][Cube_Clean[j] > 0.])
        #mean_signal = np.median(totsig[totsig > 0.])
        mean_signal=cf.get_mean_flux(Cube_Clean)
        SNRachieved = mean_signal/(sigma)
    elif  Current_Galaxy.Corruption == 'Gaussian':
        corrupt_gauss(f"{cfg.general.main_directory}{name}/",Current_Galaxy.Res_Beam,Current_Galaxy.SNR)
        mask = fits.open(f"{cfg.general.main_directory}{name}/mask.fits",uint = False, do_not_scale_image_data=True,ignore_blank = True)
        Cube = fits.open(f"{cfg.general.main_directory}{name}/Convolved_Cube.fits",uint = False, do_not_scale_image_data=True,ignore_blank = True)
        maskr = mask[0].data[:]
        sigma = (np.std(Cube[0].data[0])+np.std(Cube[0].data[-1]))/2.
        Cube_Clean = Cube[0].data
        Cube_Clean[maskr < 0.5] = 0.
        beamarea=(np.pi*abs(Cube[0].header['BMAJ']*3600.*Cube[0].header['BMIN']*3600.))/(4.*np.log(2.))
        pixperbeam=beamarea/(abs(Cube[0].header['CDELT1']*3600.)*abs(Cube[0].header['CDELT2']*3600.))
        totalsignal = np.sum(Cube_Clean)/pixperbeam
        mass = 2.36E5*Distance**2*totalsignal*Cube[0].header['CDELT3']/1000.
        #mean_signal = np.mean(Cube_Clean[maskr > 0.5])
        mean_signal = cf.get_mean_flux(Cube_Clean)
        SNRachieved = mean_signal/(sigma)
        #if we have
    else:
        print("!!!!!!!This corruption method is unknown, leaving the cube uncorrupted and unconvolved!!!!!!!!")
        # We'll create a little text file with an Overview of all the parameters
    Template['BMAJ'] = f'BMAJ= {Current_Galaxy.Res_Beam[0]}'
    Template['BMIN'] = f'BMIN= {Current_Galaxy.Res_Beam[1]}'
    Template['BPA'] = f'BPA= {Current_Galaxy.Res_Beam[2]}'
    if cfg.agc.corrupt_models:
        beam_line = f"Major axis beam = {Current_Galaxy.Res_Beam[0]}, Minor axis beam= {Current_Galaxy.Res_Beam[1]}, Beam PA = {Current_Galaxy.Res_Beam[2]}."
        corrupt_line = f"The cube was corrupted with the {Current_Galaxy.Corruption} method."
        if Current_Galaxy.Corruption == 'Casa_Sim':
            catalog_cube_name = 'Convolved_Cube_CS'
            os.system(f"mv {cfg.general.main_directory}{name}/Convolved_Cube.fits {cfg.general.main_directory}{name}/Convolved_Cube_CS.fits")
        else:
            catalog_cube_name = 'Convolved_Cube'
        os.remove(f"{cfg.general.main_directory}{name}/unconvolved_cube.fits")
    else:
        beam_line = 'This galaxy is unconvolved.'
        corrupt_line = 'This galaxy is not corrupted.'
        SNRachieved = float('NaN')
        sigma=float('NaN')
        mass=float('NaN')
        catalog_cube_name = 'unconvolved_cube'
        mean_signal,hdr,data= create_mask(f"{cfg.general.main_directory}{name}/",Current_Galaxy.Res_Beam)
        beamarea=(np.pi*abs(hdr['BMAJ']*3600.*hdr['BMIN']*3600.))/(4.*np.log(2.))
        pixperbeam=beamarea/(abs(hdr['CDELT1']*3600.)*abs(hdr['CDELT2']*3600.))

        hdr=[]
        data=[]
    with open(f"{cfg.general.main_directory}{name}/{name}-Info.txt", 'w') as overview:
        overview.write(f'''This file contains the basic parameters of this galaxy.
For the radial dependencies look at Overview.png or ModelInput.def.
Inclination = {Current_Galaxy.Inclination}.
The dispersion = {dispersion[0]:.2f}-{dispersion[1]:.2f}.
The type of galaxy = {Current_Galaxy.Mass:1e}.
PA = {Current_Galaxy.PA}.
Warp = {Current_Galaxy.Warp[0]}-{Current_Galaxy.Warp[1]}.
Which starts at {WarpStart:.2f} kpc and the 1M/pc^2 radius is {Rad_HI[0]:.2f} kpc.
Flare = {Current_Galaxy.Flare}.
Beams across the major axis = {Current_Galaxy.Beams}.
SNR Requested = {Current_Galaxy.SNR} SNR Achieved = {SNRachieved}.
Mean Signal = {mean_signal}.
Channelwidth = {Current_Galaxy.Channelwidth} and their dependency is {cfg.agc.channel_dependency}.
{beam_line}
This galaxy has {Current_Galaxy.Arms} and a {Current_Galaxy.Bar}.
It's central coordinates are RA={RAhr} DEC={DEChr} vsys={vsys:.2f} km/s.
At a Distance of {Distance:.2f} Mpc.
HI_Mass Requested {MHI:.2e} (M_solar) and an optical h {sclength:.2f} (kpc).
HI_Mass Retrieved {mass:.2e} (M_solar).
We have {pixperbeam} pix per beam.
{corrupt_line}
The final noise level is {sigma} Jy/beam.
h_z = {h_z[0]:.3f}-{h_z[-1]:.3f} (kpc).''')

    # We also want a file that contains initial estimates for all the parameters. We scramble them with gaussian variations
    cf.scrambled_initial(f"{cfg.general.main_directory}{name}/",Template)
    cf.plot_input(f"{cfg.general.main_directory}{name}/",Template, \
                Title=f'DM Mass in = {Current_Galaxy.Mass:.2e}',RHI=Rad_HI
                ,add_sbr=molecular_profile, WarpR=[WarpStart,WarpEnd],Distance=Distance )

    Template.clear()


    # and cleanup
    os.system(f"rm -f {cfg.general.main_directory}{name}/ModelInput.def")
    os.system(f"rm -f {cfg.general.main_directory}{name}/Input.fits")
    os.system(f"mv {cfg.general.main_directory}{name}/tirific.def {cfg.general.main_directory}{name}/ModelInput.def")
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

def plot_RC(set_done,Mass,Rad,Vrot,colors,max_rad,sub_ring,ax):
    '''dd the RC to the overview plot and return updated tracker'''
    if set_done[0] == 1024:
        set_done= [Mass]
        labelfont= {'family':'Times New Roman',
                 'weight':'normal',
                 'size':22}
        plt.rc('font',**labelfont)
        c = next(colors)
        plt.figure(59, figsize=(8, 8), dpi=300, facecolor='w', edgecolor='k')
        ax = plt.subplot(1, 1, 1)
        plt.plot(Rad,Vrot,c=c)
        #plt.plot(Rad,Vrot,'ko',label='M$_{\odot}$ = {:.1e}'.format(Current_Galaxy.Mass))
        plt.plot(Rad,Vrot,'o',c=c,label='M$_{\odot}$ = '+' {:.1e}'.format(Mass))
        plt.ylabel('V$_{rot}$ (km s$^{-1}$)',**labelfont)
        plt.xlabel('Radius (kpc)',**labelfont)
        ax.yaxis.set_minor_locator(AutoMinorLocator(4))
        ax.xaxis.set_minor_locator(AutoMinorLocator(4))
        for axis in ['top','bottom','left','right']:
            ax.spines[axis].set_linewidth(1.5)
        plt.tick_params(axis='both', which='minor', bottom=True,left=True,length=3)
        plt.tick_params(axis='both', which='major', labelsize=17, length=6)
        plt.tick_params(axis='both', which='both', direction = 'in', width=1.5 , bottom=True,left=True ,right =True, top=True)
        max_rad = np.max(Rad)+1

    else:
        set_done.append(Mass)
        plt.figure(59)
        c = next(colors)
        if np.max(Rad)+1 > max_rad:
            max_rad =  np.max(Rad)+1
        plt.plot(Rad,Vrot,c=c)
        #plt.plot(Rad,Vrot,'o',label=r'M$_{\odot}$ = {:.1e}'.format(Current_Galaxy.Mass),c=c)
        plt.plot(Rad,Vrot,'o',label='M$_{\odot}$ = '+' {:.1e}'.format(Mass),c=c)

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
    name=f"Mass{Current_Galaxy.Mass:.1e}-i{Current_Galaxy.Inclination}d{Current_Galaxy.Dispersion[0]}-{Current_Galaxy.Dispersion[1]}"
    name=f"{name}pa{Current_Galaxy.PA}w{Current_Galaxy.Warp[0]}-{Current_Galaxy.Warp[0]}-"
    name=f"{name}{Current_Galaxy.Flare}-ba{Current_Galaxy.Beams}SNR{Current_Galaxy.SNR}"
    name=f"{name}bm{Current_Galaxy.Res_Beam[0]}-{Current_Galaxy.Res_Beam[1]}ch{Current_Galaxy.Channelwidth}"
    name=f"{name}-{Current_Galaxy.Arms}-{Current_Galaxy.Bar}-rm{Current_Galaxy.Radial_Motions}"

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
   sort_results(All_Galaxies,Result_Galaxies,results)

PURPOSE:
    sort the results into the final output

CATEGORY:
    agc

INPUTS:
    All_Galaxies= a list with the name of all galaxies
    Result_Galaxies = Part of the results to arrange
    results= final sorted array where to plant Result_Galaxies

OPTIONAL INPUTS:

OUTPUTS:
    An organized results list

OPTIONAL OUTPUTS:

PROCEDURES CALLED:
   Unspecified

NOTE:
'''
