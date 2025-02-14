#!/usr/local/bin/ python3
# This module contains a set of functions and classes that are used in several different Python scripts in the Database.

from astropy.utils.data import download_file, clear_download_cache
from astropy.io import fits
from collections import OrderedDict  # used in Proper_Dictionary
from pyHIARD import Templates as templates
from pyHIARD.AGC.base_galaxies import Base_Galaxy
from pyHIARD.Resources import Cubes as cubes
from scipy.ndimage import gaussian_filter, rotate, zoom,map_coordinates

import copy  # Used in columndensities
import numpy as np  # Used in convertskyangle and columndensity and
import os
import pyHIARD
import re
import resource
import signal
import subprocess
import sys
import time
import traceback
import warnings

with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    import matplotlib
    matplotlib.use('pdf')
    import matplotlib.pyplot as plt
    import matplotlib.font_manager as mpl_fm

try:
    import importlib.resources as import_res
except ImportError:
    import importlib_resources as import_res
 


class SofiaFaintError(Exception):
    pass
class InputError(Exception):
    pass
class CasaInstallError(Exception):
    pass

class SofiaRunError(Exception):
    pass
# A class of ordered dictionary where keys can be inserted in at specified locations or at the end.


class Proper_Dictionary(OrderedDict):
    def insert(self, existing_key, new_key, key_value):
        done = False
        if new_key in self:
            self[new_key] = key_value
            done = True
        else:
            new_orderded_dict = self.__class__()
            for key, value in self.items():
                new_orderded_dict[key] = value
                if key == existing_key:
                    new_orderded_dict[new_key] = key_value
                    done = True
            if not done:
                new_orderded_dict[new_key] = key_value
                done = True
                print(
                    "----!!!!!!!! YOUR new key was appended at the end as you provided a non-existing key to add it after!!!!!!---------")
            self.clear()
            self.update(new_orderded_dict)

        if not done:
            print("----!!!!!!!!We were unable to add your key!!!!!!---------")
#Function to convert column densities
def ask_main_directory(current_directory,input_directory = 'NonE'):
    while not os.path.isdir(input_directory):
        if not input_directory == 'NonE':
                print(f'''The directory {input_directory} does not exist.
Please provide the correct directory''')
        input_directory =  input(\
            f'''Please provide the directory where to create the database.
(default = {current_directory}):''')
        if input_directory == '':
            input_directory = current_directory
    return input_directory

def check_main_directory(cfg):
    current_directory = os.getcwd()
    if cfg.general.main_directory != current_directory:
        correct_directory = get_bool(f'''!!! Your creation directory is not directory from which you started pyHIARD.
Are you sure you want to create the database in:
{cfg.general.main_directory}?   !!!!!!
Please type yes or no (default = yes) ''',default = True)
        if not correct_directory:
            cfg.general.main_directory = ask_main_directory(current_directory)
   
    #Check the main directory exists
    while not os.path.isdir(cfg.general.main_directory):
        cfg.general.main_directory = ask_main_directory(current_directory,\
                                        input_directory=cfg.general.main_directory )
        #if we want the full database default we check that the user wants this
    if cfg.general.main_directory[-1] != '/':
        cfg.general.main_directory = f"{cfg.general.main_directory}/"
    return cfg

def create_directory(directory, base_directory, debug=False):
    split_directory = [x for x in directory.split('/') if x]
    split_directory_clean = [x for x in directory.split('/') if x]
    split_base = [x for x in base_directory.split('/') if x]
    #First remove the base from the directory but only if the first directories are the same
    if split_directory[0] == split_base[0]:
        for dirs, dirs2 in zip(split_base, split_directory):
            if dirs == dirs2:
                split_directory_clean.remove(dirs2)
            else:
                if dirs != split_base[-1]:
                    raise InputError(
                        f"You are not arranging the directory input properly ({directory},{base_directory}).")
    for new_dir in split_directory_clean:
        if not os.path.isdir(f"{base_directory}/{new_dir}"):
            os.mkdir(f"{base_directory}/{new_dir}")
        base_directory = f"{base_directory}/{new_dir}"


create_directory.__doc__ = f'''
 NAME:
    create_directory

 PURPOSE:
    create a directory recursively if it does not exists and strip leading directories when the same fro the base directory and directory to create

 CATEGORY:
    support_functions

 INPUTS:
    directory = string with directory to be created
    base_directory = string with directory that exists and from where to start the check from

 OPTIONAL INPUTS:
    debug = False

 OUTPUTS:

 OPTIONAL OUTPUTS:
    The requested directory is created but only if it does not yet exist

 PROCEDURES CALLED:
    Unspecified

 NOTE:
'''


def check_input(cfg):
    cfg = check_main_directory(cfg)

    # if we only have a single cpu turn multiprocessing of
    if cfg.general.ncpu == 1:
        cfg.general.multiprocessing = False

    if cfg.agc.enable:
        cfg.general.tirific = find_program(cfg.general.tirific, 'TiRiFiC')

        #    cfg.general.casa = find_program(cfg.general.casa, 'CASA')
    if cfg.roc.enable:
        cfg.general.sofia2 = find_program(cfg.general.sofia2, 'SoFiA2')

    if cfg.agc.enable:

        #sets = 5  # This is the amount of base galaxies we want, i.e. the number of rotation curves
        question_base = False
        for gals in cfg.agc.base_galaxies:
            if not 1 <= gals <= 7:
                question_base = True
        if question_base:
            print(f'''Please select the base types you would like to use''')
            for i in range(1, 6):
                print(
                    f'''{i}) Galaxy {i} has the following Base parameters to vary on.''')
                print_base_galaxy(Base_Galaxy(i))

            vals = input(
                f"Or you can construct you own by selecting 7, default = 7: ")
            if vals == '':
                cfg.agc.base_galaxies = [7]
            else:
                cfg.agc.base_galaxies = [int(x) for x in re.split(
                    "\s+|\s*,\s*|\s+$", vals.strip()) if 1 <= int(x) <= 6]
                if len(cfg.agc.base_galaxies) == 0:
                    cfg.agc.base_galaxies = [7]
        while cfg.agc.corruption_method.lower() not in ['casa_sim', 'gaussian','no_corrupt', 'casa_5','tres']:
            cfg.agc.corruption_method = input(
                'Your method of corruption is not acceptable please choose from Casa_Sim, Gaussian, No_Corrupt, Tres, Casa_5 (Default = Gaussian):')
            if cfg.agc.corruption_method == '':
                cfg.agc.corruption_method = 'Gaussian'

        if cfg.agc.corruption_method.lower() == 'casa_sim' or cfg.agc.corruption_method.lower() == "c_s":
            cfg.agc.corruption_method = 'Casa_Sim'
            print(
                "You are using the casa corruption method please make sure python can access casa.")
        elif cfg.agc.corruption_method .lower() == 'gaussian' or cfg.agc.corruption_method == "" or cfg.agc.corruption_method .lower() == 'g':
            cfg.agc.corruption_method = 'Gaussian'
        elif cfg.agc.corruption_method .lower() == 'casa_5' or cfg.agc.corruption_method .lower() == "c_5":
            cfg.agc.corruption_method = 'Casa_5'
            print(
                "You are using the casa corruption method please make sure python can access casa.")
        channel_options = ['independent', 'sinusoidal', 'hanning']
        while cfg.agc.channel_dependency.lower() not in channel_options:
            cfg.agc.channel_dependency = input(f''' {cfg.agc.channel_dependency} is not a valid channel dependency
Please choose from {','.join([x for x in channel_options])}:''')

        if cfg.agc.corruption_method.lower() in ['casa_sim', 'casa_5', 'tres']:
            try:
                from casaconfig import AutoUpdatesNotAllowed
                import casatasks
            except ModuleNotFoundError:
                if pyHIARD.__casa_max__ != 'OK':
                    raise CasaInstallError(f'''Your modular casa is not installed.
Most likely because it is not available for python version {sys.version}
As such you can not run the corrupt casa method''')
                else:
                    raise CasaInstallError(f'''Your modular casa is not installed.
We do not know the reason for this but it means you cannot run the casa corrupt method.
''')
            except AutoUpdatesNotAllowed:
                raise CasaInstallError(f'''Since Casa 6.6.6 Casa wants to update its data in ~/.casa/data
that directory does not exist. Please create it or see the Casa website for other options.
Alternatively you can run pyHIARD agc.corruption_method = No_Corrupt or Gauss
''')

              

        #question_variations = False
        changes_poss = ['Inclination', 'PA', 'Beams', 'Radial_Motions', 'Flare', 'Arms',
            'Dispersion', 'Bar', 'Channelwidth', 'SNR', 'Warp', 'Mass', 'Beam_Size', 'Base']
        changes_poss_lower = [x.lower() for x in changes_poss]
        for i, variables in enumerate(cfg.agc.variables_to_vary):
            while variables.lower() not in changes_poss_lower:
                variables = input(f'''{variables} is not a legitimate variable to variations.
Please choose from {','.join([x for x in changes_poss])}
replace with:''')

            cfg.agc.variables_to_vary[i] = changes_poss[changes_poss_lower.index(
                variables.lower())]
        #We always make the base
        #if 'Base' not in cfg.agc.variables_to_vary:
        #    cfg.agc.variables_to_vary.append('Base')

        if 'inclination' in [x.lower() for x in cfg.agc.variables_to_vary]:
            for i, incs in enumerate(cfg.agc.inclination):
                while not 0. <= incs <= 90.:
                    incs = float(
                        input(f'please choose a value between 0.- 90.:'))
                cfg.agc.inclination[i] = incs
        if 'pa' in [x.lower() for x in cfg.agc.variables_to_vary]:
            for i, incs in enumerate(cfg.agc.pa):
                while not 0. <= incs <= 360.:
                    incs = float(
                        input(f'please choose a value between 0.- 360.:'))
                cfg.agc.pa[i] = incs
        #the double lists to be proper, OmegaConf will throw a error when list is not list
        list_to_check = ['warp','dispersion','beam_size']
        elements = [2,2,3]
        for var_to_check in list_to_check:
            if var_to_check in [x.lower() for x in cfg.agc.variables_to_vary]:
                current =getattr(cfg.agc,var_to_check)
                try:
                    if len(current[0]) == elements[list_to_check.index(var_to_check)]:
                        pass
                except TypeError:
                    if len(current) == elements[list_to_check.index(var_to_check)]:
                        #We assume that the user wanted a single variation
                        setattr(cfg.agc,var_to_check,[[x for x in current]])
                    else:
                        raise  InputError(f'''You wanted variations in {var_to_check} which should have {elements[list_to_check.index(var_to_check)]} per variation.
                        Your list has {len(current)} elements and we do not know what to do with that. Note that this is a double list, i.e. [[]]''')


    if cfg.roc.enable:

        path_to_resources = os.path.dirname(os.path.abspath(cubes.__file__))
        allowed_galaxies = [name for name in os.listdir(
            path_to_resources) if os.path.isdir(os.path.join(path_to_resources, name))]
        allowed_galaxies.remove('__pycache__')
        allowed_galaxies_low = [x.lower() for x in allowed_galaxies]
        for i, galaxy in enumerate(cfg.roc.base_galaxies):
            while galaxy.lower() not in allowed_galaxies_low:
                galaxy = input(f'''{galaxy} is not a legitimate galaxy to vary.
Please choose from {','.join([x for x in allowed_galaxies])}
replace with:''')
            cfg.roc.base_galaxies[i] = allowed_galaxies[allowed_galaxies_low.index(
                galaxy.lower())]
        changes_poss = ['Beams', 'SNR']
        changes_poss_lower = [x.lower() for x in changes_poss]
        for i, variables in enumerate(cfg.roc.variables_to_vary):
            while variables.lower() not in changes_poss_lower:
                variables = input(f'''{variables} is not a legitimate variable to variations.
Please choose from {','.join([x for x in changes_poss])}
replace with:''')

            cfg.roc.variables_to_vary[i] = changes_poss[changes_poss_lower.index(
                variables.lower())]
        if cfg.roc.minimum_degradation_factor < 1.1:
            print(f'''Your minimum degration factor is less than 1.1 ({cfg.roc.minimum_degradation_factor}).
That wil create a mess. We will abort this run.''')
            raise InputError(f'''Your minimum degradation factor is too small''')


        if cfg.roc.minimum_degradation_factor > cfg.roc.max_degradation_factor:
            print(f'''Your max degradation factor is less than the mimimum ({cfg.roc.max_degradation_factor} vs { cfg.roc.minimum_degradation_factor}).
We have equalized them for you. Your welcome.''')
            cfg.roc.max_degradation_factor = cfg.roc.minimum_degradation_factor
    if cfg.roc.enable and cfg.roc.delete_existing:
        to_delete = f'''rm -R {' '.join([f"{cfg.general.main_directory}{x}_*Beams_*SNR" for x in cfg.roc.base_galaxies])}'''
        print("All previous models of the requested base galaxies will be removed prior to the build. \n")

        print(f"The command {to_delete} will be run.")
        cfg.roc.delete_existing = get_bool(
            "Are you sure you want to do this? (Yes/No, default=No): ", default=False)

    return cfg


check_input.__doc__ = f'''
 NAME:
    check_input

 PURPOSE:
    Check the configuration input and that all required input is present.
    In case of missing input request the input from the user.

 CATEGORY:
    common_functions

 INPUTS:
    cfg = omega conf configuration

 OPTIONAL INPUTS:

 OUTPUTS:
    a modified cfg

 OPTIONAL OUTPUTS:

 PROCEDURES CALLED:
    Unspecified

 NOTE:
'''


def unused_columndensity(levels, systemic=100., beam=None, channel_width=1., \
                column=False, arcsquare=False, solar_mass=False):
    if beam is None:
        beam=[1., 1.]
    #set solar_mass to indicate the output should be M_solar/pc**2 or if column = True the input is
    f0 = 1.420405751786E9  # Hz rest freq
    c = 299792.458  # light speed in km / s
    pc = 3.086e+18  # parsec in cm
    solarmass = 1.98855e30  # Solar mass in kg
    mHI = 1.6737236e-27  # neutral hydrogen mass in kg

    if systemic > 10000:
        systemic = systemic/1000.
    f = f0 * (1 - (systemic / c))  # Systemic frequency
    if arcsquare:
        HIconv = 605.7383 * 1.823E18 * (2. * np.pi / (np.log(256.)))
        if column:
            # If the input is in solarmass we want to convert back to column densities
            if solar_mass:
                levels = levels*solarmass/(mHI*pc**2)
            #levels=levels/(HIconv*channel_width)
            levels = levels/(HIconv*channel_width)
        else:

            levels = HIconv*levels*channel_width
            if solar_mass:
                levels = levels/solarmass*(mHI*pc**2)
    else:
        if beam.size < 2:
            beam = [beam, beam]
        b = beam[0]*beam[1]
        if column:
            if solar_mass:
                levels = levels*solarmass/(mHI*pc**2)
            TK = levels/(1.823e18*channel_width)
            levels = TK/(((605.7383)/(b))*(f0/f)**2)
        else:
            TK = ((605.7383)/(b))*(f0/f)**2*levels
            levels = TK*(1.823e18*channel_width)
    if ~column and solar_mass:
        levels = levels*mHI*pc**2/solarmass
    return levels
        # a Function to convert the RA and DEC into hour angle (invert = False) and vice versa (default)


def convertRADEC(RA, DEC, invert=False, colon=False):

    if not invert:
        try:
            _ = (e for e in RA)
        except TypeError:
            RA = [RA]
            DEC = [DEC]
        for i in range(len(RA)):
            xpos = RA
            ypos = DEC
            xposh = int(np.floor((xpos[i]/360.)*24.))
            xposm = int(np.floor((((xpos[i]/360.)*24.)-xposh)*60.))
            xposs = (((((xpos[i]/360.)*24.)-xposh)*60.)-xposm)*60
            yposh = int(np.floor(np.absolute(ypos[i]*1.)))
            yposm = int(np.floor((((np.absolute(ypos[i]*1.))-yposh)*60.)))
            yposs = (((((np.absolute(ypos[i]*1.))-yposh)*60.)-yposm)*60)
            sign = ypos[i]/np.absolute(ypos[i])
            if colon:
                RA[i] = "{}:{}:{:2.2f}".format(xposh, xposm, xposs)
                DEC[i] = "{}:{}:{:2.2f}".format(yposh, yposm, yposs)
            else:
                RA[i] = "{}h{}m{:2.2f}".format(xposh, xposm, xposs)
                DEC[i] = "{}d{}m{:2.2f}".format(yposh, yposm, yposs)
            if sign < 0.: DEC[i] = '-'+DEC[i]
        if len(RA) == 1:
            RA = str(RA[0])
            DEC = str(DEC[0])
    else:
        if isinstance(RA, str):
            RA = [RA]
            DEC = [DEC]

        xpos = RA
        ypos = DEC

        for i in range(len(RA)):
            # first we split the numbers out
            tmp = re.split(r"[a-z,:]+", xpos[i])
            RA[i] = (float(tmp[0])+((float(tmp[1])+(float(tmp[2])/60.))/60.))*15.
            tmp = re.split(r"[a-z,:'\"]+", ypos[i])
            DEC[i] = float(np.absolute(float(tmp[0]))+((float(tmp[1])
                           + (float(tmp[2])/60.))/60.))*float(tmp[0])/np.absolute(float(tmp[0]))

        if len(RA) == 1:
            RA = float(RA[0])
            DEC = float(DEC[0])
    return RA, DEC


# function for converting kpc to arcsec and vice versa

def convertskyangle(angle, distance=1., unit='arcsec', distance_unit='Mpc', \
                        physical=False):
    try:
        _ = (e for e in angle)
    except TypeError:
        angle = [angle]

        # if physical is true default unit is kpc
    angle = np.array(angle)
    if physical and unit == 'arcsec':
        unit = 'kpc'
    if distance_unit.lower() == 'mpc':
        distance = distance * 10 ** 3
    elif distance_unit.lower() == 'kpc':
        distance = distance
    elif distance_unit.lower() == 'pc':
        distance = distance / (10 ** 3)
    else:
        print('CONVERTSKYANGLE: ' + distance_unit
              + ' is an unknown unit to convertskyangle.\n')
        print('CONVERTSKYANGLE: please use Mpc, kpc or pc.\n')
        sys.exit()
    if not physical:
        if unit.lower() == 'arcsec':
            radians = (angle / 3600.) * ((2. * np.pi) / 360.)
        elif unit.lower() == 'arcmin':
            radians = (angle / 60.) * ((2. * np.pi) / 360.)
        elif unit.lower() == 'degree':
            radians = angle * ((2. * np.pi) / 360.)
        else:
            print('CONVERTSKYANGLE: ' + unit
                  + ' is an unknown unit to convertskyangle.\n')
            print('CONVERTSKYANGLE: please use arcsec, arcmin or degree.\n')
            sys.exit()

        kpc = 2. * (distance * np.tan(radians / 2.))
    else:
        if unit.lower() == 'kpc':
            kpc = angle
        elif unit.lower() == 'mpc':
            kpc = angle * (10 ** 3)
        elif unit.lower() == 'pc':
            kpc = angle / (10 ** 3)
        else:
            print('CONVERTSKYANGLE: ' + unit
                  + ' is an unknown unit to convertskyangle.\n')
            print('CONVERTSKYANGLE: please use kpc, Mpc or pc.\n')
            sys.exit()
        radians = 2. * np.arctan(kpc / (2. * distance))
        kpc = (radians * (360. / (2. * np.pi))) * 3600.
    if len(kpc) == 1:
        kpc = float(kpc[0])
    return kpc

def select_emission(data,hdr,name,work_dir,sofia_call='sofia_call'):
    #get a mask
    Mask_Data = create_masks(data,hdr, work_dir, name, sofia_call=sofia_call)
    try:
        Vel_Units = hdr['CUNIT3'].lower()
    except KeyError:
        if hdr['CDELT3'] > 150.:
            hdr.set("CUNIT3", 'M/S', before="CTYPE3")
        else:
            hdr.set("CUNIT3", 'KM/S', before="CTYPE3")
    if hdr['CUNIT3'].lower() == 'm/s' or hdr['CDELT3'] > 150.:
        hdr['CDELT3'] = hdr['CDELT3'] / 1000.
        hdr['CUNIT3'] = 'km/s'
        hdr['CRVAL3'] = hdr['CRVAL3'] / 1000.
    # ensure we have a BPA
    if 'BPA' not in hdr:
        hdr['BPA'] = 0.
    #Apply the mask to Template
    data[Mask_Data <= 0.] = 0.
    data[np.isnan(data)] = 0.
    return data,hdr

def create_masks(data_in,hdr_in, working_dir, name, sofia_call='sofia2'):
    data = copy.deepcopy(data_in)
    hdr = copy.deepcopy(hdr_in)
    # First we smooth our template
    # We smooth this to 1.25 the input beam
    '''
    bmaj = hdr["BMAJ"]*3600.
    bmin = hdr["BMIN"]*3600.
    FWHM_conv_maj = np.sqrt((1.25 * bmaj) ** 2 - bmaj ** 2)
    FWHM_conv_min = np.sqrt((1.25 * bmin) ** 2 - bmin ** 2)
    # and in terms of pixels the sigmas
    sig_maj = (FWHM_conv_maj / np.sqrt(8 * np.log(2))) / \
               abs(hdr["CDELT1"] * 3600.)
    sig_min = (FWHM_conv_min / np.sqrt(8 * np.log(2))) / \
               abs(hdr["CDELT2"] * 3600.)
    '''
    #We replace zeros with NAN
    data[data == 0.] = float('NaN')
    # It seems that I do not smooth the mask anymore, why not?
    Tmp_Cube=data
    #Tmp_Cube = gaussian_filter(data, sigma=(0, sig_min, sig_maj), order=0)
    # Replace 0. with Nan

    #write this to the fits file
    fits.writeto(f'{working_dir}/tmp_{name}.fits',
                 Tmp_Cube, hdr, overwrite=True)
    SoFiA_Template = read_template_file('Sofia_Template.par')
    SoFiA_Template['input.data'.upper(
    )] = f'input.data = {working_dir}/tmp_{name}.fits'
    SoFiA_Template['scfind.threshold'.upper()] = 'scfind.threshold	= 3'
    SoFiA_Template['linker.minSizeZ'.upper(
    )] = f'linker.minSizeZ = {int(Tmp_Cube.shape[0]/2.)}'
    SoFiA_Template['output.filename'.upper()] = f'output.filename = tmp_{name}'
    with open(f'{working_dir}/tmp_{name}_sof.par', 'w') as file:
        file.writelines([SoFiA_Template[key] + "\n" for key in SoFiA_Template])
    run_sofia(working_dir, f'tmp_{name}_sof.par', sofia_call=sofia_call)
    #wait half a second to make sure got written
    #time.sleep(0.5)
    Mask_Outer = fits.open(f'{working_dir}/tmp_{name}_mask.fits')
    Mask_Outer[0].data[Tmp_Cube < 0.] = 0.
    #Mask_Outer[0].header['CUNIT3'] = 'KM/S'
    try:
        del Mask_Outer[0].header['HISTORY']
    except:
        pass
    #fits.writeto(f'{outdir}/Outer_{name}_mask.fits',Mask_Outer[0].data,Mask_Outer[0].header,overwrite = True)
    os.system(f"rm -f {working_dir}/tmp_{name}_mask.fits")
    # Create the inner mask
    SoFiA_Template['dilation.enable'.upper()] = 'dilation.enable	=	false'
    SoFiA_Template['scfind.threshold'.upper()] = 'scfind.threshold	= 5'
    with open(f'{working_dir}/tmp_{name}_sof.par', 'w') as file:
        file.writelines([SoFiA_Template[key] + "\n" for key in SoFiA_Template])
    run_sofia(working_dir, f'tmp_{name}_sof.par', sofia_call=sofia_call)
    #wait half a second to make sure got written
    time.sleep(0.5)
    Mask_Inner = fits.open(f'{working_dir}/tmp_{name}_mask.fits')
    #Mask_Inner[0].header['CUNIT3'] = 'KM/S'
    del Mask_Inner[0].header['HISTORY']
    #fits.writeto(f'{outdir}/Inner_{name}_mask.fits',Mask_Inner[0].data,Mask_Inner[0].header,overwrite = True)
    os.system(f"rm -f {working_dir}/tmp_{name}_mask.fits")
    os.system(f"rm -f {working_dir}/tmp_{name}.fits")
    os.system(f'rm -f {working_dir}/tmp_{name}_sof.par')
    transform = np.array(copy.deepcopy(Mask_Outer[0].data), dtype=float)
    transform[transform > 0.] = 1.
    tmp = gaussian_filter(transform, sigma=(0, 3, 3), order=0)
    transform = copy.deepcopy(tmp)
    tmp = []
    transform[Mask_Inner[0].data > 0.] = 1.
    transform = gaussian_filter(transform, sigma=(0, 1, 1), order=0)

    transform[transform < 0.1] = 0.
    del data
    del hdr
    Mask_Outer.close()
    Mask_Inner.close()

    return transform


create_masks.__doc__ = f'''
 NAME:
    create_masks

 PURPOSE:
    Create the masks for the ROC standards.

 CATEGORY:
    common_functions

 INPUTS:
    outdir = directory where to put the final masks for future use
    working_dir = directory where to run sofia
    name = Base name of the input and output files

 OPTIONAL INPUTS:
    sofia_call = command name to run sofia

 OUTPUTS:
    the cut cube is returned.

 OPTIONAL OUTPUTS:

 PROCEDURES CALLED:
    Unspecified

 NOTE:
'''


def calculate_pixel_noise(requested_noise,smoothing_sigma,tolerance=0.025):
    #firt guess of the pixel noise
    pixel_noise = requested_noise *np.mean(smoothing_sigma)*2. * np.sqrt(np.pi)
    #  we need a shape for testing. something say 100 * smoothing kernel
    test_shape= np.array([int(x*100.) for x in smoothing_sigma],dtype=int)
    rng = np.random.default_rng()
    achieved_noise = 0.
    while abs(achieved_noise - requested_noise) / requested_noise > 0.025:
        #If not the first run update our pixel noise
        if achieved_noise != 0:
            pixel_noise = pixel_noise*requested_noise/achieved_noise
        #We only need to estimate this in a single channel there is no spectral component
        #fill with gaussian values
        test_noise = rng.normal(scale=pixel_noise, size=test_shape)
        # and smooth to the final beam
        test_noise_smoothed = gaussian_filter(test_noise,sigma=smoothing_sigma,\
                                order=0)
        achieved_noise = np.std(test_noise_smoothed)
        #print(f"The current pixel noise estimate leads to {achieved_noise} mJy/beam (Requested = {requested_noise} mJy/beam).")



    return achieved_noise,pixel_noise
calculate_pixel_noise.__doc__ = f'''
 NAME:
    calculate_pixel_noise

 PURPOSE:
    calculate the noise that is required as input to random value generator to end up with the required noise after smoothing

 CATEGORY:
    common_functions

 INPUTS:
    requested_noise = the finale requested noise in smoothed images
    smoothing_sigma = the sigma of the final beam in pixels


 OPTIONAL INPUTS:
    tolerance = the accepted difference between the requested noise and the achieved noise

 OUTPUTS:
    achieved_noise = the smoothed standard deviation
    pixel_noise = the input pixel standard deviation

 OPTIONAL OUTPUTS:

 PROCEDURES CALLED:
    Unspecified

 NOTE:
'''



def cut_input_cube(file_in, sizes, debug=False):
    cube = fits.open(file_in)
    #, uint=False,do_not_scale_image_data=True, ignore_blank=True)
    hdr = cube[0].header
    for i, pos in enumerate(sizes):
        if pos[1] == -1:
            sizes[i][1] = hdr[f'NAXIS{3-i}']-1

    if hdr['NAXIS'] > 3:
        hdr['NAXIS'] = 3
        del hdr['*4']
        new_data = cube[0].data[0, sizes[0][0]:sizes[0][1],
            sizes[1][0]:sizes[1][1], sizes[2][0]:sizes[2][1]]
    else:
        new_data = cube[0].data[sizes[0][0]:sizes[0][1],
            sizes[1][0]:sizes[1][1], sizes[2][0]:sizes[2][1]]

    if f'PC01_01' in hdr:
        del hdr['PC0*_0*']
    hdr['NAXIS1'] = sizes[2][1]-sizes[2][0]
    hdr['NAXIS2'] = sizes[1][1]-sizes[1][0]
    hdr['NAXIS3'] = sizes[0][1]-sizes[0][0]

    hdr['CRPIX1'] = hdr['CRPIX1']-sizes[2][0]
    hdr['CRPIX2'] = hdr['CRPIX2']-sizes[1][0]
    hdr['CRPIX3'] = hdr['CRPIX3']-sizes[0][0]
    try:
        del hdr['HISTORY']
    except:
        pass
    if hdr['BITPIX'] < -32:
        hdr['BITPIX'] = -32
        #print(hdr['BITPIX'])
        #exit()
    elif hdr['BITPIX'] > -32:
        hdr['BITPIX'] = 32
    # Do not regid the original cubes if we do this we do this then we mess up the barollo and rc files
    cube[0].header=hdr
    cube[0].data=np.array(new_data, dtype = np.float32)
    return cube
cut_input_cube.__doc__=f'''
 NAME:
    cut_input_cube

 PURPOSE:
    Cut filename back to the size of subcube, update the header and write back to disk.

 CATEGORY:
    common_functions

 INPUTS:
    file_in = the fits file to cut
    sizes = array that contains the new size as
                [[z_min,z_max],[y_min,y_max], [x_min,x_max]]
                adhering to fits' idiotic way of reading fits files.


 OPTIONAL INPUTS:
    debug = False

 OUTPUTS:
    the cut cube is returned.

 OPTIONAL OUTPUTS:

 PROCEDURES CALLED:
    Unspecified

 NOTE:
'''

def delete_directory(dir_to_delete):
    for f in os.listdir(dir_to_delete):
        try:
            os.remove(os.path.join(dir_to_delete, f))
        except FileNotFoundError:
            pass
        except IsADirectoryError:
            delete_directory(os.path.join(dir_to_delete, f))

    try:
        os.rmdir(dir_to_delete)
    except FileNotFoundError:
        pass
delete_directory.__doc__=f'''
 NAME:
    delete_directory

 PURPOSE:
    delete a non empty directory

 CATEGORY:
    common_functions

 INPUTS:
    dir = directory to delete

 OPTIONAL INPUTS:

 OUTPUTS:

 OPTIONAL OUTPUTS:

 PROCEDURES CALLED:
    Unspecified

 NOTE:
'''


def download_cube(name, url, sizes, new_location):
    name_in=download_file(url, pkgname = 'pyHIARD', timeout = 1800)
    cube=cut_input_cube(name_in, sizes)
    fits.writeto(f'{new_location}/{name}.fits', cube[0].data, cube[0].header)
    clear_download_cache(url, pkgname = 'pyHIARD')
    return cube
download_cube.__doc__=f'''
 NAME:
    download_cube

 PURPOSE:
    Download the specified cube when not present in installation and cut down required size
    Store in installation for future use.

 CATEGORY:
    common_functions

 INPUTS:
    name = local name of the cube
    url = url of location of remote cube
    sizes = require sizes for stored cube
    new_location = location of correct directory in installation.


 OPTIONAL INPUTS:
    debug = False

 OUTPUTS:
    the cut cube is returned

 OPTIONAL OUTPUTS:

 PROCEDURES CALLED:
    Unspecified

 NOTE:
'''

def find_program(name, search):
    found=False
    while not found:
        try:
            run=subprocess.Popen(
                [name], stdout = subprocess.PIPE, stderr = subprocess.PIPE)
            run.stdout.close()
            run.stderr.close()
            os.kill(run.pid, signal.SIGKILL)
            found=True
        except Exception as e:
            traceback.print_tb(e.__traceback__)
            name=input(f'''You have indicated to use {name} for using {search} but it cannot be found.
Please provide the correct name : ''')
    return name
find_program.__doc__=f'''
 NAME:
    find_program

 PURPOSE:
    check whether a program is available for use.

 CATEGORY:
    support_functions

 INPUTS:
    name = command name of the program to run
    search = Program we are looking for
 OPTIONAL INPUTS:

 OUTPUTS:
    the correct command for running the program

 OPTIONAL OUTPUTS:

 PROCEDURES CALLED:
    Unspecified

 NOTE:
'''


def get_beam_area_in_pixels(Template_Header, beam= None):

    if beam is None:
        beam = [Template_Header["BMAJ"],Template_Header["BMIN"]]
    #  https://science.nrao.edu/facilities/vla/proposing/TBconv
    beamarea=(np.pi*abs((beam[0]*beam[1])))/(4.*np.log(2.))
    return beamarea/(abs(Template_Header['CDELT1'])*abs(Template_Header['CDELT2']))
get_beam_area_in_pixels.__doc__=f'''
NAME:
    get_beam_in_pixels(Template_Header, beam= [-1,-1.])

PURPOSE:
    calculate the area of the beam in pixels

CATEGORY:
    common_functions

INPUTS:
    name = basename of template
    beams = the requested beam in degrees, if unset it assumed that the beam is in the header.

OPTIONAL INPUTS:

OUTPUTS:
    the beam area in pixels

OPTIONAL OUTPUTS:

PROCEDURES CALLED:
   Unspecified

NOTE:
'''
# Function to input a boolean answer
def get_bool(print_str = "Please type True or False", default = True):
    invalid_input=True
    while invalid_input:
        inp=input(print_str)
        if inp == "":
            if default:
                return True
            else:
                return False
        elif inp.lower() == "true" or inp.lower() == "t" or inp.lower() == "y" or inp.lower() == "yes":
            return True
        elif inp.lower() == "false" or inp.lower() == "f" or inp.lower() == "n" or inp.lower() == "no":
            return False
        else:
            print("Error: the answer must be true/false or yes/no.")
get_bool.__doc__=f'''
NAME:
    get_bool

PURPOSE:
    get a versatile boolean response from the user

CATEGORY:
    common_functions

INPUTS:

OPTIONAL INPUTS:
    print_str = "Please type True or False"
    The question to ask the user

    default = True
    The default answer will correspond to the default

OUTPUTS:
    The corresponding boolean

OPTIONAL OUTPUTS:

PROCEDURES CALLED:
Unspecified

NOTE:
'''

def get_created_models(Catalogue, delete_existing):
    cat_exists=os.path.isfile(Catalogue)
    if not delete_existing and cat_exists:
        with open(Catalogue) as cat:
            lines=cat.readlines()
        try:
            split_line=lines[-1].split('|')
            if len(split_line) > 0:
                number_models=int(float(split_line[0])+1)
        except:
            number_models=0
    else:
        with open(Catalogue, 'w') as cat:
            cat.write('ID|Distance|Directoryname|Cubename\n')
        number_models=0
    return number_models
get_created_models.__doc__=f'''
 NAME:
    get_created_models(Catalogue,delete_exiting)

 PURPOSE:
    check whether a new catalogue should be created if not then the starting number is the end of the catalogue

 CATEGORY:
    support_functions

 INPUTS:
    name = command name of the program to run
    search = Program we are looking for
 OPTIONAL INPUTS:

 OUTPUTS:
    the correct command for running the program

 OPTIONAL OUTPUTS:

 PROCEDURES CALLED:
    Unspecified

 NOTE:
'''

def get_mask(Cube_In, factor = 5.2):
    Mask = copy.deepcopy(Cube_In)
    Top_Cube = Cube_In[Cube_In > 0.9*np.max(Cube_In)]
    Sorted_Cube = Top_Cube[np.argsort(Top_Cube)]
    Max_Flux = np.mean(Sorted_Cube[-20:])
    Mask[Cube_In > 0.005*Max_Flux] = 1.
    Mask[Cube_In < 0.005*Max_Flux] = 0.
    # If we have a lot of pixels and stuff we need smooth more
    #So we take the inverse of the initial factor with a standard of 0.75 minimum of 0.5 and a maximum of 1.25
    pix_smooth = 3.9/factor
    if pix_smooth < 0.5:
        pix_smooth=0.5
    if pix_smooth > 3:
        pix_smooth=3.
    Mask =  gaussian_filter(Mask,sigma=(0.,pix_smooth,pix_smooth),order=0)
    Mask[Mask > 0.95] = 1
    Mask[Mask < 0.05] = 0
    Mask = np.array(Mask,dtype=int)
    return Mask

def get_mean_flux(Cube_In,Mask = None):

    #As python is really the dumbest language ever invented there is different behaviour for passing np.arrays and lists
    Cube=copy.deepcopy(Cube_In)
    # To calculate the mean of all values in the mask is too sensitive to very small variations in the mask
    if Mask is not None:
        Cube[Mask < 0.5] = 0.
    # First we need to get a maximum flux value.
    Top_Cube = Cube[Cube > 0.9*np.max(Cube)]
    Sorted_Cube = Top_Cube[np.argsort(Top_Cube)]
    Max_Flux = np.mean(Sorted_Cube[-20:])
    Mean = np.mean(Cube[Cube > 0.07*Max_Flux])
    # Regridded Our new method gives 0.03289582223850146. The old gave 0.026361946919571028
    #Our new method gives 0.03307503935880944. The old gave 0.025313454549815356
    return Mean
get_mean_flux.__doc__=f'''
 NAME:
    get_created_models(Catalogue,delete_exiting)

 PURPOSE:
    check whether a new catalogue should be created if not then the starting number is the end of the catalogue

 CATEGORY:
    support_functions

 INPUTS:
    name = command name of the program to run
    search = Program we are looking for
 OPTIONAL INPUTS:

 OUTPUTS:
    the correct command for running the program

 OPTIONAL OUTPUTS:

 PROCEDURES CALLED:
    Unspecified

 NOTE:
'''
def limit_memory(maxsize):
    soft, hard=resource.getrlimit(resource.RLIMIT_AS)
    resource.setrlimit(resource.RLIMIT_AS, (maxsize, hard))

def load_text_model(filename, type = 'Tirific', Variables = None ,\
                        package_file = True):
    if Variables is None:
        Variables = ['BMIN', 'BMAJ', 'BPA', 'RMS', 'DISTANCE', 'NUR', 'RADI', \
            'VROT','Z0', 'SBR', 'INCL', 'PA', 'XPOS', 'YPOS', 'VSYS', 'SDIS', \
            'VROT_2',  'Z0_2', 'SBR_2','INCL_2', 'PA_2', 'XPOS_2', 'YPOS_2', \
            'VSYS_2', 'SDIS_2', 'CONDISP', 'CFLUX', 'CFLUX_2']

    #First check that we have a proper type
    allowed_types=['tir', 'rc', 'bar', 'tirific', 'rotcur', 'barolo', 'fat']
    while type.lower() not in allowed_types:
        type=input(f'''pyHIARD can not deal with the type {type}.
please provide one of the following types {', '.join(allowed_types)}:''')
    # Then make sure the type adheres to a single desganation
    if type.lower() in ['tir', 'tirific']:
        type='tir'
    elif type.lower() in ['rc', 'rotcur']:
        type='rc'
    elif type.lower() in ['bar', 'barolo']:
        type='bar'
    else:
        raise InputError('This should be impossible.')

    ext={'tir': 'def', 'bar': 'txt', 'rc': 'rotcur'}
    if package_file:
        python_version = sys.version_info
        if python_version[0] < 3. or \
            (python_version[0] == 3. and python_version[1] < 9.):
            model=__import__(
                f'pyHIARD.Resources.Cubes.{filename}', globals(), locals(), filename, 0)
            with import_res.open_text(model, f'{filename}.{ext[type]}') as tmp:
                unarranged = tmp.readlines()
        else:

            from importlib import import_module
            model = import_module( f'pyHIARD.Resources.Cubes.{filename}')
            with import_res.files(model).joinpath(f'{filename}.{ext[type]}').open('r') as tmp:
                unarranged = tmp.readlines()

    else:
        with open(filename, 'r') as tmp:
            unarranged=tmp.readlines()
    Values=Proper_Dictionary({})
    variable_location={}

    for key in Variables:
        Values[key]=[]
        variable_location[key]=-1
    mapping={'rc': {'RADI': 'radius', 'VSYS': 'systemic', 'VSYS_ERR': 'systemic_error',
                      'VROT': 'rotation', 'VROT_ERR': 'rotation_error', 'VRAD': 'expansion',
                      'VRAD_ERR': 'expansion_error', 'PA': 'pos.', 'PA_ERR': 'pos._error',
                      'INCL': 'incli-', 'INCL_ERR': 'incli-_error', 'XPOS': 'x-pos.',
                      'XPOS_ERR': 'x-pos._error', 'YPOS': 'y-pos.', 'YPOS_ERR': 'y-pos._error'},
               'bar': {'RADI': 'RAD(arcs)', 'VSYS': 'VSYS(km/s)', 'VSYS_ERR': 'E_VSYS',
                      'VROT': 'VROT(km/s)', 'VROT_ERR': 'E_VROT', 'VRAD': 'VRAD(km/s)',
                      'VRAD_ERR': 'E_VRAD', 'PA': 'P.A.(deg)', 'PA_ERR': 'E_PA',
                      'INCL': 'INC(deg)', 'INCL_ERR': 'E_INC', 'XPOS': 'XPOS(pix)',
                      'XPOS_ERR': 'E_XPOS', 'YPOS': 'YPOS(pix)', 'YPOS_ERR': 'E_YPOS',
                      'SDIS': 'DISP(km/s)', 'SDIS_ERR': 'E_DISP', 'Z0': 'Z0(arcs)', 'Z0_ERR': 'E_Z0'}
                      }
    for line in unarranged:
        stripped = False
        if line[0] in ['#', '!']:
            line = line[1:]
            stripped = True
        if type == 'tir':
            variable = str(line.split('=')[0].strip().upper())
            if variable in Variables:
                Values[variable] = [float(x)
                                          for x in line.split('=')[1].rsplit()]
        elif stripped:
            split_line= [x.strip() for x in line.split()]
            if len(split_line) == 0:
                continue
            if type == 'rc':
                if split_line[0] != 'radius':
                    continue
            for key in Variables:
                if key not in mapping[type]:
                    continue
                if key[-4: ] == '_ERR':
                    if type == 'bar':
                        if f'{mapping[type][key]}1' in split_line:
                            variable_location[key]= [split_line.index(f'{mapping[type][key]}1'), split_line.index(f'{mapping[type][key]}2')]
                    elif type == 'rc':
                        var= mapping[type][key].split('_')
                        if var[0] in split_line:
                            if split_line[split_line.index(var[0])+1] == var[1]:
                                variable_location[key]= split_line.index(
                                    var[0])+1
                else:
                    if mapping[type][key] in split_line:
                        variable_location[key] = split_line.index(
                            mapping[type][key])
        else:
            split_line = [float(x.strip()) for x in line.split()]
            for var in variable_location:
                if variable_location[var] == -1:
                    pass
                elif isinstance(variable_location[var], int):
                    Values[var].append(
                        float(split_line[variable_location[var]]))
                else:
                    the_av = np.mean([abs(split_line[variable_location[var][0]]), abs(split_line[variable_location[var][1]])])
                    Values[var].append(float(the_av))

    final_values = []
    for key in Variables:
        if len(Values[key]) > 0.:
            if 'RADI' in Variables:
                while len(Values[key]) < len(Values['RADI']):
                    Values[key].append(Values[key][-1])
            final_values.append(np.array(Values[key]))
        else:
            if 'RADI' in Variables:
                final_values.append(np.full(len(Values['RADI']), -1.))
            else:
                final_values.append(np.full(len(Values[0]), -1.))
    return final_values
load_text_model.__doc__ =f'''
 NAME:
    load_text_model

 PURPOSE:
    load a text file tilted ring model

 CATEGORY:
    common_functions

 INPUTS:
    name = name of the packaged galaxy or text file
    type = type of input file, options are  ['tir','rc','bar','tirific','rotcur','barolo','fat']
    Variables = variable to extract following tirific convention
    package_file = Boolean indicating whether the file is from the python package or not

 OPTIONAL INPUTS:

 OUTPUTS:
    arrays of variables in the order of variables is put in.

 OPTIONAL OUTPUTS:

 PROCEDURES CALLED:
    Unspecified

 NOTE:
    Note that rotcur and barolo are not very intelligent and hence return the
    xpos and ypos in pixel coordinates which are worthless without a corresponding cube.
    barolo errors are averaged from lower and upper error in the file.

 '''
#

def plot_input(directory, Model,add_sbr = None, Distance= 0., RHI = None \
                ,Title = 'EMPTY',WarpR=None, font_file = 'empty.ttf'):

    if add_sbr is None:
        add_sbr = [0.,0.]
    if RHI is None:
        RHI = [0.,0.,0.]
    if WarpR is None:
        WarpR=[0.,0.]
    variables_to_plot = ['SBR', 'VROT','PA','INCL','SDIS','Z0']
    plots = len(variables_to_plot)
    units = {'SBR': 'SBR (Jy km s$^{-1}$ arcsec$^{-2}$)' ,
            'VROT': 'V$_{rot}$ (km s$^{-1}$)',
            'PA': 'PA ($^{\circ}$)',
            'INCL': 'INCL ($^{\circ}$)',
            'SDIS': 'Disp. (km s$^{-1}$)',
            'Z0': 'Z0 (arcsec)'}
    plt.figure(2, figsize=(8, 12), dpi=100, facecolor='w', edgecolor='k')
    try:
        mpl_fm.fontManager.addfont(font_file)
        font_name = mpl_fm.FontProperties(fname=font_file).get_name()
    except FileNotFoundError:
        font_name = 'DejaVu Sans'
    labelfont = {'family': font_name,
                'weight': 'normal',
                'size': 18}
    plt.rc('font', **labelfont)
    radius = np.array([float(x) for x in Model['RADI'].split('=')[1].split()], dtype=float)
    bmaj = float(Model['BMAJ'].split('=')[1])

    radius_bmaj = [0]
    counter = 1.
    for i, rad in enumerate(radius):
        if rad < counter*bmaj:
            pass
        else:
            radius_bmaj.append(i)
            counter += 1.
    radius_bmaj = np.array(radius_bmaj,dtype=int)
    rad_unit = 'arcsec'
    #if we have a distance we convert all arcsec values to kpc, It is assumed all extra variables come in the right unit
    if Distance > 0.:
        radius = convertskyangle(radius, distance=Distance)
        rad_unit = 'kpc'
        units['Z0'] = 'Z0 (kpc)'
    for i, variable in enumerate(variables_to_plot):
        plt.subplot(plots, 1, plots-i)
        var_to_plot = [float(x) for x in Model[variable].split('=')[1].split()]
        while len(var_to_plot) < len(radius):
            var_to_plot.append(var_to_plot[-1])
        var_to_plot = np.array(var_to_plot,dtype=float)
        if variable == 'Z0' and Distance != 0.:
            var_to_plot = convertskyangle(var_to_plot , distance=Distance)
        plt.plot(radius, var_to_plot, 'k')
        plt.plot(radius[radius_bmaj], var_to_plot[radius_bmaj], 'ko')
        var_to_plot2 = [float(x) for x in Model[f'{variable}_2'].split('=')[1].split()]
        while len(var_to_plot2) < len(radius):
            var_to_plot2.append(var_to_plot2[-1])
        if variable == 'Z0' and Distance != 0.:
                var_to_plot2 = convertskyangle(var_to_plot2 , distance=Distance)
        if np.sum(var_to_plot2) > 0. and np.sum([x-y for x, y in zip(var_to_plot,var_to_plot2)]) != 0.:
            var_to_plot = np.array(var_to_plot2,dtype=float)
            plt.plot(radius, var_to_plot, 'r')
            plt.plot(radius[radius_bmaj], var_to_plot[radius_bmaj], 'ro')
        if i == 0.:
            plt.xlabel(f'Radius ({rad_unit})', **labelfont)
            lab_bottom = True
        else:
            lab_bottom = False
        if variable in ['SBR'] and np.sum(add_sbr) != 0:
            plt.plot(radius, add_sbr[0, :], 'b')
            plt.plot(radius[radius_bmaj], add_sbr[0, radius_bmaj], 'bo')
            plt.plot(radius, add_sbr[1, :], 'y')
            plt.plot(radius[radius_bmaj], add_sbr[1, radius_bmaj], 'yo')

        plt.ylabel(units[variable], **labelfont)
        ymin,ymax = plt.ylim()
        plt.margins(x=0., y=0.)
        if variable in ['SBR', 'VROT','SDIS','Z0'] and RHI[0] != 0:
            plt.plot([RHI[0], RHI[0]],[ymin-(ymax-ymin)*0.1,ymax+(ymax-ymin)*0.1],'b')
            if variable in ['SBR'] and RHI[1] != 0:
                plt.plot([RHI[1], RHI[1]], [ymin - (ymax - ymin)
                         * 0.1, ymax + (ymax - ymin) * 0.1], 'b--')
            if variable in ['SBR'] and RHI[2] != 0:
                plt.plot([RHI[2], RHI[2]], [ymin - (ymax - ymin)
                         * 0.1, ymax + (ymax - ymin) * 0.1], 'b--')
        if variable in ['INCL', 'PA'] and np.sum(WarpR) != 0:
            plt.plot([WarpR[1], WarpR[1]],[ymin-(ymax-ymin)*0.1,ymax+(ymax-ymin)*0.1],'g')
            if WarpR[0] != 0.:
                plt.plot([WarpR[0], WarpR[0]],[ymin-(ymax-ymin)*0.1,ymax+(ymax-ymin)*0.1],'g--')
        plt.tick_params(
        axis = 'x',  # changes apply to the x-axis
        which = 'both',  # both major and minor ticks are affected
        direction = 'in',
        bottom = True,  # ticks along the bottom edge are off
        top = False,  # ticks along the top edge are off
        labelbottom =lab_bottom)  # labels along the bottom edge are off
    if Title != 'EMPTY':
        plt.title(Title)
    plt.savefig(f"{directory}Overview_Input.png", bbox_inches ='tight')
    plt.close()
plot_input.__doc__= f'''
NAME:
    plot_input(directory,Model,add_sbr = 0, Distance= 0.,
               RHI = [0.,0.,0.] ,Title = 'EMPTY', font_file='empty.ttf')

PURPOSE:
    Plot the tirific template as output

CATEGORY:
   roc

INPUTS:
    directory = the destination of the plot
    Model = Tirific Template
    add_sbr = additional SBR profiles to plot
    Distance = Distance to have kpc instead of arcsec
    RHI = The radius corresponding to RHI [symmetric,appr,rec]
    Warp_Start = radius corresponding to the start of the warp.
    Title = title string
    font_file = location of the true type font to be used for plotting.

OPTIONAL INPUTS:

OUTPUTS:
    standardized plot for output

OPTIONAL OUTPUTS:

PROCEDURES CALLED:
   Unspecified

NOTE:
'''
def print_base_galaxy(Galaxy):
    try:
        RAhr, DEChr= convertRADEC(Galaxy.Coord[0], Galaxy.Coord[1])
    except AttributeError:
        RAhr, DEChr= 'Random', 'Random'
    print(f'''The galaxy has the central coordinates RA= {RAhr}, DEC={DEChr}
{'Inclination':15s} = {Galaxy.Inclination:<10.1f}, {'Dispersion':15s} = {Galaxy.Dispersion}
{'Mass':15s} = {Galaxy.Mass:<10.2e}, {'PA':15s} = {Galaxy.PA}
{'Beams':15s} = {Galaxy.Beams:<10.1f}, {'Warp':15s} = {Galaxy.Warp}
{'Radial_Motions':15s} = {Galaxy.Radial_Motions:<10.1f}, {'SNR':15s} = {Galaxy.SNR}
{'Channel Width':15s} = {Galaxy.Channelwidth:<10.1f}, {'Beam Size':15s} = {Galaxy.Res_Beam}
{'Arms':15s} = {Galaxy.Arms:10s}, {'Flare':15s} = {Galaxy.Flare:15s}, {'Bar':15s} = {Galaxy.Bar:10s}
''')


print_base_galaxy.__doc__=f''' NAME:
    print_base_galaxy

 PURPOSE:
    print the contents of a galaxy class in an orderly fashion

 CATEGORY:
    common_functions

 INPUTS:
    galaxy class

 OPTIONAL INPUTS:

 OUTPUTS:
    print out of the galaxies parameters.

 OPTIONAL OUTPUTS:

 PROCEDURES CALLED:
    Unspecified

 NOTE:
 '''

def scrambled_initial(directory, Model):
    variables=['VROT', 'SBR','VRAD','SDIS','PA','INCL','XPOS','YPOS','VSYS','Z0']
    profile= {}

    for var in variables:
        if var in Model:
            ini= [float(x) for x in Model[var].split('=')[1].split()]
            ini2=[float(x) for x in Model[f'{var}_2'].split('=')[1].split()]
            while len(ini) < len(ini2):
                ini.append(ini[-1])
            while len(ini) > len(ini2):
                ini2.append(ini2[-1])
            if var == 'VROT':
                profile[var]=np.mean([ini[-5:-1], ini2[-5:-1]])
            elif var == 'SBR':
                profile[var] = np.max([ini, ini2])
            else:
                profile[var]=np.mean([ini[0], ini2[0]])
        else:
            profile[var]= 0.
    rng = np.random.default_rng()
    svrot= profile['VROT']+rng.uniform(-10., 10.)
    sincl= profile['INCL']+rng.uniform(-10., 10.)
    spa=  profile['PA']+rng.uniform(-3., 3.)
    sz0= profile['Z0']
    ssbr= rng.uniform(-1*profile['SBR'], profile['SBR'])
    ssdis= profile['SDIS']+rng.uniform(-4., 4.)
    svrad=profile['VRAD']+rng.uniform(-10., 10.)
    sra=profile['XPOS']+rng.uniform(-10./3600., -10./3600.)
    sdec=profile['YPOS']+rng.uniform(-10./3600., -10./3600.)
    svsys= profile['VSYS']+rng.uniform(-4., 4.)
    with open(f"{directory}Initial_Estimates.txt", 'w') as overview:
        overview.write(f'''#This file contains the random varied initial estimates.
#{'VROT':<15s} {'INCL':<15s} {'PA':<15s} {'Z0':<15s} {'SBR':<15s} {'DISP':<15s} {'VRAD':<15s} {'RA':<15s} {'DEC':<15s}  {'VSYS':<15s}
#{'km/s':<15s} {'Degree':<15s} {'Degree':<15s} {'arcsec':<15s} {'Jy km/s/arcsec^2':<15s} {'km/s':<15s} {'km/s':<15s} {'Degree':<15s} {'Degree':<15s}  {'km/s':<15s}
{svrot:<15.2f} {sincl:<15.2f} {spa:<15.2f} {sz0:<15.3f} {ssbr:<15.7f} {ssdis:<15.2f} {svrad:<15.2f} {sra:<15.5f} {sdec:<15.5f}  {svsys:<15.2f}''')

scrambled_initial.__doc__=f'''
NAME:
   scrambled_initial(directory,Model)

PURPOSE:
    Create a file with scrambelled Initial estimates.

CATEGORY:
   common_functions

INPUTS:
    directory = the directory where to put the model
    Model = array with def template with all lines split
OPTIONAL INPUTS:

OUTPUTS:
    text file

OPTIONAL OUTPUTS:

PROCEDURES CALLED:
   Unspecified

NOTE:
'''

#Function to read simple input files that  use = as a separator between ithe required input and the values

def read_template_file(filename, package_file = True):
    if float(sys.version[:3]) < 3.9 and package_file:
        with import_res.open_text(templates, filename) as tmp:
            unarranged = tmp.readlines()
    elif package_file:
        with import_res.files(templates).joinpath(filename).open('r') as tmp:
            unarranged = tmp.readlines()
    else:
        with open(filename, 'r') as tmp:
            unarranged= tmp.readlines()
    Template_in= Proper_Dictionary({})

    # Separate the keyword names
    for tmp in unarranged:
        # python is really annoying with needing endlines. Let's strip them here and add them when writing
        Template_in[tmp.split('=', 1)[0].strip().upper()]=tmp.rstrip()
    return Template_in

def rotateCube(Cube, angle, pivot,order=1):
    padX= [int(Cube.shape[2] - pivot[0]), int(pivot[0])]
    padY= [int(Cube.shape[1] - pivot[1]), int(pivot[1])]
    imgP= np.pad(Cube, [[0, 0], padY, padX], 'constant')
    #Use nearest neighbour as it is exact enough and doesn't mess up the 0. and is a lot faster
    imgR = rotate(imgP, angle, axes =(2, 1), reshape=False,order=order)
    return imgR[:, padY[0]: -padY[1], padX[0]: -padX[1]]
rotateCube.__doc__=f'''
 NAME:
    rotateCube(Cube, angle, pivot)

 PURPOSE:
    rotate a cube in the image plane

 CATEGORY:
    common_functions

 INPUTS:
    Cube = the cube data array
    angle = the angle to rotate under
    pivot = the point around which to rotate

 OPTIONAL INPUTS:

 OUTPUTS:
    the rotated cube

 OPTIONAL OUTPUTS:

 PROCEDURES CALLED:
    Unspecified

 NOTE:
'''
def regrid_array(oldarray, Out_Shape):
    oldshape = np.array(oldarray.shape)
    newshape = np.array(Out_Shape, dtype=float)
    ratios = oldshape/newshape
        # calculate new dims
    nslices = [ slice(0,j) for j in list(newshape) ]
    #make a list with new coord
    new_coordinates = np.mgrid[nslices]
    #scale the new coordinates
    for i in range(len(ratios)):
        new_coordinates[i] *= ratios[i]
    #create our regridded array
    newarray = map_coordinates(oldarray, new_coordinates,order=1)
    if any([x != y for x,y in zip(newarray.shape,newshape)]):
        print("Something went wrong when regridding.")
    return newarray
regrid_array.__doc__ =f'''
 NAME:
regridder
 PURPOSE:
Regrid an array into a new shape through the ndimage module
 CATEGORY:
    fits_functions

 INPUTS:
    oldarray = the larger array
    newshape = the new shape that is requested

 OPTIONAL INPUTS:

 OUTPUTS:
    newarray = regridded array

 OPTIONAL OUTPUTS:

 PROCEDURES CALLED:
    scipy.ndimage.map_coordinates, np.array, np.mgrid

 NOTE:
'''

def run_sofia(working_dir, parameter_file,sofia_call = 'sofia2'):
    sfrun = subprocess.Popen([sofia_call,parameter_file],cwd =working_dir, stdout = subprocess.PIPE, stderr = subprocess.PIPE)
    sofia_run, sofia_warnings_are_annoying= sfrun.communicate()
    print(sofia_run.decode("utf-8"))
    if sfrun.returncode == 8:
        raise SofiaFaintError(
            "RUN_SOFIA:Sofia cannot find a source in the input cube")
    elif sfrun.returncode == 0:
        sofia_ok= True
    else:
        print(sofia_warnings_are_annoying.decode("utf-8"))
        raise SofiaRunError(
            "RUN_SOFIA:Sofia did not execute properly. See screen for details")
run_sofia.__doc__=f'''
 NAME:
    run_sofia

 PURPOSE:
    run sofia

 CATEGORY:
    common_functions

 INPUTS:
    working_dir = directory where to run sofia
    parameter_file = input parameters file for sofia


 OPTIONAL INPUTS:
    sofia_call = command name for running sofia

 OUTPUTS:

 OPTIONAL OUTPUTS:

 PROCEDURES CALLED:
    Unspecified

 NOTE:
'''
