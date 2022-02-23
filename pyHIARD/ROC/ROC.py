#!/usr/local/bin/ python3
#This program is a python script to create a data base of real observations at different resolutions.
# It take 6 publicly available well resolved galaxy data cubes and smooths them to 3,4,5,6,7,8,10,12,16 beams across the major axis based on the extend of the best fit model.
# The galaxies used are
from pyHIARD.Resources import Cubes as cubes
from pyHIARD.constants import c_kms
from pyHIARD import Templates as templates
from multiprocessing import Pool,get_context
from astropy.wcs import WCS
from astropy.io import fits
import copy
import importlib
import numpy as np
import os
import pyHIARD.common_functions as cf
import scipy.ndimage
import sys
import re
import warnings

try:
    import importlib.resources as import_res
except ImportError:
    # Try backported to PY<37 `importlib_resources`.
    import importlib_resources as import_res


class ProgramError(Exception):
    pass


def add_template(cfg, path_to_resources, existing_galaxies):
    confirm_directory = False
    while not confirm_directory:
        confirm_directory = cf.get_bool(
            f"We are working in the directory {cfg.general.main_directory}. Is this correct? (Yes/No, default =True)", default=True)
        if not confirm_directory:
            cfg.general.main_directory = input(
                f"Which directory would you like to work in?")
        else:
            while not os.path.isdir(cfg.general.main_directory):
                print(
                    f'The directory {cfg.general.main_directory} does not exist please provide the correct directory')
                cfg.general.main_directory = input(
                    "Please provide the directory where to create the database :")
            confirm_directory = True

    galaxy_parameters = {'Galaxy': None, 'DHIkpc': None,
        'Distance': None, 'Original_Model': None, 'RMS': None, 'MHI':  None}
    galaxy_translation = {'Galaxy': 'the name to use in the package',
                          'DHIkpc':  'the HI diameter in kpc',
                          'Distance': 'distance in Mpc',
                          'Original_Model': 'the format of the orginal model (Tirific, Barolo, Rotcur)',
                          'RMS': 'the noise in the cube in data units',
                          'MHI':  'the HI mass'}

    existing_galaxies_low = [x.lower() for x in existing_galaxies]

    for key in galaxy_parameters:
        galaxy_parameters[key] = input(
            f"Please provide {galaxy_translation[key]} of the galaxy:")
        if key == 'Galaxy':
            if galaxy_parameters[key].lower() in existing_galaxies_low:
                remove_galaxy = cf.get_bool(
                    f"The galaxy {galaxy_parameters[key]} already exists, do you want to delete it? (default = no)", default=False)
                if remove_galaxy:
                    removed = remove_template(
                        galaxy_parameters[key], path_to_resources, existing_galaxies)
                    if not removed:
                        while galaxy_parameters[key].lower() in existing_galaxies_low:
                            galaxy_parameters[key] = input(
                                f"Please provide a different name (current = {galaxy_parameters[key]}).")

                else:
                    while galaxy_parameters[key].lower() in existing_galaxies_low:
                        galaxy_parameters[key] = input(
                            f"Please provide a different name (current = {galaxy_parameters[key]}).")

        if key == 'Original_Model':
            acceptable_model = False
            while not acceptable_model:
                if galaxy_parameters[key].lower() == 'tirific':
                    galaxy_parameters[key] = 'Tir'
                    acceptable_model = True
                elif galaxy_parameters[key].lower() == 'rotcur':
                    galaxy_parameters[key] = 'RC'
                    acceptable_model = True
                elif galaxy_parameters[key].lower() == 'barolo':
                    galaxy_parameters[key] = 'Bar'
                    acceptable_model = True
                else:
                    galaxy_parameters[key] = input(
                        f"{galaxy_parameters[key]} is not yet a model pyHIARD can process please type TiRiFiC, Barolo or Rotcur as an input model: ")

            model_file = input(f'Please provide the model text file:')
            while not os.path.isfile(f"{cfg.general.main_directory}{model_file}"):
                model_file = input(
                    f'Print we can not find the file {cfg.general.main_directory}{model_file}, please provide the correct name with a path from {cfg.general.main_directory}:')

            #try:
            test = cf.load_text_model(
            f"{cfg.general.main_directory}{model_file}", package_file=False,
             type=galaxy_parameters[key], Variables=['RADI', 'VROT', 'PA', 'INCL', 'XPOS', 'YPOS',
             'VSYS', 'VROT_2', 'PA_2', 'INCL_2', 'XPOS_2', 'YPOS_2', 'VSYS_2', 'Z0', 'SDIS', 'Z0_2', 'SDIS_2', 'CONDISP', 'SBR', 'SBR_2'])
            #except:
            #    print(f"We cannot read your file, please provide a standard TiRiFiC, Barolo or RotCur model.")
            #    print("We are exiting pyHIARD")
            #    sys.exit()

    input_fits_file = input("Please provide the galaxy fits file:")
    while not os.path.isfile(f"{cfg.general.main_directory}{input_fits_file}"):
        input_fits_file = input(
            f"We cannot find {cfg.general.main_directory}{input_fits_file}. Please provide the correct path from {cfg.general.main_directory}.")

    try:
        Cube = fits.open(f"{cfg.general.main_directory}{input_fits_file}",
                         uint=False, do_not_scale_image_data=True, ignore_blank=True)
        length = Cube[0].header['NAXIS3']
    except:
        print(f"We cannot read your fits file, please provide a standard fits file with 3 axes.")
        print("We are exiting pyHIARD")
        sys.exit()

    #galaxy_parameters = {'Galaxy': 'New_Galaxy', 'DHIkpc': '9.6', 'Distance': '4.1', 'Original_Model': 'Tir', 'RMS': '0.00038', 'MHI': '0.54e9'}
    #read or template and modify it
    with import_res.open_text(templates, 'roc_galaxy_template.py') as tmp:
        module_template = tmp.readlines()

    galaxy_line = "galaxy_parameters = {"
    for key in galaxy_parameters:
        if key in ['Galaxy', 'Original_Model']:
            galaxy_line = f"{galaxy_line}'{key}': '{galaxy_parameters[key]}', "
        else:
            galaxy_line = f"{galaxy_line}'{key}': {galaxy_parameters[key]}, "
    galaxy_line = galaxy_line[:-2]+"}"

    for i, line in enumerate(module_template):
        start = line.split('=')
        if start[0].strip() == 'galaxy_parameters':
            module_template[i] = galaxy_line
        if 'Input_Name' in line:
            line = line.replace('Input_Name', galaxy_parameters['Galaxy'])
            module_template[i] = line

    #Create the new directory and place all files in there
    new_resource = os.path.join(
        path_to_resources, galaxy_parameters['Galaxy'])+'/'
    cf.create_directory(galaxy_parameters['Galaxy'], path_to_resources)

    fits.writeto(
        f'{new_resource}/{galaxy_parameters["Galaxy"]}.fits', Cube[0].data, Cube[0].header, overwrite=False)
    Cube.close()
    Mask_Inner, Mask_Outer = cf.create_masks(
        new_resource, cfg.general.main_directory, galaxy_parameters['Galaxy'], sofia_call=cfg.general.sofia2)
    ext = {'Tir': 'def', 'Bar': 'txt', 'RC': 'rotcur'}
    os.system(
        f"cp {cfg.general.main_directory}{model_file} {new_resource}{galaxy_parameters['Galaxy']}.{ext[galaxy_parameters['Original_Model']]}")
    with open(f"{new_resource}{galaxy_parameters['Galaxy']}.py", 'w') as f:
        f.writelines(module_template)


add_template.__doc__ = f'''
NAME:
    add_template

PURPOSE:
    add a galaxy template to the package

CATEGORY:
   roc

INPUTS:
    cfg=  input configuration
    path_to_resources = path to the resources cube directory
    existing_galaxies = list of existing templates

OPTIONAL INPUTS:

OUTPUTS:
    template is added inside the package

OPTIONAL OUTPUTS:

PROCEDURES CALLED:
   Unspecified

NOTE:
'''


def beam_templates(beam,Galaxy_Template):

    fits.writeto(f'Input_Template{beam}.fits',Galaxy_Template['Galaxy_Template'],Galaxy_Template['Galaxy_Template_Header'],overwrite=True)
    # first we need to calculate the shift to  apply
    print(f"We are working on {beam} beam")
    Template_Header = Galaxy_Template['Galaxy_Template_Header']
    Def_Template = Galaxy_Template['Galaxy_Model']
     # The new beam we want = the normalsize/the beams across


    newbmaj = (Galaxy_Template['DHIarcsec'])/beam

    #which means we reduce everything by a factor of the ratio of the old and new beam
    fact = newbmaj /(Template_Header['BMAJ'] * 3600.)
    # bmin changes by the same factor
    newbmin = Template_Header['BMIN'] * 3600. * fact
    # And thus the new beam area is
    beamareanew = (np.pi * abs(newbmaj * newbmin)) / (4. * np.log(2.))
    # And in pixels
    pixperbeamnew = beamareanew / \
        (abs(Template_Header['CDELT1'] * 3600.)
         * abs(Template_Header['CDELT2'] * 3600.))

    # We want to regrid to 5 pixels per beam which means
    # Our new distance is
    Distance = fact *Galaxy_Template['Distance']
    # Update the def template
    Def_Template['BMAJ'] = f'BMAJ= {newbmaj}'
    Def_Template['BMIN'] = f'BMIN= {newbmin}'
    Def_Template['BPA'] = f'BPA= {Template_Header["BPA"]}' #Neeed to sort this out
    Def_Template['DISTANCE'] = 'DISTANCE = '+str(Distance)
    # And we need to adjust the header for our new cube

    # the new vsys on the hubble flow
    # First we fix the crval vaUNIlue to the actual pixel systemic is at
    New_Systemic = Galaxy_Template['Vsys'] * fact

    Template_Header['CRPIX3'] = ( Galaxy_Template['Vsys']-Template_Header['CRVAL3']) / \
                        Template_Header['CDELT3'] + \
                            Template_Header['CRPIX3']-1
    Template_Header["CRVAL3"] = New_Systemic


    # This is actually 1+z
    New_z = np.sqrt((1 + New_Systemic / c_kms) / (1 - New_Systemic / c_kms))
    # And we want to extend our input template by an amount of pixels to account for the required smoothin
    Pix_Extend = int(np.sqrt(newbmaj ** 2)
                     / abs(Template_Header['CDELT1'] * 3600.))

    # Which means the new total flux value after adjusting the flux and the dimming factor
    Current_Flux = Galaxy_Template['Total_Flux_In']/(fact**2)
    # And our dimmed z factor relates to tolman surface brightness dimming
    Shifted_Template = copy.deepcopy(
        Galaxy_Template['Galaxy_Template'])/(fact**2)*((Galaxy_Template['z']**4)/(New_z**4))
    Current_Mean = cf.get_mean_flux(Shifted_Template)

    #And we can now update our model input
    radius= [float(x) for x in Def_Template['RADI'].split('=')[1].split()]
    conv_radi = cf.convertskyangle(
        radius, distance=Distance, physical=True)
    Def_Template['RADI'] = 'RADI = '+' '.join([str(e) for e in conv_radi])

    scaleheight =  [float(x) for x in Def_Template['Z0'].split('=')[1].split()]
    con_hz = cf.convertskyangle(
        scaleheight, distance=Distance, physical=True)
    scaleheight2 =  [float(x) for x in Def_Template['Z0_2'].split('=')[1].split()]
    con_hz2 = cf.convertskyangle(
        scaleheight2, distance=Distance, physical=True)
    if np.sum(scaleheight) != 0.:
        Def_Template['Z0'] = 'Z0 = ' + \
            " ".join(str(e) for e in con_hz if e != 0.)
        Def_Template['Z0_2'] = 'Z0_2 = ' + \
            " ".join(str(e) for e in con_hz2 if e != 0.)
    Def_Template['VSYS'] = 'VSYS = ' + str(New_Systemic)
    Def_Template['VSYS_2'] = 'VSYS_2 = ' + str(New_Systemic)
    sbr = [float(x) for x in Def_Template['SBR'].split('=')[1].split()]
    sbr2 = [float(x) for x in Def_Template['SBR_2'].split('=')[1].split()]
    if np.sum(sbr) != 0.:
        # our surface brightness is constant except for the tolman dimming
        Def_Template['SBR'] = 'SBR = ' + \
            " ".join(str(e*((Galaxy_Template['z']**4)/(New_z**4)))
                     for e in sbr if e != 0.)
        Def_Template['SBR_2'] = 'SBR_2 = ' + \
            " ".join(str(e*((Galaxy_Template['z']**4)/(New_z**4)))
                     for e in sbr2 if e != 0.)
    # To obtain our new values we need to smooth the original with the squared difference between the old beam and the new beam
    FWHM_conv_maj = np.sqrt(newbmaj ** 2 -  (Template_Header['BMAJ'] * 3600.)** 2)
    FWHM_conv_min = np.sqrt(newbmin ** 2 - (Template_Header['BMIN'] * 3600.) ** 2)
    print("Which means we need a FWHM of {} x {}".format(
        FWHM_conv_maj, FWHM_conv_min))
    # and in terms of pixels the sigmas
    sig_maj = (FWHM_conv_maj / np.sqrt(8 * np.log(2))) / \
               abs(Template_Header['CDELT2'] * 3600.)
    sig_min = (FWHM_conv_min / np.sqrt(8 * np.log(2))) / \
               abs(Template_Header['CDELT1'] * 3600.)
    # Then we want to make the cubes for the required signal to noise ratios
    sigma_new = [(newbmaj / abs(Template_Header['CDELT2']*3600.)) / (2 * np.sqrt(2 * np.log(2))),
             (newbmin / abs(Template_Header['CDELT1']*3600.)) / (2 * np.sqrt(2 * np.log(2)))]
    Ext_Template = np.zeros((
        Template_Header['NAXIS3'], Template_Header['NAXIS2']
            + 2 * Pix_Extend,
        Template_Header['NAXIS1'] + 2 * Pix_Extend))
    Ext_Template[:, Pix_Extend:Template_Header['NAXIS2'] + Pix_Extend,
        Pix_Extend:Template_Header['NAXIS1'] + Pix_Extend] = Shifted_Template
    final_clean = scipy.ndimage.gaussian_filter(
        Ext_Template, sigma=(0, sig_maj, sig_min), order=0)
    # Preserve surface brightness
    final_clean = final_clean * pixperbeamnew / Galaxy_Template['Galaxy_Beam'][2]
    # Let's make a mask from this smoothed cube to calculated the things we achieve
    Final_Mask = copy.deepcopy(final_clean)
    Final_Mask[final_clean > np.mean(final_clean[final_clean > 0.])/2.] = 1.
    Final_Mask[final_clean < np.mean(
        final_clean[final_clean > 0.]) / 2.] = 0.
    Final_Mean = cf.get_mean_flux(final_clean,Mask=Final_Mask)
    final_clean = []
    fits.writeto(f'Beam_Template{beam}.fits',Shifted_Template,Template_Header,overwrite=True)
    Template_Dictionary = {'Name': Galaxy_Template['Name'],'Beams':beam,'Galaxy_Template':Shifted_Template,\
                           'Galaxy_Template_Header':Template_Header, 'Final_Mask': Final_Mask,\
                           'Galaxy_Model':Def_Template,'Galaxy_Mask':Galaxy_Template['Galaxy_Mask'],\
                           'Galaxy_Beam':Galaxy_Template['Galaxy_Beam'],'MHI':Galaxy_Template['MHI'],\
                           'Shifted_Beam':[*sigma_new,pixperbeamnew],'Shift_Factor':fact,\
                           'Beam_Shift':[sig_maj,sig_min],'New_Beam': [newbmaj,newbmin,Template_Header["BPA"]],\
                           'Noise': Galaxy_Template['Noise'],'Mean_Flux':Final_Mean,\
                           'DHI_kpc':Galaxy_Template['DHI_kpc'],\
                           'z':New_z,'Original_z': Galaxy_Template['z'],'Distance': Distance,\
                           'Requested_SNR' :Galaxy_Template['Requested_SNR'],'Disclaimer':Galaxy_Template['Disclaimer'],
                          }
    return Template_Dictionary
beam_templates.__doc__ = f'''
NAME:
    beam_templates(cfg,beams,Galaxy_Template)

PURPOSE:
    construct all differnt beam templates

CATEGORY:
   roc

INPUTS:
    cfg=  input configuration
    path_to_resources = path to the resources cube directory
    existing_galaxies = list of existing templates

OPTIONAL INPUTS:

OUTPUTS:
    template is added inside the package

OPTIONAL OUTPUTS:

PROCEDURES CALLED:
   Unspecified

NOTE:
'''


def check_templates(name,path_to_resources,work_dir,sofia_call='sofia2'):

    #First check if the main cube is there
    file_exists = os.path.isfile(f"{path_to_resources}{name}/{name}.fits")
    if not file_exists:
        galaxy_module = importlib.import_module(f'pyHIARD.Resources.Cubes.{name}.{name}')
        Template_All=galaxy_module.get_data()
        Template_All.close()
    mask_exists =  os.path.isfile(f"{path_to_resources}{name}/{name}_mask.fits")
    if not mask_exists:
        galaxy_module = importlib.import_module(f'pyHIARD.Resources.Cubes.{name}.{name}')
        Mask = galaxy_module.get_masks(work_dir,sofia_call=sofia_call)
        mask_exists = True
    file_exists = all([file_exists,mask_exists])
    return file_exists
check_templates.__doc__= f'''
NAME:
   check_templates(name,path_to_resources)

PURPOSE:
    Check that all the templates and their masks exist, otherwise download and create them now

CATEGORY:
   roc

INPUTS:
    name = basename of template
    path_to_resources = path to the cubes directory

OPTIONAL INPUTS:

OUTPUTS:
    boolean to indicate succes

OPTIONAL OUTPUTS:

PROCEDURES CALLED:
   Unspecified

NOTE:
'''
def create_final_cube(required_noise,main_directory,Galaxy_Template):
    # Make a new directory
    dirstring = f"{Galaxy_Template['Name']}_{Galaxy_Template['Beams']:.1f}Beams_{required_noise:.1f}SNR"
    print("Processing the SNR {}".format(required_noise))
    galaxy_dir =f"{main_directory}{dirstring}/"
    galaxy_dir_exists = os.path.isdir(galaxy_dir)
    if not galaxy_dir_exists:
        os.system(f"mkdir {galaxy_dir}")
    Galaxy_Template['Disclaimer'](galaxy_dir)

    # The noise in the final cube should be
    Final_Noise = Galaxy_Template['Mean_Flux']/required_noise
    # We first find the noise in pixels that matches this
    print(f''' We are requesting a SNR of {required_noise} and the mean Flux in the shifted template is = {Galaxy_Template['Mean_Flux']}
Creating the noise cube. The Noise in the final cube should be {Final_Noise} Jy/beam.''')
    New_Noise = 0.
    #As we will be adding to the uncorrected cube we need to convert the noise back to its uncorrected value.
    Final_Noise = Final_Noise * Galaxy_Template['Galaxy_Beam'][2]/Galaxy_Template['Shifted_Beam'][2]

    # The pixelated noise estimated is
    Pix_Noise = (((Final_Noise * Galaxy_Template['Shifted_Beam'][0] * 2 * np.sqrt(np.pi)) + (
            Final_Noise * Galaxy_Template['Shifted_Beam'][1] * 2 * np.sqrt(np.pi))) / 2.)
    New_Pixel_Size = (Galaxy_Template['New_Beam'][1] / 5.) / \
                      abs(Galaxy_Template['Galaxy_Template_Header']['CDELT1'] * 3600.)
    # And we want to extend our input template by an amount of pixels to account for the required smoothin
    Pix_Extend = int(np.sqrt(Galaxy_Template['New_Beam'][0]  ** 2)
                     / abs(Galaxy_Template['Galaxy_Template_Header']['CDELT1'] * 3600.))
    # Fine tune the final noise
    while abs(New_Noise - Final_Noise) / Final_Noise > 0.025:

        # fillthe Cube with single pixel noise
        if New_Noise != 0:
            #if abs(New_Noise - Final_Noise) / Final_Noise >0.05:
            Pix_Noise = Pix_Noise / (New_Noise / Final_Noise)
            #else:
            #    Pix_Noise = Pix_Noise + (Final_Noise - New_Noise)
        Ext_Template = np.random.normal(scale=Pix_Noise, size=(
            Galaxy_Template['Galaxy_Template_Header']['NAXIS3'],
            Galaxy_Template['Galaxy_Template_Header']['NAXIS2'] + 2 * Pix_Extend,
            Galaxy_Template['Galaxy_Template_Header']['NAXIS1'] + 2 * Pix_Extend))
        Final_Template = scipy.ndimage.gaussian_filter(Ext_Template, sigma=(0,
             Galaxy_Template['Shifted_Beam'][0],
             Galaxy_Template['Shifted_Beam'][1]),order=0)
        CHANNEL1 = Final_Template[0:2, :, :]
        New_Noise = np.std(CHANNEL1[np.isfinite(CHANNEL1)])
        print("We got a noise cube with an rms of {} {} {}".format(New_Noise, Final_Noise, Pix_Noise))
    # then we constuct a noise cube at the resolution of the galaxyBPA
        #As we are going to rotatet the cube we should first extend it
    #if our beam BPA is not 0 we need to extend the template a bit more
    if Galaxy_Template['Galaxy_Template_Header']['BPA'] != 0.:
        Galaxy_Template['Galaxy_Template']= cf.rotateCube(Galaxy_Template['Galaxy_Template'],\
        (Galaxy_Template['Galaxy_Template_Header']['BPA']),\
        [Galaxy_Template['Galaxy_Template_Header']['CRPIX1'],
        Galaxy_Template['Galaxy_Template_Header']['CRPIX2']])
        if required_noise != 0.5:


            compare = copy.deepcopy(Galaxy_Template['Galaxy_Template'])*Galaxy_Template['Shifted_Beam'][2]/Galaxy_Template['Galaxy_Beam'][2]
            fits.writeto('Test.fits',compare,Galaxy_Template['Galaxy_Template_Header'],overwrite=True)
        shift = int(abs(np.sin(np.radians(Galaxy_Template['Galaxy_Template_Header']['BPA'])))*\
                    (Galaxy_Template['Galaxy_Template_Header']['NAXIS1']/2.+Pix_Extend)+3.)
        Pix_Extend= int(Pix_Extend+shift)

    Ext_Template = np.random.normal(scale=Pix_Noise, size=(
        Galaxy_Template['Galaxy_Template_Header']['NAXIS3'],
        Galaxy_Template['Galaxy_Template_Header']['NAXIS2'] + 2 * Pix_Extend,
        Galaxy_Template['Galaxy_Template_Header']['NAXIS1'] + 2 * Pix_Extend))


    Final_Template = scipy.ndimage.gaussian_filter(Ext_Template, sigma=(0,
         Galaxy_Template['Galaxy_Beam'][0],
         Galaxy_Template['Galaxy_Beam'][1]),order=0)
    # Which means that the noise we require at this resolution is
    CHANNEL1 = Final_Template[0:2, :, :]
    Temp_Noise = np.std(CHANNEL1[np.isfinite(CHANNEL1)])
    # If this noise differs significantly from the Original noise in the cube then we want to add noise to the emission part as well.

    Diff_Noise = Temp_Noise - Galaxy_Template['Noise']/\
                (Galaxy_Template['Shift_Factor']**2)*\
                ((Galaxy_Template['Original_z']**4)/(Galaxy_Template['z']**4))
    # If this difference is less than 10 % we will ignore it
    if abs(Diff_Noise) < Temp_Noise/10.:
        Diff_Noise = 0.
    print("We want the noise at native resolution to be {}. And the difference noise {}".format(Temp_Noise,Diff_Noise))
    # If the new noise is smaller than the input noise we get into trouble an we do not want to return ProgramError(f"This should not happen") as there would be a higher noise on the emission
    if Diff_Noise < 0:
        with open(f"{galaxy_dir}Why_This_Galaxy_Is_Not_There.txt",'w') as file:
            file.write(f'''Your requested noise is lower than the input noise hence the emission would be too noisy please lower SNR.
The requested noise is {Temp_Noise} and the Original noise is {Galaxy_Template['Noise']/(Galaxy_Template['Shift_Factor']**2)}.
We continue with the next SNR value.''')
        return "EMPTY"

    # If we want to add noise  we construct a new noise cube for this purpose
    if Diff_Noise > 0.:
        print("Creating the  difference noise cube. Shifted noise = {}.".format(Diff_Noise))
        Pix_Noise = ((Diff_Noise * Galaxy_Template['Galaxy_Beam'][0] * 2 * np.sqrt(np.pi)) + (
                Diff_Noise * Galaxy_Template['Galaxy_Beam'][1] * 2 * np.sqrt(np.pi))) / 2.
        New_Noise = 0.
        while abs(New_Noise - Diff_Noise) / Diff_Noise > 0.025:
            # fillthe Cube with single pixel noise
            Ext_Template = np.random.normal(scale=Pix_Noise, size=(
                Galaxy_Template['Galaxy_Template_Header']['NAXIS3'],
                Galaxy_Template['Galaxy_Template_Header']['NAXIS2'],
                Galaxy_Template['Galaxy_Template_Header']['NAXIS1']))
            Noise_Template = scipy.ndimage.gaussian_filter(Ext_Template, sigma=(0,
                Galaxy_Template['Galaxy_Beam'][0],
                Galaxy_Template['Galaxy_Beam'][1]), order=0)
            CHANNEL1 = Noise_Template[0:2, :, :]
            New_Noise = np.std(CHANNEL1[np.isfinite(CHANNEL1)])
            print("We got a difference noise cube with an rms of {}".format(New_Noise))
            Pix_Noise = Pix_Noise / (New_Noise / Diff_Noise)
        #We add this to the emission
        Galaxy_Template['Galaxy_Template'][Galaxy_Template['Galaxy_Template'] != 0.] = \
            Galaxy_Template['Galaxy_Template'][Galaxy_Template['Galaxy_Template'] != 0.] +\
            Noise_Template[Galaxy_Template['Galaxy_Template'] != 0.]
        #Current_Template[Current_Template != 0.] = Current_Template[Current_Template != 0.] + Noise_Template[Current_Template != 0.]
    # We no longer need the extended template
    #print("Finished the noise")
    Ext_Template = []
    # We make a copy of this noise the size of our template
    Noise_Template = copy.deepcopy(Final_Template[:, \
        Pix_Extend:Galaxy_Template['Galaxy_Template_Header']['NAXIS2'] + \
        Pix_Extend,Pix_Extend:Galaxy_Template['Galaxy_Template_Header']['NAXIS1']\
         + Pix_Extend])
    # Then the overlap region should be half template half new noise to avoid edges
    #Current_Template[Boundary_Mask > 0.05] = Current_Template[Boundary_Mask > 0.05]*Boundary_Mask[Boundary_Mask> 0.05]+Noise_Template[Boundary_Mask > 0.05]*(1-Boundary_Mask[Boundary_Mask > 0.05])
    Galaxy_Template['Galaxy_Template'] = Galaxy_Template['Galaxy_Template']*\
        Galaxy_Template['Galaxy_Mask']+Noise_Template *(1-Galaxy_Template['Galaxy_Mask'])
    #Current_Template = Current_Template * Boundary_Mask + Noise_Template * (1 - Boundary_Mask)
    # Finally we write the modified template into our extended template such that we can smooth it

    Final_Template[:, Pix_Extend:Galaxy_Template['Galaxy_Template_Header']['NAXIS2']\
        + Pix_Extend,Pix_Extend:Galaxy_Template['Galaxy_Template_Header']['NAXIS1'] +\
         Pix_Extend] = Galaxy_Template['Galaxy_Template']
    #print("Starting to smooth final")
    final = scipy.ndimage.gaussian_filter(Final_Template, sigma=(0, Galaxy_Template['Beam_Shift'][0], Galaxy_Template['Beam_Shift'][1]), order=0)
    #print("Finished to smooth final")

    if Galaxy_Template['Galaxy_Template_Header']['BPA'] != 0.:
        final_tmp = cf.rotateCube(final,\
        (-1.*Galaxy_Template['Galaxy_Template_Header']['BPA']),\
        [Galaxy_Template['Galaxy_Template_Header']['CRPIX1']+Pix_Extend,
        Galaxy_Template['Galaxy_Template_Header']['CRPIX2']+Pix_Extend])
        #And remove the shift
        final = copy.deepcopy(final_tmp[:, \
            shift:Galaxy_Template['Galaxy_Template_Header']['NAXIS2'] +int(2.*Pix_Extend-shift) \
            ,shift:Galaxy_Template['Galaxy_Template_Header']['NAXIS1']+int(2.*Pix_Extend-shift)])
        Pix_Extend -= shift
    # Preserve brightness temperature means to scale to the new area

    final = final * Galaxy_Template['Shifted_Beam'][2]/Galaxy_Template['Galaxy_Beam'][2]
    if required_noise != 0.5:
        fits.writeto(f"Test_ungrid_final.fits", final, Galaxy_Template['Galaxy_Template_Header'] ,
                 overwrite=True)#ANd write to our directory
    Achieved_Noise = np.std(final[0:2,:,:])
    Galaxy_Template['Galaxy_Model']['RMS'] = f"RMS = {str(Achieved_Noise)}"
    Achieved_Mean = cf.get_mean_flux(final,Mask=Galaxy_Template['Final_Mask'])
    Achieved_SNR = Achieved_Mean/Achieved_Noise
    # Regrid to 5 pix per beam
    print(" Which results in the new dimensions {} x {}".format(int(Galaxy_Template['Galaxy_Template_Header']['NAXIS2'] / New_Pixel_Size),
                                                                int(Galaxy_Template['Galaxy_Template_Header']['NAXIS1'] / New_Pixel_Size)))
    regrid = cf.regrid_array(final, Out_Shape=(
        int(Galaxy_Template['Galaxy_Template_Header']['NAXIS3']),
        int(Galaxy_Template['Galaxy_Template_Header']['NAXIS2'] / New_Pixel_Size),
        int(Galaxy_Template['Galaxy_Template_Header']['NAXIS1'] / New_Pixel_Size)))
    #also the mask Used
    regrid_mask =  cf.regrid_array(Galaxy_Template['Final_Mask'], Out_Shape=(
        int(Galaxy_Template['Galaxy_Template_Header']['NAXIS3']),
        int(Galaxy_Template['Galaxy_Template_Header']['NAXIS2'] / New_Pixel_Size),
        int(Galaxy_Template['Galaxy_Template_Header']['NAXIS1'] / New_Pixel_Size)))
    #print("Finished Regridding")
    # We have to update the header
    achieved = final.shape[1] / regrid.shape[1]
    for ext  in ['1','2']:
        #Only adjust the cdelts at the end as all the smoothing is based on the original pixel size
        Galaxy_Template['Galaxy_Template_Header'][f"CDELT{ext}"] = Galaxy_Template['Galaxy_Template_Header'][f'CDELT{ext}'] * achieved/Galaxy_Template["Shift_Factor"]
        Galaxy_Template['Galaxy_Template_Header'][f"CRPIX{ext}"] = (Galaxy_Template['Galaxy_Template_Header'][f'CRPIX{ext}']+Pix_Extend) / achieved

    Galaxy_Template['Galaxy_Template_Header']['DATAMAX'] = np.max(regrid)
    Galaxy_Template['Galaxy_Template_Header']['DATAMIN'] = np.min(regrid)

    #ANd write to our directory
    #print("Start writing")
    fits.writeto(f"{galaxy_dir}Convolved_Cube.fits", regrid, Galaxy_Template['Galaxy_Template_Header'] ,
                 overwrite=True)#ANd write to our directory
    Galaxy_Template['Galaxy_Template_Header']['DATAMAX'] = np.max(regrid_mask)
    Galaxy_Template['Galaxy_Template_Header']['DATAMIN'] = np.min(regrid_mask)
    fits.writeto(f"{galaxy_dir}mask.fits", regrid_mask, Galaxy_Template['Galaxy_Template_Header'] ,
                 overwrite=True)
    #print("Finished writing")
    # Then we also want to write some info about the galaxy
    RAdeg = float(Galaxy_Template['Galaxy_Model']['XPOS'].split('=')[1].split()[0])
    DECdeg = float(Galaxy_Template['Galaxy_Model']['YPOS'].split('=')[1].split()[0])
    RAhr, DEChr = cf.convertRADEC(RAdeg, DECdeg)
    with open(f"{galaxy_dir}{dirstring}-Info.txt", 'w') as overview:
        overview.write(f'''This file contains the basic parameters of this galaxy.
For the radial dependencies look at Overview.png or ModelInput.def.
Inclination = {Galaxy_Template['Galaxy_Model']['INCL'].split('=')[1].split()[0]}.
The dispersion = {float(Galaxy_Template['Galaxy_Model']['SDIS'].split('=')[1].split()[0]):.2f}-{float(Galaxy_Template['Galaxy_Model']['SDIS'].split('=')[1].split()[-1]):.2f}.
The type of galaxy = {Galaxy_Template['Name']}.
PA = {Galaxy_Template['Galaxy_Model']['PA'].split('=')[1].split()[0]}.
Beams across the major axis = {Galaxy_Template['Beams']}.
SNR Requested = {required_noise} SNR Achieved = {Achieved_SNR}.
Mean Signal = {Achieved_Mean}.
Channelwidth = {Galaxy_Template['Galaxy_Template_Header']['CDELT3']}.
Major axis beam = {Galaxy_Template['New_Beam'][0]} Minor axis beam= {Galaxy_Template['New_Beam'][1]}.
It's central coordinates are RA={RAhr} DEC={DEChr} vsys={float(Galaxy_Template['Galaxy_Model']['VSYS'].split('=')[1].split()[0]):.2f} km/s.
At a Distance of {Galaxy_Template['Distance']:.2f} Mpc.
HI_Mass {Galaxy_Template['MHI']:.2e} (M_solar).
The final noise level is {Achieved_Noise} Jy/beam.
h_z =  {float(Galaxy_Template['Galaxy_Model']['Z0'].split('=')[1].split()[0]):.2f}-{float(Galaxy_Template['Galaxy_Model']['Z0'].split('=')[1].split()[-1]):.2f} (arcsec).''')

    # We need to make the model input
    with open(f"{galaxy_dir}ModelInput.def", 'w') as tri:
        tri.writelines([Galaxy_Template['Galaxy_Model'][key] + "\n" for key in Galaxy_Template['Galaxy_Model']])

    # And an overview plot
    #print("Start plotting")
    # the mass leads to the radius at which we need to hit 1 M/pc^2 from Wang (2016) in kpc
    R_HI=[Galaxy_Template['DHI_kpc']/2.,0.,0.]
    #Profiles are fitted by an exponential with scale length 0.2*RHI i.e. SigHI = C*exp(-R/(0.2*RHI)) so that weay we get the Warp end at Sigma = 0.5
    Warp = [0,np.log(0.5*np.exp(-1./0.2))*-0.2*R_HI[0]]
    cf.plot_input(galaxy_dir,Galaxy_Template['Galaxy_Model']
                ,Title=f'{Galaxy_Template["Name"]} with {Galaxy_Template["Beams"]} Beams'
                ,Distance= Galaxy_Template['Distance'],RHI=R_HI,WarpR=Warp
                )
    # And a file with scrambled initial estimates
    cf.scrambled_initial(galaxy_dir,Galaxy_Template['Galaxy_Model'])

    return f"{Galaxy_Template['Distance']}|{dirstring}|Convolved_Cube\n"

create_final_cube.__doc__= f'''
NAME:
   create_final_cube(required_noise, main_directory, Galaxy_Template)

PURPOSE:
    Add the requested noise and produce the output products

CATEGORY:
   ROC

INPUTS:
    required_noise = the requested noise
    main_directory = 'The main directory to branch from'
    Galaxy_Template = dictionary as created by beam_template

OPTIONAL INPUTS:

OUTPUTS:
    Cube, ModelInput.def, Initial Estimates and Overview plot.

OPTIONAL OUTPUTS:

PROCEDURES CALLED:
   Unspecified

NOTE:
'''


def galaxy_template(name,path_to_resources,work_directory,sofia2_call):
    galaxy_module = importlib.import_module(f'pyHIARD.Resources.Cubes.{name}.{name}')
    #Get the template
    Template_All=galaxy_module.get_data()
    Template_Header = Template_All[0].header
    Template_Cube = Template_All[0].data
    Template_All.close()


    # We want to be sure our header is in km/s
    # And a corresponding model input  file
    try:
        Vel_Units = Template_Header['CUNIT3'].lower()
    except KeyError:
        if Template_Header['CDELT3'] > 150.:
            Template_Header.set("CUNIT3", 'M/S', before="CTYPE3")
        else:
            Template_Header.set("CUNIT3", 'KM/S', before="CTYPE3")

    if Template_Header['CUNIT3'].lower() == 'm/s' or Template_Header['CDELT3'] > 150.:
        Template_Header['CDELT3'] = Template_Header['CDELT3'] / 1000.
        Template_Header['CUNIT3'] = 'km/s'
        Template_Header['CRVAL3'] = Template_Header['CRVAL3'] / 1000.
    # ensure we have a BPA
    if 'BPA' not in Template_Header:
        Template_Header['BPA'] = 0.

        #We assume the cubes to be centered
    #Obtain the model
    Model_Template = get_main_template(name,Template_Header,galaxy_module)
    # get a systemic from that

    systemic =  float(Model_Template['VSYS'].split('=')[1].split()[0])
    # This means we have created the unshifted template model
    beamarea=(np.pi*abs((Template_Header["BMAJ"]*Template_Header["BMIN"])*3600.**2))/(4.*np.log(2.))
    pixperbeam=beamarea/(abs(Template_Header['CDELT1']*3600.)*abs(Template_Header['CDELT2']*3600.))
    DHIarcsec=cf.convertskyangle(galaxy_module.galaxy_parameters['DHIkpc'],distance=galaxy_module.galaxy_parameters['Distance'],physical=True)
        #The minimum degradation is if we require more beams across the degradation is skipped
    max_beams_across=(DHIarcsec)/(Template_Header["BMAJ"]*3600.)
    #Obtain our Boundary Mask

    Mask = galaxy_module.get_masks(work_directory,sofia_call=sofia2_call)
    Mask[0].header['CRVAL3'] =  Mask[0].header['CRVAL3']/1000.
    Mask[0].header['CDELT3'] =  Mask[0].header['CDELT3']/1000.
    #fits.writeto(f'{work_directory}/{name}_mask.fits',Mask[0].data,Mask[0].header,overwrite = True)
    Boundary_Mask = Mask[0].data
    Mask.close()
    #Should we create this every time? Leave for now
    #We mask all values that fall outside the outer mask
    Template_Cube[Boundary_Mask <= 0.] = 0.
    Template_Cube[np.isnan(Template_Cube)] = 0.

        #as we have now cut the cubes to a minimal size we want to make them square and add the size across
        # such that 2 beams still has a beam extra on all sides
    Add_Pix = abs(DHIarcsec/(abs(Template_Header['CDELT1']*3600.)))
    new_size = int(np.max(Template_Cube.shape[1:2])+Add_Pix)
    tmp = np.zeros([Template_Cube.shape[0],new_size,new_size])
        #Place the old cube in the extended cube correctly
    min_ext_y= int((new_size-Template_Cube.shape[1])/2)
    min_ext_x= int((new_size-Template_Cube.shape[2])/2)
    tmp[:,min_ext_y:min_ext_y+Template_Cube.shape[1],min_ext_x:min_ext_x+Template_Cube.shape[2]] = Template_Cube[:,:,:]
    #Update the header
    Template_Header['CRPIX1'] = Template_Header['CRPIX1']+min_ext_x
    Template_Header['CRPIX2'] = Template_Header['CRPIX2']+min_ext_y
    Template_Header['NAXIS1'] = new_size
    Template_Header['NAXIS2'] = new_size
    Template_Cube = copy.deepcopy(tmp)
        # do the same for the Boundary Mask
    tmp = np.zeros([Template_Cube.shape[0],new_size,new_size])
    tmp[:,min_ext_y:min_ext_y+Boundary_Mask.shape[1],min_ext_x:min_ext_x+Boundary_Mask.shape[2]] = Boundary_Mask[:,:,:]
    Boundary_Mask = copy.deepcopy(tmp)
    tmp = []
    #Now that they are square if we have a BPA we want to rotate the mask and the template
    Total_Flux_In = np.sum(Template_Cube[Template_Cube > 0.])/pixperbeam * Template_Header['CDELT3']#
    # To calculate the mean of all values in the mask is too sensitive to very small variations in the mask
    Original_Mean = cf.get_mean_flux(Template_Cube)
    print(f'''The mean as determined from the template = {Original_Mean}
The noise is {galaxy_module.galaxy_parameters["RMS"]}'''    )
        #PA is anticlockwise while rotate works clockwise so to rotate back means rotate by BPA
    #    Template_Cube =cf.rotateCube(Template_Cube,\
    #        (Template_Header['BPA']),\
    #        [Template_Header['CRPIX1'],
    #        Template_Header['CRPIX2']])
    #    Boundary_Mask =cf.rotateCube(Boundary_Mask,\
    #        (Template_Header['BPA']),\
    #        [Template_Header['CRPIX1'],
    #        Template_Header['CRPIX2']])
        #
        #shift = int(abs(np.sin(np.radians(Template_Header['BPA']))) \
        #                *(Template_Header['NAXIS1']/2.+Pix_Extend)+5)
    #fits.writeto('After_the_Rot.fits',Template_Cube,Template_Header,overwrite=True)
    #fits.writeto('After_the_Rot_Mask.fits',Boundary_Mask,Template_Header,overwrite=True)
    #print(np.mean(Template_Cube[Template_Cube > 0.]))
    #exit()
    #And then calculate the total flux in the input cube
    #Total_Flux_In = np.sum(Template_Cube[Template_Cube > 0.])/pixperbeam * Template_Header['CDELT3']# In Jy
    # and the the original SNR is the mean signal in the masked cube / by the noise
    #Original_Mean=np.mean(Template_Cube[Template_Cube > 0.])
    #This is here defined as 1
    Original_z= np.sqrt((1+systemic/c_kms)/(1-systemic/c_kms))
    sigma = [(Template_Header["BMAJ"] / abs(Template_Header['CDELT1'])) / (2 * np.sqrt(2 * np.log(2))),
             (Template_Header["BMIN"] / abs(Template_Header['CDELT2'])) / (2 * np.sqrt(2 * np.log(2)))]
    #noise = galaxy_module.galaxy_parameters["RMS"]
    #if 'Scale_Factor' in  galaxy_module.galaxy_parameters:
    #    noise= noise*galaxy_module.galaxy_parameters['Scale_Factor']
    #Template_Header['BPA'] = 0.
    Template_Dictionary = {'Name':name,'Galaxy_Template':Template_Cube,'Galaxy_Template_Header':Template_Header,\
                           'Galaxy_Model':Model_Template,'Galaxy_Mask':Boundary_Mask,\
                           'Noise': galaxy_module.galaxy_parameters["RMS"],'Mean_Flux':Original_Mean,\
                           'z':Original_z, 'Vsys':systemic,'Max_Beams_Across':max_beams_across,\
                           'Distance': galaxy_module.galaxy_parameters["Distance"],\
                           'DHIarcsec': DHIarcsec, 'Total_Flux_In': Total_Flux_In,\
                           'DHI_kpc': galaxy_module.galaxy_parameters["DHIkpc"],\
                           'Galaxy_Beam':[*sigma,pixperbeam], 'Disclaimer':galaxy_module.place_disclaimer,\
                           'MHI':galaxy_module.galaxy_parameters['MHI']}
    return Template_Dictionary
galaxy_template.__doc__=f'''
NAME:
    galaxy_templates(name,beams):

PURPOSE:
    Create all beam templates and the Model Input for a specific galaxy

CATEGORY:
   roc

INPUTS:
    name = basename of template
    beams = the requested beams

OPTIONAL INPUTS:

OUTPUTS:
    boolean to indicate succes

OPTIONAL OUTPUTS:

PROCEDURES CALLED:
   Unspecified

NOTE:
'''


def get_main_template(name,Template_Header,galaxy_module):
    galaxy_resource_path = os.path.dirname(os.path.abspath(galaxy_module.__file__))
    homogenized_template_exists = os.path.isfile(f"{galaxy_resource_path}/{name}_uniform.def")
    if homogenized_template_exists:
        Template_Model= cf.read_template_file(f"{galaxy_resource_path}/{name}_uniform.def",package_file = False)
    else:
        radius, rotation, pa,incli,xpos,ypos,systemic,rotation2, pa2,incli2,xpos2,ypos2,systemic2,\
             scaleheight, dispersion, scaleheight2, dispersion2,condisp,sbr,sbr2 = cf.load_text_model(
                name,type =galaxy_module.galaxy_parameters['Original_Model'] ,\
                Variables=['RADI','VROT','PA','INCL','XPOS','YPOS','VSYS','VROT_2',\
                'PA_2','INCL_2','XPOS_2','YPOS_2','VSYS_2','Z0','SDIS','Z0_2',\
                'SDIS_2','CONDISP', 'SBR', 'SBR_2'])


        #Homogenize the input model
        if galaxy_module.galaxy_parameters['Original_Model'] in ['Bar','RC']:
            ndisks= 1
            rotation2= rotation
            pa2=  pa
            incli2= incli
            if galaxy_module.galaxy_parameters['Original_Model'] in ['RC']:
                xpos = [x+Template_Header['CRPIX1'] for x in xpos]
                ypos = [x+Template_Header['CRPIX2'] for x in ypos]
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                input_wcs = WCS(Template_Header).celestial
            xpos,ypos = zip(*[input_wcs.wcs_pix2world(x,y,1) for x,y in zip(xpos,ypos)])
            xpos2=xpos
            ypos2=ypos

        if condisp[0] == -1.:
            condisp = [(Template_Header['CDELT3']*1.2/(2*np.sqrt(2*np.log(2.))))]
        if np.sum(dispersion) == -1.*len(dispersion):
            dispersion[:] = np.sqrt(condisp[0]**2-(Template_Header['CDELT3']*1.2/(2*np.sqrt(2*np.log(2.))))**2)
            dispersion2[:] = np.sqrt(condisp[0]**2-(Template_Header['CDELT3']*1.2/(2*np.sqrt(2*np.log(2.))))**2)
                #continue
        if xpos[0] < 0.:
            xpos = [e+360. for e in xpos]
        if xpos2[0] < 0.:
            xpos2 = [e+360. for e in xpos2]
        RAdeg = xpos[0]
        DECdeg = ypos[0]
        #All values that are in arcsec need to be converted to kpc
        radius_kpc = cf.convertskyangle([x for x in radius if x != -1.],distance=galaxy_module.galaxy_parameters['Distance'])
        if any([x != -1 for x in scaleheight]):
            scaleheight_kpc = cf.convertskyangle([x for x in scaleheight if x != -1.],distance=galaxy_module.galaxy_parameters['Distance'])
        else:
            scaleheight_kpc = [0.]
        if any([x != -1 for x in scaleheight2]):
            scaleheight2_kpc = cf.convertskyangle([x for x in scaleheight2 if x != -1.], distance=galaxy_module.galaxy_parameters['Distance'])
        else:
            scaleheight2_kpc = [0.]
        #Everything that is constant can be written to our the Template def file
        #Read the default template
        Template_Model = cf.read_template_file('Template.def')
        Template_Model['NUR'] = f'NUR = {len(radius_kpc)}'
        Template_Model['RADI'] = f'RADI= {" ".join(str(e) for e in radius_kpc)}'
        Template_Model['INCL'] = 'INCL = '+" ".join(str(e) for e in incli if e != -1.)
        Template_Model['INCL_2'] = 'INCL_2 = ' + " ".join(str(e) for e in incli2 if e != -1.)
        Template_Model['PA'] = 'PA = ' + " ".join(str(e) for e in pa if e != -1.)
        Template_Model['PA_2'] = 'PA_2 = ' + " ".join(str(e) for e in pa2 if e != -1.)
        Template_Model['VROT'] = 'VROT = ' + " ".join(str(e) for e in rotation if e != -1.)
        Template_Model['VROT_2'] = 'VROT_2 = ' + " ".join(str(e) for e in rotation2 if e != -1)
        Template_Model['XPOS'] = 'XPOS = ' + " ".join(str(e) for e in xpos if e != -1. )
        Template_Model['XPOS_2'] = 'XPOS_2 = ' + " ".join(str(e) for e in xpos2 if e != -1.)
        Template_Model['YPOS'] = 'YPOS = ' + " ".join(str(e) for e in ypos if e != -1.)
        Template_Model['YPOS_2'] = 'YPOS_2 = ' + " ".join(str(e) for e in ypos2 if e != -1.)

        if any([x != -1 for x in dispersion]):
            Template_Model['SDIS'] = 'SDIS = ' + " ".join(str(e) for e in dispersion if e != -1.)
        else:
            Template_Model['SDIS'] = 'SDIS =  0.'
        if any([x != -1 for x in dispersion2]):
            Template_Model['SDIS_2'] = 'SDIS_2 = ' + " ".join(str(e) for e in dispersion2 if e != -1.)
        else:
            Template_Model['SDIS_2'] = 'SDIS_2 = 0. '



        if any([x != -1 for x in sbr]):
            Template_Model['SBR'] = 'SBR = ' + " ".join(str(e) for e in sbr if e != -1.)
        else:
            Template_Model['SBR'] = 'SBR =  0.'
        if any([x != -1 for x in sbr2]):
            Template_Model['SBR_2'] = 'SBR_2 = ' + " ".join(str(e) for e in sbr2 if e != -1.)
        else:
            Template_Model['SBR_2'] = 'SBR_2 = 0. '


        Template_Model['Z0'] = f'Z0 = {" ".join([str(e) for e in scaleheight_kpc])}'
        Template_Model['Z0_2'] = f'Z0_2 = {" ".join([str(e) for e in scaleheight2_kpc])}'
        try:
            Template_Header.set("BMAJ", Template_Header["BMMAJ"] / 3600.,before="BMMAJ")
            Template_Header.set("BMIN", Template_Header["BMMIN"] / 3600., before="BMMIN")
        except KeyError:
            print("No BMMAJ")
        Template_Model['BMAJ'] = f'BMAJ = {Template_Header["BMAJ"]*3600.}'
        Template_Model['BMIN'] = f'BMIN = {Template_Header["BMIN"]*3600.}'
        Template_Model['BPA'] = f'BPA = {Template_Header["BPA"]}'
        Template_Model['INSET'] = f'INSET = {name}.fits'
        Template_Model['OUTSET'] = f'OUTSET = {name}.fits'
        Template_Model['LOGNAME'] = f'LOGNAME= {name}.log'
        Template_Model['TIRDEF'] = f'TIRDEF= {name}.def'
        Template_Model['CONDISP']= f'CONDISP = {condisp[0]}'
        # As we are now storing this for the main we need to transfer all other values as well
        if any([x != -1 for x in systemic]):
            Template_Model['VSYS'] = f'VSYS = {" ".join([str(e) for e in systemic])}'
        else:
            Template_Model['VSYS'] = f'VSYS = 0'
        if any([x != -1 for x in systemic2]):
            Template_Model['VSYS_2'] = f'VSYS_2 = {" ".join([str(e) for e in systemic2])}'
        else:
            Template_Model['VSYS_2'] = f'VSYS_2 = 0'
        # write this as a homogenized model so we do not have to do this again
        with open(f"{galaxy_resource_path}/{name}_uniform.def",'w') as template:
            for key in Template_Model:
                template.write(f'{Template_Model[key]} \n')
    return Template_Model
get_main_template.__doc__=f'''
NAME:
    get_main_template(name,Template_Header,galaxy_module)

PURPOSE:
    load the main template corresponding to the galaxy

CATEGORY:
   roc

INPUTS:
    name = basename of template
    Template_Header = the header of the template cube. this is required when unifomizing the templates as barollo and rotcur are pixel based
    galaxy_module= the template module

OPTIONAL INPUTS:

OUTPUTS:
    The homogenized template model

OPTIONAL OUTPUTS:

PROCEDURES CALLED:
   Unspecified

NOTE:
'''

def remove_template(galaxy_to_remove,path_to_resources,existing_galaxies):
    existing_galaxies_low = [x.lower() for x in existing_galaxies]
    to_remove = existing_galaxies[existing_galaxies_low.index(galaxy_to_remove.lower())]
    really = cf.get_bool(f"Are you certain you want to remove the template for  {to_remove}? This action can not be undone. (default = yes)",default=True)
    if really:
        print(f"We will remove the template {to_remove}.")
        cf.delete_directory(f"{path_to_resources}{to_remove}")
        if not os.path.isdir(f"{path_to_resources}{to_remove}"):
            print(f"We have succesfully removed the template {to_remove}.")
            return True
        else:
            print(f"Something went wrong. please try again or start an issue on github")
            sys.exit()
    else:
        return False
remove_template.__doc__ = f'''
NAME:
  remove_template(galaxy_to_remove,path_to_resources,existing_galaxies):

PURPOSE:
    remove a galaxy template from the package

CATEGORY:
   roc

INPUTS:
    galaxy_to_remove = Name of the template to remove
    path_to_resources = path to the resources cube directory
    existing_galaxies = list of existing templates

OPTIONAL INPUTS:

OUTPUTS:
    boolean indicating succes or abort

OPTIONAL OUTPUTS:

PROCEDURES CALLED:
   Unspecified

NOTE:
'''

def ROC(cfg,path_to_resources):

    # Let's give an over view of the database that will be created
    print(f"We will create a database with the galaxies {','.join([x for x in cfg.roc.base_galaxies])} basic sets in the directory {cfg.general.main_directory}.\n")
    print("We will vary the following parameters")

    if 'Beams' in cfg.roc.variables_to_vary:
        print("We create the model with {} beams across the major axis.\n".format(", ".join([str(e) for e in cfg.roc.beams])))
    if 'SNR' in cfg.roc.variables_to_vary:
        print("Varying the signal to noise ratio with: {}.\n".format(" ".join([str(e) for e in cfg.roc.snr])))

    Clean_Cube = True
    # If we make new models delete everything in the directory
    if cfg.roc.delete_existing:
        to_delete  =f'''rm -R {' '.join([f"{cfg.general.main_directory}{x}_*Beams_*SNR" for x in cfg.roc.base_galaxies])}'''
        os.system(to_delete)
    Catalogue=f'{cfg.general.main_directory}Output_ROC_Summary.txt'

    number_models= cf.get_created_models(Catalogue,cfg.roc.delete_existing)


    if 'Beams' not in cfg.roc.variables_to_vary:
        cfg.roc.beams= [-1]
    if 'SNR' not in cfg.roc.variables_to_vary:
        cfg.roc.snr= [-1]
    cfg.roc.beams = [float(x) for x in np.sort(cfg.roc.beams)]
    if cfg.roc.beams[0] == -1:
        cfg.roc.beams = np.roll(np.array(cfg.roc.beams),-1).tolist()
    else:
        cfg.roc.beams = list(cfg.roc.beams)
    #if we have a lot of models to create we should split up in chunks?

    Galaxies = {}
    for name in cfg.roc.base_galaxies:
        #Check whether we need to construct this
        beams_to_produce = []
        combined = []
        for nobeams in cfg.roc.beams:
            noise_to_produce = []
            for reqnoise in cfg.roc.snr:
                dirstring = f"{name}_{nobeams:.1f}Beams_{reqnoise:.1f}SNR"
                beams_to_produce.append(nobeams)
                noise_to_produce.append(reqnoise)
                galaxy_dir =f"{cfg.general.main_directory}{dirstring}/"
                galaxy_dir_exists = os.path.isdir(galaxy_dir)
                if galaxy_dir_exists:
            # Do we have a cube
                    galaxy_cube_exist = os.path.isfile(f"{galaxy_dir}Convolved_Cube.fits")

                    if galaxy_cube_exist:
                        beams_to_produce.remove(nobeams)
                        noise_to_produce.remove(reqnoise)
                        print(f"The galaxy {dirstring} galaxy appears fully produced")
            if nobeams not in beams_to_produce:
                print(f"All galaxies with {name} and {nobeams} across the major axis are  thought to be produced already")
            else:
                combined.append([nobeams,noise_to_produce])
        if len(beams_to_produce) > 0:
            Galaxies[name] = {'name':name, 'iterations': combined}



    needed_templates = [[x,path_to_resources,cfg.general.main_directory,cfg.general.sofia2] for x in Galaxies]
    if len(needed_templates) == 0:
        print(f"Seems like the ROC has nothing to do.")
        return
    #check for the templates
    with get_context("spawn").Pool(processes=cfg.general.ncpu) as pool:
        results = pool.starmap(check_templates, needed_templates)

    #Then we setup the templates for all beam iterations we wants
    with get_context("spawn").Pool(processes=cfg.general.ncpu) as pool:
        All_Galaxy_Templates = pool.starmap(galaxy_template, needed_templates)
    #clean arrays
    results = []
    needed_templates = []
    #All_Galaxy_Templates returns a list of dictionaries that should be complete for constructing all the different beam templates
    #We need to make a list with all the different beams

    different_beams= []
    for galaxy in All_Galaxy_Templates:
        key = galaxy['Name']
        beams_and_noise = [x for x in Galaxies[key]['iterations']]
        for x in beams_and_noise:
            if x[0] > galaxy['Max_Beams_Across']/1.25:
                print(f"You requested {x[0]} beams across. However this galaxy originally only has {galaxy['Max_Beams_Across']} beams across.")
                print(f"There needs to be a degradation in size for this code to work. Please just use the Original Cube for testing. Or use less than {galaxy['Max_Beams_Across']/1.25:.1f} beams across.")
                print("Continuing to the next beam sizes. If you only want different signal to noise with minimal degradation please set the beams to -1")
                continue
            elif x[0] == -1.:
                #We could not check before whether these exist
                nobeams = galaxy['Max_Beams_Across']/1.25
                noise_to_produce=[]
                for reqnoise in cfg.roc.snr:
                    if reqnoise == -1:
                        reqnoise= galaxy['Mean_Flux']/galaxy['Noise']
                    dirstring = f"{key}_{nobeams:.1f}Beams_{reqnoise:.1f}SNR"
                    noise_to_produce.append(reqnoise)
                    galaxy_dir =f"{cfg.general.main_directory}{dirstring}/"
                    galaxy_dir_exists = os.path.isdir(galaxy_dir)
                    if galaxy_dir_exists:
                # Do we have a cube
                        galaxy_cube_exist = os.path.isfile(f"{galaxy_dir}Convolved_Cube.fits")

                        if galaxy_cube_exist:
                            noise_to_produce.remove(reqnoise)
                            print(f"The galaxy {dirstring} galaxy appears fully produced")
                if len(noise_to_produce) == 0:
                    print(f"All galaxies with {name} and {nobeams} across the major axis are  thought to be produced already")
                else:
                    galaxy['Requested_SNR'] = noise_to_produce
                    different_beams.append((nobeams,galaxy))
            else:
                galaxy['Requested_SNR'] = x[1]
                if -1. in galaxy['Requested_SNR']:
                    reqnoise = galaxy['Mean_Flux']/galaxy['Noise']
                    galaxy['Requested_SNR'] = [reqnoise if x == -1. else x for x in galaxy['Requested_SNR']]
                different_beams.append((x[0],galaxy))
    All_Galaxy_Templates = []
    if len(different_beams) == 0:
        print(f"Seems like the ROC has nothing to do.")
        return
    #Now Contruct all the beam templates
    with get_context("spawn").Pool(processes=cfg.general.ncpu) as pool:
        All_Beam_Templates = pool.starmap(beam_templates,different_beams)
    different_beams = []



    # now setup a an array with the required noise input
    all_noise = []
    for galaxy in All_Beam_Templates:
        for value in galaxy['Requested_SNR']:
            all_noise.append((value,cfg.general.main_directory,galaxy))

    with get_context("spawn").Pool(processes=cfg.general.ncpu) as pool:
        results = pool.starmap(create_final_cube, all_noise)

    with open(Catalogue, 'a') as cat:
        for line in results:
            if line != 'EMPTY':
                cat.write(f"{number_models}|{line}")
                number_models += 1
    print(f"We created {number_models} models")





ROC.__doc__ = f'''
NAME:
    ROC

PURPOSE:
    Shift and degrade real galaxies.

CATEGORY:
    agc

INPUTS:
    cfg = OmegaConf Configuration file

OPTIONAL INPUTS:

OUTPUTS:
    A set of real galaxies in a setup ready for FAT fitting

OPTIONAL OUTPUTS:

PROCEDURES CALLED:
   Unspecified

NOTE:
'''
