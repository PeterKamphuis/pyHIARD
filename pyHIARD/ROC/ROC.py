#!/usr/local/bin/ python3
#This program is a python script to create a data base of real observations at different resolutions.
# It take 6 publicly available well resolved galaxy data cubes and smooths them to 3,4,5,6,7,8,10,12,16 beams across the major axis based on the extend of the best fit model.
# The galaxies used are

from pyHIARD.constants import c_kms, H_0, G_agc
from pyHIARD import Templates as templates
from multiprocessing import Pool, get_context
from astropy.wcs import WCS
from astropy.io import fits
import copy
import importlib
import numpy as np
import os
import psutil
import pyHIARD.common_functions as cf
import scipy.ndimage
import sys
import warnings

# be explicit
try:
    import importlib.resources as import_res
except ImportError:
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
    '''    
    if sys.version_info[0][:3]) < 3.9:
        with import_res.open_text(templates, 'roc_galaxy_template.py') as tmp:
            module_template = tmp.readlines()
    else:
    '''    
    with import_res.files(templates).joinpath('roc_galaxy_template.py').open('r') as tmp:
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
        f'{new_resource}/{galaxy_parameters["Galaxy"]}_Original.fits', Cube[0].data, Cube[0].header, overwrite=False)
    Cube.close()
    #Mask_Inner, Mask_Outer = cf.create_masks(
    #    new_resource, cfg.general.main_directory, galaxy_parameters['Galaxy'], sofia_call=cfg.general.sofia2)
    ext = {'Tir': 'def', 'Bar': 'txt', 'RC': 'rotcur'}
    os.system(
        f"cp {cfg.general.main_directory}{model_file} {new_resource}{galaxy_parameters['Galaxy']}.{ext[galaxy_parameters['Original_Model']]}")
    with open(f"{new_resource}{galaxy_parameters['Galaxy']}.py", 'w') as f:
        f.writelines(module_template)
    with open(f"{new_resource}__init__.py", 'w') as f:
        f.write(f'   ')

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


def beam_templates(beam_req, Galaxy_Template_In, max_degradation_factor,\
                        main_directory):
    Galaxy_Template = copy.deepcopy(Galaxy_Template_In)
    # first we need to calculate the shift to  apply
    print(f"We are processing the template {Galaxy_Template['Name']} with {beam_req} beams across the major axis.")
    #First we calculate the degradation we want

    new_beam = [(Galaxy_Template['Original_DHI_arcsec'])
                 / (3600.*beam_req)]  # let's do beams in degrees
    initial_factor = new_beam[0] / \
        Galaxy_Template['Galaxy_Template_Header']['BMAJ']
    # if this is larger than the maximal degradation we smooth and regird to safe time and memory
    #fits.writeto('Before_Regrid.fits', Galaxy_Template['Galaxy_Template_Cube'],
    #             Galaxy_Template['Galaxy_Template_Header'], overwrite=True)

    #print(
    #    f"Currently the template is {Galaxy_Template['Galaxy_Template_Cube'].size*Galaxy_Template['Galaxy_Template_Cube'].itemsize/(1024**2)} Mb")
    if initial_factor > max_degradation_factor:
        smooth_factor=initial_factor/max_degradation_factor
        extended_template_cube, extended_template_header, pixel_buffer=\
            extend_cube(Galaxy_Template['Galaxy_Template_Cube'],
                                 Galaxy_Template['Galaxy_Template_Header'],
                                 new_beam = Galaxy_Template['Galaxy_Template_Header']['BMAJ']
                                  * smooth_factor)

        template_cube, template_hdr, template_noise=smooth_and_regrid(extended_template_cube,
            extended_template_header, regrid=True, track_noise = Galaxy_Template['Original_Noise'],
            factor=smooth_factor)
        del extended_template_cube
        del extended_template_header

        del Galaxy_Template['Galaxy_Template_Cube']
        del Galaxy_Template['Galaxy_Template_Header']

        regrid_factor = new_beam[0]/(template_hdr['BMAJ'])
    else:
        #if we don't regrid the factors are the same
        smooth_factor = 1.
        regrid_factor = initial_factor
        template_cube = Galaxy_Template['Galaxy_Template_Cube']
        template_hdr = Galaxy_Template['Galaxy_Template_Header']
        template_noise = Galaxy_Template['Original_Noise']

    template_mask = cf.get_mask(template_cube,factor = initial_factor)

    new_mean = cf.get_mean_flux(template_cube,Mask=template_mask)

    # And the new bmin has the same factor
    new_beam.append(template_hdr['BMIN'] * regrid_factor)
    # calclate the pixels in the new beam area
    new_pixels_per_beam=cf.get_beam_area_in_pixels(
        template_hdr, beam = np.array(new_beam, dtype=float))


    #Start shifting the template
    # first a new distance
    new_distance=Galaxy_Template['Original_Distance']*initial_factor
    #And the corresponding systemic on the Hubble Flow +the original peculiar velocity
    peculiar_velocity=Galaxy_Template['Original_Vsys'] - \
        H_0*Galaxy_Template['Original_Distance']
    new_systemic=H_0*new_distance+peculiar_velocity
    # And the new redshift (This is actually 1+z)
    new_z=np.sqrt((1. + new_systemic / c_kms) / (1. - new_systemic / c_kms))
    #print(f"Is the z reasonable {new_z}")
    #And the shifted flux value adjusting the flux and the dimming factor
    #new_total_flux = Galaxy_Template['Original_Total_Flux_In'] / \
    #    (initial_factor**2)
    #
    new_DHI_arcsec=cf.convertskyangle(cf.convertskyangle(Galaxy_Template['Original_DHI_arcsec'],\
                    distance=Galaxy_Template['Original_Distance']),\
                    distance = new_distance, physical = True)
    #And now we want to shift the actual cube Which simply mean reducing the flux and account for Tollman dimming through the z factor
    new_template_cube=copy.deepcopy(template_cube)/(initial_factor**2) *\
            ((Galaxy_Template['Original_z']**4)/(new_z**4))
    #Calculate the noise corresponding to this

    new_noise=float(template_noise) /\
                (initial_factor**2) *\
                ((Galaxy_Template['Original_z']**4)/(new_z**4))
    new_mean=float(new_mean) /\
            (initial_factor**2) *\
            ((Galaxy_Template['Original_z']**4)/(new_z**4))

    #And update the header to the corresponding values

    template_hdr["CRVAL3"]=new_systemic
    #Now we update the TRM to new values
    new_TRM_model=copy.deepcopy(Galaxy_Template['Galaxy_TRM_Model'])
    for i, beam in enumerate(['BMAJ', 'BMIN']):
        new_TRM_model[beam]=f'{beam}= {new_beam[i]/initial_factor*3600.:.2f}'
    new_TRM_model['DISTANCE']=f'DISTANCE = {new_distance:.2f}'
    #And the radi, they should be stored in kpc in the uniform template
    radius=[float(x) for x in new_TRM_model['RADI'].split('=')[1].split()]
    conv_radi=cf.convertskyangle(
        radius, distance = new_distance, physical = True)
    new_TRM_model['RADI']='RADI = '+' '.join([str(e) for e in conv_radi])
    #same for the scale heights
    for i, var in enumerate(['', '_2']):
        scaleheight=[float(x)
                             for x in new_TRM_model[f'Z0{var}'].split('=')[1].split()]
        con_hz= cf.convertskyangle(
            scaleheight, distance=new_distance, physical=True)
        if np.sum(scaleheight) != 0.:
            new_TRM_model[f'Z0{var}'] = f'Z0{var} = {" ".join(str(e) for e in con_hz if e != 0.)}'
        else:
            new_TRM_model[f'Z0{var}'] = f'Z0{var} = 0.'
        new_TRM_model[f'VSYS{var}'] = f'VSYS{var} = {new_systemic}'
        sbr = [float(x)
                     for x in new_TRM_model[f'SBR{var}'].split('=')[1].split()]
        if np.sum(sbr) != 0.:
            new_TRM_model[f'SBR{var}'] =\
             f'SBR{var} = {" ".join(str(e*((Galaxy_Template["Original_z"]**4)/(new_z**4))) for e in sbr if e != 0.)}'
        else:
            new_TRM_model[f'SBR{var}'] = f'SBR{var} = 0.'
    # And we need to adjust the header for our new cube
    # We want to make a template at the final resolution to get some values from
    # Frist we ensure some pixel buffer after smoothing
    extended_new_template_cube, extended_new_template_hdr,pixel_buffer ,extended_mask= extend_cube(
            new_template_cube, template_hdr, new_beam=new_beam[0],Mask=template_mask,rotation_pa=template_hdr['BPA'])
    #print(f"We extended these image with {pixel_buffer} current size {extended_new_template_cube.shape} before {new_template_cube.shape}")
    template_at_final_res, hdr_at_final_res,noise_at_final_res = smooth_and_regrid(extended_new_template_cube,
        extended_new_template_hdr, factor=regrid_factor, regrid = False,track_noise=new_noise)

    # Let's make a mask from this smoothed cube to calculated the things we achieve
    Final_Mask=copy.deepcopy(template_at_final_res)
    Final_Mask[template_at_final_res > np.mean(template_at_final_res[template_at_final_res > 0.])/2.]=1
    Final_Mask[template_at_final_res < np.mean(
        template_at_final_res[template_at_final_res > 0.]) / 2.]=0
    Final_Mean=cf.get_mean_flux(template_at_final_res, Mask = Final_Mask)
    Final_Mask = np.array(Final_Mask,dtype=int)
    pxbp = cf.get_beam_area_in_pixels(hdr_at_final_res)
    Final_Total_Flux = np.sum(template_at_final_res)/pxbp *hdr_at_final_res['CDELT3']
    #Now we can check if this galaxy already exists
    if -1. in Galaxy_Template['Requested_SNR']:
        reqnoise = float(Final_Mean/noise_at_final_res)
        dirstring = f"{Galaxy_Template['Name']}_{beam_req:.1f}Beams_{reqnoise:.1f}SNR"
        galaxy_dir =f"{main_directory}{dirstring}/"
        galaxy_dir_exists = os.path.isdir(galaxy_dir)
        if galaxy_dir_exists:
    # Do we have a cube
            galaxy_cube_exist = os.path.isfile(f"{galaxy_dir}Convolved_Cube_Gauss.fits")
            if galaxy_cube_exist:
                Galaxy_Template['Requested_SNR'].remove(-1.)
                print(f"The galaxy {dirstring} galaxy appears fully produced")
    if len(Galaxy_Template['Requested_SNR']) == 0.:
        print(f"All galaxies for template {Galaxy_Template['Name']} with {beam_req} beams across the major axis are already produced.")
        return "EMPTY"


    #These values should be corrected for the increased beam
    del template_at_final_res
    del hdr_at_final_res
    #fits.writeto('Test_New_Shift.fits',extended_new_template_cube,extended_new_template_hdr,overwrite=True)
    Template_Dictionary={'Name': Galaxy_Template['Name'], 'Beams': beam_req, 'Shifted_Template_Cube': extended_new_template_cube,\
                           'Shifted_Template_Header': extended_new_template_hdr, 'Final_Mask': Final_Mask,\
                           'Pixel_Buffer': pixel_buffer,'Final_DHI_arcsec': new_DHI_arcsec,\
                           'Shifted_TRM_Model': new_TRM_model, 'Shifted_Mask': extended_mask,\
                           'M_HI': Galaxy_Template['M_HI'],'Final_Noise': noise_at_final_res,\
                           'Shift_Factor': [initial_factor, regrid_factor],\
                           'New_Beam': [*new_beam, extended_new_template_hdr["BPA"]],\
                           'Shifted_Original_Noise': new_noise, 'Final_Mean_Flux': Final_Mean,\
                           'Final_Distance': new_distance,\
                           'Requested_SNR': Galaxy_Template['Requested_SNR'], 'Disclaimer': Galaxy_Template['Disclaimer']
                          }

    del Galaxy_Template
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

def calculate_processes(templates,beams,SNR):
    available_memory = psutil.virtual_memory().total/2**30
    max_no_process = int(np.floor(available_memory/3.5)) # A single process run peaked at 7 Gb so this should suffice
    no_SNR_process = len(SNR)
    no_beam_process = 1
    no_template_process = 1
    if no_SNR_process > (max_no_process-2):
        no_SNR_process = max_no_process-2
    if (max_no_process-1) > 2.*no_SNR_process and no_SNR_process > 1:
        if len(beams)*no_SNR_process > max_no_process:
            left_process =  max_no_process-no_SNR_process-1
            while left_process > no_SNR_process:
                no_beam_process +=1
                left_process -= no_SNR_process
        else:
            no_beam_process = len(beams)
    if max_no_process > 2.*no_SNR_process*no_beam_process and no_beam_process > 1:
        if len(templates)*no_beam_process*no_SNR_process > max_no_process:
            left_process =  max_no_process-(no_SNR_process*no_beam_process)
            while left_process > no_SNR_process*no_beam_process:
                no_template_process +=1
                left_process -= no_SNR_process*no_beam_process
        else:
            no_template_process = len(templates)
    return no_template_process,no_beam_process,no_SNR_process


calculate_processes.__doc__ = f'''
NAME:
    calculate_processes(templates,beams,SNR)

PURPOSE:
    calculate the number of processes in each pool and limit the total to not go beyond the maximum memory

CATEGORY:
   roc

INPUTS:
    templates = templates to process
    beams = beams for each template
    SNR = all SNRs to process

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
    file_exists_uni = os.path.isfile(f"{path_to_resources}{name}/{name}_uniform.def")
    if not file_exists or not file_exists_uni:
        galaxy_module = importlib.import_module(f'pyHIARD.Resources.Cubes.{name}.{name}')
        Template_All=galaxy_module.get_data(work_dir,sofia_call)
        if not file_exists_uni:
            Model_Template = get_main_template(name,Template_All[0].header,galaxy_module)
        Template_All.close()

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
def create_final_cube(required_noise,main_directory,Galaxy_Template_In,font_file):
    Galaxy_Template = copy.deepcopy(Galaxy_Template_In)
    # Make a new directory
    native_noise = False
    if required_noise == -1.:
        native_noise = True
        required_noise = Galaxy_Template['Final_Mean_Flux']/Galaxy_Template['Final_Noise']

    dirstring = f"{Galaxy_Template['Name']}_{Galaxy_Template['Beams']:.1f}Beams_{required_noise:.1f}SNR"
    print(f'''{'Processing the template':>25s} = {Galaxy_Template['Name']}
{'Beams':>25s} = {Galaxy_Template['Beams']}
{'SNR':>25s} = {required_noise}''')
    galaxy_dir =f"{main_directory}{dirstring}/"
    galaxy_dir_exists = os.path.isdir(galaxy_dir)
    if not galaxy_dir_exists:
        os.system(f"mkdir {galaxy_dir}")
    #Place the disclaimer about the cube
    Galaxy_Template['Disclaimer'](galaxy_dir)
    # The noise in the final cube should be
    final_req_noise = Galaxy_Template['Final_Mean_Flux']/required_noise
    # We first find the noise in pixels that matches this
    #As we will be adding to the uncorrected cube we need to convert the noise back to its uncorrected value.
    # first we estimate what will be the noise we want for a cube with the resolution of the template
    # so we scale with the beam area new_beam/old_beam
    shifted_template_pix_beam = cf.get_beam_area_in_pixels(Galaxy_Template['Shifted_Template_Header'])
    new_beam_pix_beam = cf.get_beam_area_in_pixels(Galaxy_Template['Shifted_Template_Header'], beam=Galaxy_Template['New_Beam'][0:2])
    shifted_req_noise = final_req_noise *shifted_template_pix_beam/new_beam_pix_beam

    # The pixelated noise estimated is
    #First the smoothing in pixels required to get to our new beam
    new_pixel_sigma=(np.array(Galaxy_Template['New_Beam'][0:1],dtype=float)/ \
        np.abs(np.array([Galaxy_Template['Shifted_Template_Header'][key] for key in ['CDELT2','CDELT1'] ],dtype=float)))/\
        (2. * np.sqrt(2. * np.log(2.)))
    #And the sigma to get to the shifted template resolution
    template_pixel_sigma=(np.array([Galaxy_Template['Shifted_Template_Header'][key] for key in ['BMAJ','BMIN'] ],dtype=float)/ \
        np.abs(np.array([Galaxy_Template['Shifted_Template_Header'][key] for key in ['CDELT2','CDELT1'] ],dtype=float)))/\
        (2. * np.sqrt(2. * np.log(2.)))
    # Fine tune the pixel noise estimate to achieve the final noise
    #setup the random number generator
    achieved_shifted_noise,pixel_noise = cf.calculate_pixel_noise(shifted_req_noise,template_pixel_sigma)

    # If this noise differs significantly from the original noise in the cube then we want to add noise to the emission part as well.
    diff_noise = achieved_shifted_noise - Galaxy_Template['Shifted_Original_Noise']
    # If this difference is less than 10 % we will ignore it
    if abs(diff_noise) < achieved_shifted_noise/10. or native_noise:
        diff_noise = 0.
    #print(f'''We will get the noise at the shifted template resolution to  {achieved_shifted_noise} Jy/beam.
#This is a difference of {diff_noise} Jy/beam with the noise in the template at this shift.''')
    # If the new noise is smaller than the input noise we get into trouble an we do not want to return ProgramError(f"This should not happen") as there would be a higher noise on the emission
    if diff_noise < 0:
        with open(f"{galaxy_dir}Why_This_Galaxy_Is_Not_There.txt",'w') as file:
            file.write(f'''Your requested noise is lower than the input noise hence the emission would be too noisy please lower SNR.
The requested noise is {achieved_shifted_noise} and the original noise translated to this shift is {Galaxy_Template['Shifted_Original_Noise']}.
We continue with the next SNR value.''')
        return "EMPTY"
    #making a random number generator
    rng = np.random.default_rng()
    # calculate the noise to be added to the template
    # If we want to add noise  we construct a new noise cube for this purpose
    if diff_noise > 0.:
        achieved_diff_noise,pixel_diff_noise = cf.calculate_pixel_noise(diff_noise,template_pixel_sigma)
        #create a cube with this noise
        diff_noise_cube = rng.normal(scale=pixel_diff_noise, size=Galaxy_Template['Shifted_Template_Cube'].shape)
        # smooth it to the template resolution
        diff_noise_template = scipy.ndimage.gaussian_filter(diff_noise_cube,\
                        sigma=(0,*template_pixel_sigma),order=0)
        #We add this to the emission
        Galaxy_Template['Shifted_Template_Cube'][Galaxy_Template['Shifted_Template_Cube'] != 0.] = \
            Galaxy_Template['Shifted_Template_Cube'][Galaxy_Template['Shifted_Template_Cube'] != 0.] +\
            diff_noise_template[Galaxy_Template['Shifted_Template_Cube'] != 0.]
        del diff_noise_cube
        del diff_noise_template

    # create an extended cube with pixelize noise at the template resolution
    noise_template =  rng.normal(scale=pixel_noise, size=Galaxy_Template['Shifted_Template_Cube'].shape)
    #smooth to the template resolution
    noise_template = scipy.ndimage.gaussian_filter(noise_template,\
                        sigma=(0,*template_pixel_sigma),order=0)

    Galaxy_Template['Shifted_Template_Cube'] = Galaxy_Template['Shifted_Template_Cube']*Galaxy_Template['Shifted_Mask']+noise_template*(1.-Galaxy_Template['Shifted_Mask'])
    del noise_template

    del Galaxy_Template['Shifted_Mask']
    #
    final_cube, final_hdr,final_mask = smooth_and_regrid(Galaxy_Template['Shifted_Template_Cube'],
        Galaxy_Template['Shifted_Template_Header'], factor=Galaxy_Template['Shift_Factor'][1], Mask=Galaxy_Template['Final_Mask'],regrid=False)

    del Galaxy_Template['Shifted_Template_Cube']; del Galaxy_Template['Final_Mask']
    del Galaxy_Template['Shifted_Template_Header']
    #For accuracy it would be better to de-rotate before regridding, from an efficiency point this is better
    if final_hdr['BPA'] != 0.:
        final_cube = cf.rotateCube(final_cube,\
        (-1.*final_hdr['BPA']),\
        [final_hdr['CRPIX1'],
        final_hdr['CRPIX2']],order=1)
        final_mask = cf.rotateCube(final_mask,\
        (-1.*final_hdr['BPA']),\
        [final_hdr['CRPIX1'],
        final_hdr['CRPIX2']],order=1)
    #And regrid

        #And remove the shift
    #fits.writeto(f"{galaxy_dir}Test_final.fits", final_cube, final_hdr,
    #             overwrite=True)
    pixel_shift = int(np.sum(Galaxy_Template['Pixel_Buffer']))
    #print(f"This is the pixel shift we remove {pixel_shift} current size {final_cube.shape}")
    final_cube = copy.deepcopy(final_cube[:, \
            pixel_shift:final_hdr['NAXIS2']-pixel_shift,\
            pixel_shift:final_hdr['NAXIS1']-pixel_shift])
    final_mask=  copy.deepcopy(final_mask[:, \
            pixel_shift:final_hdr['NAXIS2']-pixel_shift,\
            pixel_shift:final_hdr['NAXIS1']-pixel_shift])

    final_hdr['NAXIS1'] = final_cube.shape[2]
    final_hdr['NAXIS2'] = final_cube.shape[1]
    final_hdr['CRPIX1'] = final_hdr['CRPIX1'] -pixel_shift
    final_hdr['CRPIX2'] = final_hdr['CRPIX2'] -pixel_shift

    #FINALLY WE REGRID
    pixel_factor = (final_hdr['BMIN'] / 5.) / abs(final_hdr['CDELT2'])
    if pixel_factor >   1.1:
        old_shape = final_cube.shape
        final_cube = cf.regrid_array(final_cube, Out_Shape=((int(final_hdr['NAXIS3']),
            int(final_hdr['NAXIS2'] /  pixel_factor),int(final_hdr['NAXIS1'] / pixel_factor))))
        final_mask = cf.regrid_array(final_mask, Out_Shape=((int(final_hdr['NAXIS3']),
            int(final_hdr['NAXIS2'] /  pixel_factor),int(final_hdr['NAXIS1'] / pixel_factor))))
        final_mask[final_mask > 0.5] = 1.
        final_mask[final_mask <= 0.5] = 0.

        for ext in ['1','2']:
            achieved_regrid = old_shape[int(ext)] / final_cube.shape[int(ext)]
            final_hdr[f'CDELT{ext}'] = final_hdr[f'CDELT{ext}']*achieved_regrid
            final_hdr[f'CRPIX{ext}'] = final_hdr[f'CRPIX{ext}']/achieved_regrid
            final_hdr[f'NAXIS{ext}'] = final_cube.shape[int(ext)]

    #Then shift the header values to the new distance
    final_hdr['CDELT2'] = final_hdr['CDELT2']/Galaxy_Template['Shift_Factor'][0]
    final_hdr['CDELT1'] = final_hdr['CDELT1']/Galaxy_Template['Shift_Factor'][0]
    final_hdr['BMAJ'] = final_hdr['BMAJ']/Galaxy_Template['Shift_Factor'][0]
    final_hdr['BMIN'] = final_hdr['BMIN']/Galaxy_Template['Shift_Factor'][0]
    #And preserve the brightness for the peam change


    Achieved_Noise = np.std(final_cube[0:2,:,:])
    Galaxy_Template['Shifted_TRM_Model']['RMS'] = f"RMS = {str(Achieved_Noise)}"
    Achieved_Mean = cf.get_mean_flux(final_cube,Mask=final_mask)
    Achieved_SNR = Achieved_Mean/Achieved_Noise


    final_hdr['DATAMAX'] = np.max(final_cube)
    final_hdr['DATAMIN'] = np.min(final_cube)
    fits.writeto(f"{galaxy_dir}Convolved_Cube_Gauss.fits", final_cube, final_hdr,
                 overwrite=True)#ANd write to our directory
    final_hdr['DATAMAX'] = np.max(final_mask)
    final_hdr['DATAMIN'] = np.min(final_mask)
    final_mask = np.array(final_mask, dtype=np.int16)
    fits.writeto(f"{galaxy_dir}mask.fits", final_mask, final_hdr,overwrite=True,output_verify='ignore')
    # This completes cube manupulation
    channel_width = final_hdr['CDELT3']
    # the mass leads to the radius at which we need to hit 1 M/pc^2 from Wang (2016) in kpc
    R_HI=[cf.convertskyangle(Galaxy_Template['Final_DHI_arcsec'],distance=\
            Galaxy_Template['Final_Distance'])/2.,0.,0.]
    V_max=float(Galaxy_Template['Shifted_TRM_Model']['VROT'].split('=')[1].split()[-1])
    DynMass=R_HI[0]*10**3*V_max**2/G_agc
    Achieved = {'SNR':Achieved_SNR, 'Mean_Signal': Achieved_Mean, 'Noise': Achieved_Noise,\
            'Channel_Width': channel_width,'Res_Beam':\
            [final_hdr['BMAJ']*3600.,final_hdr['BMIN']*3600.,final_hdr['BPA']],\
            'Pixel_Beam': cf.get_beam_area_in_pixels(final_hdr), 'Mass': DynMass,\
            'Corruption': 'Gaussian' }
    del final_mask ; del final_cube; del final_hdr
    write_overview_file(f"{galaxy_dir}{dirstring}-Info.txt", Galaxy_Template,Achieved,required_noise)
    # We need to make the model input
    with open(f"{galaxy_dir}ModelInput.def", 'w') as tri:
        tri.writelines([Galaxy_Template['Shifted_TRM_Model'][key] + "\n" for key in Galaxy_Template['Shifted_TRM_Model']])

    # And an overview plot
    #print("Start plotting")

    #Profiles are fitted by an exponential with scale length 0.2*RHI i.e. SigHI = C*exp(-R/(0.2*RHI)) so that weay we get the Warp end at Sigma = 0.5
    Warp = [0,np.log(0.5*np.exp(-1./0.2))*-0.2*R_HI[0]]
    cf.plot_input(galaxy_dir,Galaxy_Template['Shifted_TRM_Model']
                ,Title=f'{Galaxy_Template["Name"]} with {Galaxy_Template["Beams"]:.2f} Beams'
                ,Distance= Galaxy_Template['Final_Distance'],RHI=R_HI,WarpR=Warp,
                font_file = font_file )
    # And a file with scrambled initial estimates
    cf.scrambled_initial(galaxy_dir,Galaxy_Template['Shifted_TRM_Model'])
    return_line = f"{Galaxy_Template['Final_Distance']:.2f}|{dirstring}|Convolved_Cube_Gauss\n"
    del Galaxy_Template
    return return_line

create_final_cube.__doc__= f'''
NAME:
   create_final_cube(required_noise, main_directory, Galaxy_Template,font_file = 'empty.ttf')

PURPOSE:
    Add the requested noise and produce the output products

CATEGORY:
   ROC

INPUTS:
    required_noise = the requested noise
    main_directory = 'The main directory to branch from'
    Galaxy_Template = dictionary as created by beam_template
    font_file = location of the font to be used
OPTIONAL INPUTS:

OUTPUTS:
    Cube, ModelInput.def, Initial Estimates and Overview plot.

OPTIONAL OUTPUTS:    exit()

PROCEDURES CALLED:
   Unspecified

NOTE:
'''

def download_templates(galaxy_names,path_to_resources,mp = False, ncpu = 1,\
                        work_dir = './',sofia2 = 'sofia2'):
    if mp:
        #first we check that all templates are processed
        needed_templates = [[x,path_to_resources,work_dir,sofia2] for x in galaxy_names]
        if len(needed_templates) < ncpu:
            ncpu = len(needed_templates)
        with get_context("spawn").Pool(processes=ncpu) as pool:
            results = pool.starmap(check_templates, needed_templates)
        del results
    else:
        for galaxy in galaxy_names:
            result = check_templates(galaxy,path_to_resources,work_dir,sofia2)
        results = []
download_templates.__doc__= f'''
NAME:
   download_templates(galaxy_names,path_to_resources,mp = False, ncpu = 1,\
                           work_dir = './',sofia2 = 'sofia2')
PURPOSE:
    make sure that all templates in galaxy names are present. If not download them.

CATEGORY:
   ROC

INPUTS:
    galaxy_names = names to check
    path_to_resources = path to the pyHIARd installation

OPTIONAL INPUTS:

    mp = False
    Use multiprocessing

    ncpu = 1
    number of cpus to use in multiprocessing

    work_dir = './'
    loacataion where to run sofia

    sofia2 = 'sofia2'
    command to call sofia

OUTPUTS:
    Cube, ModelInput.def, Initial Estimates and Overview plot.

OPTIONAL OUTPUTS:    exit()

PROCEDURES CALLED:
   Unspecified

NOTE:
'''



def extend_cube(cube_in, hdr_in, new_beam=0., Mask=None, rotation_pa=0.):
    if Mask is None:
         Mask=[-1., -1.]
    cube = copy.deepcopy(cube_in)
    hdr = copy.deepcopy(hdr_in)
    Mask_Use = copy.deepcopy(Mask)

    #caculate the buffers
    pixel_buffer = [int(3.*new_beam/abs(hdr['CDELT2']))]
    pixel_buffer.append(int(abs(np.sin(np.radians(rotation_pa))) *\
                    ((hdr['NAXIS1']+pixel_buffer[0])/2.)+3.))

    #Create an extended cube
    extended_cube = np.zeros((hdr['NAXIS3'], hdr['NAXIS2']+2*int(np.sum(pixel_buffer)),\
        hdr['NAXIS2']+2*int(np.sum(pixel_buffer))))
    #and place the cube inside
    extended_cube[:, int(np.sum(pixel_buffer)):hdr['NAXIS2'] + int(np.sum(pixel_buffer)),
        int(np.sum(pixel_buffer)):hdr['NAXIS1'] +int(np.sum(pixel_buffer))] = cube
    #Update the header
    for ax in ['1','2']:
        hdr[f'NAXIS{ax}']= hdr[f'NAXIS{ax}']+2*int(np.sum(pixel_buffer))
        hdr[f'CRPIX{ax}']= hdr[f'CRPIX{ax}']+int(np.sum(pixel_buffer))

    if np.sum(Mask_Use) == -2.:
        return extended_cube,hdr,pixel_buffer
    else:
        extended_mask = np.zeros((hdr['NAXIS3'], Mask_Use.shape[1]+2*int(np.sum(pixel_buffer)),\
            Mask_Use.shape[2]+2*int(np.sum(pixel_buffer))))
        extended_mask[:, int(np.sum(pixel_buffer)):Mask_Use.shape[1]+int(np.sum(pixel_buffer)) ,
            int(np.sum(pixel_buffer)):Mask_Use.shape[2]  +int(np.sum(pixel_buffer))] = Mask_Use
        return extended_cube,hdr,pixel_buffer,extended_mask
extend_cube.__doc__ = f'''
NAME:
     extend_cube(cube_in, hdr_in,new_beam= '0.', Mask = [-1.,-1.])

PURPOSE:
    extend the spatial axes of a cube to create a buffer for rotation or smoothing

CATEGORY:
   roc

INPUTS:
    cube_in = the cube to be extended
    hdr_in = the corrseponding header

OPTIONAL INPUTS:
    new_beam = 0.
    the new beam to accomadate for in degree

    rotation_pa= 0.
    astronomical PA for which the rotation should be accomodated

    Mask = [-1.-1.]
    A cube that should be extended in the same way

OUTPUTS:
    Extended cube and updated header
    pixel_buffer = list of the size of the buffer added [rot_buff,smooth_buffer]

    if both the rotation and beam are set this will be a tuple

OPTIONAL OUTPUTS:

PROCEDURES CALLED:
   Unspecified

NOTE:
'''

def galaxy_template(name,path_to_resources,work_directory,sofia2_call):
    galaxy_module = importlib.import_module(f'pyHIARD.Resources.Cubes.{name}.{name}')
    #Get the template
    Template_All=galaxy_module.get_data(work_directory,sofia2_call)
    Template_Header = Template_All[0].header
    Template_Cube = Template_All[0].data
    Template_All.close()
    print(f"We are processing the template {name}.")
    print(f"We are creating the Original template.")

    #We ensure that the Template is square
    #and the new size of the spatial axes
    new_size = int(np.max(Template_Cube.shape[1:]))
    # create a new area with only zeros
    tmp = np.zeros([Template_Cube.shape[0],new_size,new_size])
    # calculate the shift on each axis
    spatial_extension = []
    for i in [1,2]:
        spatial_extension.append(int((new_size-Template_Cube.shape[i])/2.))
        #Update the header
        Template_Header[f'CRPIX{3-i:d}'] = Template_Header[f'CRPIX{3-i:d}']+spatial_extension[-1]
        Template_Header[f'NAXIS{3-i:d}'] = new_size
    # place the template in our extended cube
    tmp[:,spatial_extension[0]:spatial_extension[0]+Template_Cube.shape[1],\
          spatial_extension[1]:spatial_extension[1]+Template_Cube.shape[2]] = Template_Cube[:,:,:]
    Template_Cube = copy.deepcopy(tmp)
    #close the temporary array
    del tmp
    # If we have a BPA in the observations we want to algin the template such that the BMAJ align with the DEC axis
    if Template_Header['BPA'] != 0.:
        #BPA is anticlockwise while rotate works clockwise so to rotate back means rotate by BPA
        Template_Cube =cf.rotateCube(Template_Cube,\
            (Template_Header['BPA']),\
            [Template_Header['CRPIX1'],
            Template_Header['CRPIX2']],order=1)


    # Obtain the tilted ring model that corresponds to this observation
    Model_Template = get_main_template(name,Template_Header,galaxy_module)
    # Now set a set of values that we required further down the line
    #The systemic velocity of the galaxy
    systemic =  float(Model_Template['VSYS'].split('=')[1].split()[0])
    # We want to ensure that the cube is calibrated on vsys
    Template_Header['CRPIX3'] = (systemic-Template_Header['CRVAL3']) / \
                        Template_Header['CDELT3'] + \
                            Template_Header['CRPIX3']
    #beam area in pixels
    pixels_per_beam = cf.get_beam_area_in_pixels(Template_Header)
    # DHI in arcconds
    DHI_arcsec=cf.convertskyangle(galaxy_module.galaxy_parameters['DHIkpc'],distance=galaxy_module.galaxy_parameters['Distance'],physical=True)

    #print(f"we get originally DHI {DHI_arcsec} arcsec, {galaxy_module.galaxy_parameters['DHIkpc']} kpc, {galaxy_module.galaxy_parameters['Distance']} Mpc")
    #The number of beams in the original observation corresponding to the D_HI diameter

    max_beams_across=(DHI_arcsec)/(Template_Header["BMAJ"]*3600.)
    #!!Maybe we need to place this before the rotate
    Original_Total_Flux_In = np.sum(Template_Cube[Template_Cube > 0.])/pixels_per_beam* Template_Header['CDELT3']#
    # To calculate the mean of all values in the mask is too sensitive to very small variations in the mask
    Original_Mean = cf.get_mean_flux(Template_Cube)
    #print(f'''The mean as determined from the template = {Original_Mean}
#The noise is {galaxy_module.galaxy_parameters["RMS"]}'''    )
    Original_z= np.sqrt((1+systemic/c_kms)/(1-systemic/c_kms))

    Template_Dictionary = {'Name':name,'Galaxy_Template_Cube':Template_Cube,'Galaxy_Template_Header':Template_Header,\
                           'Galaxy_TRM_Model':Model_Template,\
                           'Original_Noise': galaxy_module.galaxy_parameters["RMS"],'Original_Mean_Flux':Original_Mean,\
                           'Original_z':Original_z, 'Original_Vsys':systemic,'Max_Beams_Across':max_beams_across,\
                           'Original_Distance': galaxy_module.galaxy_parameters["Distance"],\
                           'Original_DHI_arcsec': DHI_arcsec, 'Original_Total_Flux_In': Original_Total_Flux_In,\
                           'Disclaimer':galaxy_module.place_disclaimer,\
                           'M_HI':galaxy_module.galaxy_parameters['MHI']}
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
    #Open the galaxy's module
    galaxy_resource_path = os.path.dirname(os.path.abspath(galaxy_module.__file__))
    #Check whether we have already created a homogenized template
    homogenized_template_exists = os.path.isfile(f"{galaxy_resource_path}/{name}_uniform.def")
    if homogenized_template_exists:
        #If we have a homogenized template,  just read that stuff
        Template_Model= cf.read_template_file(f"{galaxy_resource_path}/{name}_uniform.def",package_file = False)
    else:
        #If not read the template we have
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
            pass
            #print("No BMMAJ")
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
    #if we have a lot of models to create we should split up in chunks

    All_Galaxies = {}
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
                    galaxy_cube_exist = os.path.isfile(f"{galaxy_dir}Convolved_Cube_Gauss.fits")

                    if galaxy_cube_exist:
                        beams_to_produce.remove(nobeams)
                        noise_to_produce.remove(reqnoise)
                        print(f"The galaxy {dirstring} galaxy appears fully produced")
            if nobeams not in beams_to_produce:
                print(f"All galaxies with {name} and {nobeams} across the major axis are  thought to be produced already")
            else:
                combined.append([nobeams,noise_to_produce])
        if len(beams_to_produce) > 0:
            All_Galaxies[name] = {'name':name, 'iterations': combined}

    if len(All_Galaxies) == 0:
        print(f"Seems like the ROC has nothing to do.")
        return

    templates = [x for x in All_Galaxies]
    download_templates(templates,path_to_resources, mp = cfg.general.multiprocessing,\
                        ncpu = cfg.general.ncpu,work_dir = cfg.general.main_directory,\
                        sofia2 = cfg.general.sofia2)
    if cfg.general.multiprocessing:
        no_template_processes,no_beam_processes,no_snr_processes = \
            calculate_processes(All_Galaxies,cfg.roc.beams,cfg.roc.snr)
        if  no_snr_processes == 1:
            # if we can not do multiple SNR processes even there is no point in MP
            cfg.general.multiprocessing =False


    if cfg.general.multiprocessing:
        #list with templates to still process
        results = []
        while len(templates) > 0:

            #break templates up to process
            proc_templates = [[x,path_to_resources,cfg.general.main_directory,cfg.general.sofia2] for x in templates[:no_template_processes]]
            processes = cfg.general.ncpu
            if len(proc_templates) < processes:
                processes = len(proc_templates)
            #Then we setup the templates for all beam iterations we wants
            with get_context("spawn").Pool(processes=processes) as pool:
                All_Galaxy_Templates = pool.starmap(galaxy_template, proc_templates)

            del proc_templates
        #All_Galaxy_Templates returns a list of dictionaries that should be complete for constructing all the different beam templates
        #We need to make a list with all the different beams

            different_beams= []
            for galaxy in All_Galaxy_Templates:
                key = galaxy['Name']
                beams_and_noise = [x for x in All_Galaxies[key]['iterations']]
                for x in beams_and_noise:
                    if x[0] > galaxy['Max_Beams_Across']/cfg.roc.minimum_degradation_factor:
                        print(f"You requested {x[0]} beams across. However this galaxy originally only has {galaxy['Max_Beams_Across']} beams across.")
                        print(f"There needs to be a degradation in size for this code to work. Please just use the Original Cube for testing. Or use less than {galaxy['Max_Beams_Across']/cfg.roc.minimum_degradation_factor:.1f} beams across.")
                        print("Continuing to the next beam sizes. If you only want different signal to noise with minimal degradation please set the beams to -1")
                        continue
                    elif x[0] == -1.:
                        #We could not check before whether these exist
                        nobeams = galaxy['Max_Beams_Across']/cfg.roc.minimum_degradation_factor
                        noise_to_produce=[]
                        for reqnoise in cfg.roc.snr:
                            #if reqnoise == -1:
                            #    reqnoise= galaxy['Original_Mean_Flux']/galaxy['Original_Noise']
                            dirstring = f"{key}_{nobeams:.1f}Beams_{reqnoise:.1f}SNR"
                            noise_to_produce.append(reqnoise)
                            galaxy_dir =f"{cfg.general.main_directory}{dirstring}/"
                            galaxy_dir_exists = os.path.isdir(galaxy_dir)
                            if galaxy_dir_exists:
                        # Do we have a cube
                                galaxy_cube_exist = os.path.isfile(f"{galaxy_dir}Convolved_Cube_Gauss.fits")

                                if galaxy_cube_exist:
                                    noise_to_produce.remove(reqnoise)
                                    print(f"The galaxy {dirstring} galaxy appears fully produced")
                        if len(noise_to_produce) == 0:
                            print(f"All galaxies with {key} and {nobeams} across the major axis are  thought to be produced already")
                        else:
                            galaxy['Requested_SNR'] = noise_to_produce
                            different_beams.append((nobeams,galaxy,cfg.roc.max_degradation_factor,cfg.general.main_directory))
                    else:
                        galaxy['Requested_SNR'] = x[1]
                        #if -1. in galaxy['Requested_SNR']:
                        #    reqnoise = galaxy['Original_Mean_Flux']/galaxy['Original_Noise']
                        #    galaxy['Requested_SNR'] = [reqnoise if x == -1. else x for x in galaxy['Requested_SNR']]
                        different_beams.append((x[0],galaxy,cfg.roc.max_degradation_factor,cfg.general.main_directory))
            # if we have more than one template we know that all beams and SNR fit in the pool
            if len(different_beams) == 0:
                print(f"Seems like the ROC has nothing to do on this batch of beam templates.")
                templates = templates[no_template_processes:]
                continue
            del All_Galaxy_Templates

            while len(different_beams) > 0:
                current_different_beams = different_beams[:int(no_beam_processes*no_template_processes)]
                #Now Contruct all the beam templates of the current batch
                processes = cfg.general.ncpu
                if len(current_different_beams) < processes:
                    processes = len(current_different_beams)
                with get_context("spawn").Pool(processes=processes) as pool:
                    All_Beam_Templates = pool.starmap(beam_templates,current_different_beams)
                del current_different_beams
                # now setup a an array with the required noise input
                all_noise = []
                for galaxy in All_Beam_Templates:
                    for value in galaxy['Requested_SNR']:
                        all_noise.append((value,cfg.general.main_directory,galaxy,cfg.general.font_file))

                while len(all_noise) > 0:
                    current_noise= all_noise[:int(no_snr_processes*no_beam_processes*no_template_processes)]
                    processes = cfg.general.ncpu
                    if len(current_noise) < processes:
                        processes = len(current_noise)
                    with get_context("spawn").Pool(processes=processes) as pool:
                        tmp_result = pool.starmap(create_final_cube, current_noise)
                    all_noise = all_noise[int(no_snr_processes*no_beam_processes*no_template_processes):]
                    results = results+tmp_result
                del all_noise
                different_beams = different_beams[int(no_beam_processes*no_template_processes):]

            templates = templates[no_template_processes:]
    else:
        #No MP processing
        results = []
        for key in All_Galaxies:
            galaxy_template_out = galaxy_template(key,path_to_resources,cfg.general.main_directory,cfg.general.sofia2)
            for beam_n_noise in All_Galaxies[key]['iterations']:
                if beam_n_noise[0] > galaxy_template_out['Max_Beams_Across']/cfg.roc.minimum_degradation_factor:
                    print(f"You requested {beam_n_noise[0]} beams across. However this galaxy originally only has {galaxy_template_out['Max_Beams_Across']} beams across.")
                    print(f"There needs to be a degradation in size for this code to work. Please just use the Original Cube for testing. Or use less than {galaxy_template_out['Max_Beams_Across']/cfg.roc.minimum_degradation_factor:.1f} beams across.")
                    print("Continuing to the next beam sizes. If you only want different signal to noise with minimal degradation please set the beams to -1")
                    continue
                elif beam_n_noise[0] == -1.:
                    #We could not check before whether these exist
                    nobeams = galaxy_template_out['Max_Beams_Across']/cfg.roc.minimum_degradation_factor
                    noise_to_produce=[]
                    for reqnoise in cfg.roc.snr:
                        #if reqnoise == -1:
                        #    reqnoise= galaxy['Original_Mean_Flux']/galaxy['Original_Noise']
                        dirstring = f"{key}_{nobeams:.1f}Beams_{reqnoise:.1f}SNR"
                        noise_to_produce.append(reqnoise)
                        galaxy_dir =f"{cfg.general.main_directory}{dirstring}/"
                        galaxy_dir_exists = os.path.isdir(galaxy_dir)
                        if galaxy_dir_exists:
                            # Do we have a cube
                            galaxy_cube_exist = os.path.isfile(f"{galaxy_dir}Convolved_Cube_Gauss.fits")

                            if galaxy_cube_exist:
                                noise_to_produce.remove(reqnoise)
                                print(f"The galaxy {dirstring} galaxy appears fully produced")
                    if len(noise_to_produce) == 0:
                        print(f"All galaxies with {name} and {nobeams} across the major axis are  thought to be produced already")
                        continue
                    else:
                        galaxy_template_out['Requested_SNR'] = noise_to_produce
                        beam_input = [nobeams,galaxy_template_out,cfg.roc.max_degradation_factor,cfg.general.main_directory]
                else:
                    galaxy_template_out['Requested_SNR'] =beam_n_noise[1]
                    beam_input = [beam_n_noise[0],galaxy_template_out,cfg.roc.max_degradation_factor,cfg.general.main_directory]

                beam_template = beam_templates(*beam_input)
                for req_snr in beam_template['Requested_SNR']:
                    result = create_final_cube(float(req_snr),cfg.general.main_directory,beam_template,cfg.general.font_file)
                    results.append(result)


    with open(Catalogue, 'a') as cat:
        for line in results:
            if line != 'EMPTY':
                cat.write(f"{number_models}|{line}")
                number_models += 1
    print(f"The ROC created {number_models} models")





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

def smooth_and_regrid(Cube_In,hdr_In,factor=1.5,update_header=True, Mask = None,\
                        regrid = True, track_noise = None):
    Cube=copy.deepcopy(Cube_In)
    hdr=copy.deepcopy(hdr_In)
    if Mask is not None:
        Mask_Use = copy.deepcopy(Mask)
    #first calculate the required beams
    new_beam = []
    smooth_beam = []
    sig_pix = []
    for val in ['BMAJ','BMIN']:
        new_beam.append(hdr[val]*factor)
        #Which means we need to smoothe with
        smooth_beam.append(np.sqrt(new_beam[-1] ** 2 -  (hdr[val])** 2))
        # which in pixels is a dispersion of
        if val == 'BMAJ':
            #the major axis is oriented along the declination axis
            pix_size= abs(hdr['CDELT2'])
        else:
            pix_size= abs(hdr['CDELT1'])
        sig_pix.append((smooth_beam[-1] / np.sqrt(8 * np.log(2))) / abs(pix_size))

    if track_noise != None:
        #create a channel of the cube withe the noise of the input
        original_pixel_sigma=(np.array([hdr[key] for key in ['BMAJ','BMIN'] ],dtype=float)/ \
            np.abs(np.array([hdr[key] for key in ['CDELT2','CDELT1'] ],dtype=float)))/\
            (2. * np.sqrt(2. * np.log(2.)))
        achieved_noise,pixel_noise = cf.calculate_pixel_noise(track_noise,original_pixel_sigma,tolerance=0.01)
        #create a cube with this noise
        rng = np.random.default_rng()
        noise_channel = rng.normal(scale=pixel_noise, size=Cube.shape[1:])
        noise_channel = scipy.ndimage.gaussian_filter(noise_channel, sigma=original_pixel_sigma, order=0)
        initial_noise = np.std(noise_channel)
        smoothed_noise = scipy.ndimage.gaussian_filter(noise_channel, sigma=sig_pix, order=0)
        del noise_channel

    smoothed_cube = scipy.ndimage.gaussian_filter(Cube, sigma=(0, *sig_pix), order=0)
    del Cube
    pixel_factor = (new_beam[1] / 5.) / abs(hdr['CDELT2'])
    #Preserve surface brightness
    beams  = [[*new_beam],[hdr['BMAJ'],hdr['BMIN']]]
    new_pix_perbeam = []
    beamar = []
    for beam in beams:
        beamar.append(np.pi*abs((beam[0]*beam[1])**2)/(4.*np.log(2.)))
        new_pix_perbeam.append(cf.get_beam_area_in_pixels(hdr, beam=beam))
    smoothed_cube=smoothed_cube*new_pix_perbeam[0]/new_pix_perbeam[1]
    if track_noise != None:
        smoothed_noise = smoothed_noise *new_pix_perbeam[0]/new_pix_perbeam[1]
    #now regird to 5 pixels across the minor axis
    pixel_factor = (new_beam[1] / 5.) / abs(hdr['CDELT1'])
    if regrid:
        regrid_cube = cf.regrid_array(smoothed_cube, Out_Shape=((int(hdr['NAXIS3']),
            int(hdr['NAXIS2'] /  pixel_factor),int(hdr['NAXIS1'] / pixel_factor))))
        #also the mask Used
        if track_noise != None:
            regrid_noise = cf.regrid_array(smoothed_noise, Out_Shape=((int(hdr['NAXIS2'] /  pixel_factor),int(hdr['NAXIS1'] / pixel_factor))))
            new_noise = np.std(regrid_noise)
            del regrid_noise
        if Mask is not None:
            Mask_Use =  cf.regrid_array(Mask, Out_Shape=((int(hdr['NAXIS3']),
            int(hdr['NAXIS2'] /  pixel_factor),int(hdr['NAXIS1'] / pixel_factor))))
    else:
        regrid_cube = smoothed_cube
        if track_noise != None:
            new_noise = np.std(smoothed_noise)

    # We have to update the header
    if update_header:
        hdr['BMAJ']=new_beam[0]
        hdr['BMIN']=new_beam[1]
        for ext in ['1','2']:
            achieved_regrid = smoothed_cube.shape[int(ext)] / regrid_cube.shape[int(ext)]
            hdr[f'CDELT{ext}'] = hdr[f'CDELT{ext}']*achieved_regrid
            hdr[f'CRPIX{ext}'] = hdr[f'CRPIX{ext}']/achieved_regrid
            hdr[f'NAXIS{ext}'] = regrid_cube.shape[int(ext)]
    if track_noise != None:
        del smoothed_noise
    del smoothed_cube

    return_list  = [regrid_cube,hdr]
    if Mask is not None:
        return_list.append(Mask_Use)
    if track_noise != None:
        return_list.append(track_noise*new_noise/initial_noise)

    return return_list

smooth_and_regrid.__doc__ = f'''
NAME:
    smooth_and_regrid(Cube_In,hdr_In,factor=1.5,update_header=True, Mask = ['EMPTY'])

PURPOSE:
    Shift and degrade real galaxies.

CATEGORY:
    roc

INPUTS:
    Cube_In = Cube to smooth and regrid
    hdr_In = hdr of the input cube
    factor=1.5 factor to smooth the beam by

OPTIONAL INPUTS:
    update_header=True
    booelean to indicate whether the header should be updated

    Mask = ['EMPTY']
    mask that corresponds to cube that should still correspond after regridding.


OUTPUTS:
    A set of real galaxies in a setup ready for FAT fitting

OPTIONAL OUTPUTS:

PROCEDURES CALLED:
   Unspecified

NOTE:
'''

def write_overview_file(filename,Galaxy_Template,Achieved,SNR):
    Template = Galaxy_Template['Shifted_TRM_Model']
    # Then we also want to write some info about the galaxy
    RAdeg = float(Galaxy_Template['Shifted_TRM_Model']['XPOS'].split('=')[1].split()[0])
    DECdeg = float(Galaxy_Template['Shifted_TRM_Model']['YPOS'].split('=')[1].split()[0])
    RAhr, DEChr = cf.convertRADEC(RAdeg, DECdeg)
    RHI = cf.convertskyangle(Galaxy_Template['Final_DHI_arcsec'],float(Template['DISTANCE'].split('=')[1]))/2.
    with open(filename, 'w') as overview:
            overview.write(f'''# This file contains the basic parameters of this galaxy. For the radial dependencies look at Overview.png or ModelInput.def.
#{'Variable':<14s} {'Requested':<15s} {'Achieved':15s} {'Units':<15s}
{'Inclination':<15s} {'':<15s} {float(Template['INCL'].split('=')[1].split()[0]):<15.2f} {'degree':<15s}
{'PA':<15s} {'':<15s} {float(Template['PA'].split('=')[1].split()[0]):<15.2f} {'degree':<15s}
{'Sys. Velocity':<15s} {'':<15s} {float(Template['VSYS'].split('=')[1].split()[0]):<15.2f} {'km/s':<15s}
{'RA':<15s} {'':<15s} {RAhr.strip():<15s} {'':<15s}
{'Declination':<15s} {'':<15s} {DEChr.strip():<15s} {'':<15s}
{'Dispersion':<15s} {'':<15s} {f"{float(Template['SDIS'].split('=')[1].split()[0]):.2f}-{float(Template['SDIS'].split('=')[1].split()[-1]):.2f}":<15s} {'km/s':<15s}
{'Scale height':<15s} {'':<15s} {f"{float(Template['Z0'].split('=')[1].split()[0]):.2f}-{float(Template['Z0'].split('=')[1].split()[-1]):.2f}":<15s} {'arcsec':<15s}
{'HI Radius':<15s} {'':<15s} {RHI:<15.3f} {'kpc':<15s}
{'Maj Axis':<15s} {Galaxy_Template['Beams']:<15.3f} {'':15s} {'Beams'}
{'Total Mass':<15s} {'':<15s} {Achieved['Mass']:<15.1e} {'M_solar':<15s}
{'HI Mass':<15s} {'':<15s} {Galaxy_Template['M_HI']:<15.1e} {'M_solar':<15s}
{'Chan. Width':<15s} {'':<15s} {Achieved['Channel_Width']:<15.3f} {'km/s':<15.2s}
{'SNR':<15s} {SNR:<15.3f} {Achieved['SNR']:<15.3f}
{'Mean Signal':<15s} {'':<15s} {Achieved['Mean_Signal']*1000.:<15.3f} {'mJy/beam':<15s}
{'Noise':<15s} {'':<15s} {Achieved['Noise']*1000.:<15.3f} {'mJy/beam':<15s}
{'Distance':<15s} {'':<15s} {float(Template['DISTANCE'].split('=')[1]):<15.3f} {'Mpc':<15s}
{'Maj. FWHM':<15s} {'':<15s} {Achieved['Res_Beam'][0]:<15.2f} {'arcsec':<15s}
{'Min. FWHM':<15s} {'':<15s} {Achieved['Res_Beam'][1]:<15.2f} {'arcsec':<15s}
{'Beam BPA':<15s} {'':<15s} {Achieved['Res_Beam'][2]:<15.2f} {'degree':<15s}
{'Pixel per Beam':<15s} {'':<15s} {Achieved['Pixel_Beam']:<15.2f} {'pixel':<15s}
{'Corruption':<15s} {'':<15s} {Achieved['Corruption']:<15s}''')

write_overview_file.__doc__ = f'''
NAME:
   write_overview_file

PURPOSE:
    write the overview file with requested and achieved values

CATEGORY:
   ROC

INPUTS:
   Current_Galaxy = the Current Gakaxy that is processed (The requested values)
   Template = The tirific template that is used to create the galaxy
   Achieved= Copy of Current_Galaxy that contains the values measured from the cube

OPTIONAL INPUTS:

OUTPUTS:
    file with a table with the requested and achieved values
OPTIONAL OUTPUTS:

PROCEDURES CALLED:
   Unspecified

NOTE: For now the Achieved.Flare and Achieved.Warp are merely copies of the input simply assumed to be implemented correctly
'''
