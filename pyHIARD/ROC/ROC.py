#!/usr/local/bin/ python3
#This program is a python script to create a data base of real observations at different resolutions.
# It take 6 publicly available well resolved galaxy data cubes and smooths them to 3,4,5,6,7,8,10,12,16 beams across the major axis based on the extend of the best fit model.
# The galaxies used are
import importlib
import numpy as np
import sys
import copy
import pyHIARD.common_functions as cf
import scipy.ndimage
import os
import re
from astropy.io import fits
import warnings
with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    import matplotlib
    matplotlib.use('pdf')
    import matplotlib.pyplot as plt

try:
    import importlib.resources as import_res
except ImportError:
    # Try backported to PY<37 `importlib_resources`.
    import importlib_resources as import_res
from pyHIARD.Resources import Cubes as cubes


# function to properly regrid the cube after smoothing
def Regrid_Array(Array_In, Out_Shape):
    #print("Starting Regrid")
    Array_In = np.asarray(Array_In, dtype=np.double)
    In_Shape = Array_In.shape
    if len(In_Shape) != len(Out_Shape):
        print("You are regridding to different dimensions not different sizes. This won't work")
        exit()
    #print("Obtaining Shapes")
    multipliers = []  # multipliers for the final coarsegraining
    for i in range(len(In_Shape)):
        if Out_Shape[i] < In_Shape[i]:
            multipliers.append(int(np.ceil(In_Shape[i] / Out_Shape[i])))
        else:
            multipliers.append(1)
    # shape to which to blow up
    tmp_shape = tuple([i * j for i, j in zip(Out_Shape, multipliers)])

    # stupid zoom doesn't accept the final shape. Carefully crafting the
    # multipliers to make sure that it will work.
    zoomMultipliers = np.array(tmp_shape) / np.array(In_Shape) + 0.0000001
    assert zoomMultipliers.min() >= 1
    #print("Regridding")
    # applying scipy.ndimage.zoom
    regridded = scipy.ndimage.zoom(Array_In, zoomMultipliers)
    #print("Doing some other things")
    for ind, mult in enumerate(multipliers):
        if mult != 1:
            sh = list(regridded.shape)
            assert sh[ind] % mult == 0
            newshape = sh[:ind] + [sh[ind] // mult, mult] + sh[ind + 1:]
            regridded.shape = newshape
            regridded = np.mean(regridded, axis=ind + 1)
    if regridded.shape != Out_Shape:
        print("Something went wrong when regridding.")
    return regridded



def ROC(cfg):
    #First ask for the directory to work in



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
    # If we are making new models we want to ensure this is a new file
    if cfg.roc.delete_existing:
       cat = open(Catalogue, 'w')
       cat.write('ID|Distance|Directoryname|Cubename\n')
       cat.close()
    #Input Galaxy names   DHI M 83, Circinus, UGC 7774, UGC 1281 from Wang 29016 Vizier http://vizier.u-strasbg.fr/viz-bin/VizieR?-source=J/MNRAS/460/2143
    #                     DHI NGC 3198 from Gentile https://arxiv.org/pdf/1304.4232.pdf
    #                     DHI NGC 2903 From the FAT installation check initial input
    #                         NGC 5023 from Kamphuis 2013 https://arxiv.org/pdf/1306.5312.pdf
    #                         NGC 5204 from Jozsa 2007 https://www.aanda.org/articles/aa/pdf/2007/24/aa6165-06.pdf

    Galaxies_Parameters = {'Galaxies'      :['Circinus','M_83'  ,'NGC_2903','NGC_3198','NGC_5204','NGC_5023','UGC_7774','UGC_1281'],
                   'DHIkpc'        :[ 59.52    ,58.09   ,51.9      ,60.6      ,9.6,18.6,9.76,15.46],
                   'DHIDistance'   :[ 4.2      , 4.9    ,8.7       ,13.8      ,4.1,6.6,5.5,7.9],
                   'Original_Model':['RC'      ,'RC'    ,'Tir'     ,'Tir'     ,'Tir','Tir','Tir','Tir'],
                   'RMS'           :[0.012     ,0.0025  ,0.0033    ,0.00017   ,0.00038,0.0002, 0.00021 ,   0.00043],
                   'MHI'           :[10**9.83  ,10**9.91, 3.9e9    ,1.08e10   ,0.54e9,6.1e8,10**8.51,10**8.57]
                   }
    Galaxies_In ={}
    for galaxy in cfg.roc.base_galaxies:
        for key in Galaxies_Parameters:
            if key in Galaxies_In:
                Galaxies_In[key].append(Galaxies_Parameters[key][Galaxies_Parameters['Galaxies'].index(galaxy)])
            else:
                Galaxies_In[key] = [Galaxies_Parameters[key][Galaxies_Parameters['Galaxies'].index(galaxy)]]
    # The modifications requested. The original beam sizes will be maintained

    Modifications= {}
    if 'Beams' in cfg.roc.variables_to_vary:
        Modifications['Beams_Across'] = [float(x) for x in cfg.roc.beams]
    else:
        Modifications['Beams_Across'] = [-1]
    if 'SNR' in cfg.roc.variables_to_vary:
        Modifications['SNR'] = [float(x) for x in cfg.roc.snr]
    else:
        Modifications['SNR'] = [-1]
    Template_in = cf.read_template_file('Template.def')
    H_0 = 69.6 # http://www.astro.ucla.edu/~wright/CosmoCalc.html
    c = 299792.458  # Km/s
    number_models=0.
    for i in range(len(Galaxies_In['Galaxies'])):
        name=Galaxies_In['Galaxies'][i]
        print("Assembling Galaxy {}".format(name))


        galaxy_module = importlib.import_module(f'pyHIARD.Resources.Cubes.{name}.{name}')
        Template_All=galaxy_module.get_data()

        #Template_Cube=fits.open(f"{work_dir}/Input.fits",uint = False, do_not_scale_image_data=True,ignore_blank = True)
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
            Template_Header['CUNIT3'] = 'KM/S'
            Template_Header['CRVAL3'] = Template_Header['CRVAL3'] / 1000.

        #We assume the cubes to be centered
        if Galaxies_In['Original_Model'][i] == 'RC':


            radius, width, systemic, errorsys, rotation, errorrot, expansion, errorexp,\
             pa, errorpa, incli, inclerror, xpos, xposerror, ypos, yposerror, npts, sigma = cf.read_template_RC(name)
            ndisks= 1
            dispersion = np.zeros(len(radius))
            scaleheight = np.zeros(len(radius))
            dispersion2 = np.zeros(len(radius))
            scaleheight2 =  np.zeros(len(radius))
            rotation2= rotation
            pa2=  pa
            incli2= incli
            condisp = (Template_Header['CDELT3']*1.2/(2*np.sqrt(2*np.log(2.))))
            xpos2=xpos
            ypos2=ypos
            sbr = [0.]
            sbr2 = [0.]
            #Rotcur values should be converted to degrees with
            # xpos = Template_Header["CRVAL1"]  +(xpospix[:])*(Template_Header["CDELT1"])

            dispersion[:] = 0.
            dispersion2[:] = 0.
        elif Galaxies_In['Original_Model'][i] == 'Tir':
            radius, rotation, pa,incli,xpos,ypos,systemic,rotation2, pa2,incli2,xpos2,ypos2,systemic2,\
             scaleheight, dispersion, scaleheight2, dispersion2,condisp,sbr,sbr2 = cf.load_tirific(
                name , unpack=True,Variables=['RADI','VROT','PA','INCL','XPOS','YPOS','VSYS','VROT_2','PA_2','INCL_2','XPOS_2','YPOS_2','VSYS_2','Z0','SDIS','Z0_2','SDIS_2','CONDISP', 'SBR', 'SBR_2'])

            short = np.where(rotation == 0. )[0]
            if len(short) > 0.:
                for ind in short:
                    if radius[ind] != 0.:
                        rotation[ind] = rotation2[ind]
            short = np.where(rotation2 == 0. )[0]
            if len(short) > 0.:
                for ind in short:
                    if radius[ind] != 0.:
                        rotation2[ind] = rotation[ind]
            rotation = (rotation[:]+rotation2[:])/2.
            rotation2 = rotation
            if np.sum(dispersion) == 0.:
                dispersion[:] = np.sqrt(condisp[0]**2-(Template_Header['CDELT3']*1.2/(2*np.sqrt(2*np.log(2.))))**2)
                dispersion2[:] = np.sqrt(condisp[0]**2-(Template_Header['CDELT3']*1.2/(2*np.sqrt(2*np.log(2.))))**2)
                condisp = (Template_Header['CDELT3']*1.2/(2*np.sqrt(2*np.log(2.))))

        else:
            print("This type of model is not supported")
            continue
        if xpos[0] < 0.:
            xpos = [e+360. for e in xpos]
        if xpos2[0] < 0.:
            xpos2 = [e+360. for e in xpos2]


        RAdeg = xpos[0]
        DECdeg = ypos[0]
        #All values that are in arcsec need to be converted to kpc
        radius = cf.convertskyangle(radius,distance=Galaxies_In['DHIDistance'][i])
        scaleheight = cf.convertskyangle(scaleheight,distance=Galaxies_In['DHIDistance'][i])
        scaleheight2 = cf.convertskyangle(scaleheight2, distance=Galaxies_In['DHIDistance'][i])
        #Everything that is constant can be written to the Template def file
        Template_in['NUR'] = f'NUR = {len(radius)}'
        Template_in['INCL'] = 'INCL = '+" ".join(str(e) for e in incli if e != 0.)
        Template_in['INCL_2'] = 'INCL_2 = ' + " ".join(str(e) for e in incli2 if e != 0.)
        Template_in['PA'] = 'PA = ' + " ".join(str(e) for e in pa if e != 0.)
        Template_in['PA_2'] = 'PA_2 = ' + " ".join(str(e) for e in pa2 if e != 0.)
        Template_in['VROT'] = 'VROT = ' + " ".join(str(e) for e in rotation)
        Template_in['VROT_2'] = 'VROT_2 = ' + " ".join(str(e) for e in rotation2)
        Template_in['XPOS'] = 'XPOS = ' + " ".join(str(e) for e in xpos if e != 0.)
        Template_in['XPOS_2'] = 'XPOS_2 = ' + " ".join(str(e) for e in xpos2 if e != 0.)
        Template_in['YPOS'] = 'YPOS = ' + " ".join(str(e) for e in ypos if e != 0.)
        Template_in['YPOS_2'] = 'YPOS_2 = ' + " ".join(str(e) for e in ypos2 if e != 0.)
        if np.sum(dispersion) != 0.:
            Template_in['SDIS'] = 'SDIS = ' + " ".join(str(e) for e in dispersion if e != 0.)
            Template_in['SDIS_2'] = 'SDIS_2 = ' + " ".join(str(e) for e in dispersion2 if e != 0.)
        else:
            Template_in['SDIS'] = 'SDIS =  0.'
            Template_in['SDIS_2'] = 'SDIS_2 = 0. '
        if np.sum(sbr) == 0:
            Template_in['SBR'] = 'SBR = 0.'
            Template_in['SBR_2'] = 'SBR_2 = 0.'
        if np.sum(scaleheight) == 0.:
            Template_in['Z0'] = 'Z0 = 0.'
            Template_in['Z0_2'] = 'Z0_2 = 0.'
        try:
            Template_Header.set("BMAJ", Template_Header["BMMAJ"] / 3600.,before="BMMAJ")
            Template_Header.set("BMIN", Template_Header["BMMIN"] / 3600., before="BMMIN")
        except KeyError:
            print("No BMMAJ")
        bmaj=Template_Header["BMAJ"]*3600.
        bmin=Template_Header["BMIN"]*3600.
        Template_in['BMAJ'] = 'BMAJ = '+str(bmaj)
        Template_in['BMIN'] = 'BMIN = ' + str(bmin)
        Template_in['INSET'] = 'INSET = Convolved_Cube.fits'
        Template_in['CONDISP']= 'CONDISP = '+str(condisp)
        beamarea=(np.pi*abs(bmaj*bmin))/(4.*np.log(2.))
        pixperbeam=beamarea/(abs(Template_Header['CDELT1']*3600.)*abs(Template_Header['CDELT2']*3600.))


        # We need to expand the input cube with the correct noise level for the largest increase in smoothing beam
        # minbeam=np.min(beamsreq)

        DHIarcsec=cf.convertskyangle(Galaxies_In['DHIkpc'][i],distance=Galaxies_In['DHIDistance'][i],physical=True)
        #The minimum degradation is if we require more beams across the degradation is skipped
        max_beams_across=(DHIarcsec)/(bmaj)
        #the maximum smoothing will be
        #newbmaj= (DHIarcsec)/minbeam
        # Which means we extend the cube by
        #pixextendis=int(np.sqrt(newbmaj**2-bmaj**2)/abs(Template_Header['CDELT1']*3600.))
        # if we want a clean set of cubes with only Gaussian noise we check for the existence of a mask and create that if not present to run on a cube with twice the beam size

        # if we want a cube without artifact we first smooth to twice the beamsize, Use SofIa to create a mask and cut out the emission to be replaced with gaussian noise
        Mask_Inner,Mask_Outer = galaxy_module.get_masks(cfg.general.main_directory,sofia_call=cfg.general.sofia2)



        Boundary_Mask=np.array(copy.deepcopy(Mask_Outer[0].data),dtype=float)
        Boundary_Mask[Mask_Outer[0].data > 0.] = 1.
        tmp = scipy.ndimage.gaussian_filter(Boundary_Mask, sigma=(0, 3, 3), order=0)
        Boundary_Mask = copy.deepcopy(tmp)
        tmp =[]
        Boundary_Mask[Mask_Inner[0].data > 0.] = 1.
        #We mask all values that fall outside the outer mask
        Template_Cube[Boundary_Mask < 0.05] = 0.
        Template_Cube[np.isnan(Template_Cube)] = 0.
        Mask_Inner.close()
        Mask_Outer.close()
        #as we have now cut the cubes to a minimal size we want to make them square and add the size across
        # such that 2 beams still has a beam extra on all sides
        Add_Pix = abs(DHIarcsec/(Template_Header['CDELT1']*3600.))
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
        #And then calculate the total flux in the input cube
        Total_Flux_In = np.sum(Template_Cube[Template_Cube > 0.])/pixperbeam * Template_Header['CDELT3']# In Jy
        #Get the input noise
        Original_Noise = Galaxies_In["RMS"][i]
        # and the the original SNR is the mean signal in the masked cube / by the noise
        Original_Mean=np.mean(Template_Cube[Template_Cube > 0.])
        Original_z= np.sqrt((1+systemic[0]/c)/(1-systemic[0]/c))
        #Then for each shifting the flux needs to be dimmed and the corresponding noise should be calculated
        for nobeams in Modifications['Beams_Across']:
            # First we check whether the degration is more than 1.5 beams.
            Def_Template = copy.deepcopy(Template_in)
            if nobeams > max_beams_across/1.25:
                print("You requested {} beams across. However this galaxy originally only has  {} beams across.".format(nobeams,max_beams_across))
                print("There needs to be a degradation in size for this code to work. Please just use the Original Cube for testing. Or use less than {:.1f} beams across.".format(float(max_beams_across/1.3)))
                print("Continuing to the next beam sizes. If you only want different signal to noise with minimal degradation please set the beams to -1")
                continue
            if nobeams == -1:
                nobeams = max_beams_across/1.3
            # first we need to calculate the shift to  apply
            print("We are working on {} beam".format(nobeams))
            # The new beam we want = the normalsize/the beams across
            newbmaj = (DHIarcsec) / nobeams
            #which means we reduce everything by a factor of the ratio of the old and new beam
            fact = newbmaj / bmaj
            # bmin changes by the same factor
            newbmin = bmin * fact
            # And thus the new beam area is
            beamareanew = (np.pi * abs(newbmaj * newbmin)) / (4. * np.log(2.))
            # And in pixels
            pixperbeamnew = beamareanew / (abs(Template_Header['CDELT1'] * 3600.) * abs(Template_Header['CDELT2'] * 3600.))
            # We want to regrid to 5 pixels per beam which means
            # in terms of pixels
            pixsizenew = (newbmin / 5.) / abs(Template_Header['CDELT1'] * 3600.)
            # Our new distance is
            Distance = fact * Galaxies_In['DHIDistance'][i]

            Def_Template['DISTANCE'] = 'DISTANCE = '+str(Distance)
            # And we need to adjust the header for our new cube
            hednew = copy.deepcopy(Template_Header)
            # the new vsys on the hubble flow
            # First we fix the crval vaUNIlue to the actual pixel systemic is at
            New_Systemic = systemic[0] * fact
            hednew['CRPIX3'] = (systemic[0]-Template_Header['CRVAL3'])/Template_Header['CDELT3']+Template_Header['CRPIX3']-1
            hednew["CRVAL3"] = New_Systemic
            New_z = np.sqrt((1 + New_Systemic / c) / (1 - New_Systemic / c)) #This is actually 1+z
            # And we require our new pixel size
            New_Pixel_Size = (newbmaj / 5.) / abs(Template_Header['CDELT1'] * 3600.)
            # And we want to extend our input template by an amount of pixels to account for the required smoothin
            Pix_Extend = int(np.sqrt(newbmaj ** 2) / abs(Template_Header['CDELT1'] * 3600.))

            # Which means the new total flux value after adjusting the flux and the dimming factor
            Current_Flux = Total_Flux_In/(fact**2)
            # And our dimmed z factor relates to tolman surface brightness dimming
            Shifted_Template = copy.deepcopy(Template_Cube)/(fact**2)*((Original_z**4)/(New_z**4))
            Current_Mean = np.mean(Shifted_Template[Shifted_Template >0.])
            #And we can now update our model input
            conv_radi = cf.convertskyangle(radius,distance=Distance,physical= True)
            Def_Template['RADI']='RADI = '+' '.join([str(e) for e in conv_radi])
            con_hz = cf.convertskyangle(scaleheight,distance=Distance,physical= True)
            con_hz2 = cf.convertskyangle(scaleheight2, distance=Distance, physical=True)
            if np.sum(scaleheight) != 0.:
                Def_Template['Z0'] = 'Z0 = ' + " ".join(str(e) for e in con_hz if e != 0.)
                Def_Template['Z0_2'] = 'Z0_2 = ' + " ".join(str(e) for e in con_hz2 if e != 0.)
            Def_Template['VSYS'] = 'VSYS = ' + str(New_Systemic)
            Def_Template['VSYS_2'] = 'VSYS_2 = ' + str(New_Systemic)
            if np.sum(sbr) != 0.:
                # our surface brightness is constant except for the tolman dimming
                Def_Template['SBR'] = 'SBR = ' + " ".join(str(e*((Original_z**4)/(New_z**4))) for e in sbr if e != 0.)
                Def_Template['SBR_2'] = 'SBR_2 = ' + " ".join(str(e*((Original_z**4)/(New_z**4))) for e in sbr2 if e != 0.)
            # To obtain our new values we need to smooth the original with the squared difference between the old beam and the new beam
            FWHM_conv_maj = np.sqrt(newbmaj ** 2 - bmaj ** 2)
            FWHM_conv_min = np.sqrt(newbmin ** 2 - bmin ** 2)
            print("Which means we need a FWHM of {} x {}".format(FWHM_conv_maj, FWHM_conv_min))
            # and in terms of pixels the sigmas
            sig_maj = (FWHM_conv_maj / np.sqrt(8 * np.log(2))) / abs(Template_Header['CDELT1'] * 3600.)
            sig_min = (FWHM_conv_min / np.sqrt(8 * np.log(2))) / abs(Template_Header['CDELT2'] * 3600.)
            # Then we want to make the cubes for the required signal to noise ratios
            sigma = [(Template_Header["BMAJ"] / abs(Template_Header['CDELT1'])) / (2 * np.sqrt(2 * np.log(2))),
                     (Template_Header["BMIN"] / abs(Template_Header['CDELT2'])) / (2 * np.sqrt(2 * np.log(2)))]
            sigma_new = [(newbmaj / abs(Template_Header['CDELT1']*3600.)) / (2 * np.sqrt(2 * np.log(2))),
                     (newbmin / abs(Template_Header['CDELT2']*3600.)) / (2 * np.sqrt(2 * np.log(2)))]
            Ext_Template = np.zeros((
                Template_Header['NAXIS3'], Template_Header['NAXIS2'] + 2 * Pix_Extend,
                Template_Header['NAXIS1'] + 2 * Pix_Extend))
            Ext_Template[:, Pix_Extend:Template_Header['NAXIS2'] + Pix_Extend,Pix_Extend:Template_Header['NAXIS1'] + Pix_Extend] = Shifted_Template
            final_clean = scipy.ndimage.gaussian_filter(Ext_Template, sigma=(0, sig_min, sig_maj), order=0)
            # Preserve surface brightness
            final_clean = final_clean* pixperbeamnew / pixperbeam
            Final_Mean = np.mean(final_clean[final_clean > np.mean(final_clean[final_clean > 0.])/2.])
            # Let's make a mask from this smoothed cube to calculated the things we achieve
            Final_Mask =  final_clean
            Final_Mask[final_clean > np.mean(final_clean[final_clean > 0.])/2.] = 1
            Final_Mask[final_clean < np.mean(final_clean[final_clean > 0.]) / 2.] = 0.
            final_clean = []
            for reqnoise in Modifications['SNR']:
                Current_Template =copy.deepcopy(Shifted_Template)
                # To keep the same SNR as in AGC we want SNR to be mean(in)/noise
                if reqnoise == -1:
                    reqnoise = Original_Mean/Original_Noise

                print("Processing the SNR {}".format(reqnoise))
                # The noise in the final cube should be
                Final_Noise = Final_Mean/reqnoise
                # We first find the noise in pixels that matches this
                print("Creating the noise cube. The Noise in the final cube should be {} Jy/beam.".format(Final_Noise))
                New_Noise = 0.
                #As we will be adding to the uncorrected cube we need to convert the noise back to its uncorrected value.
                Final_Noise = Final_Noise * pixperbeam / pixperbeamnew

                # The pixelated noise estimated is
                Pix_Noise = (((Final_Noise * sigma_new[0] * 2 * np.sqrt(np.pi)) + (
                        Final_Noise * sigma_new[1] * 2 * np.sqrt(np.pi))) / 2.)
                # Fine tune the final noise
                while abs(New_Noise - Final_Noise) / Final_Noise > 0.025:

                    # fillthe Cube with single pixel noise
                    if New_Noise != 0:
                        #if abs(New_Noise - Final_Noise) / Final_Noise >0.05:
                        Pix_Noise = Pix_Noise / (New_Noise / Final_Noise)
                        #else:
                        #    Pix_Noise = Pix_Noise + (Final_Noise - New_Noise)
                    Ext_Template = np.random.normal(scale=Pix_Noise, size=(
                        Template_Header['NAXIS3'], Template_Header['NAXIS2'] + 2 * Pix_Extend,
                        Template_Header['NAXIS1'] + 2 * Pix_Extend))
                    Final_Template = scipy.ndimage.gaussian_filter(Ext_Template, sigma=(0, sigma_new[1], sigma_new[0]),
                                                                   order=0)
                    CHANNEL1 = Final_Template[0:2, :, :]
                    New_Noise = np.std(CHANNEL1[np.isfinite(CHANNEL1)])
                    print("We got a noise cube with an rms of {} {} {}".format(New_Noise, Final_Noise, Pix_Noise))
                # then we constuct a noise cube at the resolution of the galaxy
                Ext_Template = np.random.normal(scale=Pix_Noise, size=(
                    Template_Header['NAXIS3'], Template_Header['NAXIS2'] + 2 * Pix_Extend,
                    Template_Header['NAXIS1'] + 2 * Pix_Extend))
                Final_Template = scipy.ndimage.gaussian_filter(Ext_Template, sigma=(0, sigma[1], sigma[0]),
                                                               order=0)
                # Which means that the noise we require at this resolution is
                CHANNEL1 = Final_Template[0:2, :, :]
                Temp_Noise = np.std(CHANNEL1[np.isfinite(CHANNEL1)])
                # If this noise differs significantly from the Original noise in the cube then we want to add noise to the emission part as well.
                Diff_Noise = Temp_Noise - Original_Noise/(fact**2)*((Original_z**4)/(New_z**4))
                # If this difference is less than 10 % we will ignore it
                if abs(Diff_Noise) < Temp_Noise/10.:
                    Diff_Noise = 0.
                print("We want the noise at native resolution to be {}. And the difference noise {}".format(Temp_Noise,Diff_Noise))
                # If the new noise is smaller than the input noise we get into trouble an we do not want to continue as there would be a higher noise on the emission
                if Diff_Noise < 0:
                    print("Your requested noise is lower than the input noise hence the emission would be too noisy please lower SNR")
                    print("The requested noise is {} and the Original noise is {}".format(Temp_Noise,Original_Noise/(fact**2)))
                    print("We continue with the next SNR value")
                    continue

                # If we want to add noise  we construct a new noise cube for this purpose
                if Diff_Noise > 0.:
                    print("Creating the the difference noise cube. Shifted noise = {}.".format(Diff_Noise))
                    Pix_Noise = ((Diff_Noise * sigma[0] * 2 * np.sqrt(np.pi)) + (
                            Diff_Noise * sigma[1] * 2 * np.sqrt(np.pi))) / 2.
                    New_Noise = 0.
                    while abs(New_Noise - Diff_Noise) / Diff_Noise > 0.025:
                        # fillthe Cube with single pixel noise
                        Ext_Template = np.random.normal(scale=Pix_Noise, size=(
                            Template_Header['NAXIS3'], Template_Header['NAXIS2'],
                            Template_Header['NAXIS1']))
                        Noise_Template = scipy.ndimage.gaussian_filter(Ext_Template, sigma=(0, sigma[1], sigma[0]), order=0)
                        CHANNEL1 = Noise_Template[0:2, :, :]
                        New_Noise = np.std(CHANNEL1[np.isfinite(CHANNEL1)])
                        print("We got a difference noise cube with an rms of {}".format(New_Noise))
                        Pix_Noise = Pix_Noise / (New_Noise / Diff_Noise)
                    #We add this to the emission
                    Current_Template[Current_Template != 0.] = Current_Template[Current_Template != 0.] + Noise_Template[Current_Template != 0.]
                # We no longer need the extended template
                #print("Finished the noise")
                Ext_Template = []
                # We make a copy of this noise the size of our template
                Noise_Template = copy.deepcopy(Final_Template[:, Pix_Extend:Template_Header['NAXIS2'] + Pix_Extend,Pix_Extend:Template_Header['NAXIS1'] + Pix_Extend])
                # Then the overlap region should be half template half new noise to avoid edges
                #Current_Template[Boundary_Mask > 0.05] = Current_Template[Boundary_Mask > 0.05]*Boundary_Mask[Boundary_Mask> 0.05]+Noise_Template[Boundary_Mask > 0.05]*(1-Boundary_Mask[Boundary_Mask > 0.05])
                Current_Template = Current_Template * Boundary_Mask + Noise_Template * (1 - Boundary_Mask)
                # Finally we write the modified template into our extended template such that we can smooth it

                Final_Template[:, Pix_Extend:Template_Header['NAXIS2'] + Pix_Extend,Pix_Extend:Template_Header['NAXIS1'] + Pix_Extend] = Current_Template
                #print("Starting to smooth final")
                final = scipy.ndimage.gaussian_filter(Final_Template, sigma=(0, sig_min, sig_maj), order=0)
                #print("Finished to smooth final")
                # Preserve brightness temperature means to scale to the new area
                final = final * pixperbeamnew / pixperbeam
                Achieved_SNR = np.mean(final[Final_Mask > 0])/np.std(final[0:2,:,:])
                Achieved_Noise = np.std(final[0:2,:,:])
                Def_Template['RMS'] = 'RMS = ' + str(Achieved_Noise)
                Achieved_Mean = np.mean(final[Final_Mask > 0])
                # Regrid to 5 pix per beam
                print(" Which results in the new dimensions {} x {}".format(int(Template_Header['NAXIS2'] / pixsizenew),
                                                                            int(Template_Header['NAXIS1'] / pixsizenew)))
                regrid = Regrid_Array(final, Out_Shape=(
                    int(Template_Header['NAXIS3']), int(Template_Header['NAXIS2'] / pixsizenew),
                    int(Template_Header['NAXIS1'] / pixsizenew)))
                #also the mask Used
                regrid_mask = Regrid_Array(Final_Mask, Out_Shape=(
                    int(Template_Header['NAXIS3']), int(Template_Header['NAXIS2'] / pixsizenew),
                    int(Template_Header['NAXIS1'] / pixsizenew)))
                #print("Finished Regridding")
                # We have to update the header
                achieved = final.shape[1] / regrid.shape[1]
                hednew["CDELT1"] = Template_Header['CDELT1'] * achieved/fact
                hednew["CDELT2"] = Template_Header['CDELT2'] * achieved/fact
                hednew["CRPIX1"] = (Template_Header['CRPIX1']+Pix_Extend) / achieved
                hednew["CRPIX2"] = (Template_Header['CRPIX2']+Pix_Extend) / achieved
                # Make a new directory
                dirstring = "{}_{:.1f}Beams_{:.1f}SNR".format(name,nobeams,reqnoise)
                galaxy_dir =f"{cfg.general.main_directory}{dirstring}/"
                galaxy_dir_exists = os.path.isdir(galaxy_dir)
                if not galaxy_dir_exists:
                    os.system(f"mkdir {galaxy_dir}")
                else:
                    # Do we have a cube
                    galaxy_cube_exist = os.path.isfile(f"{galaxy_dir}Convolved_Cube.fits")

                    if galaxy_cube_exist:
                        print("This galaxy appears fully produced")
                        continue
                    else:
                        print("The directory was made but there is no full cube avalaible")
                        print("Reproducing the galaxy. Be aware of Double Table entries")
                galaxy_module.place_disclaimer(f"{galaxy_dir}")
                #ANd write to our directory
                #print("Start writing")
                fits.writeto(f"{galaxy_dir}Convolved_Cube.fits", regrid, hednew,
                             overwrite=True)#ANd write to our directory
                fits.writeto(f"{galaxy_dir}mask.fits", regrid_mask, hednew,
                             overwrite=True)
                #print("Finished writing")
                # Then we also want to write some info about the galaxy
                overview = open(f"{galaxy_dir}{dirstring}-Info.txt", 'w')
                overview.write("This file contains the basic parameters of this galaxy\n")
                overview.write("For the radial dependencies look at Overview.png or ModelInput.def\n")
                overview.write("Inclination = {}\n".format(incli[0]))
                overview.write("The dispersion = {:.2f}-{:.2f}\n".format(dispersion[0], dispersion[-1]))
                overview.write("The type of galaxy = {}\n".format(name))
                overview.write("PA = {}\n".format(pa[0]))
                #overview.write("Warp = {}-{}\n".format(Current_Galaxy.Warp[0], Current_Galaxy.Warp[1]))
                #overview.write(
                #    "Which starts at {:.2f} kpc and the 1M/pc^2 radius is {:.2f} kpc \n".format(WarpStart, Rad_HI))
                #overview.write("Flare = {}\n".format(Current_Galaxy.Flare))
                overview.write("Beams across the major axis = {}\n".format(nobeams))
                overview.write("SNR Requested = {} SNR Achieved = {}  \n".format(reqnoise, Achieved_SNR))
                overview.write("Mean Signal = {}  \n".format(Achieved_Mean))
                overview.write("Channelwidth = {}\n".format(hednew['CDELT3']))
                overview.write("Major axis beam = {} Minor axis beam= {}\n".format(bmaj,
                                                                               bmin))
                RAhr, DEChr = cf.convertRADEC(RAdeg, DECdeg)
                #overview.write("This galaxy has {} and a {}\n".format(Current_Galaxy.Arms, Current_Galaxy.Bar))
                overview.write("It's central coordinates are RA={} DEC={} vsys={:.2f} km/s\n".format(RAhr, DEChr, New_Systemic))
                overview.write("At a Distance of {:.2f} Mpc \n".format(Distance))
                #overview.write("HI_Mass Requested {:.2e} (M_solar) and an optical h {:.2f} (kpc)\n".format(MHI, sclength))
                overview.write("HI_Mass {:.2e} (M_solar) \n".format(Galaxies_In['MHI'][i]))
                #overview.write("We have {} pix per beam \n".format(pixperbeam))
                overview.write("The cube was corrupted with the {} method \n".format('Gaussian'))
                overview.write("The final noise level is {} Jy/beam \n".format(Achieved_Noise))
                overview.write("h_z = {:.3f}-{:.3f} (kpc)".format(scaleheight[0], scaleheight[-1]))
                overview.close()
                # We need to make the model input
                with open(f"{galaxy_dir}ModelInput.def", 'w') as tri:
                    tri.writelines([Def_Template[key] + "\n" for key in Def_Template])
            
                # And an overview plot
                #print("Start plotting")
                plt.figure(2, figsize=(8, 8), dpi=100, facecolor='w', edgecolor='k')
                plt.subplot(4, 1, 4)
                labelfont = {'family': 'Times New Roman',
                             'weight': 'normal',
                             'size': 22}
                plt.rc('font', **labelfont)
                plt.xlabel('Radius (arcsec)', **labelfont)
                plt.plot(conv_radi, rotation, 'k')
                plt.plot(conv_radi, rotation, 'ko')
                ymin, ymax = plt.ylim()

                plt.ylabel('V$_{rot}$ (km s$^{-1}$)', **labelfont)

                plt.margins(x=0., y=0.)
                labelfont = {'family': 'Times New Roman',
                       'weight': 'normal',
                         'size': 18}
                plt.rc('font', **labelfont)
                plt.tick_params(
                axis='x',  # changes apply to the x-axis
                which='both',  # both major and minor ticks are affected
                direction='in',
                bottom=True,  # ticks along the bottom edge are off
                top=False,  # ticks along the top edge are off
                labelbottom=True)  # labels along the bottom edge are off
                plt.subplot(4, 1, 3)
                plt.plot(conv_radi, pa, 'k')
                plt.plot(conv_radi, pa, 'ko')
                if np.sum(pa2) != 0.:
                    plt.plot(conv_radi, pa2, 'r')
                    plt.plot(conv_radi, pa2, 'ro')
                plt.tick_params(
                    axis='x',  # changes apply to the x-axis
                    which='both',  # both major and minor ticks are affected
                    direction='in',
                    bottom=True,  # ticks along the bottom edge are off
                    top=False,  # ticks along the top edge are off
                    labelbottom=False)  # labels along the bottom edge are off
                plt.ylabel('PA ($^{\circ}$)', **labelfont)
                plt.subplot(4, 1, 2)
                plt.plot(conv_radi, incli, 'k')
                plt.plot(conv_radi, incli, 'ko')
                if np.sum(incli2) != 0.:
                    plt.plot(conv_radi, incli2, 'r')
                    plt.plot(conv_radi, incli2, 'ro')
                plt.tick_params(
                    axis='x',  # changes apply to the x-axis
                    which='both',  # both major and minor ticks are affected
                    direction='in',
                    bottom=True,  # ticks along the bottom edge are off
                    top=False,  # ticks along the top edge are off
                    labelbottom=False)  # labels along the bottom edge are off
                plt.ylabel('INCL ($^{\circ}$)', **labelfont)
                plt.subplot(4, 1, 1)
                plt.plot(conv_radi, dispersion, 'k')
                plt.plot(conv_radi, dispersion, 'ko')
                if np.sum(dispersion2) != 0.:
                    plt.plot(conv_radi, dispersion2, 'r')
                    plt.plot(conv_radi, dispersion2, 'ro')
                plt.tick_params(
                    axis='x',  # changes apply to the x-axis
                    which='both',  # both major and minor ticks are affected
                    direction='in',
                    bottom=True,  # ticks along the bottom edge are off
                    top=False,  # ticks along the top edge are off
                    labelbottom=False)  # labels along the bottom edge are off
                plt.ylabel('Disp. (km s$^{-1}$)', **labelfont)
                plt.savefig(f"{galaxy_dir}Overview_Input.png", bbox_inches='tight')
                plt.close()
                #print("Finisshed Plotting")
                #And finally we need to add an entry to our catalogue
                cat = open(Catalogue, 'a')
                cat.write('{:d}|{:.2f}|{}|Convolved_Cube\n'.format(int(number_models), Distance, dirstring))
                cat.close()
                # And a file with scrambled initial estimates
                overview = open(f"{galaxy_dir}Initial_Estimates.txt", 'w')
                overview.write("#This file contains the initial estimates \n")

                SBRprof = [sbr,sbr2]
                overview.write(
                    "#{:<15} {:<15} {:<15} {:<15} {:<15} {:<15} {:<15} {:<15} {:<15}  {:<15}\n ".format('VROT', 'INCL',
                                                                                                        'PA', 'Z0', 'SBR',
                                                                                                        'DISP', 'VRAD',
                                                                                                        'RA', 'DEC',
                                                                                                        'VSYS'))
                overview.write(
                    "#{:<15} {:<15} {:<15} {:<15} {:<15} {:<15} {:<15} {:<15} {:<15}  {:<15}\n ".format('(km/s)', 'Degree',
                                                                                                        'Degree', 'kpc',
                                                                                                        'Jy km/s/arcsec^2',
                                                                                                        'km/s', 'km/s',
                                                                                                        'Degree', 'Degree',
                                                                                                        'km/s'))
                overview.write(
                    "{:<16.2f} {:<15.2f} {:<15.2f} {:<15.3f} {:<15.7f} {:<15.2f} {:<15.2f} {:<15.5f} {:<15.5f}  {:<15.2f}\n ".format(
                        np.mean(rotation[-5:-1]) + np.random.randn(1)[0] * 10.,
                        incli[0] + np.random.randn(1)[0] * 10.,
                        pa[0] + np.random.randn(1)[0] * 3.,
                        con_hz[0],
                        np.mean(0.) + np.random.randn(1)[0] * np.max(SBRprof) / 10,
                        np.mean(dispersion) + np.random.randn(1)[0] * 4.,
                        0.,
                        np.mean(RAdeg) + np.random.randn(1)[0] * 10. / 3600.,
                        np.mean(DECdeg) + np.random.randn(1)[0] * 10. / 3600.,
                        float(New_Systemic) + np.random.randn(1)[0] * 4.))

                overview.close()
                #print(f"Finished {dirstring}")
                number_models += 1.

if __name__ == '__main__':
    ROC()
