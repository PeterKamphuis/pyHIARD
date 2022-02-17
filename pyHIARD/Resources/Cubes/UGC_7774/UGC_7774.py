

from pyHIARD.common_functions import download_cube,create_masks
from astropy.io import fits
import os

galaxy_parameters = {'Galaxy': 'UGC_7774', 'DHIkpc': 9.76 ,'Distance': 5.5, 'Original_Model': 'Tir', 'RMS': 0.00021 , 'MHI': 10**8.51  }


def get_data():
    '''Download the data for this galaxy and prepare the cube for usage'''
    succes= False
    outdir = os.path.dirname(os.path.abspath(__file__))
    try:
        Cube = fits.open(f"{outdir}/UGC_7774.fits",uint = False, do_not_scale_image_data=True,ignore_blank = True)
    except FileNotFoundError:
        url = 'https://zenodo.org/record/3715549/files/UGC7774-HR-cube.fits'
        name = 'UGC_7774'
        sizes=[[22,83],[200,310],[165,350]]
        Cube = download_cube(name,url,sizes,outdir)

    #place_disclaimer(dir_to_place)
    return Cube
get_data.__doc__=f'''
NAME:
   get_data

PURPOSE:
   Download the data for this galaxy and prepare the cube for usage

CATEGORY:
   agc

INPUTS:

OPTIONAL INPUTS:


OUTPUTS:
   Cube = Cube in astropy format

OPTIONAL OUTPUTS:

PROCEDURES CALLED:
   Unspecified

NOTE:
'''

def get_masks(dir_to_place,sofia_call='sofia2'):
    '''Get or create the masks for a galaxy'''
    name = 'UGC_7774'
    outdir = os.path.dirname(os.path.abspath(__file__))
    mask_exists = os.path.isfile(f"{outdir}/{name}_mask.fits")
    if mask_exists:
        Mask = fits.open(f"{outdir}/{name}_mask.fits", uint=False,
                               do_not_scale_image_data=True, ignore_blank=True)
    else:
        Mask = create_masks(outdir,dir_to_place,name,sofia_call=sofia_call)
        fits.writeto(f'{outdir}/{name}_mask.fits',Mask[0].data,Mask[0].header,overwrite = True)
    return Mask
get_masks.__doc__=f'''
NAME:
   get_masks

PURPOSE:
   Get or create the masks for a galaxy

CATEGORY:
   agc

INPUTS:
    dir_to_place = The directory where the galaxy is to be created.

OPTIONAL INPUTS:
    sofia_call = command name for sofia

OUTPUTS:
   Mask = the Blurring Mask for the template

OPTIONAL OUTPUTS:

PROCEDURES CALLED:
   Unspecified

NOTE:
'''

def place_disclaimer(dir_to_place):
    disclaimer = '''----- UGC 7774 -----

    is part of the HALOGAS Survey (Heald et al. 2011)
    If you make use of these data products in any publication or presentation, we kindly ask you to cite the following paper(s):

    Primary HALOGAS reference: Heald et al. (2011)
    and to include the following acknowledgement:

    This research made use of data from WSRT HALOGAS-DR1. The Westerbork Synthesis Radio Telescope is operated by ASTRON (Netherlands Institute for Radio Astronomy) with support from the Netherlands Foundation for Scientific Research NWO.

    Please download the original as the HR cube from
    https://www.astron.nl/halogas/data.php

    the model is from Kamphuis 2008 Chapter 5
    http://adsabs.harvard.edu/abs/2008PhDT........21K
    http://adsabs.harvard.edu/abs/2001ASPC..240..451V

    Original Name UGC7774-HR-cube.fits
'''
    with open(f'{dir_to_place}/ACKNOWLEDGE_HALOGAS.txt', 'w') as file:
        file.writelines(disclaimer)
place_disclaimer.__doc__=f'''
NAME:
   place_disclaimer

PURPOSE:
   Place a disclaimer about the source, models and acknowledgements in the directory.

CATEGORY:
   agc

INPUTS:
    dir_to_place = The directory where the galaxy is cretated.

OPTIONAL INPUTS:


OUTPUTS:

OPTIONAL OUTPUTS:

PROCEDURES CALLED:
   Unspecified

NOTE:
'''
