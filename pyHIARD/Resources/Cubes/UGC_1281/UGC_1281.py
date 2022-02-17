

from pyHIARD.common_functions import download_cube,create_masks
from astropy.io import fits
import os
import copy

galaxy_parameters = {'Galaxy': 'UGC_1281', 'DHIkpc': 15.46  ,'Distance': 7.9, 'Original_Model': 'Tir', 'RMS':  0.00043 , 'MHI': 10**8.57  }


class PackageError(Exception):
    pass
def get_data():
    '''Download the data for this galaxy and prepare the cube for usage'''
    succes= False
    outdir = os.path.dirname(os.path.abspath(__file__))
    try:
        Cube = fits.open(f"{outdir}/UGC_1281.fits",uint = False, do_not_scale_image_data=True,ignore_blank = True)
    except FileNotFoundError:
        raise PackageError(f'''This file should been downloaded with the install.
Please list an issue on the Github.''')
        url = ''
        name = 'UGC_1281'
        sizes=[[12,57],[180,350],[190,330]]
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
    name = 'UGC_1281'
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
    disclaimer = '''---- UGC 1281 -----

    Is observed and modelled with TiRiFiC by P. Kamphuis.
    If you use this cube for scientific analysis please acknowledge it through a reference to Kamphuis et al. 2011
    http://adsabs.harvard.edu/abs/2011MNRAS.414.3444K
'''
    with open(f'{dir_to_place}/ACKNOWLEDGE.txt', 'w') as file:
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
