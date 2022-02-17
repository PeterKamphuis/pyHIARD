

from pyHIARD.common_functions import download_cube,create_masks
from astropy.io import fits
import os


galaxy_parameters = {'Galaxy': 'ESO_223_G009', 'DHIkpc': 21.72 ,'Distance': 6.49, 'Original_Model': 'RC', 'RMS': 0.001781 , 'MHI': 1.02e+09   }


def get_data():
    '''Download the data for this galaxy and prepare the cube for usage'''
    succes= False
    outdir = os.path.dirname(os.path.abspath(__file__))
    try:
        Cube = fits.open(f"{outdir}/ESO_223_G009.fits",uint = False, do_not_scale_image_data=True,ignore_blank = True)
    except FileNotFoundError:
        url = 'https://www.atnf.csiro.au/research/LVHIS/data/LVHIS-cubes/LVHIS071.na.icln.fits'
        name = 'ESO_223_G009'
        sizes=[[8,44],[320,720],[320,720]]
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
    name = 'ESO_223_G009'
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
    disclaimer = '''---- ESO 223-G009 -----
    The cube used here comes from the LVHIS Survey. If you use this cube for scientific analysis please acknowledge it by citing the LVHIS overview paper (Koribalski et al. 2018, http://adsabs.harvard.edu/doi/10.1093/mnras/sty479) - Thank you.

    Please download the original cube from

    http://www.atnf.csiro.au/research/LVHIS/LVHIS-database.html

    LVHIS 071 na cube

    Original Cube name LVHIS071.na.icln.fits

    The Model is created by S.-H. Oh with ROTCUR and presented in Kamphuis et al. 2015 and Oh. et al 2018
    http://adsabs.harvard.edu/abs/2015MNRAS.452.3139K
    http://adsabs.harvard.edu/abs/2018MNRAS.473.3256O
'''
    with open(f'{dir_to_place}/ACKNOWLEDGE_LVHIS.txt', 'w') as file:
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
