

from pyHIARD.common_functions import download_cube,select_emission
from astropy.io import fits
import os

galaxy_parameters = {'Galaxy': 'NGC_5204', 'DHIkpc': 9.6   ,'Distance': 4.1, 'Original_Model': 'Tir', 'RMS': 0.00038  , 'MHI': 0.54e9  }


class PackageError(Exception):
    pass
def get_data(work_dir,sofia_call='sofia2'):
    '''Download the data for this galaxy and prepare the cube for usage'''
    succes= False
    outdir = os.path.dirname(os.path.abspath(__file__))
    try:
        Cube = fits.open(f"{outdir}/NGC_5204.fits",uint = False, do_not_scale_image_data=True,ignore_blank = True)
    except FileNotFoundError:
        url = 'https://github.com/PeterKamphuis/pyHIARD/raw/refs/heads/main/pyHIARD/Resources/Cubes/NGC_5204/NGC_5204_Original.fits'
        name = 'NGC_5204'
        #sizes=[[12,58],[15,180],[20,180]] 
        sizes = [[0,-1],[0,-1],[0,-1]]
        try:
            Cube = fits.open(f"{outdir}/{name}_Original.fits",uint = False, do_not_scale_image_data=True,ignore_blank = True)
        except:
            #raise PackageError(f'''This file should been downloaded with the install.
#Please list an issue on the Github.''')
            Cube = download_cube(f'{name}_Original',url,sizes,outdir)
        Clean_Cube,hdr = select_emission(Cube[0].data,Cube[0].header,name,work_dir,sofia_call=sofia_call)
        fits.writeto(f"{outdir}/{name}.fits",Clean_Cube,hdr,overwrite = False)
        Cube[0].data=Clean_Cube
        Cube[0].header=hdr
        del Clean_Cube
        del hdr
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

def place_disclaimer(dir_to_place):
    disclaimer = '''------ NGC 5204 -----
    Is observed and modelled with TiRiFiC by P. G.I.G. Jozsa.
    If you use this cube for scientific analysis please acknowledge it through a reference to Jozsa et al. 2007
    http://adsabs.harvard.edu/abs/2007A%26A...468..903J

    Original name n5204_lo.fits

    --------------------------------------------------------------
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
