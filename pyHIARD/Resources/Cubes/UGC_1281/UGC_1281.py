

from pyHIARD.common_functions import download_cube,create_masks
from astropy.io import fits
import os

class PackageError(Exception):
    pass
def get_data():
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

def get_masks(dir_to_place,sofia_call='sofia2'):
    name = 'UGC_1281'
    outdir = os.path.dirname(os.path.abspath(__file__))
    try:
        for type in ['inner','outer']:
            if type == 'inner':
                Mask_Inner = fits.open(f"{outdir}/Inner_{name}_mask.fits", uint=False,
                               do_not_scale_image_data=True, ignore_blank=True)
            else:
                Mask_Outer = fits.open(f"{outdir}/Outer_{name}_mask.fits", uint=False,
                               do_not_scale_image_data=True, ignore_blank=True)
    except FileNotFoundError:
        Mask_Inner, Mask_Outer = create_masks(outdir,dir_to_place,name,sofia_call=sofia_call)
    return Mask_Inner,Mask_Outer

def place_disclaimer(dir_to_place):
    disclaimer = '''---- UGC 1281 -----

    Is observed and modelled with TiRiFiC by P. Kamphuis.
    If you use this cube for scientific analysis please acknowledge it through a reference to Kamphuis et al. 2011
    http://adsabs.harvard.edu/abs/2011MNRAS.414.3444K
'''
    with open(f'{dir_to_place}/ACKNOWLEDGE.txt', 'w') as file:
        file.writelines(disclaimer)
