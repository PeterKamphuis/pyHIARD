

from pyHIARD.common_functions import download_cube,create_masks
from astropy.io import fits
import os

def get_data():
    succes= False
    outdir = os.path.dirname(os.path.abspath(__file__))
    try:
        Cube = fits.open(f"{outdir}/ESO_223_G009.fits",uint = False, do_not_scale_image_data=True,ignore_blank = True)
    except FileNotFoundError:
        url = 'https://www.atnf.csiro.au/research/LVHIS/data/LVHIS-cubes/LVHIS071.na.icln.fits'
        name = 'ESO_223_G009'
        sizes=[[0,-1],[0,-1],[0,-1]]
        Cube = download_cube(name,url,sizes,outdir)

    #place_disclaimer(dir_to_place)
    return Cube

def get_masks(dir_to_place,sofia_call='sofia2'):
    name = 'ESO_223_G009'
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
