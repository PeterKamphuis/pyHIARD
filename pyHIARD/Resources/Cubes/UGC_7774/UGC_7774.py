

from pyHIARD.common_functions import download_cube,create_masks
from astropy.io import fits
import os

def get_data():
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

def get_masks(dir_to_place,sofia_call='sofia2'):
    name = 'UGC_7774'
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
