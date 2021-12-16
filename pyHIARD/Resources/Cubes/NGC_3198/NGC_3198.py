

from pyHIARD.common_functions import download_cube,create_masks
from astropy.io import fits
import os

def get_data():
    succes= False
    outdir = os.path.dirname(os.path.abspath(__file__))
    try:
        Cube = fits.open(f"{outdir}/NGC_3198.fits",uint = False, do_not_scale_image_data=True,ignore_blank = True)
    except FileNotFoundError:
        url = 'https://zenodo.org/record/3715549/files/NGC3198-HR-cube.fits'
        name = 'NGC_3198'
        sizes=[[1,92],[320,730],[380,670]]
        Cube = download_cube(name,url,sizes,outdir)

    #place_disclaimer(dir_to_place)
    return Cube

def get_masks(dir_to_place,sofia_call='sofia2'):
    name = 'NGC_3198'
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
    disclaimer = '''----- NGC 3198 -----
    is part of the HALOGAS Survey (Heald et al. 2011)
    If you make use of these data products in any publication or presentation, we kindly ask you to cite the following paper(s):

    Primary HALOGAS reference: Heald et al. (2011)
    and to include the following acknowledgement:

    This research made use of data from WSRT HALOGAS-DR1. The Westerbork Synthesis Radio Telescope is operated by ASTRON (Netherlands Institute for Radio Astronomy) with support from the Netherlands Foundation for Scientific Research NWO.

    Original Cube name NGC3198-HR-Cube.fits

    Please download the HR cube from
    https://www.astron.nl/halogas/data.php

    The model is created by G. Gentile and published in
    "HALOGAS: Extraplanar gas in NGC 3198" Gentile, G.; Jozsa, G. I. G.; Serra, P.; Heald, G. H.; de Blok, W. J. G.; Fraternali, F.; Patterson, M. T.; Walterbos, R. A. M.; Oosterloo, T. 2013, A&A, 554, 125

    An additional model is by E. de Blok with ROTCUR and presented in the paper
    http://adsabs.harvard.edu/abs/2008AJ....136.2648D
'''
    with open(f'{dir_to_place}/ACKNOWLEDGE_HALOGAS.txt', 'w') as file:
        file.writelines(disclaimer)
