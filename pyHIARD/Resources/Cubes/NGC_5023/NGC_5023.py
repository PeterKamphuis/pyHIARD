

from pyHIARD.common_functions import download_cube,create_masks
from astropy.io import fits
import os

def get_data():
    succes= False
    outdir = os.path.dirname(os.path.abspath(__file__))
    try:
        Cube = fits.open(f"{outdir}/NGC_5023.fits",uint = False, do_not_scale_image_data=True,ignore_blank = True)
    except FileNotFoundError:
        url = 'https://zenodo.org/record/3715549/files/NGC5023-HR-cube.fits'
        name = 'NGC_5023'
        sizes=[[7,62],[420,600],[430,570]]
        Cube = download_cube(name,url,sizes,outdir)

    #place_disclaimer(dir_to_place)
    return Cube

def get_masks(dir_to_place,sofia_call='sofia2'):
    name = 'NGC_5023'
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
    disclaimer = '''----- NGC 5023 -----

    is part of the HALOGAS Survey (Heald et al. 2011)
    If you make use of these data products in any publication or presentation, we kindly ask you to cite the following paper(s):

    Primary HALOGAS reference: Heald et al. (2011)
    and to include the following acknowledgement:

    This research made use of data from WSRT HALOGAS-DR1. The Westerbork Synthesis Radio Telescope is operated by ASTRON (Netherlands Institute for Radio Astronomy) with support from the Netherlands Foundation for Scientific Research NWO.

    Please download the original HR cube from
    https://www.astron.nl/halogas/data.php

    Original Cube Name  NGC5023-HR-cube.fits

    This model is created by P.Kamphuis and published in
    HALOGAS observations of NGC 5023 and UGC 2082: Modeling of non-cylindrically symmetric gas distributions in edge-on galaxies" Kamphuis, P.; Rand, R. J.; Józsa, G. I. G.; Zschaechner, L. K.; Heald, G. H.; Patterson, M. T.; Gentile, G.; Walterbos, R. A. M.; Serra, P.; & de Blok, W. J. G. 2013, MNRAS, 434, 2069

'''
    with open(f'{dir_to_place}/ACKNOWLEDGE_HALOGAS.txt', 'w') as file:
        file.writelines(disclaimer)
