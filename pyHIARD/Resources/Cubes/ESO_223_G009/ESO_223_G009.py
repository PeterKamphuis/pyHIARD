

from pyHIARD.common_functions import download_cube,create_masks,select_emission
from astropy.io import fits
import os


galaxy_parameters = {'Galaxy': 'ESO_223_G009', 'DHIkpc': 21.72 ,'Distance': 6.49, 'Original_Model': 'RC', 'RMS': 0.001781 , 'MHI': 1.02e+09   }
class PackageError(Exception):
    pass

def get_data(work_dir,sofia_call='sofia2'):
    '''Download the data for this galaxy and prepare the cube for usage'''
    succes= False
    outdir = os.path.dirname(os.path.abspath(__file__))
    try:
        Cube = fits.open(f"{outdir}/ESO_223_G009.fits",uint = False, do_not_scale_image_data=True,ignore_blank = True)
    except FileNotFoundError:
        #url = 'https://www.atnf.csiro.au/research/LVHIS/data/LVHIS-cubes/LVHIS071.na.icln.fits
        url = 'https://github.com/PeterKamphuis/pyHIARD/raw/refs/heads/main/pyHIARD/Resources/Cubes/ESO_223_G009/ESO_223_G009_Original.fits'
        #url = ''
        name = 'ESO_223_G009'
        #sizes=[[8,44],[320,720],[320,720]]
        sizes = [[0,-1],[0,-1],[0,-1]]
        try:
            Cube = fits.open(f"{outdir}/{name}_Original.fits",uint = False, do_not_scale_image_data=True,ignore_blank = True)
        except:
                        #Cause LVHIS is somehow not available anymore
            #raise PackageError(f'''This file should been downloaded with the install.
#Please list an issue on the Github.''')
            Cube = download_cube(f'{name}_Original',url,sizes,outdir)
        Clean_Cube,hdr = select_emission(Cube[0].data,Cube[0].header,name,work_dir,sofia_call=sofia_call)
        fits.writeto(f"{outdir}/{name}.fits",Clean_Cube,hdr,overwrite = False)
        Cube[0].data=Clean_Cube
        Cube[0].header=hdr
        #if url != '':
        #    os.system(f"rm -f {outdir}/{name}_Original.fits")
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
