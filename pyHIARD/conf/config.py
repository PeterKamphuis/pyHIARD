from dataclasses import dataclass, field
#from multiprocessing import cpu_count
import psutil
from omegaconf import MISSING
from typing import List, Optional
import os
import sys


import pyHIARD
import pyHIARD.common_functions as cf

from omegaconf import OmegaConf, MissingMandatoryValue
from pyHIARD.AGC.AGC import AGC
from pyHIARD.ROC.ROC import ROC, add_template, remove_template,download_templates

from pyHIARD.AGC.base_galaxies import Base_Galaxy
from pyHIARD.Resources import Cubes as cubes


#The default total database currently make 229 galaxies
@dataclass
class AGC:
    #The default database currently make 149 galaxies
    enable: bool = True
    delete_existing: bool = False  # Delete all models already existing in the directory?
    # 1-6 6 creates user galaxy by asking questions
    base_galaxies: List = field(
        default_factory=lambda: [1, 2, 3, 4, 5, 6])
    inhomogenous: bool = True  # Add homogenieties?
    symmetric: bool = False  # Keep galaxies symmetric
    # options are Casa_Sim, Gaussian, No_Corrupt, Tres, Casa_5
    corruption_method: str = 'Tres'
    # The channel dependency
    # 'Options are independent, sinusoidal, hanning
    channel_dependency: str = 'independent'
    retain_unconvolved_model: bool = False
    #Produce and retain the simobserve graphics. This will crash the code when running in screen as it needs to connect to the local host
    sim_observe_graphics: bool = False
    variables_to_vary: List = field(default_factory=lambda: ['Base', 'Inclination', 'Beams', 'Radial_Motions',
                                         'Flare', 'Arms', 'Bar', 'Mass', 'Channelwidth', 'SNR', 'Warp', 'Mass', 'Beam_Size'])
    # Each base is created with the variations in the following parameters if they are listed to be varied.
    inclination: List= field(default_factory=lambda: [
                                     15., 20., 30., 50., 70., 80., 88., 90.])
    pa: List = field(default_factory=lambda: [0., 360.])
    warp: List = field(default_factory=lambda: [
                              [0.15, 0.05], [0.05, 0.2]])
    radial_motions: List = field(default_factory=lambda: [-10., -20.])
    dispersion: List = field(default_factory=lambda: [[30.0, 8.0]])
    #The flare, arms and bar will be swapped when incuded in the swap lisr
    beams:  List = field(default_factory=lambda: [
                                2., 4., 5., 6., 8., 10., 12.])
    # Beam across the major axis. This also set the distance as the size in kpc will be determined by Wang 2016 from the SBR profile
    snr: List = field(default_factory=lambda: [0.5, 1., 3., 5.])
    # These  are average signal to noise ratios
    channelwidth: List = field(default_factory=lambda: [2.])

    beam_size: List = field(default_factory=lambda: [[5., 5., 0.]])
    #Resolution of the beam in arcsec
    masses:  List = field(default_factory=lambda: [2.5e11])


@dataclass
class ROC:
    enable: bool = True
    add_template: bool = False
    remove_template: bool = False
    download_templates: bool = False 
    delete_existing: bool = False
    base_galaxies: List = field(default_factory=lambda: [
                                     'M_83', 'Circinus', 'NGC_5023', 'NGC_2903', 'NGC_3198', 'NGC_5204', 'UGC_1281', 'UGC_7774', 'ESO_223_G009'])
    variables_to_vary: List = field(
        default_factory=lambda: ['Beams', 'SNR'])
    beams: List = field(default_factory=lambda: [2., 4., 6., 8., -1.])
    # Beam across the major axis. This also set the distance as the size in kpc
    #will be determined by Wang 2016 from the SBR profile. -1 means the maximum possible for the ROC.
    # These  are average signal to noise ratios
    # The output template needs to be at least this much smaller (the exact number of beam is achieved with -1)
    minimum_degradation_factor: float = 1.25
    max_degradation_factor: float = 1.6
    # This is the maximum difference between the shifted template and the final template. If the shift factor is bigger we will smooth and regrid the template before the noise calculations.
    #A lower number  is for speed and low memory usage, a high number provides beter integration between the noise and the template. If lower than minimum_degradation_factor it will be set to the minimum_degradation_factor
    snr: List = field(default_factory=lambda: [0.5, 1., 3.])


@dataclass
class General:
    try:
        ncpu: int = len(psutil.Process().cpu_affinity())
    except AttributeError:
        ncpu: int = psutil.cpu_count()
    main_directory: str = os.getcwd()
    tirific: str = "tirific"  # Command to call tirific
    sofia2: str = "sofia"  # Command to call sofia 2
    multiprocessing: bool = True
    debug: bool = False
    font_file: str = "/usr/share/fonts/truetype/msttcorefonts/Times_New_Roman.ttf"

@dataclass
class Config:
    print_examples: bool = False
    print_bases: bool = False
    configuration_file: Optional[str] = None
    general: General = field(default_factory=General)
    agc: AGC = field(default_factory=AGC)
    roc: ROC = field(default_factory=ROC)
def help_message():
    '''The message to be printed'''
    help_message = '''
    Use pyHIARD in this way:

        pyHIARD

    will produce the standard database in the current working directory.


        pyHIARD configuration_file=pyHIARD_defaults.yml

    where configuration_file specifies a yaml file with specific settings
    such as the catalog requirements. If a required setting is left out of the yml or pyHIARD is ran without a configuration file. The code will ask about seetings without defaults.

        pyHIARD -h

    prints this message.

        pyHIARD print_examples=True

    Prints a example yml file in with the defaults in current working directory.

    '''


    return help_message

def mask_parameters(cfg_in,parameters):
    '''Mask parameters in the config file'''
    mask = []
    for key in cfg_in.__dict__['_content']:
        if not key in parameters:
            mask.append(key)
    cfg_out = OmegaConf.masked_copy(cfg_in ,\
                            mask)
    return cfg_out

def print_bases(roc_galaxies):
    '''Print the base galaxies'''
    print(f"The AGC has the following base galaxies to create variations on.")
    print(f"These are selected through a integer list in agc.base_galaxies from 1 - 5.")

    for i in range(1, 6):
        print(
            f'''{i}) Galaxy {i} has the following Base parameters to vary on.''')
        cf.print_base_galaxy(Base_Galaxy(i))

    print(f"Or you can select 6 to create your own base (1 per run)")
    print("----------------------------------------------")
    print("The ROC can create noise and size variation on the galaxies:")
    print(f"{', '.join([x for x in roc_galaxies])}")
    sys.exit()



def setup_conf():
    '''Setup and check the config'''
    #get the line input
    argv = sys.argv[1:]
    #If request version print and exit
    if '-v' in argv or '--version' in argv:
        print(f"This is version {pyHIARD.__version__} of the program.")
        sys.exit()
    #If request help print and exit
    if '-h' in argv or '--help' in argv:
        print(help_message())
        sys.exit()

    # Load the default
    cfg = OmegaConf.structured(Config)
    # Load command line arguments
    inputconf = OmegaConf.from_cli(argv)
     #combine default and line, if we do not combine here then the parameters are unset if not given
    cfg_input = OmegaConf.merge(cfg, inputconf)
    # if print_examples  print the defaults
    if cfg_input.print_examples:
        # mask the parameters that only work from command line
        cfg_masked = mask_parameters(cfg ,\
                            ['print_examples','print_bases','configuration_file'])
        with open('pyHIARD_defaults.yml', 'w') as default_write:
            default_write.write(OmegaConf.to_yaml(cfg_masked))

        print(f'''We have printed the file pyHIARD_defaults.yml in {os.getcwd()}.
''')
        sys.exit()
    #obtain the path to resources
    path_to_resources = os.path.dirname(os.path.abspath(cubes.__file__))+'/'
    #obtain the roc galaxies
    roc_galaxies = [name for name in os.listdir(
        path_to_resources) if os.path.isdir(os.path.join(path_to_resources, name))]
    roc_galaxies.remove('__pycache__')
    #print the info for the bases if required
    if cfg_input.print_bases:
        print_bases(roc_galaxies)
    #read the configuration file if requested
    if cfg_input.configuration_file:
        succes = False
        while not succes:
            try:
                yaml_config = OmegaConf.load(cfg_input.configuration_file)
        #merge yml file with defaults
                cfg = OmegaConf.merge(cfg, yaml_config)
                succes = True
            except FileNotFoundError:
                cfg_input.configuration_file = input(f'''
You have provided a config file ({cfg_input.configuration_file}) but it can't be found.
If you want to provide a config file please give the correct name.
Else press CTRL-C to abort.
configuration_file = ''')
    #overwrite with the command line input
    cfg = OmegaConf.merge(cfg, inputconf)
    # check the input
    cfg = cf.check_input(cfg)
    if cfg.roc.add_template:
        add_template(cfg, path_to_resources, roc_galaxies)
        sys.exit()
    if cfg.roc.remove_template:
        succes = False
        while not succes:
            galaxy = input(f'''Which of the following galaxies do you wat to remove?
{', '.join(roc_galaxies)} (Type exit to quit): ''')
            if galaxy.lower() == 'exit':
                succes = True
            else:
                succes = remove_template(
                    galaxy, path_to_resources, roc_galaxies)
        sys.exit()
    if cfg.roc.download_templates:
        download_templates(roc_galaxies,path_to_resources, mp = cfg.general.multiprocessing,\
                            ncpu = cfg.general.ncpu,work_dir = cfg.general.main_directory,\
                            sofia2 = cfg.general.sofia2)

        sys.exit()


    return cfg,path_to_resources



