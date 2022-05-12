#!/usr/local/bin/ python3
# This is the main program to easily create the real and false database

import os
import pyHIARD
import pyHIARD.common_functions as cf
import sys

from omegaconf import OmegaConf, MissingMandatoryValue
from pyHIARD.AGC.AGC import AGC
from pyHIARD.ROC.ROC import ROC, add_template, remove_template
from pyHIARD.conf.config import Config
from pyHIARD.AGC.base_galaxies import Base_Galaxy
from pyHIARD.Resources import Cubes as cubes

def main(argv):
    #Get default settings
    if '-v' in argv or '--version' in argv:
        print(f"This is version {pyHIARD.__version__} of the program.")
        sys.exit()

    help_message = '''
    Use pyHIARD in this way:

        pyHIARD

    will produce the standard database in the current working directory.


        pyHIARD configuration_file=FAT_Input.yml

    where configuration_file specifies a yaml file with specific settings
    such as the catalog requirements. If a required setting is left out of the yml or pyHIARD is ran without a configuration file. The code will ask about seetings without defaults.

        pyHIARD -h

    prints this message.

        pyHIARD print_examples=True

    Prints a example yml file in with the dafaults in current working directory.

    '''

    if '-h' in argv or '--help' in argv:
        print(help_message)
        sys.exit()

    cfg = OmegaConf.structured(Config)

    inputconf = OmegaConf.from_cli(argv)
    cfg_input = OmegaConf.merge(cfg, inputconf)
    if cfg_input.print_examples:
        with open('pyHIARD_defaults.yml', 'w') as default_write:
            default_write.write(OmegaConf.to_yaml(cfg))

        print(f'''We have printed the file pyHIARD_defaults.yml in {os.getcwd()}.
''')
        sys.exit()
    path_to_resources = os.path.dirname(os.path.abspath(cubes.__file__))+'/'
    roc_galaxies = [name for name in os.listdir(
        path_to_resources) if os.path.isdir(os.path.join(path_to_resources, name))]
    roc_galaxies.remove('__pycache__')
    if cfg_input.print_bases:
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

    cfg = OmegaConf.merge(cfg, inputconf)
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
    if cfg.agc.enable:
        AGC(cfg)

    if cfg.roc.enable:
        ROC(cfg,path_to_resources)
    Catalogue = f'{cfg.general.main_directory}/Output_Summary.txt'
    # If we are making new models we want to ensure this is a new file
    with open(Catalogue, 'w') as cat:
        cat.write('ID|Distance|Directoryname|Cubename\n')
        if cfg.agc.enable:
            with open(
                f'{cfg.general.main_directory}/Output_AGC_Summary.txt', 'r') as cat_AGC:
                AGCline = cat_AGC.readlines()
                for line in AGCline:
                    tmp = line.split('|')
                    try:
                        lastnum = int(tmp[0])
                    except:
                        continue
                    cat.write(line)
                lastnum += 1
        else:
            lastnum = 0.
        if cfg.roc.enable:
            with open(
                f'{cfg.general.main_directory}/Output_ROC_Summary.txt', 'r') as cat_ROC:
                ROCline = cat_ROC.readlines()
                for line in ROCline:
                    tmp = line.split('|')
                    try:
                        l = int(tmp[0])
                    except:
                        continue
                    cat.write(
                        str(int(lastnum+int(tmp[0])))+'|'+tmp[1]+'|'+tmp[2]+'|'+tmp[3])
                lastroc = int(tmp[0])+1
        else:
            lastroc = 0
    print(f"In this pyHIARD run we have created {lastnum} AGC models and {lastroc} ROC cubes. In total this makes this database {int(lastnum+lastroc)} models big")
