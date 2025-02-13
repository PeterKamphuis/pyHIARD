#!/usr/local/bin/ python3
# This is the main program to easily create the real and false database

from pyHIARD.AGC.AGC import AGC
from pyHIARD.ROC.ROC import ROC
from pyHIARD.conf.config import setup_conf


def main_with_trace():
    from viztracer import VizTracer
    with VizTracer(output_file="pyHIARD_Viztracer.json",min_duration=1000) as tracer:
        main()

def main():
    cfg,path_to_resources = setup_conf()
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
