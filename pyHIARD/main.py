#!/usr/local/bin/ python3
# This is the main program to easily create the real and false database

import pyHIARD.common_functions as cf
import pyHIARD.AGC.AGC as AGC
import pyHIARD.ROC.ROC as ROC
import os


def main():
    start_dir = os.getcwd()
    work_in_current = cf.get_bool(f"Do you want want to create the Database in {start_dir}? (Yes/No, default = Yes ) : ")
    if work_in_current:
        work_dir_in = start_dir
    else:
        work_dir_in = input("Please provide the directory where to create the database :")  #
        while not os.path.isdir(work_dir_in):
            print("That is not a valid directory please try again :")
            work_dir_in = input("Please provide the directory where to create the database :")
    run_default = cf.get_bool("Do you want want to create the default data base for both real and artificial galaxies? (Yes/No, default = Yes ) : ")
    if not run_default:
        run_AGC = cf.get_bool("Do you want want to create artificial galaxies? (Yes/No, default = Yes ) : ")
        if run_AGC:
            run_AGC_default = cf.get_bool("Should it be the default data base for the artificial galaxies? (Yes/No, default = Yes ) : ")

        run_ROC = cf.get_bool("Do you want want to manipulate the real observations? (Yes/No, default = Yes ) : ")
        if run_ROC:
            run_ROC_default = cf.get_bool("Do you want want to create the default data base for the real galaxies?  (Yes/No, default = Yes ): ")

        if run_AGC and not run_ROC_default:
            tmp = input("Please acknowledge that you are aware that after the artificial galaxies are created you will have to answer some questions. ")
    else:
        run_AGC = True
        run_ROC = True
        run_AGC_default = True
        run_ROC_default = True
    if run_AGC:
        AGC.AGC(work_dir=work_dir_in,running_default=run_AGC_default)
        os.chdir(start_dir)
    if run_ROC:
        ROC.ROC(work_dir=work_dir_in,running_default=run_ROC_default)
    Catalogue = work_dir_in + '/Output_Summary.txt'
    # If we are making new models we want to ensure this is a new file
    cat = open(Catalogue, 'w')
    cat.write('number|Distance|Directoryname|Cubename\n')
    cat_AGC = open(work_dir_in+'/Output_AGC_Summary.txt','r')
    AGCline = cat_AGC.readlines()
    for line in AGCline:
        tmp = line.split('|')
        try:
            lastnum=int(tmp[0])
        except:
            continue
        cat.write(line)
    lastnum += 1
    cat_AGC.close()
    cat_ROC = open(work_dir_in+'/Output_ROC_Summary.txt','r')
    ROCline = cat_ROC.readlines()
    for line in ROCline:
        tmp = line.split('|')
        try:
            l=int(tmp[0])
        except:
            continue
        cat.write(str(int(lastnum+int(tmp[0])))+'|'+tmp[1]+'|'+tmp[2]+'|'+tmp[3])
    cat_ROC.close()
    cat.close()

    cat.close()
main()
