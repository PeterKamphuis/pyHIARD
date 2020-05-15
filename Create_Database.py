#!/usr/local/bin/ python3
# This is the main program to easily create the real and false database

import common_functions as cf
import AGC
import ROC
import os


def main():
    work_dir_in = input("Please provide the directory where to create the database :")  #
    start_dir = os.getcwd()
    while not os.path.isdir(work_dir_in):
        print("That is not a valid directory please try again :")
        work_dir_in = input("Please provide the directory where to create the database :")
    run_default = cf.get_bool("Do you want want to create the default data base for both real and artificial galaxies? (Yes/No, default = Yes ) : ")
    if not run_default:
        run_AGC_default = cf.get_bool("Do you want want to create the default data base for the artificial galaxies? (Yes/No, default = Yes ) : ")
        run_ROC_default = cf.get_bool("Do you want want to create the default data base for the real galaxies? \n "
                                      "Please note that if you answer no you will have to answer some questions after the Artificial galaxies are created (Yes/No, default = Yes ) : ")
    else:
        run_AGC_default = True
        run_ROC_default = True
    AGC.AGC(work_dir=work_dir_in,running_default=run_AGC_default)
    os.chdir(start_dir)
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
