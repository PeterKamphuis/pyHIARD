#!/usr/local/bin/ python3
#Different classes used in the database
import re

# Function to input a boolean answer
def get_bool(print_str="Please type True or False",default=True):
    invalid_input = True
    while invalid_input:
        inp = input(print_str)
        if inp == "":
            if default:
                return True
            else:
                return False
        elif inp.lower() == "true" or inp.lower() == "t" or inp.lower() == "y" or inp.lower() == "yes":
            return True
        elif inp.lower() == "false" or inp.lower() == "f" or inp.lower() == "n" or inp.lower() == "no":
            return False
        else:
            print("Error: the answer must be true/false or yes/no.")
#Function for loading the variables of a tirific def file into a set of variables to be used
#This class defines a set of Base galaxy parameters
class Base_Galaxy:
    def __init__(self, num):
        if num == 1:
            self.Inclination = 60.
            self.Dispersion = [14.,7.5] #[13., 7.5]
            self.Mass= 2.5e12 #1e12
            self.PA = 35.
            self.Warp = [0.,0.]
            self.Flare = "No_Flare" #"Flare"
            self.Beams= 12. #18 #16.
            self.SNR= 8.0
            self.Channelwidth = 4.
            self.Res_Beam = [20.,20.,0.]
            self.Arms = "Arms"
            self.Bar = "No_Bar"
            self.Radial_Motions= 0.
        elif num == 2:
             # options are inclination, PA, flare, warp, beams, SNR, Channelwidth, Res_Beam, Arms, Bar, Radial_Motions
            self.Inclination = 55.
            self.Dispersion = [13.,9.] # [9., 8.]
            self.PA = 45.
            self.Warp = [0.,0.] #[0.03, 0.111] # in radians.  in Theta and phi
            self.Flare = "No_Flare"
            self.Beams= 14. #17 #16
            self.SNR= 8.
            self.Channelwidth = 4.
            self.Res_Beam = [15.,15.,0.]
            self.Arms = "Arms"
            self.Bar = "Bar"
            self.Radial_Motions= 0.
            self.Mass = 7.5e11 #5e11 # in km/s
        elif num == 3:
            self.Inclination = 65.
            self.Dispersion = [8., 8.]
            self.PA = 145.
            self.Warp = [0.05, 0.025] # in radians.
            self.Flare = "Flare"
            self.Beams= 16. #16
            self.SNR= 7.
            self.Channelwidth = 5.4
            self.Res_Beam = [25.,25.,0.]
            self.Arms = "No_Arms"
            self.Bar = "No_Bar"
            self.Radial_Motions= 0.
            self.Mass= 2.5e11
        elif num == 4:
            self.Inclination = 48.
            self.Dispersion = [13., 7.5]
            self.PA = 115.
            self.Warp = [0.07, 0.15] # in radians.
            self.Flare = "No_Flare"
            self.Beams= 15.
            self.SNR= 6.
            self.Channelwidth = 8.
            self.Res_Beam = [15.,10.,0.]
            self.Arms = "No_Arms"
            self.Bar = "Bar"
            self.Radial_Motions= 0.
            self.Mass= 7.5e10
            #self.Mass= 1e10
        elif num == 5:
            self.Inclination = 42.
            self.Dispersion = [15., 12.]
            self.PA = 115.
            self.Warp = [0.1, 0.07] # in radians.
            self.Flare = "Flare"
            self.Beams= 14.
            self.SNR= 8.
            self.Channelwidth = 6.
            self.Res_Beam = [12.5,10.,17.]
            self.Arms = "No_Arms"
            self.Bar = "No_Bar"
            self.Radial_Motions= 0.
            self.Mass= 2.5e10
        elif num == 6:
            self.Inclination = 60.
            self.Dispersion = [14.,7.5] #[13., 7.5]
            self.Mass= 7.5e9 #1e12
            self.PA = 35.
            self.Warp = [0.,0.]
            self.Flare = "No_Flare" #"Flare"
            self.Beams= 12 #18 #16.
            self.SNR= 8
            self.Channelwidth = 4.
            self.Res_Beam = [20.,20.,0.]
            self.Arms = "Arms"
            self.Bar = "No_Bar"
            self.Radial_Motions= 0.
        else:
            #Self construct a base by asking
            single = ['Inclination', 'PA', 'Beams','SNR','Channelwidth','Radial_Motions','Mass']
            double = ['Warp','Dispersion']
            triple = ['Beam_resolution']
            bool = ['Flare','Arms','Bar']
            for parameter in single:
                val = ''
                while not isinstance(val,float):
                    val = float(input(f"please provide value for {parameter} : "))
                self.__dict__[parameter] = val
            for parameter in double:
                val = ''
                while len(val) != 2:
                    vals = input(f"Please provide two and only two values for {parameter}: ")
                    tmp = re.split("\s+|\s*,\s*|\s+$",vals.strip())
                    try:
                        val = [float(tmp[0]),float(tmp[1])]
                    except:
                        val = ''
                self.__dict__[parameter] = val
            for parameter in triple:
                val=''
                while len(val) != 3:
                    vals = input(f"Please provide 1, 2 or 3 values for the {parameter}: ")
                    tmp = re.split("\s+|\s*,\s*|\s+$",vals.strip())
                    if isinstance(tmp,str):
                        tmp=[tmp]
                    try:
                        if len(tmp) == 1:
                            val = [float(tmp[0]),float(tmp[0]), 0.]
                        elif len(tmp) == 2:
                            val = [float(tmp[0]),float(tmp[1]), 0.]
                        elif len(tmp) == 3:
                            val = [float(tmp[0]),float(tmp[1]), float(tmp[2])]
                    except:
                        continue
                self.__dict__['Res_Beam'] = val
            for parameter in bool:
                add = get_bool(f"Do you want to add {parameter}? (Yes/No, default=No): ",default=False)
                if add:
                    self.__dict__[parameter] = f"{parameter}"
                else:
                    self.__dict__[parameter] = f"No_{parameter}"
