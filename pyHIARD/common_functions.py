#!/usr/local/bin/ python3
# This module contains a set of functions and classes that are used in several different Python scripts in the Database.


from collections import OrderedDict #used in Proper_Dictionary
import numpy as np # Used in convertskyangle and columndensity and
import copy # Used in columndensities
from pyHIARD import Templates as templates
from pyHIARD.Resources import Known_Models as KM
try:
    import importlib.resources as import_res
except ImportError:
    # Try backported to PY<37 `importlib_resources`.
    import importlib_resources as import_res
# A class of ordered dictionary where keys can be inserted in at specified locations or at the end.
class Proper_Dictionary(OrderedDict):
    def insert(self, existing_key, new_key, key_value):
        done = False
        if new_key in self:
            self[new_key] = key_value
            done = True
        else:
            new_orderded_dict = self.__class__()
            for key, value in self.items():
                new_orderded_dict[key] = value
                if key == existing_key:
                    new_orderded_dict[new_key] = key_value
                    done = True
            if not done:
                new_orderded_dict[new_key] = key_value
                done = True
                print(
                    "----!!!!!!!! YOUR new key was appended at the end as you provided a non-existing key to add it after!!!!!!---------")
            self.clear()
            self.update(new_orderded_dict)

        if not done:
            print("----!!!!!!!!We were unable to add your key!!!!!!---------")
#Function to convert column densities

def columndensity(levels,systemic = 100.,beam=[1.,1.],channel_width=1.,column= False,arcsquare=False,solar_mass =False):
    #set solar_mass to indicate the output should be M_solar/pc**2 or if column = True the input is
    f0 = 1.420405751786E9 #Hz rest freq
    c = 299792.458 # light speed in km / s
    pc = 3.086e+18 #parsec in cm
    solarmass = 1.98855e30 #Solar mass in kg
    mHI = 1.6737236e-27 #neutral hydrogen mass in kg

    if systemic > 10000:
        systemic = systemic/1000.
    f = f0 * (1 - (vsys / c)) #Systemic frequency
    if arcsquare:
        HIconv = 605.7383 * 1.823E18 * (2. *np.pi / (np.log(256.)))
        if column:
            # If the input is in solarmass we want to convert back to column densities
            if solar_mass:
                levels=levels*solarmass/(mHI*pc**2)
            #levels=levels/(HIconv*channel_width)
            levels = levels/(HIconv*channel_width)
        else:

            levels = HIconv*levels*channel_width
            if solar_mass:
                levels = levels/solarmass*(mHI*pc**2)
    else:
        if beam.size <2:
            beam= [beam,beam]
        b=beam[0]*beam[1]
        if column:
            if solar_mass:
                levels=levels*solarmass/(mHI*pc**2)
            TK = levels/(1.823e18*channel_width)
            levels = TK/(((605.7383)/(b))*(f0/f)**2)
        else:
            TK=((605.7383)/(b))*(f0/f)**2*levels
            levels = TK*(1.823e18*channel_width)
    if ~column and solar_mass:
        levels = levels*mHI*pc**2/solarmass
    return levels
        # a Function to convert the RA and DEC into hour angle (invert = False) and vice versa (default)
def convertRADEC(RA,DEC,invert=False, colon=False):

    if not invert:
        try:
            _ = (e for e in RA)
        except TypeError:
            RA= [RA]
            DEC =[DEC]
        for i in range(len(RA)):
            xpos=RA
            ypos=DEC
            xposh=int(np.floor((xpos[i]/360.)*24.))
            xposm=int(np.floor((((xpos[i]/360.)*24.)-xposh)*60.))
            xposs=(((((xpos[i]/360.)*24.)-xposh)*60.)-xposm)*60
            yposh=int(np.floor(np.absolute(ypos[i]*1.)))
            yposm=int(np.floor((((np.absolute(ypos[i]*1.))-yposh)*60.)))
            yposs=(((((np.absolute(ypos[i]*1.))-yposh)*60.)-yposm)*60)
            sign=ypos[i]/np.absolute(ypos[i])
            if colon:
                RA[i]="{}:{}:{:2.2f}".format(xposh,xposm,xposs)
                DEC[i]="{}:{}:{:2.2f}".format(yposh,yposm,yposs)
            else:
                RA[i]="{}h{}m{:2.2f}".format(xposh,xposm,xposs)
                DEC[i]="{}d{}m{:2.2f}".format(yposh,yposm,yposs)
            if sign < 0.: DEC[i]='-'+DEC[i]
        if len(RA) == 1:
            RA = str(RA[0])
            DEC = str(DEC[0])
    else:
        if isinstance(RA,str):
            RA=[RA]
            DEC=[DEC]

        xpos=RA
        ypos=DEC

        for i in range(len(RA)):
            # first we split the numbers out
            tmp = re.split(r"[a-z,:]+",xpos[i])
            RA[i]=(float(tmp[0])+((float(tmp[1])+(float(tmp[2])/60.))/60.))*15.
            tmp = re.split(r"[a-z,:'\"]+",ypos[i])
            DEC[i]=float(np.absolute(float(tmp[0]))+((float(tmp[1])+(float(tmp[2])/60.))/60.))*float(tmp[0])/np.absolute(float(tmp[0]))

        if len(RA) == 1:
            RA= float(RA[0])
            DEC = float(DEC[0])
    return RA,DEC


# function for converting kpc to arcsec and vice versa

def convertskyangle(angle, distance=1., unit='arcsec', distance_unit='Mpc', physical=False):
    try:
        _ = (e for e in angle)
    except TypeError:
        angle = [angle]

        # if physical is true default unit is kpc
    angle = np.array(angle)
    if physical and unit == 'arcsec':
        unit = 'kpc'
    if distance_unit.lower() == 'mpc':
        distance = distance * 10 ** 3
    elif distance_unit.lower() == 'kpc':
        distance = distance
    elif distance_unit.lower() == 'pc':
        distance = distance / (10 ** 3)
    else:
        print('CONVERTSKYANGLE: ' + distance_unit + ' is an unknown unit to convertskyangle.\n')
        print('CONVERTSKYANGLE: please use Mpc, kpc or pc.\n')
        sys.exit()
    if not physical:
        if unit.lower() == 'arcsec':
            radians = (angle / 3600.) * ((2. * np.pi) / 360.)
        elif unit.lower() == 'arcmin':
            radians = (angle / 60.) * ((2. * np.pi) / 360.)
        elif unit.lower() == 'degree':
            radians = angle * ((2. * np.pi) / 360.)
        else:
            print('CONVERTSKYANGLE: ' + unit + ' is an unknown unit to convertskyangle.\n')
            print('CONVERTSKYANGLE: please use arcsec, arcmin or degree.\n')
            sys.exit()

        kpc = 2. * (distance * np.tan(radians / 2.))
    else:
        if unit.lower() == 'kpc':
            kpc = angle
        elif unit.lower() == 'mpc':
            kpc = angle / (10 ** 3)
        elif unit.lower() == 'pc':
            kpc = angle * (10 ** 3)
        else:
            print('CONVERTSKYANGLE: ' + unit + ' is an unknown unit to convertskyangle.\n')
            print('CONVERTSKYANGLE: please use kpc, Mpc or pc.\n')
            sys.exit()
        radians = 2. * np.arctan(kpc / (2. * distance))
        kpc = (radians * (360. / (2. * np.pi))) * 3600.
    if len(kpc) == 1:
        kpc = float(kpc[0])
    return kpc

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
def load_tirific(name,Variables = ['BMIN','BMAJ','BPA','RMS','DISTANCE','NUR','RADI','VROT',
                 'Z0', 'SBR', 'INCL','PA','XPOS','YPOS','VSYS','SDIS','VROT_2',  'Z0_2','SBR_2',
                 'INCL_2','PA_2','XPOS_2','YPOS_2','VSYS_2','SDIS_2','CONDISP','CFLUX','CFLUX_2'],
                 unpack = True ):

    model = __import__(f'pyHIARD.Resources.Known_Models.{name}', globals(), locals(), name,0)
    with import_res.open_text(model, f'{name}.def') as tmp:
        unarranged = tmp.readlines()


    Variables = np.array([e.upper() for e in Variables],dtype=str)
    numrings = [int(e.split('=')[1].strip()) for e in unarranged if e.split('=')[0].strip().upper() == 'NUR']
    outputarray=np.zeros((numrings[0],len(Variables)),dtype=float)
    
    # Separate the keyword names
    for line in unarranged:
        var_concerned = str(line.split('=')[0].strip().upper())
        if len(var_concerned) < 1:
            var_concerned = 'xxx'
        varpos = np.where(Variables == var_concerned)[0]
        if varpos.size > 0:
            tmp =  np.array(line.split('=')[1].rsplit(),dtype=float)
            outputarray[0:len(tmp),int(varpos)] = tmp[0:len(tmp)]
        else:
            if var_concerned[0] == '#':
                varpos = np.where(var_concerned[2:] == Variables)[0]
                if varpos.size > 0:
                    tmp = np.array(line.split('=')[1].rsplit(),dtype=float)
                    outputarray[0:len(tmp),int(varpos)] = tmp[:]
    if unpack:
        return (*outputarray.T,)
    else:
        return outputarray


def read_template_RC(name,type= 'RC'):
    #temp = __import__('spam.ham', globals(), locals(), ['eggs', 'sausage'], 0)
    model = __import__(f'pyHIARD.Resources.Known_Models.{name}', globals(), locals(), name,0)
    with import_res.open_text(model, f'{name}.rotcur') as tmp:
        unarranged = tmp.readlines()

    Template_in = Proper_Dictionary({})
    counter = 0.
    # Separate the keyword names
    values=[]
    for tmp in unarranged:
        if tmp[0] != '!':
            range = tmp.split()
            if counter == 0:
                for item in range:
                    values.append([float(item)])
                counter = 1.
            else:
                for i,item in enumerate(range):
                    values[i].append(float(item))
    return np.array(values)



def read_template_file(filename):
    with import_res.open_text(templates, filename) as tmp:
        unarranged = tmp.readlines()
    Template_in = Proper_Dictionary({})

    # Separate the keyword names
    for tmp in unarranged:
        # python is really annoying with needing endlines. Let's strip them here and add them when writing
        Template_in[tmp.split('=',1)[0].strip().upper()]=tmp.rstrip()
    return Template_in

#Function to read simple input files that  use = as a separator between ithe required input and the values
def read_input_file(filename):


    tmpfile = open(filename, 'r')
    File = Proper_Dictionary({})
    unarranged = tmpfile.readlines()
    # Separate the keyword names
    for tmp in unarranged:
        # python is really annoying with needing endlines. Let's strip them here and add them when writing
        File[tmp.split('=', 1)[0].strip().upper()] = tmp.rstrip()
    return File
