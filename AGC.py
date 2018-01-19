#!/usr/local/bin/ python3
#This program is a pyhton script to create a data base of artificial galaxies.
#It first creates a model with tirific at a high resolution and then runs it through casa to get obtain a realistic observation.
# once a numerical list is set in length we can convert it to a numpy array in order to do operations faster.
# first we import numpy
import numpy as np
import warnings
import subprocess
import array
import scipy
import sys
import copy
from collections import OrderedDict
import scipy.ndimage as ndimage
import os
import math
import re
from scipy import interpolate
from scipy import integrate
from scipy.integrate import quad
from astropy.io import fits
import matplotlib.pyplot as plt
import matplotlib
from matplotlib.ticker import AutoMinorLocator
import matplotlib.patches as patches
#del matplotlib.font_manager.weight_dict['roman']
#matplotlib.font_manager._rebuild()

#Modules
#This class defines a set of Base galaxy parameters
class Base_Galaxy:
    def __init__(self, num):
        if num == 0:
            self.inclination = 50.
            self.dispersion = [13., 7.5]
            self.rc_name = 'Massive'
            self.rc_radii = [0.,0.88,2.22,3.11,4.00,4.89,5.78,6.67,7.56,8.44,9.33,10.22,11.11,12.00,12.89,13.78,14.67,15.56,16.44, 17.] # in kpc
            self.rc_speed = [0.,234.5,191.8,211.7,223.3,222.5,223.6,224.2,226.1,226.4,226.6,226.6,223.9,220.1,217.6,217.2,215.6,210.1,207.6,207.6] # in km/s
            self.PA = 160.
            self.warp = [0.05, 0.15]
            self.flare = "Flared"
            self.beams= 32
            self.SNR= 16
            self.Channelwidth = 8.
            self.Coord = [3.56406541347E+01,4.23492081422E+01]
            self.Res_Beam = [20.,20.]
            self.Arms = "Arms"
            self.Bar = "No_Bar"
            self.Radial_Motions= 0.         
        elif num == 1:
             # options are inclination, PA, flare, warp, beams, SNR, Channelwidth, Res_Beam, Arms, Bar, Radial_Motions
            self.inclination = 50.
            self.dispersion = [8., 8.]
            self.PA = 45
            self.warp = [0.03, 0.111] # in radians.  in Theta and phi
            self.flare = "No_Flare"
            self.beams= 20
            self.SNR= 8.
            self.Channelwidth = 4.
            self.Res_Beam = [15.,15.]
            self.Arms = "Arms"
            self.Bar = "No_Bar"
            self.Radial_Motions= 0.
            # RA and DEC in degrees, Only used in casa sim 
            self.Coord = [50.81166937467079,57.76644335595375]
            self.rc_name = 'Intermediate'
            self.rc_radii = [0.00000,0.269944,0.539889,0.809833,1.07978,1.34972,1.61967,1.88961,2.15955,2.42950,2.69944,2.96997,3.24049,3.50811,3.77864,4.04916,4.31969,4.59022,4.85783,5.12836,5.39889,5.66941,5.93994,6.20756,6.47808,6.74861,7.01913,7.28966,7.55728,7.82780,8.09833,8.36886,8.63938,8.90700,9.17753,9.44805,9.71858,9.98910,10.2567,10.5272,10.7978,11.0683,11.3388,11.6064,11.8770,12.1475,12.4180,12.6885,12.9562,13.2267,13.4972,13.7677,14.0383,14.3059,14.5764,14.8469,15.1175,15.3880,15.6556,15.9261,16.1967,16.4672,16.7377,17.0053,17.2759,17.5464,17.8169,18.0874,18.3551,18.6256,18.8961,19.1666,19.4372,19.7048,19.9753,20.2458,20.5164,20.7869,21.0545,21.3250,21.5956,21.8661,22.1366,22.4042,22.6748,22.9453,23.2158,23.4863,23.7540,24.0245,24.2950,24.5655,24.8361,25.1037,25.3742,25.6447,25.9153,26.1858,26.4534,26.7239,26.9945] # in kpc
            self.rc_speed = [0.00100000,17.4200,24.7000,30.2400,34.9000,39.0300,42.7200,46.0900,49.2400,52.1600,54.9200,57.5100,59.9700,62.3300,64.5700,66.7100,68.7800,70.7500,72.6600,74.5000,76.2800,78.0000,79.6700,81.2800,82.8500,84.3600,85.8300,87.2700,88.6600,90.0100,91.3400,92.6200,93.8700,95.0900,96.2800,97.4300,98.5700,99.6600,100.700,101.800,102.800,103.800,104.800,105.700,106.700,107.600,108.400,109.300,110.200,111.000,111.800,112.600,113.300,114.100,114.800,115.500,116.200,116.900,117.600,118.200,118.900,119.500,120.100,120.700,121.300,121.900,122.400,123.000,123.500,124.100,124.600,125.100,125.600,126.100,126.600,127.000,127.500,128.000,128.400,128.800,129.300,129.700,130.100,130.500,130.900,131.300,131.700,132.100,132.500,132.900,133.200,133.600,134.000,134.300,134.700,135.100,135.400,135.800,136.200,136.600,137.100] # in km/s
        elif num == 2:
            self.inclination = 70.
            self.dispersion = [10., 10.]
            self.PA = 145
            self.warp = [0.05, 0.025] # in radians.
            self.flare = "Flared"
            self.beams= 16
            self.SNR= 4.
            self.Channelwidth = 4.
            self.Coord = [50.81166937467079,57.76644335595375]
            self.Res_Beam = [15.,15.]
            self.Arms = "Arms"
            self.Bar = "No_Bar"
            self.Radial_Motions= 0.
            self.rc_name = 'Dwarf'
            self.rc_radii = [0.00000,0.116937,0.233874,0.350811,0.467748,0.584685,0.701622,0.818559,0.935497,1.05243,1.16937,1.28631,1.40324,1.52018,1.63712,1.75406,1.87099,1.98793,2.10487,2.22180,2.33874,2.45568,2.57262,2.68955,2.80649,2.92343,3.04036,3.15730,3.27424,3.39117,3.50811,3.62505,3.74199,3.85892,3.97586,4.09280,4.20973,4.32667,4.44361,4.56055,4.67748,4.79442,4.91136,5.02829,5.14523,5.26217,5.37911,5.49604,5.61298,5.72992,5.84685,5.96379,6.08073,6.19767,6.31460,6.43154,6.54848,6.66541,6.78235,6.89929,7.01622,7.13316,7.25010,7.36704,7.48397,7.60091,7.71785,7.83478,7.95172,8.06866,8.18560,8.30253,8.41947,8.53641,8.65334,8.77028,8.88722,9.00416,9.12109,9.23803,9.35497,9.47190,9.58884,9.70578,9.82272,9.93965,10.0566,10.1735,10.2905,10.4074,10.5243,10.6413,10.7582,10.8751,10.9921,11.1090,11.2260,11.3429,11.4598,11.5768,11.6937] # in kpc
            self.rc_speed = [0.0010000000,1.3521181,2.7061880,4.0629544,5.4086103,6.7455344,8.0742311,9.3950300,10.704199,12.002402,13.280368,14.550604,15.812838,17.047775,18.273697,19.489975,20.671967,21.839514,22.995636,24.140942,25.264465,26.370249,27.450041,28.517067,29.571920,30.596693,31.600136,32.588028,33.557526,34.503735,35.433346,36.342438,37.236427,38.115829,38.961163,39.789677,40.604469,41.394070,42.160137,42.912304,43.651333,44.362362,45.056591,45.732300,46.395130,47.044331,47.663960,48.266659,48.856632,49.425865,49.975082,50.511177,51.028709,51.537907,52.035069,52.506584,52.961025,53.407246,53.845123,54.261375,54.667686,55.064053,55.436180,55.799145,56.154839,56.499504,56.835819,57.159958,57.467716,57.768078,58.060429,58.341400,58.614807,58.878239,59.136856,59.388737,59.630547,59.861759,60.089817,60.315182,60.525566,60.734787,60.944874,61.134556,61.322639,61.511303,61.694168,61.877083,62.060680,62.222191,62.381886,62.541080,62.704517,62.864269,63.019512,63.176380,63.330254,63.480808,63.647068,63.824081,64.011932] # in km/s
           
        else :
            print("There are not {} bases".format(num))
            sys.exit()
# First we define some useful routines for converting and copying
# A routine to copy disks
def copy_disk(olddisk,newdisk,Template):
    start = 0
    startlast = 0.
    if int(Template["NDISKS"].split('=')[1]) < newdisk:
        Template["NDISKS"] = "NDISKS = {:d}".format(newdisk)
    copkeys ="Empty"
    for key in Template.keys():
        if start == 2:
            if copkeys == "Empty":
                copkeys = [key]
            else:
                copkeys.append(key)         
        if key == 'RADI':
            start += 1
            if olddisk == 1:
                start += 1
        if '_' in key:
            key_ext = key.split('_')
            if key_ext[1] == str(olddisk) and start == 1:
                start += 1
                if olddisk > 1:
                    copkeys = [key]
            elif (key_ext[1] != str(olddisk))  and start == 2:
                del copkeys[-1]
                start += 1
            if key_ext[1] == str(newdisk-1) and startlast == 0:
                startlast = 1
            if key_ext[1] == str(newdisk-1) and startlast == 1:
                last = key
            if key_ext[1] != str(newdisk-1) and startlast == 1:
                startlast +=1
        if key == 'CONDISP' and (start == 2 or startlast ==1):
            if start == 2:
                del copkeys[-1]
            startlast += 1
            start += 1
            
    for key in reversed(copkeys):
        name = key.split('_')[0]
        Template.insert(last,name+"_{:d}".format(newdisk),name+"_{:d} =".format(newdisk)+Template[key].split('=')[1])
    Template.insert("CFLUX_{:d}".format(newdisk-1),"CFLUX_{:d}".format(newdisk),"CFLUX_{:d} =".format(newdisk)+Template["CFLUX_{:d}".format(newdisk-1)].split('=')[1])    
    if '_' in copkeys[-1]:
        return copkeys[-1].split('_')[0]+"_{:d}".format(newdisk)
    else:
        return copkeys[-1]+"_{:d}".format(newdisk)
            
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
            if (sign < 0.): DEC[i]='-'+DEC[i]
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
            tmp = re.split(r"[a-z,:,',\"]+",ypos[i])
            DEC[i]=float(np.absolute(float(tmp[0]))+((float(tmp[1])+(float(tmp[2])/60.))/60.))*float(tmp[0])/np.absolute(float(tmp[0]))

        if len(RA) == 1:
            RA= float(RA[0])
            DEC = float(DEC[0])             
    return RA,DEC


# Obtaining the derivative at any give point of a function
def derivative(x,func):
    h=x/1000.
    der=(func(x+h)-func(x-h))/(2*h)
    return der


# function for converting kpc to arcsec and vice versa
def convertskyangle(angle,distance=1, unit = 'arcsec', distance_unit = 'Mpc', physical = False):
    try:
        _ = (e for e in angle)
    except TypeError:
        angle = [angle]
        

        # if physical is true default unit is kpc
    angle=np.array(angle)
    if physical and unit == 'arcsec':
        unit = 'kpc'
    if distance_unit.lower() == 'mpc' : distance=distance*10**3
    elif distance_unit.lower() == 'kpc' : distance=distance
    elif distance_unit.lower() == 'pc' : distance=distance/(10**3)
    else:
        print('CONVERTSKYANGLE: '+distance_unit+' is an unknown unit to convertskyangle.\n')
        print('CONVERTSKYANGLE: please use Mpc, kpc or pc.\n')    
        sys.exit()
    if not physical:
        if unit.lower() == 'arcsec' : radians=(angle/3600.)*((2.*np.pi)/360.) 
        elif unit.lower() == 'arcmin' : radians=(angle/60.)*((2.*np.pi)/360.)  
        elif unit.lower() == 'degree' : radians=(angle)*((2.*np.pi)/360.)
        else:
            print('CONVERTSKYANGLE: '+unit+' is an unknown unit to convertskyangle.\n')
            print('CONVERTSKYANGLE: please use arcsec, arcmin or degree.\n')    
            sys.exit()
      
        kpc=2.*(distance*np.tan(radians/2.))
    else:
        if unit.lower() == 'kpc' : kpc = angle
        elif unit.lower() == 'mpc' : kpc = angle/(10**3)
        elif unit.lower() == 'pc' : kpc = angle*(10**3)
        else:
            print('CONVERTSKYANGLE: '+unit+' is an unknown unit to convertskyangle.\n')
            print('CONVERTSKYANGLE: please use kpc, Mpc or pc.\n')    
            sys.exit()
        radians=2.*np.arctan((kpc)/(2.*distance)) 
        kpc=(radians*(360./(2.*np.pi)))*3600. 
    if len(kpc) == 1:
        kpc=float(kpc[0])
    return kpc       

# Code for creating a proper dictionary instead of this python unordered nonsense.
# This also includes an insert class.
class MyOrderedDict(OrderedDict):
     def insert(self,existing_key,new_key,key_value):
        done = False 
        if new_key in self:
            self[new_key]=key_value
            done = True
        else:
            new_orderded_dict=self.__class__()
            for key, value in self.items():
                new_orderded_dict[key]=value
                if key==existing_key:
                    new_orderded_dict[new_key]=key_value
                    done = True
            if not done:
                new_orderded_dict[new_key]=key_value
                done = True
                print("----!!!!!!!! YOUR new key was appended at the end as you provided a non-existing key to add it after!!!!!!---------")
            self.clear()
            self.update(new_orderded_dict)
            
        if not done:
            print("----!!!!!!!!We were unable to add your key!!!!!!---------")


# Then the functions that create the different components of the galaxy


# function to get the surface brightness profile based on the last element in the rotation curve            
def build_sbr_prof(num,beams,Template):
    # First we use the rotational velocity to get a magnitude
    # Avanti used 	V_Ropt = (200*(l**0.41))/(0.80 +0.49*np.log10(l) +(0.75*np.exp(-0.4*l))/(0.47 + 2.25*(l**0.4)))**0.5 (Persic & Sallucci)
    # to get a B-band magnitude but mass picking are in K-band
    # Ponomareva 2017 has
    # M^T,b,i_[3.6] = (−9.52 ± 0.32) × log(2V_flat) + 3.3 ± 0.8
    # So can we find a 3.6 scaling to HI
    # Mag_36 = -9.52*np.log10(2*Base_Galaxy(num).rc_speed[-1])+3.3
    # in luminosity L_solar
    # log(LT ,b,i[3.6] ) = (3.7 ± 0.11) × log(2Vflat) + 1.3 ± 0.3,
    #L_36 = 10**(3.7*np.log10(2.*Base_Galaxy(num).rc_speed[-1])+1.3)
    MK=-9.77*np.log10(2*Base_Galaxy(num).rc_speed[-1])+1.44
    LK=10**((MK-3.29)/(-2.5))
    # We convert that to a stellar mass
    # McGaugh 2014  M_*/Lk = 0.6
    m_star = 0.6*LK
   
    # and convert stellar mass to an HI Mass following van Driel (2016)
    # log(MH I/M?) = −0.43 log (M?) + 3.75.
    MHI = m_star*10**(-0.43*np.log10(m_star)+3.75)
    # the mass leads to the radius at which we need to hit 1 M/pc^2 from Wang (2016) in kpc
    HIrad=10**(0.506*np.log10(MHI)-3.293)/2.
    # from these we use the the prescription of Martinsson to get an HI profile.
    # As done by A. Gogate in initialparam.py
    #First we set the radii at with a 5 elements per rings out to 1.5 times HIrad 
    Rad= np.arange(0,1.5*HIrad, 1.5*HIrad/(5*beams))	# in kpc
    # get the index of HIrad
    Hiradindex= np.where((Rad > HIrad - HIrad/(5*beams)) & (Rad < HIrad + HIrad/(5*beams)))[0][0]
    a = (Rad)*10**3 				#in parsecs
    Rhi_p = HIrad*10**3
    print("This the HI Radius in kpc {}".format(HIrad))
	#std deviation or dispersion of gaussian
    s1 = 0.36*Rhi_p
	#gaussian 1
    #Sig1 = np.exp(-(a + 0.4*Rhi_p)**2/(2*(s1)**2))
    # hole in the center should be amplified because of the conversion to H_2
    # with a turn over at ~ 120 V_max where it should start adding back
    # the H_2 follows an exponential distribution similar to the stellar disk (ref?)
    # The scale length is h=Vmax^2/(0.88**2*np.pi*G*sig0) (Freeman 1970)
    # G=6.674×10−20 #km^3⋅kg^−1⋅s^−2 pc=3.086e+13 km solarmass=1.98855e30 kg 
    # let's transform pc M_sol^-1 km^2 s^-2
    G= 4.30058*10**-3
    # central surface density ~ 20*the Dynamical mass/(2*np.pi*Rhi_p)
    # velocity at RadHI
    velkpc = interpolate.interpolate.interp1d(Base_Galaxy(num).rc_radii,Base_Galaxy(num).rc_speed,fill_value="extrapolate")
    Dynmass=Rhi_p*velkpc(HIrad)**2/G
    # assume the central density is 10 times the average density? Why?
    sig0 = Dynmass/(np.pi*Rhi_p**2)*10.     #in msol/pc^2. There is no good reason for the multiplication of 10 except that it gives correct scale lengths
    #To be coherent with the flaring the central density should be
    # sig0 = Dynmass/((4./3.)*np.pi*Rhi_p**3)*2*h_z
    # However it should be a surface density so         
    print("This the Dynamical Mass {} and central density {} and the velocity {}".format(Dynmass,sig0,velkpc(HIrad))) 
    # then from Graham rederived relation in cal_scl.py 
    h_r=(-4.13422991-0.31576291*MK)*1000. 
    # We assume that at 120 km/s v_max the H_2 disk imprint on the HI disk is adequately described by the prescription of Martinsson.
    # Lower some of the disk is added back to the HI higher the central hole is more pronounced
    I_cen=((Base_Galaxy(num).rc_speed[-1]/120.)**0.5-1)
    print("This the scale length {} and central brightness {}".format(h_r,I_cen))
    # So our molecular profile is
    Exp=I_cen*np.exp(-a/h_r)
    newHI=MHI-integrate.simps((2*np.pi*a)*Exp, a)
    print("Well {} {}".format(MHI,newHI))
    HIrad=10**(0.506*np.log10(newHI)-3.293)/2.
    Rhi_p=HIrad*10**3
	#gaussian2 From Martinsson 2015
    Sig2 = np.exp(-(a - 0.4*Rhi_p)**2/(2*(s1)**2))
	#total 
    Sigma = Sig2-Exp
    Sigma[Sigma < 0.] = 0.					#for negative sigma max, does not include negative values
    #scale Sigma such that it is one at HI rad
    new = 1/Sigma[Hiradindex]
    Sigma = new * Sigma
    # get the HI Mass in the profile
    OutHIMass = integrate.simps((2*np.pi*a)*Sigma, a)
    # And check that it matches the required HI mas within 5%
    counter = 1 
    while np.absolute(MHI-OutHIMass) > MHI/20.:
        # if not rescale sigma1
        if MHI-OutHIMass > 0:
            s1 = (0.36-(0.005*counter))*Rhi_p
        else:
            s1 = (0.36+(0.005*counter))*Rhi_p
        # and recalculate the profile    
        Sig2 = np.exp(-(a - 0.4*Rhi_p)**2/(2*(s1)**2))
        Sigma = Sig2-Exp
        Sigma[Sigma < 0.] = 0.
        new = 1/Sigma[Hiradindex]
        Sigma = new * Sigma
        OutHIMass=integrate.simps((2*np.pi*a)*Sigma, a)
        counter += 1

	#final HI radial distribution by renormalisation
    print("Is this going correct?")
    print(Sigma[Hiradindex],MHI,OutHIMass)
    S = Sigma*(1.24756e+20)
    # Where A. Gogate contribution stops
    # S is column densities but tirific takes Jy * km/s/arcsec^2 so
    conv_column_arsec=605.7383*1.823E18*(2.*np.pi/(np.log(256.))) #S(mJy/beam)*conv_column_arcsec=N_HI
    sbr_prof = S/(conv_column_arsec*1000.)
    #Let's write these to the Template immediately
    # The surface brightness profile, which is still symmetric
    Template["SBR"]="SBR = "+" ".join(str(e) for e in sbr_prof)
    Template["SBR_2"]="SBR_2 = "+" ".join(str(e) for e in sbr_prof)
    # And we plot the  profile to an overview plot
    overview.plot(Rad,sbr_prof)
    overview.plot(Rad,Exp*new*1.24756e+20/(conv_column_arsec*1000.))
    plt.plot([HIrad,HIrad],[0.,np.max(sbr_prof)])  
    plt.xlabel('Radius (kpc)')
    plt.ylabel('SBR (Jy km/s arcsec^-2)')
    return sbr_prof,Rad,h_r/1000.,OutHIMass, HIrad
   

# Thi function creates the flarin g which is based on the rotation curve and dispersion according Puche et al. 1992 
def create_flare(Radii,velocity,dispersion,flare,Template,Max_Rad,distance=1.):
    #make sure we have arrays so we can do numerical operations
    Radii = np.array(Radii)
    velocity = np.array(velocity)
    # first we need the dispersion
    # if there is no change from inside to outside it is easy 
    if dispersion[0] == dispersion[1]:
       disp=np.full(len(Radii),dispersion[0])
    else:
        # else we use a tangent change with the center at halfway 
       disp=-1*np.arctan((Radii-np.mean(Max_Rad/2.))/(Radii[-1]/10.))/(np.pi)*np.absolute(dispersion[0]-dispersion[1])+np.mean(dispersion)
    # We recalculate the Dynamical Mass and densities
    # G=6.674×10−20 #km^3⋅kg^−1⋅s^−2 pc=3.086e+13 km solarmass=1.98855e30 kg 
    # let's transform pc M_sol^-1 km^2 s^-2
    G= 4.30058*10**-3 # radii are in kpc

    # Dynamical Mass
    
    Dynmass=(Radii*10**3)*velocity**2/G
    # This is in a Volume of
    Volume=(4./3.)*np.pi*Radii**3
    Dynmass[0]=1.
    Volume[0]=1.
    Density=Dynmass/Volume*2.  #  In M_solar/kpc^3 the two comes from Puche 1992 but seems random    
    # Then we check wether we want a flare or not
    G2=G/(3.086e+13**2) #pc^3 M_sol^-1 s^-2
    halfint=int((len(Radii[Radii < Max_Rad])+10)/2.)
    if flare.lower() == 'flared':
        flare =disp/((4.*np.pi*G2*Density/1000.**3)**0.5*3.086e+16) # in kpc
        flare[:halfint-10] = flare[halfint-10]
        fact=np.arange(1/20,1,1./20)
        flare[halfint-10:halfint+9] = (1-fact)*flare[halfint-10]+fact*flare[halfint-10:halfint+9]
    elif flare.lower() == 'no_flare':       
        flare = np.full(len(Radii),disp[halfint]/((4.*np.pi*G2*Density[halfint]/1000.**3)**0.5*3.086e+16))
    else:
        print("{} is not an option for the flare. Choose Flared or No_Flare".format(flare))
        sys.exit()
    
    flare[0]=flare[1]    
    plt.figure(2)
    plt.subplot(6,1,2)
    plt.plot(Rad,disp)   
    plt.ylabel('Disp. (km/s)')
    plt.xticks([])
    print("The dispersion runs from {:5.2f} to {:5.2f} km/s".format(disp[0],disp[-1]))
    plt.subplot(6,1,1)
    plt.plot(Rad,flare)   
    plt.ylabel('Scale height (kpc)')
    plt.title(Base_Galaxy(base).rc_name)
    plt.xticks([])
    print("The flare runs from {:10.6f} to {:10.6f} kpc".format(flare[0],flare[-1]) )
    # convert the scale heights to arcsec
    h_z_arcsec = convertskyangle(flare,distance=Distance,physical = True)
    # and write both to the Template
    Template["SDIS"]="SDIS = "+" ".join(str(e) for e in disp)
    Template["SDIS_2"]="SDIS_2 = "+" ".join(str(e) for e in disp)
    # The scaleheight
    Template["Z0"]="Z0 = "+" ".join(str(e) for e in h_z_arcsec)
    Template["Z0_2"]="Z0_2 = "+" ".join(str(e) for e in h_z_arcsec)
    return flare,disp



# A function for varying the PA and inclination as function of radius
def create_warp(Radii,
                PA,
                inclination,
                warp_change,
                warp_radii,
                Template,
                disk=1):
    Radii =np.array(Radii)         
    if np.sum(warp_change) != 0:
        # First we need to obtain the vector that constitutes the inner area
        #it runs exactly counter to inclination
        inclination=90-inclination
        # For this the PA has to be between 0-90
        mult=np.floor(PA/90.)
        inPA= PA-mult*90.
        #avoid singularities
        if inPA == 0.: inPA = 0.000001
        if inclination == 0: inclination = 0.000001
        # define the angular momentum vector of the plane and the outer most ring    
        theta=np.arctan(np.tan(inclination*(np.pi/180.))*np.tan(inPA*(np.pi/180.)))
        phi=np.arctan(np.tan(inPA*(np.pi/180.))/np.sin(theta))
        # and the maximum values at Rad_HI
        thetamax=theta+warp_change[0]
        phimax=phi+warp_change[1]
                # indices of the warp start and the Rad_HI
        start_index = np.sum(Radii < warp_radii[0])
        end_index = np.sum(Radii < warp_radii[1])
        # step size of theta and phi
        # As we will increase the step size triangular we need the total number of point in the sequence
        warprange = end_index-start_index
        increasetheta=(thetamax-theta)/(0.5*warprange*(warprange+1))
        increasephi=(phimax-phi)/(0.5*warprange*(warprange+1))
        # calculate theta
        thetarings = np.array(np.full(len(Radii),theta))
        index_array=np.arange(start_index,len(Radii))-start_index
        thetarings[start_index:] = theta+0.5*index_array*(index_array+1)*increasetheta
        #calculate phi
        phirings = np.array(np.full(len(Radii),phi))
        phirings[start_index:] = phi+0.5*(index_array)*(index_array+1)*increasephi
        # return to PA
        PA= np.arctan(np.sin(thetarings)*np.tan(phirings))*(360./(2*np.pi))+mult*90 
        # Pa boundary adjustments
        PA[phirings > 0.5*np.pi]=PA[phirings > 0.5*np.pi]+180.
        #PA[np.where(PA >360.)]=PA[np.where(PA >360.)]-360
        # return inclination
        inc=90-np.arctan(1./(np.cos(thetarings)*np.tan(phirings)))*(360./(2*np.pi))
        # return inclination boundary adjustements
        inc[np.where(inc > 90.)] = 180 - inc[np.where(inc > 90.)]
        inc[np.where(inc < 0.)] = -1 * inc[np.where(inc < 0.)]
         # return a correct quadrant phirings
        phirings=phirings+mult/2.*np.pi
    else:
        #if there is no warp then all values are the same
        PA =np.full(len(Radii), PA)
        inc =np.full(len(Radii), inclination)
        phirings=np.full(len(Radii),0)
        thetarings=np.full(len(Radii),0)
    
    warnings.filterwarnings("ignore")
    plt.subplot(6,1,4)
    plt.plot(Rad,PA)
    plt.ylabel('PA (deg)')
    plt.xticks([])
    plt.subplot(6,1,3)
    plt.plot(Rad,inc)
    plt.ylabel('Inc. (deg)')
    plt.xticks([])
      
    #write to our template file
    # let's see if we can retrace intrinsic phi with the formula's from Peters    
    theta_test = np.arctan((np.sin(PA*np.pi/180.)*np.sin(inc*np.pi/180.)-np.cos(PA*np.pi/180.)*np.sin(inc*np.pi/180.))/(np.cos(PA*np.pi/180.)*np.cos(inc*np.pi/180.)-np.sin(PA*np.pi/180.)*np.cos(inc*np.pi/180.)))*180./np.pi  
    phirings=phirings*180./np.pi
    thetarings=thetarings*180./np.pi
    #print("phi {} and theta {}".format((phirings-phirings[0]),(thetarings-thetarings[0])))
    # This seems to work mostly but not at some extremes exactly for some reason 
    angle_adjust=(PA[0]-PA)*np.cos(inc*np.pi/180.)
    #print("Adjust Angle {}".format(angle_adjust))
    if disk == 1:
        Template["INCL"]="INCL = "+" ".join(str(e) for e in inc)
        Template["PA"]="PA = "+" ".join(str(e) for e in PA)
        try:
            phase =  Template["AZ1P"].split('=')[1]
        except KeyError:
            phase =  0.
            Template.insert("AZ1W","AZ1P","AZ1P = "+" ".join(str(phase+(e)) for e in angle_adjust))
    else:
        Template["INCL_{:d}".format(disk)]="INCL_{:d} =".format(disk)+" ".join(str(e) for e in inc)
        Template["PA_{:d}".format(disk)]="PA_{:d} =".format(disk)+" ".join(str(e) for e in PA)
        try:
            phase =  float(Template["AZ1P_{:d}".format(disk)].split('=')[1])
        except KeyError:
            phase =  0.
        Template.insert("AZ1W_{:d}".format(disk),"AZ1P_{:d}".format(disk),"AZ1P_{:d} = 180.".format(disk)) 
        Template.insert("AZ1W_{:d}".format(disk),"AZ1P_{:d}".format(disk),"AZ1P_{:d} =".format(disk)+" ".join(str(phase+(e)) for e in angle_adjust))
    return PA,inc,phirings


    
def create_arms(velocity,Radii,disk_brightness,Template, disk=1,WarpStart=-1,Bar="No_Bar"):
    if WarpStart == -1: WarpStart = Rad[-1]
    max_vrot=np.max(velocity)
    max_rad=Rad[-1]
    #The pattern speed at a given radius is vrot/radii
    V_Rot = interpolate.interpolate.interp1d(Radii,velocity,fill_value="extrapolate")
    # The radius of co-ration can be approximated by the extend of the visible disk ~ Warpstart (Roberts et al. 1975)
    Omega_CR=V_Rot(WarpStart)/WarpStart
    # From this we can estimate the inner and outer Lindblad resonances (Eq 38 Dobbs & Baba 2014)
     #The epicyclic frequency k^2=R*dOmega^2/dR+4*Omega^2
    # f'(x) = (f(x+h)-f(x-h))/2h
    h=WarpStart/1000.
    derive=(V_Rot(float(WarpStart+h))**2/(WarpStart+h)**2-V_Rot(float(WarpStart-h))**2/(WarpStart-h)**2)/(2*h)
    k_CR = (WarpStart * derive+4*(Omega_CR)**2)**0.5
    # So the ILR =
    if Bar == "Barred":
        num_arms=2
    else:
        #other wise 2 arms when a bulge is present and 4 arms when not        
        if Radii[np.where(max_vrot == velocity)[0]] < np.mean(Radii):
            num_arms=2
        else:
            num_arms=4   

    
    LLR = Omega_CR-k_CR/num_arms
    ULR = Omega_CR+k_CR/num_arms
    Radii[0]=0.1
    om= interpolate.interpolate.interp1d(Radii,velocity/Radii,fill_value="extrapolate")
    Radii[0]=0.
    r_cur=Radii[1]
    while om(r_cur) > ULR and r_cur < max_rad:
        r_cur += 0.1
    ILR = 0.75*r_cur
    r_cur= Radii[1]
    while om(r_cur) > LLR and r_cur < max_rad:
        r_cur += 0.1
    OLR = 0.75*r_cur
    
    
  
    # From Seigar et al. 2006 we get the relation between shear (S) and pitch angle 
    S = 0.5*(1-Radii/velocity*derivative(Radii,V_Rot))
    pitch2= 64.25-73.24*S 
    #As we assume a constant pitch angle we will take the average between ILR and OLR as the pitch angle
    tmp = np.where((Radii > ILR) & (Radii < OLR))[0]
    pitch = np.sum(pitch2[tmp])/len(tmp)
    print("This is the average pitch angle {}".format(pitch))
    #Using Kennicut's prescription.This prescription incorrectly states cot(P) instead of tan(P) see Davis et. al 2012
    # The arms start at the inner Lindblad Resonance and hence the phase is 0 there
    phase=np.log(Radii/ILR)/np.tan(pitch*np.pi/180.)*180/np.pi
    phase[0]=phase[1]
    #How many arms do we make
    # with bar it is always a grand design
  
    # we take a brightness in the arms 1/no_arms the brightness of the disk
    brightness=1./np.sqrt(num_arms)*disk_brightness
    # and only between the resonances
    index= np.where((Radii < ILR) | (Radii > OLR))[0]
    brightness[index]=0.
    # with a ten ring transition
    brightness[tmp[0]:tmp[0]+10]=brightness[tmp[0]:tmp[0]+10]*(1-1/np.arange(1,11))
    brightness[tmp[-1]-10:tmp[-1]]=brightness[tmp[-1]-10:tmp[-1]]*(1/np.arange(1,11))
    # For the width we take 15% of the full circle but scaled to the total galaxy size a good description would be swell.
    width=0.15*2.*np.pi*Radii
    if WarpStart < 10.:
        width= width*10./WarpStart    
    ndisk=int(Template["NDISKS"].split('=',1)[1])
    ndisk+=1
    last_add = copy_disk(disk,ndisk,Template)
    Template["SBR_{:d}".format(ndisk)] = "SBR_{:d} = 0.".format(ndisk)
    # To simulate strems towrds the arms we spin this up by 10 km/s
    Template["VROT_{:d}".format(ndisk)] = "VROT_{:d} = ".format(ndisk)+" ".join(str(e+20.) for e in velocity) 
    phaseshift=360./num_arms
    #we offset the phase by 37 degrees
    for i in range(num_arms):
        Template.insert(last_add,"GA{:d}A_{:d}".format(i+1,ndisk),"GA{:d}A_{:d} =".format(i+1,ndisk)+" ".join(str(e) for e in brightness))
        Template.insert("GA{:d}A_{:d}".format(i+1,ndisk),"GA{:d}P_{:d}".format(i+1,ndisk),"GA{:d}P_{:d} =".format(i+1,ndisk)+" ".join(str(e+i*phaseshift+37) for e in phase))
        Template.insert("GA{:d}P_{:d}".format(i+1,ndisk),"GA{:d}D_{:d}".format(i+1,ndisk),"GA{:d}D_{:d} =".format(i+1,ndisk)+" ".join(str(e) for e in width))
    # and we add radial motions of 20 km/s
    Template.insert("VROT_{:d}".format(ndisk),"VRAD_{:d}".format(ndisk),"VRAD_{:d} = -40.".format(ndisk))

    return phase,brightness,width
# Routine to create a bar
def create_bar(velocity,Radii,disk_brightness,Template, disk=1,WarpStart=-1):
    if WarpStart == -1: WarpStart = Rad[-1]
    bar_width = 1+0.05*WarpStart # kpc + 5% of the total optical disk    
    max_vrot=np.max(velocity)
    #The pattern speed at a given radius is vrot/radii
    V_Rot = interpolate.interpolate.interp1d(Radii,velocity,fill_value="extrapolate")
    # The radius of co-ration can be approximated by the extend of the visible disk ~ Warpstart (Roberts et al. 1975)
    Omega_CR=V_Rot(WarpStart)/WarpStart
    # From this we can estimate the inner and outer Lindblad resonances (Eq 38 Dobbs & Baba 2012)
     #The epicyclic frequency k^2=R*d/drOmega^2+4*Omega^2
    # f'(x) = (f(x+h)-f(x-h))/2h
    h=WarpStart/1000.
    derive=(V_Rot(float(WarpStart+h))**2/(WarpStart+h)**2-V_Rot(float(WarpStart-h))**2/(WarpStart-h)**2)/(2*h)
    k_CR = (WarpStart * derive+4*(Omega_CR)**2)**0.5
    # So the ILR =
    LLR = Omega_CR-k_CR/2.
    ULR = Omega_CR+k_CR/2.
    Radii[0]=0.1
    om= interpolate.interpolate.interp1d(Radii,velocity/Radii,fill_value="extrapolate")
    Radii[0]=0.
    r_cur=Radii[1]
    while om(r_cur) > ULR and r_cur < max_rad:
        r_cur += 0.1
    ILR = 0.75*r_cur
    r_cur= Radii[1]
    while om(r_cur) > LLR and r_cur < max_rad:
        r_cur += 0.1
    OLR = 0.75*r_cur
    # We set the full brightness to the maximum of the disk
    brightness= np.zeros(len(disk_brightness))
    brightness[np.where(Radii < ILR)[0]]=np.max(disk_brightness)
    # The width has to be 180 when R < width and 180*width/(pi*r)
    width =  np.zeros(len(disk_brightness))
    width[Radii <= bar_width] =180.
    width[Radii > bar_width]=360./np.pi*np.arcsin(bar_width/Radii[Radii > bar_width]) #the angle made up of the radius and width *2.
    # Get the number of disks present
    ndisk=int(Template["NDISKS"].split('=',1)[1])
    ndisk+=1
    print("We are adding disk no {:d}".format(ndisk))
    # we also need streaming motions
    vrad_bar=np.zeros(len(disk_brightness))
    vrad_bar[:]=-50.
    #copy the input disk
    last_add = copy_disk(disk,ndisk,Template)
    #We offset by 37 deg.
    Template.insert("VSYS_{:d}".format(ndisk), "AZ1P_{:d}".format(ndisk), "AZ1P_{:d} = 37.".format(ndisk))
    Template.insert("AZ1P_{:d}".format(ndisk), "AZ1W_{:d}".format(ndisk), "AZ1W_{:d} = ".format(ndisk)+" ".join(str(e) for e in width))   
    Template.insert("AZ1W_{:d}".format(ndisk), "AZ2P_{:d}".format(ndisk), "AZ2P_{:d} = 217.".format(ndisk))
    Template.insert("AZ2P_{:d}".format(ndisk), "AZ2W_{:d}".format(ndisk), "AZ2W_{:d} = ".format(ndisk)+" ".join(str(e) for e in width))
    # And we add streaming motions to the bar km/s
    Template.insert("VROT_{:d}".format(ndisk), "VRAD_{:d}".format(ndisk), "VRAD_{:d} = ".format(ndisk)+" ".join(str(e) for e in vrad_bar))
    Template["SBR_{:d}".format(ndisk)]="SBR_{:d} =".format(ndisk)+" ".join(str(e) for e in brightness)
    return ILR
  
# afunction for creating inhomogeneities acros the disk intrinsicially
def create_inhomogeneity(Template,disks=1):
    print("wtf {} {}".format(disks,len(disks)))
    if len(disks) == 0:
        disks=[disks]
    print("wtf {} {}".format(disks,len(disks)))    
    for i in range(len(disks)):
        print(i,disks[i])
        ndisk=int(Template["NDISKS"].split('=',1)[1])
        ndisk+=1
        last_add = copy_disk(disks[i],ndisk,Template)
        sbr=np.array(Template["SBR_{:d}".format(ndisk)].split('=',1)[1].strip().split(" "),dtype=np.float64)
        rad=np.array(Template["RADI"].split('=',1)[1].strip().split(" "),dtype=np.float64)
        rad=rad*4.*np.pi/max(rad)
        newprof=np.sin(rad)*sbr*2.5
        nextprof=-1*newprof
        req_flux=abs(np.sum(newprof*np.pi*rad*2)/100.)
        Template["SBR_{:d}".format(ndisk)]="SBR_{:d} = ".format(ndisk)+" ".join(str(e) for e in newprof)
        Template["CFLUX_{:d}".format(ndisk)]="CFLUX_{:d} = {}" .format(ndisk,req_flux)
        ndisk=ndisk+1
        last_add = copy_disk(disks[i],ndisk,Template)
        Template["SBR_{:d}".format(ndisk)]="SBR_{:d} = ".format(ndisk)+" ".join(str(e) for e in nextprof)
        Template["CFLUX_{:d}".format(ndisk)]="CFLUX_{:d} = {}" .format(ndisk,req_flux)                                                                             
    print("We will create inhomogeneities on top of disk(s) ="+" ".join(str(e) for e in [disks]))


# This is a routine to corrupt the cube with some standard gaussian noise    
def corrupt_gauss(work_dir,beam,noise):
    # First open the model cube
    dummy = fits.open(work_dir+'/unconvolved_cube.fits',uint = False, do_not_scale_image_data=True,ignore_blank = True)
    dummy[0].header.append('RESTFRQ')
    dummy[0].header['RESTFRQ'] = 1.420405752E+09
    # Make a moment zero
    mom=np.sum(dummy[0].data,axis=0)
    fits.writeto(work_dir+'/mom0_un.fits',mom,dummy[0].header, overwrite = True)
  # In order to create a cleaning mask we smooth to twice the beam size and cut at 1e-5 Jy/beam
    sig=abs(beam[0]/(dummy[0].header['CDELT1']*3600.))
    sigmin=abs(beam[1]/(dummy[0].header['CDELT2']*3600.))
    reals = dummy[0].data[dummy[0].data > 1e-6]
    # Remember that python is a piece of shiit so z,y,x
    cutoff = (np.mean(reals)/2.)
    print("This is the cutoff {}".format(cutoff))
    smooth = ndimage.gaussian_filter(dummy[0].data, sigma=(0,sigmin, sig), order=0)
    smooth[smooth < cutoff]=0
    smooth[smooth > cutoff]=1
    fits.writeto(work_dir+'/mask.fits',smooth,dummy[0].header, overwrite = True)
    
    # Calculate the sigma's from the required beam size
    sigma=[(beam[0]/abs(dummy[0].header['CDELT1']*3600.))/(2*np.sqrt(2*np.log(2))),(beam[1]/(dummy[0].header['CDELT2']*3600.))/(2*np.sqrt(2*np.log(2)))]
    # Calculate the area of the beam in arcsec
    beamarea=(np.pi*abs(beam[0]*beam[1]))/(4.*np.log(2.))
    # Convert arcsec to pixels
    pixperbeam=beamarea/(abs(dummy[0].header['CDELT1']*3600.)*abs(dummy[0].header['CDELT2']*3600.))
    # Make an the size of the model with random values distributed as the noise
    # first we need to scale to the proper units where does the division by 4 come from? 
    noisescl=noise*sigma[0]*2*np.sqrt(np.pi)/4.
    cuberms = np.random.normal(scale=noisescl,size=np.shape(dummy[0].data))
    # combine the two cubes    
    noisedcube=dummy[0].data+cuberms
    # Smooth to the requred resolution
    print("The beam in pixels {} x {}".format(sigma[0]*(2.*np.sqrt(2*np.log(2))),sigma[1]*(2.*np.sqrt(2*np.log(2)))))
    final = ndimage.gaussian_filter(noisedcube, sigma=(0,sigma[1], sigma[0]), order=0)
    # to preserve brightness temperature this shoul be multiplied with the increas in beam area
    # which is the same as the amount of pixels in the beam as we go from 1 pixel area to an area the size of the beam
    final=final*pixperbeam
    # And write this final cube to the directory
    dummy[0].header['BMAJ']=beam[0]/3600.
    dummy[0].header['BMIN']=beam[1]/3600.
    fits.writeto(work_dir+'/Convolved_Cube.fits',final,dummy[0].header, overwrite = True)
     
# this is a routine to use the casa task simobserve to corrupt the observations.    
def corrupt_casa(work_dir,beam,Template_Casa,noise):
    #casa takes issue with the tirific header.
    dummy = fits.open(work_dir+'/unconvolved_cube.fits',uint = False, do_not_scale_image_data=True,ignore_blank = True)
    dummy[0].header.append('RESTFRQ')
    dummy[0].header['RESTFRQ'] = 1.420405752E+09
    fits.writeto(work_dir+'/unconvolved_cube.fits',dummy[0].data,dummy[0].header, overwrite = True)
    # Calculate the required noise in the end cube
    sigma=[(beam[0]/abs(dummy[0].header['CDELT1']*3600.))/(2*np.sqrt(2*np.log(2))),(beam[1]/(dummy[0].header['CDELT2']*3600.))/(2*np.sqrt(2*np.log(2)))]
    noisescl=noise*sigma[0]*2*np.sqrt(np.pi)/4.

    
    # We need to know the location in the sky dec > 30 we use WSRT 30 > dec > -30 --> VLA -30 > dec  --> Atca
    
    # Full synthesis 12 hrs leads to noise levels 12: 0.000701111 ,24:  0.00065377   ,48:0.000499886 ,96:0.000436001 ,300:0.000213708 ,600:0.000178363   ,900:0.000149931 1200:0.000119707 2400:9.85791e-05 , 4800: 8.18873e-05 # Based on a simulation with no taper
    noisin = [0.000701111,0.00065377 ,0.000499886,0.000436001,0.000213708,0.000178363,0.000149931,0.000119707 ,9.85791e-05 , 8.18873e-05]
    timein = [12.,24,48,96,300,600,900,1200,2400,4800]
    a1,a2 =np.polyfit(noisin,timein,1)
    totaltime = a1*noise+a2
    if totaltime < 12:
        totaltime = 12.
    RA,DEC =convertRADEC(dummy[0].header['CRVAL1'],dummy[0].header['CRVAL2'])
    tri = open(work_dir+'/pntings.txt', 'w') 
    tri.writelines("#Epoch     RA          DEC      TIME(optional) \n ")
    tri.writelines("J2000     {}          {}      {} ".format(RA,DEC,str(int(12*3600.))))
    
    tri.close()
    
    if  dummy[0].header['CRVAL2'] > 90:
        Template_Casa['simobserve_antennalist']="antennalist = 'WSRT.cfg'  #  interferometer antenna position file"
        Template_Casa['imhead_hdvalue'] = "hdvalue = 'WSRT'"
        Template_Casa['tclean_vis'] = "vis = 'simulated/simulated.WSRT.noisy.ms'"
    elif dummy[0].header['CRVAL2'] > -30:
        Template_Casa['simobserve_antennalist']="antennalist = 'vla.b.cfg'  #  interferometer antenna position file"
        Template_Casa['imhead_hdvalue'] = "hdvalue = 'VLA'"
        Template_Casa['tclean_vis'] = "vis = 'simulated/simulated.vla.b.noisy.ms'"
    else:
        Template_Casa['simobserve_antennalist']="antennalist = 'atca_6c.cfg'  #  interferometer antenna position file"
        Template_Casa['imhead_hdvalue'] = "hdvalue = 'ATCA'"
        Template_Casa['tclean_vis'] = "vis = 'simulated/simulated.atca_6c.noisy.ms'"
    #let's assure the same cube size
    # We want 2000 integrations
    Template_Casa['simobserve_integration'] = "integration    = '{:d}s'".format(int(totaltime*3600/2000.))
    Template_Casa['simobserve_totaltime'] = "totaltime    = '{:d}'".format(int(totaltime/12.))  
    Template_Casa['tclean_cell']= "cell = ['{}arcsec','{}arcsec']".format(abs(dummy[0].header['CDELT1']*3600.),abs(dummy[0].header['CDELT2']*3600.))
    Template_Casa['tclean_imsize'] = "imsize=[{:d},{:d}]".format(int(abs(dummy[0].header['NAXIS1'])),int(abs(dummy[0].header['NAXIS2'])))
    Template_Casa['tclean_scales'] = "scales=[0,{:d},{:d}]".format(int(2.*beam[0]/abs(dummy[0].header['CDELT1']*3600.)),int(5.*beam[1]/abs(dummy[0].header['CDELT1']*3600.)))
    Template_Casa['tclean_threshold'] = "threshold = '{}Jy/beam'".format(noise/2.)
    Template_Casa['tclean_uvtaper'] = "uvtaper = ['{}arcsec','{}arcsec']".format(beam[0],beam[1])
     
    # In order to create a cleaning mask we smooth to twice the beam size and cut at 1e-5 Jy/beam
    sig=abs(beam[0]/(dummy[0].header['CDELT1']*3600.))
    sigmin=abs(beam[1]/(dummy[0].header['CDELT2']*3600.))
    reals = dummy[0].data[dummy[0].data > 1e-6]
    # Remember that python is a piece of shiit so z,y,x
    cutoff = (np.mean(reals)/2.)
    print("This is the cutoff {}".format(cutoff))
    smooth = ndimage.gaussian_filter(dummy[0].data, sigma=(0,sigmin, sig), order=0)
    smooth[smooth < cutoff]=0
    smooth[smooth > cutoff]=1
    fits.writeto(work_dir+'/mask.fits',smooth,dummy[0].header, overwrite = True)
    tri = open(work_dir+'/run_casa.py', 'w') 
    tri.writelines([Template_Casa[key]+"\n" for key in Template_Casa])
    tri.close()
    
    os.chdir(work_dir)
    bla = subprocess.call(['casa','--nologger','--nogui','-c','run_casa.py'])
    dummy = fits.open(work_dir+'/Convolved_Cube.fits',uint = False, do_not_scale_image_data=True,ignore_blank = True)
    # We cut 20 pixels around the edges
    newsize = int(np.shape(dummy[0].data)[2]-60.)
    print(newsize)
    newdummy=np.zeros((np.shape(dummy[0].data)[0],newsize,newsize))
    newdummy[:,:,:]=dummy[0].data[:,int(np.floor(dummy[0].header['CRPIX2']-newsize/2.)):int(np.floor(dummy[0].header['CRPIX2']+newsize/2.)),int(np.floor(dummy[0].header['CRPIX1']-newsize/2.)):int(np.floor(dummy[0].header['CRPIX1']+newsize/2.))]
    dummy[0].header['NAXIS1']=newsize+1
    dummy[0].header['NAXIS2']=newsize+1
              
    dummy[0].header['CRPIX1']=dummy[0].header['CRPIX1']-np.floor(dummy[0].header['CRPIX1']-newsize/2.)
    dummy[0].header['CRPIX2']=dummy[0].header['CRPIX2']-np.floor(dummy[0].header['CRPIX2']-newsize/2.)
    fits.writeto(work_dir+'/Convolved_Cube.fits',newdummy,dummy[0].header, overwrite = True)
    # Also the mask then
    dummy = fits.open(work_dir+'/mask.fits',uint = False, do_not_scale_image_data=True,ignore_blank = True)
    # We cut 20 pixels around the edges
    newsize = int(np.shape(dummy[0].data)[2]-60.)
    print(newsize)
    newdummy=np.zeros((np.shape(dummy[0].data)[0],newsize,newsize))
    newdummy[:,:,:]=dummy[0].data[:,int(np.floor(dummy[0].header['CRPIX2']-newsize/2.)):int(np.floor(dummy[0].header['CRPIX2']+newsize/2.)),int(np.floor(dummy[0].header['CRPIX1']-newsize/2.)):int(np.floor(dummy[0].header['CRPIX1']+newsize/2.))]
    dummy[0].header['NAXIS1']=newsize+1
    dummy[0].header['NAXIS2']=newsize+1
              
    dummy[0].header['CRPIX1']=dummy[0].header['CRPIX1']-np.floor(dummy[0].header['CRPIX1']-newsize/2.)
    dummy[0].header['CRPIX2']=dummy[0].header['CRPIX2']-np.floor(dummy[0].header['CRPIX2']-newsize/2.)
    fits.writeto(work_dir+'/mask.fits',newdummy,dummy[0].header, overwrite = True)

sets = 3  # This is the amount of base galaxies we want, i.e. the number of rotation curves 
# do we wanr inhomogeneities
inhomogeneity = True
# How do we want to corrupt the model, For now the options are Gaussian or Casa_Sim
corruption = 'Casa_Sim'
corruption = 'Gaussian'

changes = ['inclination']
#'warp','flare','beams','SNR','Channelwidth','Res_Beam'] # parameters that should vary
# options are inclination, PA, flare, warp, beams, SNR, Channelwidth, Res_Beam, Arms, Bar, Radial_Motions
# The requested variations in parameters
# The indices for these will run with i_variable
inclination= [10,15,20,25,30,40,80,85,88,90]
inclination=[40.]
warp=[[0.15,0.05],[0.05,0.2]] 
flare="Flared" # Flared, No_Flare
beams=[4,6,7,8,9,10,11,12,16,32] # Beam across the major axis. This also set the distance as the size in kpc will be determined by Wang 2016 from the SBR profile
SNR=[2,4,8,16] # These  are average signal to noise ratios
Channelwidth=[2.,8.] 
Res_Beam=[[15,15]] # Resolution of the beam in arcsec
Arms = "No_Arms"   #Either "No_Arms or Arms"
Bar = "No_Bar"   #Either "No_Bar or Bar"
Radial_Motions = [0] # in km/s
# The directory to work in

#work_dir='/Users/Peter/Database/Artificial/'
work_dir='/home/peter/Database/Artificial/'
#Parameter to force a new set of models being made in this directory
makenewmodels =True
# do we want symmetric models or not?
symmetric = False
# If we make new models delete everything in the directory
if makenewmodels: 
    os.system('rm -R '+work_dir+'/*')


# Let's just make 1 catalogue with everything we need and adjust
# the fitting program to have this one catalogue as input
    
Catalogue=work_dir+'/Output_Summary.txt'
# If we are making new models we want to ensure this is a new file
if makenewmodels:
   cat = open(Catalogue, 'w')
   cat.write('number|Distance|Directoryname|Cubename\n')
   cat.close()



#Copy a fits file from the WHISP data base to use as template if it is not ther3 yet
# Check for the existence of a template fits file
templatethere= os.path.isfile(work_dir+'/Input.fits')
templatethere =False
# If it doesn't exist copy it into the directory
if not templatethere:
    os.system('cp Input.fits '+work_dir+'/Input.fits')
   # let's make sure it has BMAJ, BMIN and BPA
    dummy = fits.open(work_dir+'/Input.fits',uint = False, do_not_scale_image_data=True,ignore_blank = True)
    sizex = 500.
    sizey = 500.
                                #first we do a check wether the BMAJ
                                #is written correctly for the python
                                #program. We set a bunch of generic header values.
    dummy[0].header['BMAJ'] = 0.
    dummy[0].header['BMIN'] = 0.
    dummy[0].header['BPA'] = 0.
    dummy[0].header['CRPIX1'] = 200.
    dummy[0].header['CRPIX2'] = 200.
    dummy[0].header['CRPIX3'] = 60.
    dummy[0].header['CDELT1'] = -4./3600.
    dummy[0].header['CDELT2'] = 4./3600.
    dummy[0].header['CDELT3'] = 4.
    dummy[0].header['CUNIT3'] = 'M/S'
    dummy[0].header['CRVAL3'] = 500.
    dummy[0].header['CTYPE1'] = 'RA---SIN'
    dummy[0].header['CTYPE2'] = 'DEC--SIN'
    dummy[0].header['CTYPE3'] = 'VELO-HEL'
    dummy[0].header['BITPIX'] = -32
    dummy[0].header['NAXIS'] = 3
    dummy[0].header['NAXIS1'] = 500
    dummy[0].header['NAXIS2'] = 500
    dummy[0].header['NAXIS3'] = 120.
    del  dummy[0].header['CROTA1']
    del  dummy[0].header['CROTA2']
    del  dummy[0].header['DRVAL3']
    del  dummy[0].header['DTYPE3']
    del  dummy[0].header['DUNIT3']
    del  dummy[0].header['HISTORY']
    del  dummy[0].header['BMMAJ']
    del  dummy[0].header['BMMIN']
    del  dummy[0].header['BMPA']
    
    # make the cube a typical size
    dummy2 = np.zeros((120,500,500))
    dummy2[60,250,250] = 5
    fits.writeto(work_dir+'/Input.fits',dummy2,dummy[0].header, output_verify='silentfix+ignore', overwrite = True)
    dummy.close()


# All this went well
templatethere= os.path.isfile(work_dir+'/Input.fits')
if templatethere:
    print(" The Input Template is found and all is well")
else:
    print(" The Input Template is NOT found !!! ABORTING")
    sys.exit()


#Read a template def file, We typically read them from the FAT installation

#tmp = open('/Users/Peter/GitHub/FAT-GDL-Beta/Support/2ndfit.def','r')
tmp = open('Template.def','r')
Template_in = MyOrderedDict({})
unarranged = tmp.readlines()
# Separate the keyword names
for tmp in unarranged:
    # python is really annoying with needing endlines. Let's strip them here and add them when writing
    Template_in[tmp.split('=',1)[0].strip().upper()]=tmp.rstrip()

#If we want to corrupt in the casa way we'd need to read the file corruption file
if corruption =='Casa_Sim':
    tmp = open('Template_Casa.py','r')
    Template_Casa_In = MyOrderedDict({})
    unarranged = tmp.readlines()
    # Separate the keyword names
    current_task="pre_task"
    counter=0
    tasks=["pre_task"]
    task_counter = 0
    deccounter = 0
    decin= np.arange(43.951946,84,5.)    
    for tmp in unarranged:
    # python is really annoying with needing endlines. Let's strip them here and add them when writing
        if re.split('\(|\)',tmp)[0] == 'default':
            task_counter = 0
            for_task = 0
            current_task=re.split('\(|\)',tmp)[1]
            while current_task in tasks:
                for_task += 1
                current_task=re.split('\(|\)',tmp)[1]+"_{:d}".format(for_task)
            tasks.append(current_task)    
        if len(tmp.split('=',1)) >1:
            variable=tmp.split('=',1)[0].strip().lower()
            same_var = 0
            while current_task+'_'+variable in Template_Casa_In:
                same_var += 1
                variable=tmp.split('=',1)[0].strip().lower()+"_{:d}".format(same_var)
        elif len(tmp.split('#',1)) >1:
            counter+=1
            variable = "comment_{:d}".format(counter)
        else:
            task_counter +=1
            variable = "run_{:d}".format(task_counter)       
        Template_Casa_In[current_task+'_'+variable]=tmp.rstrip()
 
H_0 = 69.6 # http://www.astro.ucla.edu/~wright/CosmoCalc.html
# start a loop over the various base galaxies
number_models = 0.
set_done= [1024]
colors=['k','b','r','c','m','y','w']
for base in range(sets):

    
    # From here we go into a loop to adjust variables over the bases
    for ix in range(len(changes)):
        if changes[ix] == 'inclination':numloops=len(inclination)
        elif changes[ix] == 'PA': numloops=len(PA)
        elif changes[ix] == 'flare': numloops=len(flare[0,:])
        elif changes[ix] == 'warp': numloops=len(warp[0,:])
        elif changes[ix] == 'beams': numloops=len(beams)
        elif changes[ix] == 'SNR': numloops=len(SNR)
        elif changes[ix] == 'Channelwidth': numloops=len(Channelwidth)
        elif changes[ix] == 'Res_Beam': numloops=len(Res_Beam[0,:])
        elif changes[ix] == 'Arms': numloops=len(Arms)
        elif changes[ix] == 'Bar': numloops=len(Bar)
        elif changes[ix] == 'Radial_Motions': numloops=len(Radial_Motions)
        else:print("This is not a supported parameter")
        for jx in range (numloops):
            Current_Galaxy = Base_Galaxy(base)
            number_models += 1
            if changes[ix] == 'inclination':Current_Galaxy.inclination = inclination[jx]
            elif changes[ix] == 'PA': Current_Galaxy.PA = PA[jx]
            elif changes[ix] == 'flare': Current_Galaxy.flare = flare[jx]
            elif changes[ix] == 'warp': Current_Galaxy.warp = [warp[0,jx],warp[1,jx]]
            elif changes[ix] == 'beams':Current_Galaxy.beams = beams[jx]
            elif changes[ix] == 'SNR': Current_Galaxy.SNR = SNR[jx]
            elif changes[ix] == 'Channelwidth': Current_Galaxy.Channelwidth = Channelwidth[jx]
            elif changes[ix] == 'Res_Beam': Current_Galaxy.Res_Beam = [Res_Beam[0,jx],Res_Beam[1,jx]]
            elif changes[ix] == 'Arms': Current_Galaxy.Arms = Arms[jx]
            elif changes[ix] == 'Bar': Current_Galaxy.Bar = Bar[jx]
            elif changes[ix] == 'Radial_Motions': Current_Galaxy.Radial_Motions = Radial_Motions[jx]
            else:print("This is not a supported parameter")
            Current_Galaxy.Res_Beam[0:1] = np.sort(Current_Galaxy.Res_Beam[0:1])
            if len(Current_Galaxy.Res_Beam) == 2: Current_Galaxy.Res_Beam.append(0) 
            print("This is the parameter {}. And this is the input {}".format(changes[ix],Current_Galaxy.inclination))         
            # Build a name and a directory where to stor the specific output
            name="{}-i{}d{}-{}pa{}w{}-{}-{}-ba{}SNR{}bm{}-{}ch{}-{}-{}-rm{}".format(Current_Galaxy.rc_name,Current_Galaxy.inclination,Current_Galaxy.dispersion[0],Current_Galaxy.dispersion[1],Current_Galaxy.PA,Current_Galaxy.warp[0],Current_Galaxy.warp[1],Current_Galaxy.flare,Current_Galaxy.beams,Current_Galaxy.SNR,Current_Galaxy.Res_Beam[0],Current_Galaxy.Res_Beam[1],Current_Galaxy.Channelwidth,Current_Galaxy.Arms,Current_Galaxy.Bar,Current_Galaxy.Radial_Motions)
            print("{} is the name of the current galaxy".format(name))
            #We want to make a figure with all the Rotation curves
            if base not in set_done:
                if set_done[0] == 1024:
                    set_done= [base]
                    labelfont= {'family':'Times New Roman',
                             'weight':'normal',
                             'size':22}
                    plt.rc('font',**labelfont)    
                    plt.figure(59, figsize=(8, 8), dpi=300, facecolor='w', edgecolor='k')
                    ax = plt.subplot(1, 1, 1)
                    plt.plot(Current_Galaxy.rc_radii,Current_Galaxy.rc_speed,'k') 
                    plt.plot(Current_Galaxy.rc_radii,Current_Galaxy.rc_speed,'ko',label=Current_Galaxy.rc_name)   
                    plt.ylabel('V$_{rot}$ (km s$^{-1}$)',**labelfont)
                    plt.xlabel('Radius (kpc)',**labelfont)
                    ax.yaxis.set_minor_locator(AutoMinorLocator(4))
                    ax.xaxis.set_minor_locator(AutoMinorLocator(4))                    
                    for axis in ['top','bottom','left','right']:
                        ax.spines[axis].set_linewidth(1.5)
                    plt.tick_params(axis='both', which='minor', bottom='on',left='on',length=3)
                    plt.tick_params(axis='both', which='major', labelsize=17, length=6)
                    plt.tick_params(axis='both', which='both', direction = 'in', width=1.5 , bottom='on',left='on' ,right ='on', top='on')
                    max_rad = np.max(Current_Galaxy.rc_radii)+1
                    names_used=[Current_Galaxy.rc_name]
                    colors_used=[colors[base]]
                else: 
                    set_done.append(base)
                    plt.figure(59)
                    if np.max(Current_Galaxy.rc_radii)+1 > max_rad:
                        max_rad =  np.max(Current_Galaxy.rc_radii)+1
                    plt.plot(Current_Galaxy.rc_radii,Current_Galaxy.rc_speed,colors[base]) 
                    plt.plot(Current_Galaxy.rc_radii,Current_Galaxy.rc_speed,colors[base]+'o',label=Current_Galaxy.rc_name)
                    names_used.append(Current_Galaxy.rc_name)
                    colors_used.append(colors[base])
                    
                    
            # Make a dirctory
            # Check for the existence of the directory
            constructstring="mkdir "+work_dir+name
            
            galaxy_dir= os.path.isdir(work_dir+name)
            if not galaxy_dir:
                os.system(constructstring)
            else:
                # Do we have a cube
                galaxy_cube_exist = os.path.isfile(work_dir+name+"/Cube.fits")
                if galaxy_cube_exist:
                    print("This galaxy appears fully produced")
                    continue
                else:
                    print("The directory was made but there is no full cube avalaible")
                    print("Reproducing the galaxy. Be aware of Double Table entries")
            # Then we copy the original Template def
            Template=copy.deepcopy(Template_in)
            # We also open a figure to plot all info on 
            plt.figure(2, figsize=(8, 12), dpi=100, facecolor='w', edgecolor='k')
            overview = plt.subplot(6,1,6)
            
             # First set the beam to 0.            
            Template["BMAJ"]= "BMAJ = 0."
            Template["BMIN"]= "BMIN = 0."
            Template["BPA"]= "BPA = 0."
            #Then we need to build the Surface Brightnes profile    
            SBRprof,Rad,sclength,MHI,Rad_HI = build_sbr_prof(base,Current_Galaxy.beams,Template) #Column densities,Raii in kpc, Opt_scalelength in kpc, HI mass in M_solar 
            #We need central coordinates the vsys will come from the required distance and hubble flow. The RA and dec should not matter hance it will be the only random component in the code as we do want to test variations of them
            Sky_Size = np.radians(Current_Galaxy.Res_Beam[0]*Current_Galaxy.beams/3600.)
            Distance = (Rad_HI/(np.tan(Sky_Size/2.)))/1000.
            print("The Distance is {:5.2f} Mpc".format(Distance))
            vsys = Distance*H_0
            if corruption == 'Gaussian':           
                RAdeg=np.random.uniform()*360
                DECdeg=(np.arccos(2*np.random.uniform()-1)*(360./(2.*np.pi)))-90
            else:
                RAdeg = np.random.uniform()*360
                DECdeg=-60
                while DECdeg < -20.:
                    DECdeg = (np.arccos(2*np.random.uniform()-1)*(360./(2.*np.pi)))-90
            RAhr,DEChr= convertRADEC(RAdeg,DECdeg)
            print("It's central coordinates are RA={} DEC={} vsys={} km/s".format(RAhr,DEChr,vsys))
            # Write them to the template
            Template["XPOS"]="XPOS= {}".format(RAdeg)
            Template["XPOS_2"]="XPOS_2= {}".format(RAdeg)
            Template["YPOS"]="YPOS= {}".format(DECdeg)
            Template["YPOS_2"]="YPOS_2= {}".format(DECdeg)
            Template["VSYS"]="VSYS= {}".format(vsys)
            Template["VSYS_2"]="VSYS_2= {}".format(vsys)
            # The Distance, this is not required but can be usefull
            Template["DISTANCE"]="DISTANCE= {}".format(Distance)            
            # With the distance we can also convert our radii
             # we need our radii in the corresponding arcsec
            Rad_arcsec = convertskyangle(Rad,distance=Distance,physical = True)
             # then the number of rings to the total number of rings
            Template["NUR"]="NUR= {}".format(len(Rad))
            # The radii in arcsec
            Template["RADI"]="RADI = "+" ".join(str(e) for e in Rad_arcsec)
            # Then we need to get a starting radius for the warp.
            # the warp should start at the edge of the optical radius which is the HI scale length/0.6
            # which are about ~ 4 * h_r
            #WarpStart=4*sclength
            WarpStart = 4.*sclength
            # we need the rotation speed to correspond to our radii
            V_Rot = interpolate.interpolate.interp1d(Current_Galaxy.rc_radii,Current_Galaxy.rc_speed,fill_value="extrapolate")
            # Write it to the Template
            Template["VROT"]="VROT = "+" ".join(str(e) for e in V_Rot(Rad))
            Template["VROT_2"]="VROT_2 = "+" ".join(str(e) for e in V_Rot(Rad))
            #and the Overview Figure
            plt.figure(2)
            plt.subplot(6,1,5)
            plt.plot(Rad,V_Rot(Rad))
            plt.plot([Rad_HI,Rad_HI],[0.,np.max(V_Rot(Rad))])   
            plt.ylabel('V_Rot (km/s)')
            plt.xticks([])
            # We need a scale height and dispersion for each ring. They are coupled and hence they are both created in create_flare
            h_z,dispersion=create_flare(Rad,V_Rot(Rad),Current_Galaxy.dispersion,Current_Galaxy.flare,Template,Rad_HI,distance=Distance)

            # Finally we need to set the warping
            PA,inc,phirings = create_warp(Rad,Current_Galaxy.PA,Current_Galaxy.inclination,Current_Galaxy.warp,[WarpStart,Rad_HI],Template)
            # we want an assymetric warp so we redo the PA and inc but swap the variation
            PA_2,inc_2,phirings_2 = create_warp(Rad,Current_Galaxy.PA,Current_Galaxy.inclination,[Current_Galaxy.warp[0]-Current_Galaxy.warp[1],Current_Galaxy.warp[1]/2.+Current_Galaxy.warp[0]/2.],[WarpStart,Rad_HI],Template,disk=2)
            
            # This comes from FAT. If I remember correctly this is the sine response to the channels *1.2/(2*SQRT(2*ALOG(2.))))
            # However in our input we want independent channels
            Template["CONDISP"]="CONDISP = {}".format(Current_Galaxy.Channelwidth)
            # We need to set the input and output cube
            Template["INSET"]="INSET = ../Input.fits"
            Template["OUTSET"]="OUTSET = unconvolved_cube.fits"
            #Some tirific varaiables
            Template["LOOPS"]="LOOPS = 0 "
            # We need models with about 3 million particles but not more as it takes too long
            TotFlux=MHI/(2.36e5*(Distance)**2)
            Template["CFLUX"]="CFLUX = {}".format(TotFlux/5e6)
            Template["CFLUX_2"]="CFLUX_2 = {}".format(TotFlux/5e6)
            Template.insert("RMS","NDISKS","NDISKS = 2")
            Template["TIRDEF"]="TIRDEF = ModelInput.def"
            Template["GR_DEVICE"]="GR_DEVICE = "
            Template["GR_CONT"]="GR_CONT = "
            Template.insert("INSET","ACTION","ACTION = 1")
            Template["PROGRESSLOG"]="PROGRESSLOG = "
            Template["LOGNAME"]= "LOGNAME= "
            #-------------------------------This finishes the basic disk the following are optional components----------------------

            
            # The possible arms
            if Current_Galaxy.Arms == 'Arms':
                phase,arm_brightness,arm_width = create_arms(V_Rot(Rad),Rad,SBRprof,Template,WarpStart=WarpStart, Bar=Current_Galaxy.Bar)
                phase,arm_brightness,arm_width = create_arms(V_Rot(Rad),Rad,SBRprof,Template,disk=2,WarpStart=WarpStart, Bar=Current_Galaxy.Bar)
            #Template["SBR_2"]="SBR_2=0."
            # A possible Bar     
            if Current_Galaxy.Bar == 'Barred':
                bar_length = create_bar(V_Rot(Rad),Rad,SBRprof,Template,WarpStart=WarpStart)
            # and possible inhomogeneities
            if inhomogeneity:
                inhomogeneity_amp = create_inhomogeneity(Template,disks=[1,2])
            # we write the def files        
            tri = open(work_dir+'/'+name+'/tirific.def', 'w') 
            tri.writelines([Template[key]+"\n" for key in Template])
            tri.close()
            # So we need to  modify the input file to the correct coordinates else we'll get an empty cube
            dummy = fits.open(work_dir+'/Input.fits',uint = False, do_not_scale_image_data=True,ignore_blank = True)
            #first we do a check wether the BMAJ
            #is written correctly for the python
            #program. We set a bunch of generic header values.
            #First we make sure that it fits
           
            size = 2.*Rad_arcsec[-1] +(Current_Galaxy.Res_Beam[0])
            pix_size = (Rad_arcsec[1]-Rad_arcsec[0])
            if corruption == 'Casa_Sim':
                size += 60*pix_size
            required_pixels=int(np.ceil(size/pix_size))
            vel_max=2.5*np.max([V_Rot(Rad)*np.sin(inc*np.pi/180.)+2*dispersion,V_Rot(Rad)*np.sin(inc_2*np.pi/180.)+2*dispersion])
            velpix=int(np.ceil(vel_max/Current_Galaxy.Channelwidth) )                   
            dummy[0].header['CRPIX1'] = np.floor(required_pixels/2.)
            dummy[0].header['CRPIX2'] = np.floor(required_pixels/2.)
            dummy[0].header['CRPIX3'] = np.floor(velpix/2.)
            # Stupid astropy doesn't account for minus in header of cdelt1 then cdelt has different precision
            tmp=int(pix_size/3600.*1e15)/1e15
            dummy[0].header['CDELT1'] = -1*tmp
            dummy[0].header['CDELT2'] = tmp
            dummy[0].header['CDELT3'] = Current_Galaxy.Channelwidth*1000.
            dummy[0].header['CRVAL1'] = RAdeg
            dummy[0].header['CRVAL2'] = DECdeg
            dummy[0].header['CRVAL3'] = vsys*1000.
            dummy[0].header['NAXIS1'] = required_pixels
            dummy[0].header['NAXIS2'] = required_pixels
            dummy[0].header['NAXIS3'] = velpix
            dummy[0].header['BMAJ'] = 0.
            dummy[0].header['BMIN'] = 0.
         
            # make the cube a typical size
            dummy2 = np.zeros((velpix,required_pixels,required_pixels),dtype= np.float32)
            dummy2[int(np.floor(velpix/2.)),int(np.floor(required_pixels/2.)),int(np.floor(required_pixels/2.))] = 5
            fits.writeto(work_dir+'/Input.fits',dummy2,dummy[0].header, output_verify='silentfix+ignore', overwrite = True)
            dummy.close()
            os.chdir(work_dir+'/'+name)
            os.system("tirific deffile=tirific.def")
            
            plt.savefig('Overview.png', bbox_inches='tight')
            plt.savefig('Overview.ps', bbox_inches='tight')
            plt.close()
            os.chdir(work_dir)
            Template.clear()
          
            # Now we want to corrupt this cube with some realistic noise
            # For this we first want to get the noise we want in terms of Jansky per beam
            # we will define the SNR as the mean(Intensity)/noiselevel hence noise =mean(In)/SNR
            unconvolved = fits.open(work_dir+name+'/unconvolved_cube.fits',uint = False, do_not_scale_image_data=True,ignore_blank = True)
            if inhomogeneity:
                unconvolved[0].data[unconvolved[0].data < 0] = 0.
                fits.writeto(work_dir+name+'/unconvolved_cube.fits',unconvolved[0].data,unconvolved[0].header, output_verify='silentfix+ignore', overwrite = True)
            # In order to corrupt we need to know the average signal.
            # we do this taking the mean in each chaneel above a tenth of the max and then take the mean of that profile
            
            totsig=np.zeros(len(unconvolved[0].data[:]))
            print("This is the z-axis length {}".format(len(unconvolved[0].data[:])))
            for j in range(len(totsig)):
                if  len(unconvolved[0].data[j][unconvolved[0].data[j] > np.max(unconvolved[0].data)/10.]) > 0:
                    totsig[j]=np.mean(unconvolved[0].data[j][unconvolved[0].data[j] > np.max(unconvolved[0].data)/10.])
            mean_signal = np.mean(totsig[totsig > 0.])        
            mean_signal_avg = np.mean(unconvolved[0].data[unconvolved[0].data > np.max(unconvolved[0].data)/10.])
            noise=mean_signal/Current_Galaxy.SNR
            print("This our mean signal {} and noise {} in Jy/pixel".format(mean_signal,noise))
            #Then we need to corrupt and concolve our model
            if corruption == 'Casa_Sim':
                Template_Casa=copy.deepcopy(Template_Casa_In)
                corrupt_casa(work_dir+name+'/',Current_Galaxy.Res_Beam,Template_Casa,noise)
                os.chdir(work_dir)
                mask = fits.open(work_dir+name+'/mask.fits',uint = False, do_not_scale_image_data=True,ignore_blank = True)
                Cube = fits.open(work_dir+name+'/Convolved_Cube.fits',uint = False, do_not_scale_image_data=True,ignore_blank = True)
                maskr = mask[0].data[1:]
                sigma = (np.std(Cube[0].data[0])+np.std(Cube[0].data[-1]))/2.
                Cube_Clean = Cube[0].data
                Cube_Clean[maskr < 0.5] = 0.
                beamarea=(np.pi*abs(Cube[0].header['BMAJ']*3600.*Cube[0].header['BMIN']*3600.))/(4.*np.log10(2.))
                pixperbeam=beamarea/(abs(Cube[0].header['CDELT1']*3600.)*abs(Cube[0].header['CDELT2']*3600.))
                totalsignal = np.sum(Cube_Clean)/pixperbeam
                mass = 2.36E5*Distance**2*totalsignal*Cube[0].header['CDELT3']/1000.
                totsig=np.zeros(len(Cube_Clean[:]))
                for j in range(len(totsig)):
                    if  len(Cube_Clean[j][Cube_Clean[j] > 0.]) > 0:
                        totsig[j]=np.mean(Cube_Clean[j][Cube_Clean[j] > 0.])
                mean_signal = np.mean(totsig[totsig > 0.]) 
                SNR = mean_signal/(sigma)
                
            elif corruption == 'Gaussian':
                corrupt_gauss(work_dir+name+'/',Current_Galaxy.Res_Beam,noise)
                mask = fits.open(work_dir+name+'/mask.fits',uint = False, do_not_scale_image_data=True,ignore_blank = True)
                Cube = fits.open(work_dir+name+'/Convolved_Cube.fits',uint = False, do_not_scale_image_data=True,ignore_blank = True)
                maskr = mask[0].data[:]
                sigma = (np.std(Cube[0].data[0])+np.std(Cube[0].data[-1]))/2.
                Cube_Clean = Cube[0].data
                Cube_Clean[maskr < 0.5] = 0.
                beamarea=(np.pi*abs(Cube[0].header['BMAJ']*3600.*Cube[0].header['BMIN']*3600.))/(4.*np.log10(2.))
                pixperbeam=beamarea/(abs(Cube[0].header['CDELT1']*3600.)*abs(Cube[0].header['CDELT2']*3600.))
                totalsignal = np.sum(Cube_Clean)/pixperbeam
                mass = 2.36E5*Distance**2*totalsignal*Cube[0].header['CDELT3']/1000.
                totsig=np.zeros(len(Cube_Clean[:]))
                for j in range(len(totsig)):
                    if  len(Cube_Clean[j][Cube_Clean[j] > 0.]) > 0:
                        totsig[j]=np.mean(Cube_Clean[j][Cube_Clean[j] > 0.])
                mean_signal = np.mean(totsig[totsig > 0.]) 
                SNR = mean_signal/(sigma)
            else:
                print("!!!!!!!This corruption method is unknown, leaving the cube uncorrupted and unconvolved!!!!!!!!")
                # We'll create a little text file with an Overview of all the parameters
            overview = open(work_dir+name+'/'+name+'-Info.txt', 'w')                           
            overview.write("This file contains the basic parameters of this galaxy\n")
            overview.write("For the radial dependencies look at Overview.png or ModelInput.def\n")
            overview.write("Inclination = {}\n".format(Current_Galaxy.inclination))
            overview.write("The dispersion = {:.2f}-{:.2f}\n".format(dispersion[0],dispersion[-1]))
            overview.write("The type of galaxy = {}\n".format(Current_Galaxy.rc_name))
            overview.write("PA = {}\n".format(Current_Galaxy.PA))
            overview.write("Warp = {}-{}\n".format(Current_Galaxy.warp[0],Current_Galaxy.warp[1]))
            overview.write("Which starts at {:.2f} kpc and the 1M/pc^2 radius is {:.2f} kpc \n".format(WarpStart,Rad_HI))
            overview.write("Flare = {}\n".format(Current_Galaxy.flare))
            overview.write("Beams across the major axis = {}\n".format(Current_Galaxy.beams))
            overview.write("SNR Requested = {} SNR Achieved = {}  \n".format(Current_Galaxy.SNR,SNR))
            overview.write("Channelwidth = {}\n".format(Current_Galaxy.Channelwidth))
            overview.write("Major axis beam = {} Minor axis beam= {}\n".format(Current_Galaxy.Res_Beam[0],Current_Galaxy.Res_Beam[1]))
            overview.write("This galaxy has {} and a {}\n".format(Current_Galaxy.Arms,Current_Galaxy.Bar))
            overview.write("It's central coordinates are RA={} DEC={} vsys={:.2f} km/s\n".format(RAhr,DEChr,vsys))
            overview.write("At a Distance of {:.2f} Mpc \n".format(Distance))
            overview.write("HI_Mass Requested {:.2e} (M_solar) and an optical h {:.2f} (kpc)\n".format(MHI,sclength))
            overview.write("HI_Mass Retrieved {:.2e} (M_solar) \n".format(mass))
            overview.write("The cube was corrupted with the {} method \n".format(corruption))
            overview.write("The final noise level is {} Jy/beam \n".format(sigma))            
            overview.write("h_z = {:.3f}-{:.3f} (kpc)".format(h_z[0],h_z[-1]))
            overview.close()    
print("We created {} models".format(number_models))
plt.figure(59)
ax.set_ylim(ymin=0)
ax.set_xlim(xmin=0, xmax=max_rad)
plt.legend(loc='upper right',fontsize=12)
plt.savefig('Rotation_Curves.ps', bbox_inches='tight')
plt.close()
