#!/usr/local/bin/ python3
# -*- coding: utf-8 -*-
#This program is a pyhton script to create a data base of artificial galaxies.
#It first creates a model with tirific at a high resolution and then runs it through casa to get obtain a realistic observation.
# once a numerical list is set in length we can convert it to a numpy array in order to do operations faster.
# first we import numpy

import numpy as np
import common_functions as cf
import warnings
import subprocess
import sys
import copy
import scipy.ndimage as ndimage
import os
import re
from scipy import interpolate
from scipy import integrate
from astropy.io import fits
import matplotlib.pyplot as plt
from matplotlib.ticker import AutoMinorLocator

#This class defines a set of Base galaxy parameters
class Base_Galaxy:
    def __init__(self, num):
        if num == 0:
            self.Inclination = 60.
            self.Dispersion = [14.,7.5] #[13., 7.5]
            self.Mass= 2.5e12 #1e12
            self.PA = 35.
#            self.Warp = [0.1, 0.15]
            self.Warp = [0.,0.]
            self.Flare = "No_Flare" #"Flare"
            self.Beams= 12 #18 #16.
            self.SNR= 8
            self.Channelwidth = 4.
            self.Coord = [3.56406541347E+01,4.23492081422E+01]
            self.Res_Beam = [20.,20.]
            self.Arms = "Arms"
            self.Bar = "No_Bar"
            self.Radial_Motions= 0.
            self.RC_Radii= [0.]
            self.RC_Vrot= [0.]
        elif num == 1:
             # options are inclination, PA, flare, warp, beams, SNR, Channelwidth, Res_Beam, Arms, Bar, Radial_Motions
            self.Inclination = 55.
            self.Dispersion = [13,9.] # [9., 8.]
            self.PA = 45
            self.Warp = [0.,0.] #[0.03, 0.111] # in radians.  in Theta and phi
            self.Flare = "No_Flare"
            self.Beams= 14 #17 #16
            self.SNR= 8.
            self.Channelwidth = 4.
            self.Res_Beam = [15.,15.]
            self.Arms = "Arms"
            self.Bar = "No_Bar"
            self.Radial_Motions= 0.
            # RA and DEC in degrees, Only used in casa sim
            self.Coord = [50.81166937467079,57.76644335595375]
            self.Mass = 2.5e10 #5e11 # in km/s
            self.RC_Radii= [0.]
            self.RC_Vrot= [0.]
        elif num == 2:
            self.Inclination = 65.
            self.Dispersion = [8., 8.]
            self.PA = 145
            self.Warp = [0.05, 0.025] # in radians.
            self.Flare = "Flared"
            self.Beams= 16 #16
            self.SNR= 8.
            self.Channelwidth = 4.
            self.Coord = [50.81166937467079,57.76644335595375]
            self.Res_Beam = [25.,25.]
            self.Arms = "No_Arms"
            self.Bar = "No_Bar"
            self.Radial_Motions= 0.
            self.Mass= 5e11
            self.RC_Radii= [0.]
            self.RC_Vrot= [0.]
        elif num == 3:
            self.Inclination = 48.
            self.Dispersion = [13., 7.5]
            self.PA = 115
            self.Warp = [0.07, 0.15] # in radians.
            self.Flare = "No_Flare"
            self.Beams= 15
            self.SNR= 8.
            self.Channelwidth = 4.
            self.Coord = [50.81166937467079,57.76644335595375]
            self.Res_Beam = [10.,10.]
            self.Arms = "No_Arms"
            self.Bar = "No_Bar"
            self.Radial_Motions= 0.
            self.Mass= 2.5e11
            self.RC_Radii= [0.]
            self.RC_Vrot= [0.]
        elif num == 4:
            self.Inclination = 42.
            self.Dispersion = [15., 12.]
            self.PA = 115
            self.Warp = [0.1, 0.07] # in radians.
            self.Flare = "Flared"
            self.Beams= 14
            self.SNR= 8.
            self.Channelwidth = 4.
            self.Coord = [50.81166937467079,57.76644335595375]
            self.Res_Beam = [10.,10.]
            self.Arms = "No_Arms"
            self.Bar = "No_Bar"
            self.Radial_Motions= 0.
            self.Mass= 5e10
            self.RC_Radii= [0.]
            self.RC_Vrot= [0.]
        else :
            print("There are not {} bases".format(num))
            sys.exit()
# First we define some useful routines for converting and copying
# A routine to copy disks
def copy_disk(olddisk,newdisk):
    start = 0
    startlast = 0.
    if int(Template["NDISKS"].split('=')[1]) < newdisk:
        Template["NDISKS"] = "NDISKS = {:d}".format(newdisk)
    copkeys ="Empty"
    last = 'RADI'
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
        var_name = key.split('_')[0]
        Template.insert(last,var_name+"_{:d}".format(newdisk),var_name+"_{:d} =".format(newdisk)+Template[key].split('=')[1])
    Template.insert("CFLUX_{:d}".format(newdisk-1),"CFLUX_{:d}".format(newdisk),"CFLUX_{:d} =".format(newdisk)+Template["CFLUX_{:d}".format(newdisk-1)].split('=')[1])
    if '_' in copkeys[-1]:
        return copkeys[-1].split('_')[0]+"_{:d}".format(newdisk)
    else:
        return copkeys[-1]+"_{:d}".format(newdisk)

# Obtaining the derivative at any give point of a function
def derivative(x,func):
    h=x/1000.
    der=(func(x+h)-func(x-h))/(2*h)
    return der





# Then the functions that create the different components of the galaxy



# function to get the surface brightness profile based on the last element in the rotation curve
def build_sbr_prof(Current_Galaxy):
    # First we use the mass to build a rotation curve and sbr profile
    # Avanti used 	V_Ropt = (200*(l**0.41))/(0.80 +0.49*np.log10(l) +(0.75*np.exp(-0.4*l))/(0.47 + 2.25*(l**0.4)))**0.5 (Persic & Sallucci)
    # We will not use the URC as it get's silly at high mass and remains showing flat parts at low mass
    # We will use the parameterisation of Courteau 1997
    # G=6.674 x 10^neg20 #km^3⋅kg^neg1⋅s^neg2 pc=3.086e+13 km solarmass=1.98855e30 kg
    # transform (pc x km^2)/(s^2 x solarmass)
    G= 4.30058*10**-3
    # About 18.5% is cosmic baryon fraction from planck (Omb=0.048, OmDM = 0.2589). In a galactic halo this is roughly 70%-90% of the cosmic mean (See Crain et al. 2006, Pezzulli 2019)
    bary_frac = 0.1854
    diff=0.9
    counter = 0
    while abs(diff) > 0.1:
    # if above 10**10 then -10**9 to account for gas
        if Current_Galaxy.Mass > 10**10:
            #m_star = bary_frac* Current_Galaxy.Mass-(bary_frac*Current_Galaxy.Mass-10**9)*10**(-0.43*np.log10((bary_frac*Current_Galaxy.Mass-10**9))+3.75)
            m_star=bary_frac* Current_Galaxy.Mass-10**9
        else:
            # else half baryonic in stars
            #m_star  = bary_frac*Current_Galaxy.Mass-(bary_frac*Current_Galaxy.Mass/2.)*10**(-0.43*np.log10((bary_frac*Current_Galaxy.Mass/2.))+3.75)
            m_star=bary_frac* Current_Galaxy.Mass/2.
        # and convert stellar mass to an HI Mass following van Driel (2016)
        # log(MH I/M?) = −0.43 log (M?) + 3.75.
        m_star_prev=1
        while abs(m_star_prev-m_star)/m_star > 0.1:
            MHI = m_star*10**(-0.43*np.log10(m_star)+3.75)
            m_star_prev=copy.deepcopy(m_star)
            m_star=bary_frac* Current_Galaxy.Mass-MHI*1.4
        print("This is MHI {:.2e} and mstar {:.2e}".format(MHI,m_star,m_star_prev))
        # the mass leads to the radius at which we need to hit 1 M/pc^2 from Wang (2016) in kpc
        HIrad=10**(0.506*np.log10(MHI)-3.293)/2.
        # let's calculate the final mass we obtain
        masses = (MHI+m_star)/bary_frac
        # McGaugh 2014  M_*/Lk = 0.6
        # to get a B-band magnitude but mass picking are in K-band
        # Ponomareva 2017 has
        # M^T,b,i_[3.6] = (−9.52 ± 0.32) × log(2V_flat) + 3.3 ± 0.8
        #    So can we find a 3.6 scaling to HI
        # Mag_36 = -9.52*np.log10(2*Base_Galaxy(num).rc_speed[-1])+3.3
        # in luminosity L_solar
        # log(LT ,b,i[3.6] ) = (3.7 ± 0.11) × log(2Vflat) + 1.3 ± 0.3,
        #L_36 = 10**(3.7*np.log10(2.*Base_Galaxy(num).rc_speed[-1])+1.3)
        MK = np.log10(m_star/0.6)*-2.5+3.29
        v_c=10**((MK-1.44)/-9.77)/2.
        #  Let's check that this roughly gives our DM mass at the virial radius in pc
        r200= (Current_Galaxy.Mass*G/(100.* H_0**2)*1e12)**(1./3.)
        v200square=Current_Galaxy.Mass*G/r200
        xr=HIrad/(r200/10**3)
        c200=10**(0.905-0.101*np.log10(Current_Galaxy.Mass/(10**12*100./H_0))) #From A. Dutton 2014
        NFWvflat=np.sqrt(v200square*(1./xr)*((np.log(1.+c200*xr)-(c200*xr)/(1+c200*xr))/(np.log(1+c200)-c200/(1+c200))))
        v_star=np.sqrt(m_star*G/(HIrad*10**3))
        v_HI=np.sqrt(MHI*1.4*G/(HIrad*10**3))
        v_c_tot=np.sqrt(NFWvflat**2+v_star**2+v_HI**2)
        DynMassNFW=HIrad*10**3*v_c_tot**2/G

        # The core radius should be smaller for more massive galaxies with a minimum of 1.5 kpc see transition function .py for how these parameters run
       # Made to match the sparcs data set
        height=1.5
        center=235
        disp=30.
        z=(v_c-center)/disp
        R_0 =height*np.exp(-1*z**2/2.)*(np.arctan((v_c-center+disp)/2.)+0.4)*0.7+(HIrad-HIrad/10.)/(v_c/60.)**2.+HIrad/20.
        # The behaviour at large radii should be to grow slowly above 120 km/s with a maximum of 0.5
        # above 212 v_c it should decline again
        beta = np.arctan((v_c-120.)/60.)*0.6
        if beta < 0: beta = 0.
        # The turnover sharpness should increase while beta < 0.45 above that it should decline again
        height=0.9
        center=212
        disp=100.
        z=(v_c-center)/disp
        gamma_tmp =height*np.exp(-1*z**2/2.)+1.4
        gamma_tmpagain=-0.2*(np.arctan((v_c-300)/30.)-0.5)
        gamma=gamma_tmp+gamma_tmpagain
        x=R_0/HIrad
        vfinal=v_c*(1+x)**beta*(1+x**gamma)**(-1./gamma)

        #diff=(Current_Galaxy.Mass*(3*xr)-DynMass)/DynMass
        print("The NFW velocity at R_HI = {}  and the retrieved curve velocity at R_HI = {} diff {}".format(v_c_tot,vfinal,v_c_tot-vfinal))

        diff=(vfinal-v_c_tot)/vfinal
        print("The current fractional difference = {} and the baryon fraction = {}".format(diff,bary_frac))

        if abs(diff) > 2.:
            diffch = 2.
        else:
            diffch = abs(diff)
        if diff < 0.:
            bary_frac = bary_frac+0.05*diffch
        else:
            bary_frac = bary_frac-0.05*diffch
        if bary_frac < 0.01:
            bary_frac= 0.01
        counter += 1
        if bary_frac >  0.1854*1.1 or bary_frac < 0.1854*0.2:
            counter=101
        if counter > 100:
            v_c=(v_c+v_c_tot)/2.
               # The core radius should be smaller for more massive galaxies with a minimum of 1.5 kpc see transition function .py for how these parameters run
            # Made to match the sparcs data set
            height=1.5
            center=235
            disp=30.
            z=(v_c-center)/disp
            R_0 =height*np.exp(-1*z**2/2.)*(np.arctan((v_c-center+disp)/2.)+0.4)*0.7+(HIrad-HIrad/10.)/(v_c/60.)**2.+HIrad/20.
            # The behaviour at large radii should be to grow slowly above 120 km/s with a maximum of 0.5
            # above 212 v_c it should decline again
            beta = np.arctan((v_c-120.)/60.)*0.6
            if beta < 0: beta = 0.
            # The turnover sharpness should increase while beta < 0.45 above that it should decline again
            height=0.9
            center=212
            disp=100.
            z=(v_c-center)/disp
            gamma_tmp =height*np.exp(-1*z**2/2.)+1.4
            gamma_tmpagain=-0.2*(np.arctan((v_c-300)/30.)-0.5)
            gamma=gamma_tmp+gamma_tmpagain
            break
    DynMass=HIrad*10**3*v_c**2/G
    print("The input Mass = {:.2e}  and the retrieved NFW Dynamical mass = {:.2e} and Dynamical Mass based on v_c = {:.2e}".format(Current_Galaxy.Mass,DynMassNFW,DynMass))
    print("The current fractional difference = {:.5f} and the baryon fraction = {:.5f}".format(diff,bary_frac))
        #First we set the radii at with a 5 elements per rings out to 1.5 times HIrad. However we always need at least 25 rings for transition purposes
    if Current_Galaxy.Beams < 5:
        sub_ring = int(np.ceil(25./Current_Galaxy.Beams))
    else:
        sub_ring = 5
    Rad= np.arange(0,1.5*HIrad, 1.5*HIrad/(sub_ring*Current_Galaxy.Beams))
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        x = R_0/Rad
    x[0]=1e-8
    # Finaaly the courteau presciption
    vrot = v_c*(1+x)**beta*(1+x**gamma)**(-1./gamma)
    vrot[0] = 0.
    # from these we use the the prescription of Martinsson to get an HI profile.
    # As done by A. Gogate in initialparam.py
    # First we set the radii at with a 5 elements per rings out to 1.5 times HIrad
    a = (Rad) * 10 ** 3  # in parsecs
    fracHIrad = (HIrad - np.floor(HIrad))
    # hole in the center should be amplified because of the conversion to H_2
    # with a turn over at ~ 120 V_max where it should start adding back
    # the H_2 follows an exponential distribution similar to the stellar disk (ref?)
    # The scale length is h=Vmax^2/(0.88**2*np.pi*G*sig0) (Freeman 1970)

    # The scale length for the molecular disk is then from Graham rederived relation in cal_scl.py
    h_r = (-4.13422991 - 0.31576291 * MK) * 1000.
    if fracHIrad < 0.5:
        Rhi_p = np.array([(HIrad - HIrad / 15.) * 10 ** 3, (HIrad + HIrad / 15.) * 10 ** 3])
        # Rhi_p2 = HIrad + HIrad / 20. * 10 ** 3
        h_r = np.array([h_r + h_r / 5, h_r - h_r / 5])
        # h_r2 = h_r - h_r / 20.
    else:
        Rhi_p = np.array([(HIrad + HIrad / 15.) * 10 ** 3, (HIrad - HIrad / 15.) * 10 ** 3])
        # Rhi_p2 = HIrad - HIrad / 20. * 10 ** 3
        h_r = np.array([h_r - h_r / 5., h_r + h_r / 5.])
        # h_r2 = h_r + h_r / 20.
    Hiradindex = [0, 0]
    Hiradindex[0] = np.where((Rad > (Rhi_p[0] - Rhi_p[0] / (sub_ring * Current_Galaxy.Beams)) / 10 ** 3) & (
            Rad < (Rhi_p[0] + Rhi_p[0] / (sub_ring * Current_Galaxy.Beams)) / 10 ** 3))[0][0]
    Hiradindex[1] = np.where((Rad > (Rhi_p[1] - Rhi_p[1] / (sub_ring * Current_Galaxy.Beams)) / 10 ** 3) & (
            Rad < (Rhi_p[1] + Rhi_p[1] / (sub_ring * Current_Galaxy.Beams)) / 10 ** 3))[0][0]

    # print("This the HI Radius in kpc {}".format(HIrad))
    # std deviation or dispersion of gaussian
    s1 = np.zeros(2)
    s1 = 0.36 * Rhi_p
    # s2 = 0.36 * Rhi_p2
    # We assume that at 120 km/s v_max the H_2 disk imprint on the HI disk is adequately described by the prescription of Martinsson.
    # Lower some of the disk is added back to the HI higher the central hole is more pronounced
    I_cen = ((v_c / 120.) ** 0.5 - 1)
    print("This the scale length {} and central brightness {}".format(h_r, I_cen))
    # So our molecular profile is
    Exp = np.zeros([len(a), 2])
    Exp[:, 0] = I_cen * np.exp(-a / h_r[0])
    Exp[:, 1] = I_cen * np.exp(-a / h_r[1])


    # gaussian2 From Martinsson 2015
    Sig2 = np.zeros([len(a), 2])
    Sig2[:, 0] = np.exp(-(a - 0.4 * Rhi_p[0]) ** 2 / (2 * (s1[0]) ** 2))
    Sig2[:, 1] = np.exp(-(a - 0.4 * Rhi_p[1]) ** 2 / (2 * (s1[1]) ** 2))
    # total
    Sigma = np.zeros([len(a), 2])
    Sigma = Sig2 - Exp
    Sigma[Sigma < 0.] = 0.  # for negative sigma max, does not include negative values
    # scale Sigma such that it is one at HI rad
    new = np.zeros(2)
    new[0] = 1. / Sigma[Hiradindex[0], 0]
    new[1] = 1. / Sigma[Hiradindex[1], 1]
    Sigma[:, 0] = new[0] * Sigma[:, 0]
    Sigma[:, 1] = new[1] * Sigma[:, 1]
    Sigma[0:4,0]=(Sigma[0:4,0]+Sigma[0:4,1])/2.
    Sigma[0:4, 1] = Sigma[0:4, 0]
    # get the HI Mass in the profile
    OutHIMass = integrate.simps((np.pi * a) * Sigma[:, 0], a) + integrate.simps((np.pi * a) * Sigma[:, 1], a)
    # And check that it matches the required HI mas within 5%
    counter = 1
    while np.absolute(MHI - OutHIMass) > MHI / 50.:
        # if not rescale sigma1
        if MHI - OutHIMass > 0:
            s1 = (0.36 - (0.0025 * counter)) * Rhi_p
        else:
            s1 = (0.36 + (0.0025 * counter)) * Rhi_p
        # and recalculate the profile
        Sig2[:, 0] = np.exp(-(a - 0.4 * Rhi_p[0]) ** 2 / (2 * (s1[0]) ** 2))
        Sig2[:, 1] = np.exp(-(a - 0.4 * Rhi_p[1]) ** 2 / (2 * (s1[1]) ** 2))
        Sigma = Sig2 - Exp
        Sigma[Sigma < 0.] = 0.
        new[0] = 1. / Sigma[Hiradindex[0], 0]
        new[1] = 1. / Sigma[Hiradindex[1], 1]
        Sigma[:, 0] = new[0] * Sigma[:, 0]
        Sigma[:, 1] = new[1] * Sigma[:, 1]
        Sigma[0:4,0]=(Sigma[0:4,0]+Sigma[0:4,1])/2.
        Sigma[0:4, 1]=Sigma[0:4,0]
        OutHIMass = integrate.simps((np.pi * a) * Sigma[:, 0], a) + integrate.simps((np.pi * a) * Sigma[:, 1], a)
        counter += 1
    print("The requested HI mass = {:.2e} and the retrieved HI mass = {:.2e}".format(MHI, OutHIMass))
    # final HI radial distribution by renormalisation
    S = Sigma * (1.24756e+20)
    # Where A. Gogate contribution stops
    # S is column densities but tirific takes Jy * km/s/arcsec^2 so
    conv_column_arsec = 605.7383 * 1.823E18 * (2. * np.pi / (np.log(256.)))  # S(mJy/beam)*conv_column_arcsec=N_HI
    sbr_prof = S / (conv_column_arsec * 1000.)
    # Let's write these to the Template immediately
    # The surface brightness profile, which is still symmetric
    Template["SBR"] = "SBR = " + " ".join(str(e) for e in sbr_prof[:, 0])
    Template["SBR_2"] = "SBR_2 = " + " ".join(str(e) for e in sbr_prof[:, 1])
    # And we plot the  profile to an overview plot
    overview.plot(Rad,sbr_prof[:,0],'k')
    overview.plot(Rad[0::sub_ring],sbr_prof[0::sub_ring,0],'ko')
    overview.plot(Rad, sbr_prof[:, 1], 'r')
    overview.plot(Rad[0::sub_ring], sbr_prof[0::sub_ring, 1], 'ro')
    overview.plot(Rad,Exp[:,0]*new[0]*1.24756e+20/(conv_column_arsec*1000.),'b')
    overview.plot(Rad[0::sub_ring],Exp[0::sub_ring,0]*new[0]*1.24756e+20/(conv_column_arsec*1000.),'bo')
    overview.plot(Rad, Exp[:,1] * new[1] * 1.24756e+20 / (conv_column_arsec * 1000.), 'y')
    overview.plot(Rad[0::sub_ring], Exp[0::sub_ring,1] * new[1] * 1.24756e+20 / (conv_column_arsec * 1000.), 'yo')
    ymin,ymax=plt.ylim()
    plt.margins(x=0., y=0.)
    labelfont= {'family':'Times New Roman',
                'weight':'normal',
                'size':18}
    plt.rc('font',**labelfont)
    plt.plot([HIrad,HIrad],[ymin-(ymax-ymin)*0.1,ymax+(ymax-ymin)*0.1],'b')
    plt.plot([Rhi_p[0] / 10 ** 3, Rhi_p[0] / 10 ** 3], [ymin - (ymax - ymin) * 0.1, ymax + (ymax - ymin) * 0.1], 'b--')
    plt.plot([Rhi_p[1] / 10 ** 3, Rhi_p[1] / 10 ** 3], [ymin - (ymax - ymin) * 0.1, ymax + (ymax - ymin) * 0.1], 'b--')

    #0.,np.max(sbr_prof)],'b')
    plt.xlabel('Radius (kpc)',**labelfont)
    plt.ylabel('SBR (Jy km s$^{-1}$ arcsec$^-2$)',**labelfont)
    # as the rest on of the code is based on a single SBR profile let's average the value
    HIrad = np.mean(Rhi_p) / 10 ** 3
    print("The average HI radius = {} from {} in the appraoching and {} in the receding side".format(HIrad,
                                                                                                     Rhi_p[0] / 10 ** 3,
                                                                                                     Rhi_p[
                                                                                                         1] / 10 ** 3))
    h_r = np.mean(h_r)
    tmp = np.zeros(len(sbr_prof[:, 0]))
    for i in range(len(sbr_prof)):
        tmp[i] = np.mean(sbr_prof[i, :])
    return tmp,Rad,h_r/1000.,OutHIMass, HIrad,vrot,sub_ring


# Thi function creates the flarin g which is based on the rotation curve and dispersion according Puche et al. 1992
def create_flare(Radii,velocity,dispersion,flare,Max_Rad,sub_ring,distance=1.):
    #make sure we have arrays so we can do numerical operations
    Radii = np.array(Radii)
    velocity = np.array(velocity)
    # first we need the dispersion
    # if there is no change from inside to outside it is easy
    if dispersion[0] == dispersion[1]:
       disp=np.full(len(Radii),dispersion[0])
    else:
        # else we use a tangent change with the center at halfway
       disp=-1*np.arctan((Radii-np.mean(Max_Rad/2.))/(Radii[-1]/10.))/np.pi*np.absolute(dispersion[0]-dispersion[1])+np.mean(dispersion)
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
        fact=np.arange(1/21,1,1./21)
        flare[halfint-10:halfint+10] = (1-fact)*flare[halfint-10]+fact*flare[halfint-10:halfint+10]
    elif flare.lower() == 'no_flare':
        flare = np.full(len(Radii),disp[halfint]/((4.*np.pi*G2*Density[halfint]/1000.**3)**0.5*3.086e+16))
    else:
        print("{} is not an option for the flare. Choose Flared or No_Flare".format(flare))
        sys.exit()

    flare[0]=flare[1]
    plt.figure(2)
    plt.subplot(6,1,2)
    labelfont= {'family':'Times New Roman',
                'weight':'normal',
                'size':18}
    plt.rc('font',**labelfont)
    plt.plot(Radii,disp,'k')
    plt.plot(Radii[0::sub_ring],disp[0::sub_ring],'ko')
    ymin,ymax=plt.ylim()
    plt.plot([Max_Rad,Max_Rad],[ymin-(ymax-ymin)*0.1,ymax+(ymax-ymin)*0.1],'b')
    plt.margins(x=0., y=0.0)
    plt.tick_params(
        axis='x',          # changes apply to the x-axis
        which='both',      # both major and minor ticks are affected
        direction = 'in',
        bottom=True,      # ticks along the bottom edge are off
        top=False,         # ticks along the top edge are off
        labelbottom=False) # labels along the bottom edge are off
    plt.ylabel('Disp. (km s$^{-1}$)',**labelfont)
    print("The dispersion runs from {:5.2f} to {:5.2f} km/s".format(disp[0],disp[-1]))
    plt.subplot(6,1,1)
    plt.plot(Radii,flare,'k')
    plt.plot(Radii[0::sub_ring],flare[0::sub_ring],'ko')
    ymin,ymax=plt.ylim()
    plt.plot([Max_Rad,Max_Rad],[ymin-(ymax-ymin)*0.1,ymax+(ymax-ymin)*0.1],'b')
    plt.margins(x=0., y=0.0)
    plt.tick_params(
        axis='x',          # changes apply to the x-axis
        which='both',      # both major and minor ticks are affected
        direction = 'in',
        bottom=True,      # ticks along the bottom edge are off
        top=False,         # ticks along the top edge are off
        labelbottom=False) # labels along the bottom edge are off
    plt.ylabel('Scale height (kpc)',labelpad=30,**labelfont)
    #plt.xticks([])
    print("The flare runs from {:10.6f} to {:10.6f} kpc".format(flare[0],flare[-1]) )
    # convert the scale heights to arcsec
    h_z_arcsec = cf.convertskyangle(flare,distance=distance,physical = True)
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
                sub_ring,
                disk=1):
    Radii =np.array(Radii)
    if ((np.sum(warp_change) != 0) and (warp_radii[0] < warp_radii[1])).all():
        # First we need to obtain the vector that constitutes the inner area
        #it runs exactly counter to inclination
        inclination=90-inclination
        # For this the PA has to be between 0-90
        mult=np.floor(PA/90.)
        inPA= PA-mult*90.
        #avoid singularities
        if inPA == 0.:
            inPA = 0.001
        if inclination < 0.001 :
            inclination = 0.001
        # define the angular momentum vector of the plane and the outer most ring
        theta=np.arctan(np.tan(inclination*(np.pi/180.))*np.tan(inPA*(np.pi/180.)))
        phi=np.arctan(np.tan(inPA*(np.pi/180.))/np.sin(theta))
        # and the maximum values at Rad_HI
        thetamax=theta+warp_change[0]
        phimax=phi+warp_change[1]
                # indices of the warp start and the Rad_HI
        start_index = int(np.sum(Radii < warp_radii[0]))
        end_index = int(np.sum(Radii < warp_radii[1]))
        # step size of theta and phi
        # As we will increase the step size triangular we need the total number of point in the sequence
        warprange = end_index-start_index
        if warprange < 2:
            thetamax=theta
            phimax=phi
            warprange=1
        increasetheta=(thetamax-theta)/(0.5*warprange*(warprange+1))
        increasephi=(phimax-phi)/(0.5*warprange*(warprange+1))
        #print(warprange,thetamax,phimax,Radii[1]-Radii[0],warp_radii[1]-warp_radii[0])
        # calculate theta
        thetarings = np.array(np.full(len(Radii),theta))
        index_array=np.arange(start_index,len(Radii))-start_index
        thetarings[start_index:] = theta+0.5*index_array*(index_array+1)*increasetheta
        #calculate phi
        phirings = np.array(np.full(len(Radii),phi))
        phirings[start_index:] = phi+0.5*index_array*(index_array+1)*increasephi
        # return to PA

        if (phirings[0] < np.pi/2.) and (phirings[-1] > np.pi/2) and (inclination < 5.):
            PA= np.arctan(np.sin(thetarings)*np.tan(phirings-np.pi/2.))*(360./(2*np.pi))+mult*90+np.arctan(np.sin(theta)*np.tan(phi))*(360./(2*np.pi))
        else:
            PA= np.arctan(np.sin(thetarings)*np.tan(phirings))*(360./(2*np.pi))+mult*90
            PA[phirings > 0.5*np.pi]=PA[phirings > 0.5*np.pi]+180.

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

    if disk == 1:
        c='k'
    else:
        c='r'
        warnings.filterwarnings("ignore")
    plt.subplot(6,1,4)
    labelfont= {'family':'Times New Roman',
                'weight':'normal',
                'size':18}
    plt.rc('font',**labelfont)
    plt.plot(Radii,PA,c)
    plt.plot(Radii[0::sub_ring],PA[0::sub_ring],c+'o')
    ymin,ymax=plt.ylim()
    plt.plot([warp_radii[1],warp_radii[1]],[ymin-(ymax-ymin)*0.1,ymax+(ymax-ymin)*0.1],'g')
    plt.margins(x=0., y=0.)
    plt.tick_params(
        axis='x',          # changes apply to the x-axis
        which='both',      # both major and minor ticks are affected
        direction = 'in',
        bottom=True,      # ticks along the bottom edge are off
        top=False,         # ticks along the top edge are off
        labelbottom=False) # labels along the bottom edge are off
    plt.ylabel('PA ($^{\circ}$)',**labelfont)
    plt.subplot(6,1,3)
    plt.plot(Radii,inc,c)
    plt.plot(Radii[0::sub_ring],inc[0::sub_ring],c+'o')
    ymin,ymax=plt.ylim()
    plt.plot([warp_radii[1],warp_radii[1]],[ymin-(ymax-ymin)*0.1,ymax+(ymax-ymin)*0.1],'g')
    plt.margins(x=0., y=0.0)
    plt.tick_params(
        axis='x',          # changes apply to the x-axis
        which='both',      # both major and minor ticks are affected
        direction = 'in',
        bottom=True,      # ticks along the bottom edge are off
        top=False,         # ticks along the top edge are off
        labelbottom=False) # labels along the bottom edge are off
    plt.ylabel('Inc. ($^{\circ}$)',**labelfont)

    #write to our template file
    # let's see if we can retrace intrinsic phi with the formula's from Peters
    #theta_test = np.arctan((np.sin(PA*np.pi/180.)*np.sin(inc*np.pi/180.)-np.cos(PA*np.pi/180.)*np.sin(inc*np.pi/180.))/(np.cos(PA*np.pi/180.)*np.cos(inc*np.pi/180.)-np.sin(PA*np.pi/180.)*np.cos(inc*np.pi/180.)))*180./np.pi
    phirings=phirings*180./np.pi
    # thetarings=thetarings*180./np.pi
    # This seems to work mostly but not at some extremes exactly for some reason
    # According to josh tan is required, and yes that makes it work at large angles as well.
    angle_adjust=np.tan((PA[0]-PA)*np.cos(inc*np.pi/180.)*np.pi/180)*180/np.pi
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


def create_arms(velocity,Radii,disk_brightness, disk=1,WarpStart=-1,Bar="No_Bar"):
    if WarpStart == -1: WarpStart = Radii[-1]
    max_vrot=np.max(velocity)
    max_rad=Radii[-1]
    #The pattern speed at a given radius is vrot/radii
    V_Rot = interpolate.interpolate.interp1d(Radii,velocity,fill_value="extrapolate")
    # The radius of co-ration can be approximated by the extend of the visible disk ~ Warpstart (Roberts et al. 1975)
    Omega_CR=V_Rot(WarpStart)/WarpStart
    # From this we can estimate the inner and outer Lindblad resonances (Eq 38 Dobbs & Baba 2014)
     #The epicyclic frequency k^2=R*dOmega^2/dR+4*Omega^2
    # f'(x) = (f(x+h)-f(x-h))/2h
    h=WarpStart/1000.
    derive=(V_Rot(float(WarpStart+h))**2/(WarpStart+h)**2-V_Rot(float(WarpStart-h))**2/(WarpStart-h)**2)/(2*h)
    k_CR = (WarpStart * derive+4*Omega_CR**2)**0.5
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
    last_add = copy_disk(disk,ndisk)
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
    if WarpStart == -1: WarpStart = Radii[-1]
    max_rad = Radii[-1]
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
    k_CR = (WarpStart * derive+4*Omega_CR**2)**0.5
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
    last_add = copy_disk(disk,ndisk)
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
def create_inhomogeneity(mass,SNR,disks=1.):
    if len(disks) == 0:
        disks=[disks]

    for i in range(len(disks)):
        ndisk=int(Template["NDISKS"].split('=',1)[1])
        ndisk+=1
        last_add = copy_disk(disks[i],ndisk)
        sbr=np.array(Template["SBR_{:d}".format(ndisk)].split('=',1)[1].strip().split(" "),dtype=np.float64)
        rad=np.array(Template["RADI"].split('=',1)[1].strip().split(" "),dtype=np.float64)
        rad=rad*4.*np.pi/max(rad)
        # This part is tricky as we do not want negative emission but cutting out the negatives would result in added flux
        #Hence these profiles should always be smaller than the sbr
        newprof=np.sin(rad)*sbr*np.log10(mass)/11.*8./SNR*0.9
        print("We are modifying by factor {} for mass {} and SNR {}".format(np.log10(mass)/11.*8./SNR, np.log10(mass),SNR))
        nextprof=-1*newprof

        req_flux=abs(np.sum(newprof*np.pi*rad*2)/200.)
        if req_flux == 0:
            req_flux = 1e-5
        Template["SBR_{:d}".format(ndisk)]="SBR_{:d} = ".format(ndisk)+" ".join(str(e) for e in newprof)
        Template["CFLUX_{:d}".format(ndisk)]="CFLUX_{:d} = {}" .format(ndisk,req_flux)
        ndisk=ndisk+1
        last_add = copy_disk(disks[i],ndisk)
        Template["SBR_{:d}".forprint(mat(ndisk)]="SBR_{:d} = ".format(ndisk)+" ".join(str(e) for e in nextprof)
        Template["CFLUX_{:d}".format(ndisk)]="CFLUX_{:d} = {}" .format(ndisk,req_flux)
    print("We will create inhomogeneities on top of disk(s) ="+" ".join(str(e) for e in [disks]))


# This is a routine to corrupt the cube with some standard gaussian noise
def corrupt_gauss(work_dir,beam,SNR):
    # First open the model cube
    dummy = fits.open(work_dir+'/unconvolved_cube.fits',uint = False, do_not_scale_image_data=True,ignore_blank = True)
    # Add to the hedear some info
    dummy[0].header.append('RESTFRQ')
    dummy[0].header['RESTFRQ'] = 1.420405752E+09
    # downgrade the velocity resolution to the one we want
    tmp_cube=np.zeros([int(len(dummy[0].data[:,0,0])/3),len(dummy[0].data[0,:,0]),len(dummy[0].data[0,0,:])])
    for n in range(1,len(dummy[0].data[:,0,0]),3):
        tmp_cube[int((n-1)/3),:,:] = np.mean(dummy[0].data[n-1:n+2,:,:],axis=0)
    dummy[0].data = tmp_cube
    print(dummy[0].data.shape)
    # reset the header
    dummy[0].header['CDELT3'] = dummy[0].header['CDELT3']*3.
    dummy[0].header['CRPIX3'] = dummy[0].header['CRPIX3']/3.
    dummy[0].header['NAXIS3'] = int(dummy[0].header['NAXIS3']/3.)
    # In order to create a cleaning mask we smooth to the beam size and cut at 1e-5 Jy/beam
    #   Calculate the sigma's from the required beam size
    sigma=[(beam[0]/abs(dummy[0].header['CDELT1']*3600.))/(2*np.sqrt(2*np.log(2))),(beam[1]/abs(dummy[0].header['CDELT2']*3600.))/(2*np.sqrt(2*np.log(2)))]
    #
    reals = dummy[0].data[dummy[0].data > 1e-5]
    # Remember that python is a piece of shiit so z,y,x
    cutoff = (np.mean(reals)/2.)
    # correct the cutoff for the smoothing
    cutoff=cutoff/(0.5*(2*np.sqrt(np.pi)*sigma[0]+sigma[1]*np.sqrt(np.pi)*2.))
    # smooth the image
    smooth = ndimage.gaussian_filter(dummy[0].data, sigma=(0,sigma[1], sigma[0]), order=0)

    # We want the mean signal in the smoothed cube
    # !!!!! This one is not correct for conserved surface brightness temperature but we want it to estimate the noise per pixel!!!!!!!!
    mean_signal= np.mean(smooth[smooth > cutoff])
    # let's check that the minimum is less than the noise.
    #if np.min(smooth) < mean_signal/SNR*-3:
    #    print(np.min(smooth),mean_signal/SNR*-3)
    #    exit()
    print("This is the minimum signal {} and comparable noise {}".format( np.min(smooth),mean_signal/SNR*-3))
    # Then create the mask
    smooth[smooth < cutoff]=0
    smooth[smooth > cutoff]=1
    # Write mask
    fits.writeto(work_dir+'/mask.fits',smooth,dummy[0].header, overwrite = True)
    # Calculate the area of the beam in arcsec
    beamarea=(np.pi*abs(beam[0]*beam[1]))/(4.*np.log(2.))
    # Convert arcsec to pixels
    pixperbeam=beamarea/(abs(dummy[0].header['CDELT1']*3600.)*abs(dummy[0].header['CDELT2']*3600.))
    # Make an the size of the model with random values distributed as the noise
    # The noise we want is the mean signal divided by the signal to noise ratio
    # that value needs to be deconvolved so from https://en.wikipedia.org/wiki/Gaussian_blur
    # The formula is uncited but works
    noisescl=(mean_signal/SNR*sigma[0]*2*np.sqrt(np.pi))
    print("This our mean signal {} and noise {} in Jy/pixel. ".format(mean_signal,noisescl))
    cuberms = np.random.normal(scale=noisescl,size=np.shape(dummy[0].data))
    # combine the two cubes
    noisedcube=dummy[0].data+cuberms
    # Smooth to the requred resolution
    print("The beam in pixels {} x {}".format(sigma[0]*(2.*np.sqrt(2*np.log(2))),sigma[1]*(2.*np.sqrt(2*np.log(2)))))
    final = ndimage.gaussian_filter(noisedcube, sigma=(0,sigma[1], sigma[0]), order=0)
    # to preserve brightness temperature this should be multiplied with the increase in beam area
    # which is the same as the amount of pixels in the beam as we go from 1 pixel area to an area the size of the beam which is assuming two gaussians so we need to correct with a difference factor of the area of the. Last factor is to correct for the fact that the unconvolved cube hassquare pixels not a circular beam.
    final=final*pixperbeam
    # And write this final cube to the directory
    dummy[0].header['BMAJ']=beam[0]/3600.
    dummy[0].header['BMIN']=beam[1]/3600.
    fits.writeto(work_dir+'/Convolved_Cube.fits',final,dummy[0].header, overwrite = True)

# this is a routine to use the casa task simobserve to corrupt the observations.
def corrupt_casa(work_dir,beam,Template_Casa,SNR):
    #!!!!!!!!!This needs better testing
    #casa takes issue with the tirific header.
    dummy = fits.open(work_dir+'/unconvolved_cube.fits',uint = False, do_not_scale_image_data=True,ignore_blank = True)
    dummy[0].header.append('RESTFRQ')
    dummy[0].header['RESTFRQ'] = 1.420405752E+09
    # downgrade the velocity resolution to the one we want
    tmp_cube = np.zeros(
        [int(len(dummy[0].data[:, 0, 0]) / 3), len(dummy[0].data[0, :, 0]), len(dummy[0].data[0, 0, :])])
    for n in range(1, len(dummy[0].data[:, 0, 0]), 3):
        tmp_cube[int((n - 1) / 3), :, :] = np.mean(dummy[0].data[n - 1:n + 2, :, :], axis=0)
    dummy[0].data = tmp_cube
    # reset the header
    dummy[0].header['CDELT3'] = dummy[0].header['CDELT3'] * 3.
    dummy[0].header['CRPIX3'] = dummy[0].header['CRPIX3'] / 3.
    dummy[0].header['NAXIS3'] = int(dummy[0].header['NAXIS3'] / 3.)
    fits.writeto(work_dir+'/unconvolved_cube.fits',dummy[0].data,dummy[0].header, overwrite = True)
    # Calculate the required noise in the end cube
    sigma=[(beam[0]/abs(dummy[0].header['CDELT1']*3600.))/(2*np.sqrt(2*np.log(2))),(beam[1]/abs(dummy[0].header['CDELT2']*3600.))/(2*np.sqrt(2*np.log(2)))]
    #
    reals = dummy[0].data[dummy[0].data > 1e-5]
    # Remember that python is a piece of shiit so z,y,x
    cutoff = (np.mean(reals)/2.)
    # correct the cutoff for the smoothing
    cutoff=cutoff/(0.5*(2*np.sqrt(np.pi)*sigma[0]+sigma[1]*np.sqrt(np.pi)*2.))
    # smooth the image
    smooth = ndimage.gaussian_filter(dummy[0].data, sigma=(0,sigma[1], sigma[0]), order=0)
    # We want the mean signal in the smoothed cube
    # !!!!! This one is not correct for conserved surface brightness temperature but we want it to estimate the noise per pixel!!!!!!!!
    mean_signal= np.mean(smooth[smooth > cutoff])
    # Then create the mask
    smooth[smooth < cutoff]=0
    smooth[smooth > cutoff]=1
    # Write mask
    fits.writeto(work_dir+'/mask.fits',smooth,dummy[0].header, overwrite = True)

    # In order to corrupt we need to know the average signal.
    # we do this taking the mean in each chaneel above a tenth of the max and then take the mean of that profile
    # This is the noise in the final cube
    noise=mean_signal/SNR
    #*(abs(dummy[0].header['CDELT1']*3600.)*abs(dummy[0].header['CDELT2']*3600.))/1.1330900354567984


    # We need to know the location in the sky dec > 30 we use WSRT 30 > dec > -30 --> VLA -30 > dec  --> Atca

    # Full synthesis 12 hrs leads to noise levels 12: 0.000701111 ,24:  0.00065377   ,48:0.000499886 ,96:0.000436001 ,300:0.000213708 ,600:0.000178363   ,900:0.000149931 1200:0.000119707 2400:9.85791e-05 , 4800: 8.18873e-05 # Based on a simulation with no taper
    noisin = [0.000701111,0.00065377 ,0.000499886,0.000436001,0.000213708,0.000178363,0.000149931,0.000119707 ,9.85791e-05 , 8.18873e-05]
    timein = [12.,24,48,96,300,600,900,1200,2400,4800]
    a1,a2 =np.polyfit(noisin,timein,1)
    totaltime = a1*noise+a2
    if totaltime < 12:
        totaltime = 12.
    RA,DEC =cf.convertRADEC(dummy[0].header['CRVAL1'],dummy[0].header['CRVAL2'])
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

    tri = open(work_dir+'/run_casa.py', 'w')
    tri.writelines([Template_Casa[key]+"\n" for key in Template_Casa])
    tri.close()

    os.chdir(work_dir)
    bla = subprocess.call(['casa','--nologger','--nogui','-c','run_casa.py'])
    dummy = fits.open(work_dir+'/Convolved_Cube.fits',uint = False, do_not_scale_image_data=True,ignore_blank = True)
    # We cut 20 pixels around the edges
    newsize = int(np.shape(dummy[0].data)[2]-60.)
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
    newdummy=np.zeros((np.shape(dummy[0].data)[0],newsize,newsize))
    newdummy[:,:,:]=dummy[0].data[:,int(np.floor(dummy[0].header['CRPIX2']-newsize/2.)):int(np.floor(dummy[0].header['CRPIX2']+newsize/2.)),int(np.floor(dummy[0].header['CRPIX1']-newsize/2.)):int(np.floor(dummy[0].header['CRPIX1']+newsize/2.))]
    dummy[0].header['NAXIS1']=newsize+1
    dummy[0].header['NAXIS2']=newsize+1

    dummy[0].header['CRPIX1']=dummy[0].header['CRPIX1']-np.floor(dummy[0].header['CRPIX1']-newsize/2.)
    dummy[0].header['CRPIX2']=dummy[0].header['CRPIX2']-np.floor(dummy[0].header['CRPIX2']-newsize/2.)
    fits.writeto(work_dir+'/mask.fits',newdummy,dummy[0].header, overwrite = True)



#------------------------------!!!!!!!!!!!!!!!!!!!!!!!!!!!!!Start of the main program!!!!!!!!!!!!!!!!!!!!!----------------------
def AGC(work_dir='',running_default = 'Not_Set'):

    Mass = []
    #First ask for the directory to work in
    if work_dir == '':
        work_dir = input("Please provide the directory where to create the database :")  #
    while not os.path.isdir(work_dir):
        print("That is not a valid directory please try again :")
        work_dir = input("Please provide the directory where to create the database :")
    # Then do we want to set all values indivdually or not
    if work_dir[-1] != '/':
        work_dir = work_dir+'/'
    if running_default == 'Not_Set':
        run_default = cf.get_bool("Do you want want to create the default Artificial data base? (Yes/No, default = Yes ) : ")
    else:
        run_default = bool(running_default)

    if not run_default:
       #sets = 5  # This is the amount of base galaxies we want, i.e. the number of rotation curves
        sets = int(input("How many base galaxies do you want to use? (1-5, default= 5): ") or 5)
        if sets > 5 or sets < 1:
            print("{} is not a valid number of sets".format(sets))
            exit()
        # do we wanr inhomogeneities
        #inhomogeneity = True
        inhomogeneity = cf.get_bool("Do you want the galaxies to have inhomogeneities? (True/False, default = True) : ")
        makenewmodels = cf.get_bool("Do you want to erase all existing models in the working directory? (Yes/No, default=No): ",default=False)

        # How do we want to corrupt the model, For now the options are Gaussian or Casa_Sim
        # corruption = 'Casa_Sim'
        corruption = 'Gaussian'
        #corruption = 'Casa_5' # This uses casa for every 5th model
        corruption = input("How do you want to corrupt the model (Casa_Sim, Gaussian, Casa_5, default = Gaussian) :")
        if corruption.lower() == 'casa_sim' or corruption.lower() == "c_s":
            corruption = 'Casa_Sim'
            print("You are using the casa corruption method please make sure python can access casa.")
        elif  corruption.lower() == 'gaussian' or corruption == "" or corruption.lower() == 'g' :
            corruption = 'Gaussian'
        elif corruption.lower() == 'casa_5' or corruption.lower() == "c_5":
            corruption = 'Casa_5'
            print("You are using the casa corruption method please make sure python can access casa.")
        else:
            print("We can not decipher your requested method. We will use the default.")
            corruption = 'Gaussian'
        symmetric = cf.get_bool("Do you want symmetric galaxies? (Yes/No, default = No):", default=False)
        changes_poss = ['Inclination', 'PA','Beams','Radial_Motions','Flare','Arms','Bar','Channelwidth','SNR','Warp','Mass','Res_Beam']
        changes = ['Base']
        for opts in changes_poss:
            inc_current_opt = cf.get_bool("Do you want to include a variation in {} (Yes/No, default = No): ".format(opts),default=False)
            if inc_current_opt:
                changes.append(opts)
                if opts == 'Inclination':
                    vals = input("please provide the input parameters to vary {} over. : ".format(opts))
                    Inclination = [float(x) for x in re.split("\s+|\s*,\s*|\s+$",vals.strip())]
                elif opts == 'PA':
                    vals = input("please provide the input parameters to vary {} over. : ".format(opts))
                    PA = [float(x) for x in re.split("\s+|\s*,\s*|\s+$",vals.strip())]
                elif opts == 'Flare':
                    Flare = []
                    do_flare = cf.get_bool("Do you want to add a Flare when not present in the Base? (Yes/No, default = yes): ")
                    if do_flare: Flare.append('Flared')
                    do_not_flare = cf.get_bool("Do you want to remove the Flare when  present in the Base? (Yes/No, default = yes): ")
                    if do_not_flare: Flare.append('No_Flare')
                elif opts == 'Warp':
                    Warp = []
                    vals = input("Please provide the variation in theta and phi: ")
                    initial = re.split("\s+|\s*,\s*|\s+$",vals.strip())
                    while len(initial) != 2:
                        vals = input("Please provide two and only two values: ")
                        initial = re.split("\s+|\s*,\s*|\s+$",vals.strip())
                    Warp.append([float(initial[0]),float(initial[1])])
                    another = cf.get_bool("Do you want to add another set of variation? (Yes/No, default=No): ",default=False)
                    while another:
                        vals = input("Please provide the variation in theta and phi : ")
                        initial = re.split("\s+|\s*,\s*|\s+$",vals.strip())
                        while len(initial) != 2:
                             vals = input("Please provide two and only two values: ")
                             initial = re.split("\s+|\s*,\s*|\s+$",vals.strip())
                        Warp.append([float(initial[0]),float(initial[1])])
                        another = cf.get_bool("Do you want to add another set of variation? (Yes/No, default=No): ",default=False)
                elif opts == 'Beams':
                    vals = input("please provide the input parameters to vary {} over. : ".format(opts))
                    Beams =  [float(x) for x in re.split("\s+|\s*,\s*|\s+$",vals.strip())]
                elif opts == 'SNR':
                    vals = input("please provide the input parameters to vary {} over. : ".format(opts))
                    SNR =  [float(x) for x in re.split("\s+|\s*,\s*|\s+$",vals.strip())]
                elif opts == 'Channelwidth':
                    vals = input("please provide the input parameters to vary {} over. : ".format(opts))
                    Channelwidth = [float(x) for x in re.split("\s+|\s*,\s*|\s+$",vals.strip())]
                elif opts == 'Res_Beam':
                    Res_Beam = []
                    vals = input("Please provide the major and minor beam axis: ")
                    initial = re.split("\s+|\s*,\s*|\s+$",vals.strip())
                    while len(initial) != 2:
                             vals = input("Please provide two and only two values: ")
                             initial = re.split("\s+|\s*,\s*|\s+$",vals.strip())
                    Res_Beam.append([float(initial[0]),float(initial[1])])
                    another = cf.get_bool("Do you want to add another set of variation? (Yes/No, default=No): ",default=False)
                    while another:
                        vals = input("Please provide tthe major and minor beam axis: ")
                        initial = re.split("\s+|\s*,\s*|\s+$",vals.strip())
                        while len(initial) != 2:
                             vals = input("Please provide two and only two values: ")
                             initial = re.split("\s+|\s*,\s*|\s+$",vals.strip())
                        Res_Beam.append([float(initial[0]),float(initial[1])])
                        another = cf.get_bool("Do you want to add another set of variation? (Yes/No, default=No): ",default=False)
                elif opts == 'Arms':
                    Arms = []
                    do_flare = cf.get_bool("Do you want to add Arms when not present in the Base? (Yes/No, default = yes): ")
                    if do_flare: Arms.append('Arms')
                    do_not_flare = cf.get_bool("Do you want to remove the Arms when  present in the Base? (Yes/No, default = yes): ")
                    if do_not_flare: Arms.append('No_Arms')
                elif opts == 'Bar':
                    Bar = []
                    do_flare = cf.get_bool("Do you want to add a Bar when not present in the Base? (Yes/No, default = yes):")
                    if do_flare: Bar.append('Bar')
                    do_not_flare = cf.get_bool("Do you want to remove the Bar when  present in the Base? (Yes/No, default = yes):")
                    if do_not_flare: Bar.append('No_Bar')
                elif opts == 'Radial_Motions':
                    vals = input("please provide the input parameters to vary {} over: ".format(opts))
                    Radial_Motions = [float(x) for x in re.split("\s+|\s*,\s*|\s+$",vals.strip())]
                elif opts == 'Mass':
                    vals = input("please provide the input parameters to vary {} over: ".format(opts))
                    Mass =  [float(x) for x in re.split("\s+|\s*,\s*|\s+$",vals.strip())]
                else:
                    print("This is not a supported parameter")
    else:


        #work_dir='/home/peter/Database/Test/'
        while work_dir == '':
            print("There is no default")
            work_dir = input("Please provide the directory where to create the database :")
        sets = 5
        inhomogeneity = True
        corruption = 'Gaussian'
        changes = ['Base','Inclination','Beams','Radial_Motions','Flare','Arms','Bar','Mass','Channelwidth','SNR','Warp','Mass','Res_Beam']
        Inclination= [15,20,30,50,70,80,88,90]
        Warp=[[0.15,0.05],[0.05,0.2]]
        Flare=["Flared","No_Flare"] # Flared, No_Flare
        Beams=[2,4,6,7,8,10,12] # Beam across the major axis. This also set the distance as the size in kpc will be determined by Wang 2016 from the SBR profile
        SNR=[1,3,5] # These  are average signal to noise ratios
        Channelwidth=[2.,8.]
        Res_Beam=[[5.,5.]] # Resolution of the beam in arcsec
        Arms = ["Arms","No_Arms"]   #Either "No_Arms or Arms"
        Bar = ["Bar","No_Bar"]   #Either "No_Bar or Bar"
        Radial_Motions = [-5.,-10.] # in km/s
        Mass = [2.5e11]
        #Parameter to force a new set of models being made in this directory
        # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        # please note that if set to true at the start of the script everything in work_dir wil be deleted
        # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        makenewmodels = True
        # do we want symmetric models or not?
        symmetric = False
    # Let's give an over view of the database that will be created
    print("We will create a database with {} basic sets in the directory {}.\n".format(sets,work_dir))
    if makenewmodels:
        print("All previous models will be removed prior to the build. \n")
    else:
        print("We will retain previously build models. \n")
    if inhomogeneity:
        print("We will use {} as corruption method and inhomogeneities will be added.\n".format(corruption))
    else:
        print("We will use {} as corruption method.\n".format(corruption))
    print("We will vary the following parameters")
    if 'Inclination' in changes:
        print("We vary the inclination with the following values: {}.\n".format(" ".join([str(e) for e in Inclination])))
    if 'PA' in changes:
        print("We vary the PA with the following values: {}.\n".format(" ".join([str(e) for e in PA])))
    if 'Beams' in changes:
        print("We create the model with {} beams across the major axis.\n".format(", ".join([str(e) for e in Beams])))
    if 'Radial_Motions' in changes:
        print("Inject radial motions with speeds of {} km/s.\n".format(" ".join([str(e) for e in Radial_Motions])))
    if 'Flare' in changes:
        print("Varying the scale height with: {}.\n".format(" ".join([str(e) for e in Flare])))
    if 'Arms' in changes:
        print("Varying the arms with: {}.\n".format(" ".join([str(e) for e in Arms])))
    if 'Bar' in changes:
        print("Varying the bar with: {}.\n".format(" ".join([str(e) for e in Bar])))
    if 'Mass' in changes:
        print("We add the following masses to each base set: {}.\n".format(" ".join(["{:10.2e}".format(float(e)) for e in Mass])))
    if 'Channelwidth' in changes:
        print("Varying the channel width with: {} km/s.\n".format(" ".join([str(e) for e in Channelwidth])))
    if 'SNR' in changes:
        print("Varying the signal to noise ratio with: {}.\n".format(" ".join([str(e) for e in SNR])))
    if 'Warp' in changes:
        print("Varying the theta angle of the angular momentum vector with: {}.\n".format(" ".join([str(e) for e in Warp[0][:]])))
        print("Varying the phi angle of the angular momentum vector with: {}.\n".format(" ".join([str(e) for e in Warp[1][:]])))
    if 'Res_Beam' in changes:
        print("Varying the beam size with: {}.\n".format(" ".join([str(e) for e in Res_Beam])))


    # If we make new models delete everything in the directory
    if makenewmodels:
        os.system('rm -R '+work_dir+'/Mass*rm*')


    # Let's just make 1 catalogue with everything we need and adjust
    # the fitting program to have this one catalogue as input

    Catalogue=work_dir+'/Output_AGC_Summary.txt'
    # If we are making new models we want to ensure this is a new file
    if makenewmodels:
       cat = open(Catalogue, 'w')
       cat.write('number|Distance|Directoryname|Cubename\n')
       cat.close()



    #Copy a fits file from the WHISP data base to use as template if it is not ther3 yet
    # Check for the existence of a template fits file
    templatethere= os.path.isfile(work_dir+'/Input.fits')
    #templatethere =False
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
    Template_in = cf.Proper_Dictionary({})
    unarranged = tmp.readlines()
    # Separate the keyword names
    for tmp in unarranged:
        # python is really annoying with needing endlines. Let's strip them here and add them when writing
        Template_in[tmp.split('=',1)[0].strip().upper()]=tmp.rstrip()

    #If we want to corrupt in the casa way we'd need to read the file corruption file
    if (corruption == 'Casa_Sim') or (corruption == 'Casa_5') :
        tmp = open('Template_Casa.py','r')
        Template_Casa_In = cf.Proper_Dictionary({})
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
    global H_0
    H_0 = 69.6 # http://www.astro.ucla.edu/~wright/CosmoCalc.html
    # start a loop over the various base galaxies
    number_models = 0.
    set_done= [1024]


    colors=iter(plt.cm.rainbow(np.linspace(0,1,sets-1+len(Mass))))
    for base in range(sets):
        # From here we go into a loop to adjust variables over the bases
        for ix in range(len(changes)):
            if changes[ix] == 'Inclination':numloops=len(Inclination)
            elif changes[ix] == 'PA': numloops=len(PA)
            elif changes[ix] == 'Flare': numloops=len(Flare)
            elif changes[ix] == 'Warp': numloops=len(Warp)
            elif changes[ix] == 'Beams': numloops=len(Beams)
            elif changes[ix] == 'SNR': numloops=len(SNR)
            elif changes[ix] == 'Channelwidth': numloops=len(Channelwidth)
            elif changes[ix] == 'Res_Beam': numloops=len(Res_Beam)
            elif changes[ix] == 'Arms': numloops=len(Arms)
            elif changes[ix] == 'Bar': numloops=len(Bar)
            elif changes[ix] == 'Radial_Motions': numloops=len(Radial_Motions)
            elif changes[ix] == 'Mass': numloops=len(Mass)
            elif changes[ix] == 'Base': numloops=1
            else:
                print("This is not a supported parameter")
                exit()
            for jx in range (numloops):
                Current_Galaxy = Base_Galaxy(base)

                if changes[ix] == 'Inclination':Current_Galaxy.Inclination = Inclination[jx]
                elif changes[ix] == 'PA': Current_Galaxy.PA = PA[jx]
                elif changes[ix] == 'Flare': Current_Galaxy.Flare = Flare[jx]
                elif changes[ix] == 'Warp': Current_Galaxy.Warp = [Warp[jx][0],Warp[jx][1]]
                elif changes[ix] == 'Beams':Current_Galaxy.Beams = Beams[jx]
                elif changes[ix] == 'SNR': Current_Galaxy.SNR = SNR[jx]
                elif changes[ix] == 'Channelwidth': Current_Galaxy.Channelwidth = Channelwidth[jx]
                elif changes[ix] == 'Res_Beam': Current_Galaxy.Res_Beam = [Res_Beam[jx][0],Res_Beam[jx][1]]
                elif changes[ix] == 'Arms': Current_Galaxy.Arms = Arms[jx]
                elif changes[ix] == 'Bar': Current_Galaxy.Bar = Bar[jx]
                elif changes[ix] == 'Radial_Motions': Current_Galaxy.Radial_Motions = Radial_Motions[jx]
                elif changes[ix] == 'Mass': Current_Galaxy.Mass = Mass[jx]
                elif changes[ix] == 'Base': Current_Galaxy = Base_Galaxy(base)
                else:print("This is not a supported parameter")
                Current_Galaxy.Res_Beam[0:1] = np.sort(Current_Galaxy.Res_Beam[0:1])
                if len(Current_Galaxy.Res_Beam) == 2: Current_Galaxy.Res_Beam.append(0)
                print("This is the parameter {}. And this is the input {}".format(changes[ix],Current_Galaxy.Inclination))
                # Build a name and a directory where to stor the specific output

                print("This is the Current Mass = {:.1e}".format(Current_Galaxy.Mass))
                name="Mass{:.1e}-i{}d{}-{}pa{}w{}-{}-{}-ba{}SNR{}bm{}-{}ch{}-{}-{}-rm{}".format(Current_Galaxy.Mass,Current_Galaxy.Inclination,Current_Galaxy.Dispersion[0],Current_Galaxy.Dispersion[1],Current_Galaxy.PA,Current_Galaxy.Warp[0],Current_Galaxy.Warp[1],Current_Galaxy.Flare,Current_Galaxy.Beams,Current_Galaxy.SNR,Current_Galaxy.Res_Beam[0],Current_Galaxy.Res_Beam[1],Current_Galaxy.Channelwidth,Current_Galaxy.Arms,Current_Galaxy.Bar,Current_Galaxy.Radial_Motions)
                print("{} is the name of the current galaxy".format(name))

                # Make a dirctory
                # Check for the existence of the directory
                constructstring="mkdir "+work_dir+'/'+name
                checkdir=False
                galaxy_dir= os.path.isdir(work_dir+'/'+name)
                if not galaxy_dir:
                    os.system(constructstring)
                else:
                    # Do we have a cube
                    galaxy_cube_exist = os.path.isfile(work_dir+name+"/Convolved_Cube.fits")
                    if galaxy_cube_exist:
                        print("This galaxy appears fully produced")
                        checkdir = True
                        continue
                    else:
                        print("The directory was made but there is no full cube avalaible")
                        print("Reproducing the galaxy. Be aware of Double Table entries")
                        print("This is too dangerous. Breaking the code.")
                        exit()
                if checkdir:
                    print('This should be impossible')
                    exit()
                # Then we copy the original Template def
                number_models += 1
                global Template
                Template=copy.deepcopy(Template_in)
                # We also open a figure to plot all info on
                plt.figure(2, figsize=(8, 12), dpi=100, facecolor='w', edgecolor='k')
                global overview
                overview = plt.subplot(6,1,6)

                 # First set the beam to 0.
                Template["BMAJ"]= "BMAJ = 0."
                Template["BMIN"]= "BMIN = 0."
                Template["BPA"]= "BPA = 0."
                labelfont= {'family':'Times New Roman',
                            'weight':'normal',
                            'size':22}
                plt.rc('font',**labelfont)
                #Then we need to build the Surface Brightnes profile
                SBRprof,Rad,sclength,MHI,Rad_HI,Vrot,sub_ring = build_sbr_prof(Current_Galaxy) #Column densities,Raii in kpc, Opt_scalelength in kpc, HI mass in M_solar
                            #We want to make a figure with all the Rotation curves
                if Current_Galaxy.Mass not in set_done:
                    if set_done[0] == 1024:
                        set_done= [Current_Galaxy.Mass]
                        labelfont= {'family':'Times New Roman',
                                 'weight':'normal',
                                 'size':22}
                        plt.rc('font',**labelfont)
                        plt.figure(59, figsize=(8, 8), dpi=300, facecolor='w', edgecolor='k')
                        ax = plt.subplot(1, 1, 1)
                        plt.plot(Rad,Vrot,'k')
                        #plt.plot(Rad,Vrot,'ko',label='M$_{\odot}$ = {:.1e}'.format(Current_Galaxy.Mass))
                        plt.plot(Rad,Vrot,'ko',label='M$_{\odot}$ = '+' {:.1e}'.format(Current_Galaxy.Mass))
                        plt.ylabel('V$_{rot}$ (km s$^{-1}$)',**labelfont)
                        plt.xlabel('Radius (kpc)',**labelfont)
                        ax.yaxis.set_minor_locator(AutoMinorLocator(4))
                        ax.xaxis.set_minor_locator(AutoMinorLocator(4))
                        for axis in ['top','bottom','left','right']:
                            ax.spines[axis].set_linewidth(1.5)
                        plt.tick_params(axis='both', which='minor', bottom=True,left=True,length=3)
                        plt.tick_params(axis='both', which='major', labelsize=17, length=6)
                        plt.tick_params(axis='both', which='both', direction = 'in', width=1.5 , bottom=True,left=True ,right =True, top=True)
                        max_rad = np.max(Rad)+1
                        names_used=[Current_Galaxy.Mass]
                        colors_used=[0.]
                    else:
                        set_done.append(Current_Galaxy.Mass)
                        plt.figure(59)
                        c=next(colors)
                        if np.max(Rad)+1 > max_rad:
                            max_rad =  np.max(Rad)+1
                        plt.plot(Rad,Vrot,c=c)
                        #plt.plot(Rad,Vrot,'o',label=r'M$_{\odot}$ = {:.1e}'.format(Current_Galaxy.Mass),c=c)
                        plt.plot(Rad,Vrot,'o',label='M$_{\odot}$ = '+' {:.1e}'.format(Current_Galaxy.Mass),c=c)
                        names_used.append(Current_Galaxy.Mass)
                        colors_used.append(c)


                #We need central coordinates the vsys will come from the required distance and hubble flow. The RA and dec should not matter hance it will be the only random component in the code as we do want to test variations of them
                Sky_Size = np.radians(Current_Galaxy.Res_Beam[0]*Current_Galaxy.Beams/3600.)
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
                RAhr,DEChr= cf.convertRADEC(RAdeg,DECdeg)
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
                Rad_arcsec = cf.convertskyangle(Rad,distance=Distance,physical = True)
                 # then the number of rings to the total number of rings
                Template["NUR"]="NUR= {}".format(len(Rad))
                # The radii in arcsec
                Template["RADI"]="RADI = "+" ".join(str(e) for e in Rad_arcsec)
                # Then we need to get a starting radius for the warp.
                # the warp should start at the edge of the optical radius which is the HI scale length/0.6
                # which are about ~ 4 * h_r
                WarpStart = 4.*sclength
                WarpEnd=Rad[np.where(SBRprof >= 4.98534620064e-05)[0][-1]]
                # Write it to the Template
                Template["VROT"]="VROT = "+" ".join(str(e) for e in Vrot)
                Template["VROT_2"]="VROT_2 = "+" ".join(str(e) for e in Vrot)
                #and the Overview Figure
                plt.figure(2)
                plt.subplot(6,1,5)
                plt.plot(Rad,Vrot,'k')
                plt.plot(Rad[0::sub_ring],Vrot[0::sub_ring],'ko')
                ymin,ymax = plt.ylim()
                plt.margins(x=0., y=0.)
                labelfont= {'family':'Times New Roman',
                            'weight':'normal',
                            'size':18}
                plt.rc('font',**labelfont)
                plt.tick_params(
                    axis='x',          # changes apply to the x-axis
                    which='both',      # both major and minor ticks are affected
                    direction = 'in',
                    bottom=True,      # ticks along the bottom edge are off
                    top=False,         # ticks along the top edge are off
                    labelbottom=False) # labels along the bottom edge are off
                plt.plot([Rad_HI,Rad_HI],[ymin-(ymax-ymin)*0.1,ymax+(ymax-ymin)*0.1],c='b')
                plt.ylabel('V$_{rot}$ (km s$^{-1}$)',labelpad=20,**labelfont)
                #plt.xticks([])
                # We need a scale height and dispersion for each ring. They are coupled and hence they are both created in create_flare
                h_z,dispersion=create_flare(Rad,Vrot,Current_Galaxy.Dispersion,Current_Galaxy.Flare,Rad_HI,sub_ring,distance=Distance)

                # Finally we need to set the warping
                PA,inc,phirings = create_warp(Rad,Current_Galaxy.PA,Current_Galaxy.Inclination,Current_Galaxy.Warp,[WarpStart,WarpEnd],sub_ring)
                if symmetric:
                    PA_2,inc_2,phirings_2 = create_warp(Rad,Current_Galaxy.PA,Current_Galaxy.Inclination,Current_Galaxy.Warp,[WarpStart,WarpEnd],sub_ring,disk=2)
                else:
                    # we want an assymetric warp so we redo the PA and inc but swap the variation
                    PA_2,inc_2,phirings_2 = create_warp(Rad,Current_Galaxy.PA,Current_Galaxy.Inclination,[Current_Galaxy.Warp[0]-Current_Galaxy.Warp[1],Current_Galaxy.Warp[1]/2.+Current_Galaxy.Warp[0]/2.],[WarpStart,WarpEnd],sub_ring,disk=2)

                #If we want Radial Motions then they need to be inserted
                if Current_Galaxy.Radial_Motions != 0.:
                      Template.insert("VROT","VRAD","VRAD = {}".format(Current_Galaxy.Radial_Motions))
                      Template.insert("VROT_2","VRAD_2","VRAD_2 = {}".format(Current_Galaxy.Radial_Motions))

                # This comes from FAT. If I remember correctly this is the sine response to the channels *1.2/(2*SQRT(2*ALOG(2.))))
                # However in our input we want independent channels which means we should set this to 0.
                # Template["CONDISP"]="CONDISP = 0."
                Template["CONDISP"]="CONDISP = {}".format(Current_Galaxy.Channelwidth*1.2/(2*np.sqrt(2*np.log(2.))))
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
                plt.subplot(6,1,1)
                plt.title("DM Mass in = {:.2e}".format(Current_Galaxy.Mass))
                #-------------------------------This finishes the basic disk the following are optional components----------------------


                # The possible arms
                if Current_Galaxy.Arms == 'Arms':
                    phase,arm_brightness,arm_width = create_arms(Vrot,Rad,SBRprof,WarpStart=WarpStart, Bar=Current_Galaxy.Bar)
                    phase,arm_brightness,arm_width = create_arms(Vrot,Rad,SBRprof,disk=2,WarpStart=WarpStart, Bar=Current_Galaxy.Bar)
                # A possible Bar
                if Current_Galaxy.Bar == 'Barred':
                    bar_length = create_bar(Vrot,Rad,SBRprof,Template,WarpStart=WarpStart)
                # and possible inhomogeneities
                if inhomogeneity:
                    inhomogeneity_amp = create_inhomogeneity(MHI,Current_Galaxy.SNR,disks=[1,2])
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
                if  Current_Galaxy.Beams < 6.:
                    size = 2.*Rad_arcsec[-1] +(6.*Current_Galaxy.Res_Beam[0])
                pix_size = (Rad_arcsec[1]-Rad_arcsec[0])
                if (corruption == 'Casa_5' and (int(number_models/5.) == number_models/5.)) or (corruption == 'Casa_Sim'):
                    size += 60*pix_size
                required_pixels=int(np.ceil(size/pix_size))
                vel_max=2.5*np.max([Vrot*np.sin(inc*np.pi/180.)+4*dispersion,Vrot*np.sin(inc_2*np.pi/180.)+4*dispersion])
                velpix=int(np.ceil(vel_max/Current_Galaxy.Channelwidth)*3 )
                dummy[0].header['CRPIX1'] = np.floor(required_pixels/2.)
                dummy[0].header['CRPIX2'] = np.floor(required_pixels/2.)
                dummy[0].header['CRPIX3'] = np.floor(velpix/2.)
                # Stupid astropy doesn't account for minus in header of cdelt1 then cdelt has different precision
                tmp=int(pix_size/3600.*1e15)/1e15
                dummy[0].header['CDELT1'] = -1*tmp
                dummy[0].header['CDELT2'] = tmp
                dummy[0].header['CDELT3'] = Current_Galaxy.Channelwidth*1000./3.
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
                plt.savefig('Overview_Input.png', bbox_inches='tight')

                plt.close()
                os.chdir(work_dir)
                Template.clear()

                # Now we want to corrupt this cube with some realistic noise
                # For this we first want to get the noise we want in terms of Jansky per beam
                # we will define the SNR as the mean(Intensity)/noiselevel hence noise =mean(In)/SNR

                if (corruption == 'Casa_5' and (int(number_models/5.) == number_models/5.)) or (corruption == 'Casa_Sim'):
                    Template_Casa=copy.deepcopy(Template_Casa_In)
                    corrupt_casa(work_dir+name+'/',Current_Galaxy.Res_Beam,Template_Casa,SNR)
                    os.chdir(work_dir)
                    mask = fits.open(work_dir+name+'/mask.fits',uint = False, do_not_scale_image_data=True,ignore_blank = True)
                    Cube = fits.open(work_dir+name+'/Convolved_Cube.fits',uint = False, do_not_scale_image_data=True,ignore_blank = True)
                    maskr = mask[0].data[1:]
                    sigma = (np.std(Cube[0].data[0])+np.std(Cube[0].data[-1]))/2.
                    Cube_Clean = Cube[0].data
                    Cube_Clean[maskr < 0.5] = 0.
                    beamarea=(np.pi*abs(Cube[0].header['BMAJ']*3600.*Cube[0].header['BMIN']*3600.))/(4.*np.log(2.))
                    pixperbeam=beamarea/(abs(Cube[0].header['CDELT1']*3600.)*abs(Cube[0].header['CDELT2']*3600.))
                    totalsignal = np.sum(Cube_Clean)/pixperbeam
                    mass = 2.36E5*Distance**2*totalsignal*Cube[0].header['CDELT3']/1000.
                    totsig=np.zeros(len(Cube_Clean[:]))
                    for j in range(len(totsig)):
                        if  len(Cube_Clean[j][Cube_Clean[j] > 0.]) > 0:
                            totsig[j]=np.mean(Cube_Clean[j][Cube_Clean[j] > 0.])
                    mean_signal = np.median(totsig[totsig > 0.])
                    SNRachieved = mean_signal/(sigma)

                elif (corruption == 'Gaussian' or corruption == 'Casa_5'):
                    corrupt_gauss(work_dir+name+'/',Current_Galaxy.Res_Beam,Current_Galaxy.SNR)
                    mask = fits.open(work_dir+name+'/mask.fits',uint = False, do_not_scale_image_data=True,ignore_blank = True)
                    Cube = fits.open(work_dir+name+'/Convolved_Cube.fits',uint = False, do_not_scale_image_data=True,ignore_blank = True)
                    maskr = mask[0].data[:]
                    sigma = (np.std(Cube[0].data[0])+np.std(Cube[0].data[-1]))/2.
                    Cube_Clean = Cube[0].data
                    Cube_Clean[maskr < 0.5] = 0.
                    beamarea=(np.pi*abs(Cube[0].header['BMAJ']*3600.*Cube[0].header['BMIN']*3600.))/(4.*np.log(2.))


                    pixperbeam=beamarea/(abs(Cube[0].header['CDELT1']*3600.)*abs(Cube[0].header['CDELT2']*3600.))
                    totalsignal = np.sum(Cube_Clean)/pixperbeam
                    mass = 2.36E5*Distance**2*totalsignal*Cube[0].header['CDELT3']/1000.
                    mean_signal = np.mean(Cube_Clean[maskr > 0.5])
                    SNRachieved = mean_signal/(sigma)
                else:
                    print("!!!!!!!This corruption method is unknown, leaving the cube uncorrupted and unconvolved!!!!!!!!")
                    # We'll create a little text file with an Overview of all the parameters
                overview = open(work_dir+name+'/'+name+'-Info.txt', 'w')
                overview.write("This file contains the basic parameters of this galaxy\n")
                overview.write("For the radial dependencies look at Overview.png or ModelInput.def\n")
                overview.write("Inclination = {}\n".format(Current_Galaxy.Inclination))
                overview.write("The dispersion = {:.2f}-{:.2f}\n".format(dispersion[0],dispersion[-1]))
                overview.write("The type of galaxy = {:1e}\n".format(Current_Galaxy.Mass))
                overview.write("PA = {}\n".format(Current_Galaxy.PA))
                overview.write("Warp = {}-{}\n".format(Current_Galaxy.Warp[0],Current_Galaxy.Warp[1]))
                overview.write("Which starts at {:.2f} kpc and the 1M/pc^2 radius is {:.2f} kpc \n".format(WarpStart,Rad_HI))
                overview.write("Flare = {}\n".format(Current_Galaxy.Flare))
                overview.write("Beams across the major axis = {}\n".format(Current_Galaxy.Beams))
                overview.write("SNR Requested = {} SNR Achieved = {}  \n".format(Current_Galaxy.SNR,SNRachieved))
                overview.write("Mean Signal = {}  \n".format(mean_signal))
                overview.write("Channelwidth = {}\n".format(Current_Galaxy.Channelwidth))
                overview.write("Major axis beam = {} Minor axis beam= {}\n".format(Current_Galaxy.Res_Beam[0],Current_Galaxy.Res_Beam[1]))
                overview.write("This galaxy has {} and a {}\n".format(Current_Galaxy.Arms,Current_Galaxy.Bar))
                overview.write("It's central coordinates are RA={} DEC={} vsys={:.2f} km/s\n".format(RAhr,DEChr,vsys))
                overview.write("At a Distance of {:.2f} Mpc \n".format(Distance))
                overview.write("HI_Mass Requested {:.2e} (M_solar) and an optical h {:.2f} (kpc)\n".format(MHI,sclength))
                overview.write("HI_Mass Retrieved {:.2e} (M_solar) \n".format(mass))
                overview.write("We have {} pix per beam \n".format(pixperbeam))
                overview.write("The cube was corrupted with the {} method \n".format(corruption))
                overview.write("The final noise level is {} Jy/beam \n".format(sigma))
                overview.write("h_z = {:.3f}-{:.3f} (kpc)".format(h_z[0],h_z[-1]))
                overview.close()
                cat = open(Catalogue, 'a')
                cat.write('{:d}|{:.2f}|{}|Convolved_Cube\n'.format(int(number_models),Distance,name))
                cat.close()
                # We also want a file that contains initial estimates for all the parameters. We scramble them with gaussian variations
                overview = open(work_dir+name+'/Initial_Estimates.txt', 'w')
                overview.write("#This file contains the initial estimates \n")
                scale=h_z[0]+np.random.randn(1)[0]*0.1+0.1
                if scale < 0.05: scale = 0.05
                overview.write("#{:<15} {:<15} {:<15} {:<15} {:<15} {:<15} {:<15} {:<15} {:<15}  {:<15}\n ".format('VROT','INCL','PA','Z0','SBR','DISP','VRAD','RA','DEC','VSYS'))
                overview.write("#{:<15} {:<15} {:<15} {:<15} {:<15} {:<15} {:<15} {:<15} {:<15}  {:<15}\n ".format('(km/s)','Degree','Degree','kpc','Jy km/s/arcsec^2','km/s','km/s','Degree','Degree','km/s'))
                overview.write("{:<16.2f} {:<15.2f} {:<15.2f} {:<15.3f} {:<15.7f} {:<15.2f} {:<15.2f} {:<15.5f} {:<15.5f}  {:<15.2f}\n ".format(np.mean(Vrot[-5:-1])+np.random.randn(1)[0]*10.,
                                                                                                                  Current_Galaxy.Inclination+np.random.randn(1)[0]*10.,
                                                                                                                  Current_Galaxy.PA+np.random.randn(1)[0]*3.,
                                                                                                                  scale,
                                                                                                                  np.mean(SBRprof)+np.random.randn(1)[0]*np.max(SBRprof)/10,
                                                                                                                  np.mean(dispersion)+np.random.randn(1)[0]*4.,
                                                                                                                  0.,
                                                                                                                  RAdeg+np.random.randn(1)[0]*10./3600.,
                                                                                                                  DECdeg+np.random.randn(1)[0]*10./3600.,
                                                                                                                  float(vsys)+np.random.randn(1)[0]*4.))

                overview.close()

                # and cleanup
                os.system("rm -f "+work_dir+name+'/ModelInput.def')
                os.system("mv "+work_dir+name+'/tirific.def '+work_dir+name+'/ModelInput.def')

    print("We created {} models".format(number_models))
    plt.figure(59)
    ax.set_ylim(ymin=0)
    ax.set_xlim(xmin=0, xmax=max_rad)
    plt.legend(loc='lower right',fontsize=12)
    plt.savefig('Rotation_Curves.ps', bbox_inches='tight')
    plt.close()
if __name__ == '__main__':
    AGC()
