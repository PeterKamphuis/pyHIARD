from dataclasses import dataclass, field
from multiprocessing import cpu_count
import omegaconf
from omegaconf import MISSING
from typing import List, Optional
from datetime import datetime
import os

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
    # The channel dependency
    # 'Options are independent, sinusoidal, hanning
    channel_dependency: str = 'sinusoidal'
    retain_unconvolved_model: bool = False


@dataclass
class ROC:
    enable: bool = True
    add_template: bool = False
    remove_template: bool = False
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
    #A lower number  is for speed and low memory usage, a high number provides beter integrateion between the noise and the template. If lower than minimum_degradation_factor it will be set to the minimum_degradation_factor
    snr: List = field(default_factory=lambda: [0.5, 1., 3.])


@dataclass
class General:
    ncpu: int = cpu_count()-1
    main_directory: str = os.getcwd()
    tirific: str = "tirific"  # Command to call tirific
    sofia2: str = "sofia2"  # Command to call sofia 2
    multiprocessing: bool = True


@dataclass
class Config:
    print_examples: bool = False
    print_bases: bool = False
    configuration_file: Optional[str] = None
    general: General = General()
    agc: AGC = AGC()
    roc: ROC = ROC()
