from dataclasses import dataclass, field
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
    base_galaxies: List[int] = field(
        default_factory=lambda: [1, 2, 3, 4, 5, 6])
    inhomogenous: bool = True  # Add homogenieties?
    symmetric: bool = False  # Keep galaxies symmetric
    corruption_method: str = 'Gaussian'  # options are Casa_Sim, Gaussian, Casa_5
    variables_to_vary: List[str] = field(default_factory=lambda: ['Inclination', 'Beams', 'Radial_Motions',
                                         'Flare', 'Arms', 'Bar', 'Mass', 'Channelwidth', 'SNR', 'Warp', 'Mass', 'Beam_Resolution'])
    # Each base is created with the variations in the following parameters if they are listed to be varied.
    inclination: List[float] = field(default_factory=lambda: [
                                     15., 20., 30., 50., 70., 80., 88., 90.])
    pa: List[float] = field(default_factory=lambda: [0., 360.])
    warp: List[float] = field(default_factory=lambda: [
                              [0.15, 0.05], [0.05, 0.2]])
    radial_motions: List[float] = field(default_factory=lambda: [-10. , -20.])
    #The flare, arms and bar will be swapped when incuded in the swap lisr
    beams:  List[float] = field(default_factory=lambda: [
                                2., 4., 6., 7., 8., 10., 12.])
    # Beam across the major axis. This also set the distance as the size in kpc will be determined by Wang 2016 from the SBR profile
    snr: List[float] = field(default_factory=lambda: [0.5, 1., 3., 5.])
    # These  are average signal to noise ratios
    channelwidth: List[float] = field(default_factory=lambda: [2., 8.])

    beam_size: List[float] = field(default_factory=lambda: [[5., 5.]])
    #Resolution of the beam in arcsec
    masses:  List[float] = field(default_factory=lambda: [2.5e11])
    # The channel dependency
    # 'Options are independent, sinusoidal, hanning
    channel_dependency: str = 'sinusoidal'
    corrupt_models: bool = True


@dataclass
class ROC:
    enable: bool = True
    add_template: bool = False
    remove_template: bool = False
    delete_existing: bool = False
    base_galaxies: List[str] = field(default_factory=lambda: [
                                     'M_83', 'Circinus', 'NGC_5023', 'NGC_2903', 'NGC_3198', 'NGC_5204', 'UGC_1281', 'UGC_7774', 'ESO_223_G009'])
    variables_to_vary: List[str] = field(
        default_factory=lambda: ['Beams', 'SNR'])
    beams: List[float] = field(default_factory=lambda: [2., 4., 6., 8., -1.])
    # Beam across the major axis. This also set the distance as the size in kpc
    #will be determined by Wang 2016 from the SBR profile. -1 means the maximum possible for the ROC.
    # These  are average signal to noise ratios
    snr: List[float] = field(default_factory=lambda: [1., 3.])


@dataclass
class General:
    ncpu: int = 6
    main_directory: str = os.getcwd()
    tirific: str = "tirific"  # Command to call tirific
    sofia2: str = "sofia2"  # Command to call sofia 2
    casa: str = "casa"  # Command to call sofia 2


@dataclass
class Config:
    print_examples: bool = False
    print_bases: bool = False
    configuration_file: Optional[str] = None
    general: General = General()
    agc: AGC = AGC()
    roc: ROC = ROC()
