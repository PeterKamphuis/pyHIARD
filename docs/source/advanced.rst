Setting your fit preferences through a yaml file.
=================================

Introduction
--------

Within pyHIARD it is possible to completely setup your own artificial input galaxy and to build variations on top of that. In principle any type of artificial galaxy can be build based on a set of input parameters.
In principle it is also possible to add any galaxy to the real galaxies although for now any addition on the available galaxies requires hard coding the package.
The settings of pyHIARD are maintained through the use of omegaconf package which means that pyHIARD can be ran with its defaults (standard database creation for comparison to other codes using pyHIARD for verification) or input can be modified in a yml file or through the command line.
In this page we explain how the database can be adapted by providing a yaml input parameters for pyHIARD. An example file with all the settings can be printed by runnning 'pyHIARD print_examples=True' after installing pyHIARD. This command prints an example yaml file with all the defaults listed.
Below we explain what each section in the yaml file does.

The yaml file is divided  in three different sections and has three independent keywords. All keywords are optional.
All these options can also be called directly from the command line when calling pyHIARD. For example whether or not to delete the existing artificial galaxies can easily be adapted by calling 'pyHIARD agc.delete_existing=True'. In the case of a list the option has to be bracketed in apostrophes, i.e. 'pyHIARD "agc.base_galaxies=[2,4]"'.

Individual Keywords
 --------
*No specifier*

**print_examples**:

  *bool, optional, default = False*

  Print an example input yaml file and an example catalogue.

**print_bases**:

    *bool, optional, default = False*

    Print the pyHIARD standard included galaxies.


**configuration_file**:

    *str, optional, default = None*

    configuration input file


The Artificial Galaxy Catalogue (AGC) keywords
--------
*Specified with agc*

**enable**:

    *bool, optional, default = True*

    Boolean to indicate whether artificial galaxies should be made (True) or not (False)

**delete_existing**

    *bool, optional, default = False*

    If true the database checks whether the base galaxies already has existing variations of the base galaxies and removes all these directories.
    If false, the code checks whether existing directories already have the final convolved cube and skips the galaxy if so.

**base_galaxies**:

    *integer List, optional, default = [1,2,3,4,5]*

    An integer list that selectes the base for the artificial galaxies. pyHIARD has 5 standard base galaxies for the AGC. These can be listed by running pyHIARD print_bases=True.
    Selecting 6 will allow the user to create their own base galaxy through answering questions about the input parameters.

**inhomogenous**:

    *bool, optional, default = True*

    If set to True random variation on top of the cylindrically symmetric disks will be created.

**symmetric**:

    *bool, optional, default = False*

    If set to true the approaching and receding side of the galaxies will differ in their warping and surface brightness profiles

**corruption_method**:

    *str, optinal, default = 'Gaussian'*

    How the artificial galaxies are corrupted. The options are Gaussian, Casa_Sim, and Casa_5.

      -Gaussian: Random noise with a gaussian distribution is added.

      -Casa_Sim: Casa's simulation method is used to invert the initial artificial
                 galaxy to the uv-plane and then converted back to the image plane and cleaned

      -Casa_5: As the simulation method is is expensive this option allows the user to only simulate every fifth artificial galaxy.

**variables_to_vary**:

    *str List, optional, default  = ['Inclination','Beams','Radial_Motions','Flare','Arms','Bar','Mass','Channelwidth','SNR','Warp','Mass','Beam_Resolution']*

    A list of the variables that should be varied within each base. If no variations on the base are required set this to ['Base'].
    The Arms, Bar and Flare will be swapped when present in this list. For the other variables values should be set in the yml file or the defaults will be used.

**masses**:

    *float List, optional, default = [2.5e11]*

    List of variations of the base mass. The mass determines the rotation curve, surface brightness profile and scale height.

**inclination**:

    *float List, optional, default = [15.,20.,30.,50.,70.,80.,88.,90.]*

    List of variations of the base inclination

**warp**:

    *float List, optional, default = [[0.15,0.05],[0.05,0.2]]*

    List of variations of the base warp angles

**radial_motons**:

    *float List, optional, default = [-5.,-10.]*

    List of variations of the radial motions in the galaxy. Negative values indicate inflows positive values indicate outflows.


**beams**:

    *float List, optional, default = [2.,4.,6.,7.,8.,10.,12.]*

    List of variations of the base beams across the major axis. Note that these also set v_sys through distance on a pure Hubble flow as the size in kpc will be determined by the relation presented in Wang et al 2016 from the SBR profile.

**snr**:

    *float List, optional, default = [1.,3.,5.]*

    List of variations of the base signal to noise ration. This is the average SNR over the whole galaxy.

**channelwidth**:

    *float List, optional, default = [2.,8.]*

    List of variations of the base channel width in km/s.

**beam_size**:

    *float List, optional, default = [[5.,5.]]*

    List of variations of the resolution of the synthesized beam in arcsec


The Real Observations Catalogue (AGC) keywords
--------
*Specified with roc*

**enable**:

    *bool, optional, default = True*

    Boolean to indicate whether observed and shifted galaxies should be made (True) or not (False)

**delete_existing**

    *bool, optional, default = False*

    If true the database checks whether the base galaxies already has existing variations of the base galaxies and removes all these directories.
    If false, the code checks whether existing directories already have the final convolved cube and skips the galaxy if so.

**base_galaxies**:

    *integer List, optional, default = ['M_83','Circinus','NGC_5023','NGC_2903','NGC_3198','NGC_5204','UGC_1281','UGC_7774']*

    List of base galaxies to vary. The default contains all possible options. Only the galaxies NGC_5204 and UGC_1281 come with the python distribution.
    Other options are downloaded from their respective survey websites and thus require an internet connection the first time they are used.

**variables_to_vary**:

    *str List, optional, default  = ['Beams','SNR']*

    A list of the variables that should be varied within each base. If no variations on the base are required set this to ['Beams'] and set roc.beams: [-1].

**beams**:

    *float List, optional, default = [2.,4.,6.,8.,-1.]*

    List of variations of the base beams across the major axis. Note that these also set v_sys through distance on a pure Hubble flow.
    -1 indicates the largest galaxy possible. This is slightly smaller than the input source extend to ensure proper blending with the noise.
    The maximums are M 83 = 21.9, Circinus = 18.2, NGC 5023 = 23.4, NGC 2903 = 16.5, NGC 3198 = 36.9, NGC 5204 = 15.0, UGC 1281 = 12.4, UGC 7774 = 13.8
    If the number of beams requested is larger than this the galaxy is skipped.

**snr**:

    *float List, optional, default = [1.,3.]*

    List of variations of the base signal to noise ration. This is the average SNR over the whole galaxy. -1 indicates the original ratio.


General settings
--------
*Specified with general*

**main_directory**:

  *str, optional, default = './'*

  The main directory where the create the database. The default is the directory where pyHIARD is started

**tirific**:

  *str, optional, default = tirific*

  Command used to call tirific from the python subprocess

**sofia2**:

  *str, optional, default = sofia2*

  Command to call sofia 2 from the python subprocess

  **tirific**:

    *str, optional, default = tirific*

    Command used to call tirific from the python subprocess

**casa**:

  *str, optional, default = casa*

  Command to call casa from the python subprocess. Note that aliases will likely not be seen by python.

**ncpu**:

  *int, optional, default = 6*

  Option to select the number of cpus. pyHIARD for not is not parallelized so this is defunct.