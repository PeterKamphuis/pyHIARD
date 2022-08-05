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


General Settings
--------
*Specified with general*

  **main_directory**:

    *str, optional, default = './'*

    The main directory where to create the database. The default is the directory where pyHIARD is started

  **tirific**:

    *str, optional, default = tirific*

    Command used to call tirific from the python subprocess

  **sofia2**:

    *str, optional, default = sofia2*

    Command to call sofia 2 from the python subprocess

  **tirific**:

      *str, optional, default = tirific*

      Command used to call tirific from the python subprocess

  **ncpu**:

    *int, optional, default = no of available cores -1*

    Option to select the number of cpus. If multiprocessing is set to False this parameters is defunct.

  **multiprocessing**

    *bool, optional, default = True*

    pyHIARD alows for multiprocessing through the multiprocessing module of python. This works like a charm for the AGC with Gaussian corruption.
    In case of casa corruption or the ROC the system might freeze in limited amounts of RAM. In this case please turn of the multiprocessing to create the data base.

  **font_file**

    *str, optional, default =  "/usr/share/fonts/truetype/msttcorefonts/Times_New_Roman.ttf"*

    As fonts are a nightmare in matplotlib one can set the location of their preferred font for the plotting.
    On Ubuntu the default can be obtained by installing apt-get install ttf-mscorefonts-installer. This should point to the actual font, if the file is not fond we will fall back to DejaVu Sans.


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

      *str, optinal, default = 'Tres'*

      How the artificial galaxies are corrupted. The options are No_Corrupt, Gaussian, Casa_Sim, Tres and Casa_5.

        -No_Corrupt: Do not add any noise to the model simply smooth to the required beam. Cubes names will end in '_UC' to indicate uncorrupted.

        -Gaussian: Random noise with a gaussian distribution is added. Cubes names will end in '_Gauss'

        -Casa_Sim: CASA's simulation method is used to invert the initial artificial galaxy to the uv-plane and then converted back to the image plane and cleaned. Cubes names will end in '_CS'.

        -Tres: Mix the the three corruption methods. If set pyHIARD produces three Gaussian corrupted galaxies, one uncorrupted, and one Casa corrupted galaxy.

        -Casa_5: As CASA's simulation method is is expensive this option allows the user to only simulate every fifth artificial galaxy.

  **channel_dependency**:

    *str, optional, default = sinusoidal*

    How the channels of the input cubes overlap. Possible options are independent, sinusoidal, hanning.

  **retain_unconvolved_model**

    *bool, optional, default = False*

    retain a version of the unconvolved model, if you are only interested in these models the best way to run pyHIARD is to set the corruption method to No_Corrupt

  **sim_observe_graphics**

    *bool,optional, default = False*

    Option to have sim observe produce it's graphical output. CASA's simobserve produces graphics through the front end Tcl, this means the code needs to be attached to localhost and turning this on would crash the code when running in, e.g. screen.

  **variables_to_vary**:

      *str List, optional, default  = ['Inclination','Beams','Radial_Motions','Flare','Arms','Bar','Mass','Channelwidth','SNR','Warp','Mass','Beam_Size']*

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




The Real Observations Catalogue (ROC) keywords
--------
*Specified with roc*

  **enable**:

      *bool, optional, default = True*

      Boolean to indicate whether observed and shifted galaxies should be made (True) or not (False)

  **add_template**:

      *bool, optional, default = False*

      Boolean to allow a new galaxy template to the ROC. If set to True all other options are ignored.

  **remove_template**:

      *bool, optional, default = False*

      Boolean to the removal of unwanted galaxy templates in the ROC. If set to True all other options are ignored.

  **delete_existing**

      *bool, optional, default = False*

      If true the database checks whether the base galaxies already has existing variations of the base galaxies and removes all these directories.
      If false, the code checks whether existing directories already have the final convolved cube and skips the galaxy if so.

  **base_galaxies**:

      *integer List, optional, default = ['M_83','Circinus','NGC_5023','NGC_2903','NGC_3198','NGC_5204','UGC_1281','UGC_7774','ESO_223_G009']*

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

  **minimum_degradation_factor**

      *float, optional, default = 1.25*

      The templates in the ROC have to be smoothed with a gaussian in order to have a smooth connection between the template and the artificial noise. The minimum degradation factor control how much the original template at least has to be shrank.
      This sets the maximum size in the output of the ROC.

  **max_degradation_factor**

      *float, optional, default = 1.6*

      The ROC can take a lot of memory this factor sets the size at which pixel resolution the beam template is maintained. At large degradadtion this can save a lot of memory.
