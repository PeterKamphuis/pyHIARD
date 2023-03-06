The Configuration File
=====

Introduction
----

pyHIARD uses OmegaConf (https://github.com/omry/omegaconf) to handle the input settings. pyHIARD can be ran with default settings or settings from a yml configuration file or with command line input. For batch fitting one needs to set an input catalogue with input for each galaxy.

In a run pyHIARD first checks the defaults, then the configuration yaml and finally the command line input. This mean that if a value is set in all three input methods the one from the command line is used.

The yml file has three main sections agc, roc, general that contain several parameters that can adjust how pyHIARD runs. All parameters are described in in the section Advanced Settings <advanced.rst> one by one. Here we give a more general overview of setting up a yml configuration file.

Individual Keywords
----

PyHIARD uses three individual keywords  that are not under any section and can be called directly from the command line.

These keywords and their defaults are:

  print_examples: False

  print_bases: False

  configuration_file: None

A configuration file with all default values can be printed by

  pyHIARD print_examples=True

and the values and names of all the base galaxies available in the code can be printed with

  pyHIARD print_bases=True

To keep a universal testing database we recommend running pyHIARD with the default settings. However if you are looking for a specific orientation or effect you can set the creation parameters in the configuration file.
This is easiest done through a yml file which can be provided to pyHIARD by:

  pyHIARD configuration_file=HIARD_Input.yml

The different sections in the yml file are explained in general terms below and all individual keywords can be found in Advanced Settings <advanced.rst>.

The Artificial Galaxies Catalogue (AGC) Section
----
This section controls how the completely artificial galaxies are created. The advantage here is that we know exactly the structure of the input galaxies.
The intrinsic shape and size are almost all derived through empirical relations from the mass. It allows for different corruption methods.


The Real Observations Catalogue (ROC) Section
----
This section controls how a set of real galaxies are corrupted and shifted. The advantage here is that the input galaxies are completely realistic.
The possible base galaxies are  M 83, Circinus, NGC 5023, NGC 2903, NGC 3198, NGC 5204, UGC 1281, UGC 7774. Note that only the galaxies NGC_5204 and UGC_1281 come with the python distribution.
Other options are downloaded from their respective survey websites and thus require an internet connection the first time they are used.
If one wants to test how their own model fitting would degregrade under smoothing or increased noise one easily add their own model (FAT, Tirific, RotCur or Barolo) and fits file as an template through the add_template option.



The General section
----
This section controls several setting that are used throughout the code
