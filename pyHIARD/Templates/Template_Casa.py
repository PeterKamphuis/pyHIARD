import os
print "-----This script is for putting the model in the visibilities-----"
#First import the cube
os.system('rm -Rf in_cube')
default(importfits)
fitsimage = 'unconvolved_cube.fits'  # Name of input image FITS file
imagename = 'in_cube'  # Name of output CASA image
whichrep = 0  # If fits image has multiple coordinate
whichhdu = -1  # If its file contains multiple images,
zeroblanks = True  # Set blanked pixels to zero (not NaN)
overwrite = True  # Overwrite pre-existing imagename
defaultaxes = False  # Add the default 4D coordinate axes
defaultaxesvalues = []  # List
beam = []  # List of values to be used to define
importfits()
# And import a cleaning mask
os.system('rm -Rf smask.image')
default(importfits)
fitsimage = 'mask.fits'  # Name of input image FITS file
imagename = 'smask.image'  # Name of output CASA image
whichrep = 0  # If fits image has multiple coordinate
whichhdu = -1  # If its file contains multiple images,
zeroblanks = True  # Set blanked pixels to zero (not NaN)
overwrite = True  # Overwrite pre-existing imagename
defaultaxes = True  # Add the default 4D coordinate axes
defaultaxesvalues = ['12h42m40.9s', '+14d17m45s', '1.42064Ghz', 'I']  # List
beam = []  # List of values to be used to define
importfits()
#Make sure that the header values play nice with casa
default(imhead)
imagename = 'smask.image'  # Name of the input image
mode = 'put'  # Mode of operation: "add", "del",
hdkey = 'telescope'
hdvalue = 'WSRT'
imhead()
hdkey = 'date-obs'
hdvalue = '2017/12/24/12:00:00'
imhead()
mode = 'list'
imhead()
# Create a mask to use in the cleaning
ia.open('smask.image')
ia.calcmask(mask='smask.image > 0.5', name='mymask1')
ia.close()
os.system('rm -Rf mask.image')
default(makemask)
mode = 'expand'
inpimage = 'smask.image'
inpmask = 'smask.image:mymask1'
output = 'mask.image'
overwrite = True
inp(makemask)
makemask()
# Create the simulated galaxies
os.system('rm -Rf simulated')
default(simobserve)
#  simobserve :: visibility simulation task
project = 'simulated'  # root prefix for output file names
skymodel = 'in_cube'  # model image to observe
complist = ''  # componentlist to observe
setpointings = False
ptgfile = 'pntings.txt'
integration = '900s'  # integration (sampling) time
direction = ''  # "J2000 19h00m00 -40d00m00" or "" to center on model
mapsize = ['', '']  # angular size of map or "" to cover model
maptype = 'square'  # hexagonal, square (raster), ALMA, etc
# spacing in between pointings or "0.25PB" or "" for ALMA default INT=lambda/D/sqrt(3), SD=lambda/D/3
pointingspacing = ''
# observation mode to simulate [int(interferometer)|sd(singledish)|""(none)]
obsmode = 'int'
antennalist = 'WSRT.cfg'  # interferometer antenna position file
refdate = '2017/12/25'  # date of observation - not critical unless concatting simulations
# hour angle of observation center e.g. "-3:00:00", "5h", "-4.5" (a number without units will be
hourangle = 'transit'
#   interpreted as hours), or "transit"
totaltime = '12h'  # total time of observation or number of repetitions
caldirection = ''  # pt source calibrator [experimental]
calflux = '1Jy'
thermalnoise = 'tsys-atm'  # add thermal noise: [tsys-atm|tsys-manual|""]
t_sky = 260
tau0 = 0.01
user_pwv = 0.  # Precipitable Water Vapor in mm
t_ground = 283.0  # ambient temperature
seed = 11111  # random number seed
leakage = 0.0  # cross polarization (interferometer only)
# display graphics at each stage to [screen|file|both|none]. Have to be off to be able to run in screen.
graphics = 'none'
verbose = False
overwrite = True  # overwrite files starting with $project
inp(simobserve)
simobserve()
# And invert and clean the visibilities using the pre-determined mask
default(tclean)
vis = 'simulated/simulated.WSRT.noisy.ms'
usemask = 'user'
restart = False
os.system('rm -Rf Final_Cube*')
imagename = 'Final_Cube'
niter = 1000
threshold = '1e-4Jy/beam'
mask = 'mask.image'
imsize = [256, 256]
smallscalebias = 0.6
cell = ['3arcsec', '3arcsec']
scales = [0, 10, 25]
uvtaper = ['20arcsec', '20arcsec', '0']
pblimit = -1.0
pbmask = 0.0
restoringbeam = 'common'
chanchunks = -1
overwrite = False
specmode = 'cube'
start = 1
outframe = 'LSRK'
width = 1
restfreq = '1.420405752GHz'
veltype = 'radio'
uvrange = ''
field = '0'
spw = ''
weighting = 'briggs'
gridder = 'wproject'
deconvolver = 'multiscale'
wprojplanes = 128
robust = -2.0
interactive = False
inp(tclean)
tclean()
# export the final Cube to a fits file
default(exportfits)
imagename = 'Final_Cube.image'  # Name of input CASA image
fitsimage = 'Convolved_Cube.fits'  # Name of output image FITS file
velocity = True  # Use velocity (rather than frequency) as spectral axis
optical = False  # Use the optical (rather than radio) velocity convention
bitpix = -32  # Bits per pixel
# Minimum pixel value (if minpix > maxpix, value is automatically determined)
minpix = 0
# Maximum pixel value (if minpix > maxpix, value is automatically determined)
maxpix = -1
overwrite = True  # Overwrite pre-existing imagename
dropstokes = False  # Drop the Stokes axis?
stokeslast = True  # Put Stokes axis last in header?
history = True  # Write history to the FITS image?
dropdeg = True  # Drop all degenerate axes (e.g. Stokes and/or Frequency)?
exportfits()
# And clean up. If you want to keep the casa products just comment out this line
os.system('mkdir Casa_Log')
os.system('mv simulated/*.png ./Casa_Log/')
os.system('mv *.last ./Casa_Log/')
os.system('mv run_casa.py ./Casa_Log/')
os.system('mv pntings.txt ./Casa_Log/')
os.system('rm -Rf in_cube smask.image mask.image simulated Final_Cube.* casa*.log')
