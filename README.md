# pyHIARD
A project to create a comprehensive database of Artificial and real galaxies that can be used for testing the various tilted ring fitting software that is available.


Requires working versions of tirific,casa and sofia2 installed as sofia2 for full functionality

tirific is required for the Artificial branch and can be installed trough kern-suite https://kernsuite.info/ for ubuntu
The Real Observed Catalogue should work without any external installations as long as no new galaxies are added.


Please see the wiki pages for discussion and set up.

Short installation explanation
create a virtual enviroment:

`python -m venv my_DB_venv`

activate the environment:

`source my_DB_venv/bin/activate`

install pyHIARD:

`pip install /path to pyHIARD code`

cd to the directory where you want the database to exist.
run pyHIARD by typing:

`pyHIARD`

answer the questions about what you want to create.

Note: Degrading large real galaxies to very small ones takes some time due to the convolution and regridding of the cubes.
