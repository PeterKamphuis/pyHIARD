# -*- coding: future_fstrings -*-
try:
    from importlib.metadata import version
    from importlib.util import find_spec
except ImportError:
    # Try backported to PY<37 `importlib_resources`.
    # For Py<3.9 files is not available
    from importlib_metadata import version 

import os
import subprocess

def report_version():
    # Distutils standard  way to do version numbering
    try:
        __version__ = version("pyHIARD")
    except:
        __version__ = "dev"
    # perhaps we are in a github with tags; in that case return describe
    path = os.path.dirname(os.path.abspath(__file__))
    try:
        # work round possible unavailability of git -C
        result = subprocess.check_output(
            'cd %s; git describe --tags' % path, shell=True, stderr=subprocess.STDOUT).rstrip().decode()
    except subprocess.CalledProcessError:
        result = None
    if result != None and 'fatal' not in result:
        # will succeed if tags exist
        return result
    else:
        # perhaps we are in a github without tags? Cook something up if so
        try:
            result = subprocess.check_output(
                'cd %s; git rev-parse --short HEAD' % path, shell=True, stderr=subprocess.STDOUT).rstrip().decode()
        except subprocess.CalledProcessError:
            result = None
        if result != None and 'fatal' not in result:
            return __version__+'-'+result
        else:
            # we are probably in an installed version
            return __version__


__version__ = report_version()

try:
    if find_spec('casatasks') is not None:
        __casa_max__ =  'OK'
    else:
        __casa_max__ =  'Smaller'
except:
    try:
        import casatasks
        __casa_max__ =  'OK'
    except:
        __casa_max__ =  'Smaller'

