from os import system

import sys
if sys.version_info.major > 2:
    from subprocess import STDOUT, run
else:
    from subprocess32 import STDOUT, run

def thermal(noscat=None, nphot_therm=None, nphot_scat=None, setthreads=1, \
        inclfreefree=None, nofreefree=None, inclgascont=None, nogascont=None, \
        verbose=True, timelimit=7200, nice=None, thermpath='thermal/', radmc3d_cmd='radmc3d'):

    if nice != None:
        command="nice -{0:d} {1:s} mctherm".format(nice, radmc3d_cmd)
    else:
        command="{} mctherm".format(radmc3d_cmd)

    if (noscat == True):
        command += " noscat"
    if (nphot_therm != None):
        command += " nphot_therm {0:d}".format(nphot_therm)
    if (setthreads != 1):
        command += " setthreads {0:d}".format(setthreads)
    if (inclfreefree == True):
        command += " inclfreefree"
    if (nofreefree == True):
        command += " nofreefree"
    if (inclgascont == True):
        command += " inclgascont"
    if (nogascont == True):
        command += " nogascont"

    if not verbose:
        f = open(thermpath + "radmc3d.out","w")
        output = run(command.split(" "), cwd=thermpath, stdout=f, stderr=f, timeout=timelimit)
        f.close()
    else:
        output = run(command.split(" "), cwd=thermpath, stderr=STDOUT, timeout=timelimit)

def localfield(nphot_mono=None, verbose=True, timelimit=7200, thermpath='thermal/', radmc3d_cmd='radmc3d'):

    command="{} mcmono".format(radmc3d_cmd)

    if nphot_mono != None:
        command += " nphot_mono {0:d}".format(nphot_mono)

    if not verbose:
        f = open(thermpath + "radmc3d.out","w")
        output = run(command.split(" "), cwd=thermpath, stdout=f, stderr=f, timeout=timelimit)
        f.close()
    else:
        output = run(command.split(" "), cwd=thermpath, stderr=STDOUT, timeout=timelimit)

def image(npix=None, lambda_micron=None, iline=None, incl=None, verbose=True, timelimit=7200, thermpath='thermal/', radmc3d_cmd='radmc3d'):
    command="{} image".format(radmc3d_cmd)
    if npix != None and lambda_micron != None and incl != None:
        command += "npix {0:d} lambda {1:.6f} incl {2:.1f}".format(npix, lambda_micron, incl)


    if iline == None and lambda_micron == None:
        print('Error: must provide either iline or lambda_micron')
        return

    elif iline != None and lambda_micron == None:
        command += " iline {0:d}".format(iline)
    elif iline == None and lambda_micron != None:
        command += " lambda {0:.6f}".format(lambda_micron)

    if npix != None and lambda_micron != None and incl != None:
        command += " npix {0:d} lambda {1:.6f} incl {2:.1f}".format(npix, lambda_micron, incl)

    if not verbose:
        f = open(thermpath + "radmc3d.out","w")
        output = run(command.split(" "), cwd=thermpath, stdout=f, stderr=f, timeout=timelimit)
        f.close()
    else:
        output = run(command.split(" "), cwd=thermpath, stderr=STDOUT, timeout=timelimit)
