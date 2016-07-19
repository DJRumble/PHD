"""Takes input file, applies a kappa command, and returns a specific parameter"""
import numpy as np
import os

kapdir = '/star/bin/kappa'

def PARGET(file,parameter,kappa):
    """Function takes a file, opens it in either STATS or NDFTRACE and pulls out a specific parameter, returns it and bins the tempoary file."""
    cmd = '%s/%s ndf=%s QUIET'%(kapdir,kappa,file)
    os.system(cmd)
    cmd = '%s/parget %s %s > parameter.txt'%(kapdir,parameter,kappa)
    os.system(cmd)
    a = np.loadtxt('parameter.txt',dtype='string')
    os.remove('parameter.txt')
    return a

if __name__ == "__main__":
    PARGET(file,parameter,kappa)
