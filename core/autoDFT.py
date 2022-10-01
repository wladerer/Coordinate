from pexpect import popen_spawn
import yaml 
from turbomoleio.input.define import DefineRunner
import os

def define():
    """
        Takes xyz file and parameter.yaml file and creates all requisite files for TURBOMOLE calculation
    """
    
    # os.system(f"x2t {file} > coord")

    with open("/home/wladerer/github/Coordinate/utils/parameters.yaml", "r") as fh:
        dp = yaml.load(fh, Loader=yaml.SafeLoader)
        dr = DefineRunner(parameters=dp)
        dr.run_full()
    
def fixturbo():
    """
        fixes a bug often seen in older versions of turbomole where the orbital shift does not update
    """
    os.system(f'sed -i "s/scforbitalshift  closedshell=.05/scforbitalshift  closedshell=.3 /" control')

define()