from pexpect import popen_spawn
import yaml 
from turbomoleio.input.define import DefineRunner
import os

def define(yaml_file: str ="parameters.yaml", convert_xyz: bool=False, xyz_file: str = "*.xyz"):
    """
        Takes coord file and parameter.yaml file and creates all requisite files for TURBOMOLE calculation

        yaml_file: filename for the yaml file that contains the user's specified TURBOMOLE parameters
        convert_xyz: a boolean to specify if there is an existing xyz file to be converted to TURBOMOLE coord format
        xyz_file: chooses the xyz file present in the current directory as default, but can be changed to an alternative file

        Output: alpha, basis, auxbasis, control, coord (iff convert_xyz == True) 
    """

    if convert_xyz: 
        
        os.system(f"x2t {xyz_file} > coord")

    with open(yaml_file, "r") as fh:
        dp = yaml.load(fh, Loader=yaml.SafeLoader)
        dr = DefineRunner(parameters=dp)
        dr.run_full()
    
def fixturbo():
    """
        fixes a bug often seen in older versions of turbomole where the orbital shift does not update
    """
    os.system(f'sed -i "s/scforbitalshift  closedshell=.05/scforbitalshift  closedshell=.3 /" control')

