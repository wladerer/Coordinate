import sys
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
        
        os.system(f"x2t {xyz_file} > coord") #use the builtin TURBOMOLE command 'x2t' to convert coord

    with open(yaml_file, "r") as fh:
        dp = yaml.load(fh, Loader=yaml.SafeLoader)
        dr = DefineRunner(parameters=dp)
        dr.run_full()
    

def isConverged(dir: str =".") -> bool:

    if os.path.exists("GEO_OPT_FAILED"):
        print("This geometry optimization failed, please consider the parameters and orginal geometry given at the beginning of the optimization \n Your \"coord\" file is likely to be incorrect \n Would you like to continue? [y/n] \n : ")
        return False
    
    if os.path.exists("GEO_OPT_RUNNING"):
        print("Your calculation timed out, consider restarting your calculation")
        return False
    
    else:
        print("Converged")
        return True


def fixturbo():
    """
        fixes a bug often seen in older versions of turbomole where the orbital shift does not update
    """
    os.system(f'sed -i "s/scforbitalshift  closedshell=.05/scforbitalshift  closedshell=.3 /" control')


def nextCycle():

    if isConverged == False:
        print("It is inadvisable to continue, breaking now")
        quit()
    
    os.system("mkdir next_step ; cd next_step ; cp ../coord . ; cp ../*.pbs ; cp ../parameters.yaml . ; python3 ~/.execs/autoDFT.py define")



#Allows you to call function in terminal 
if __name__ == '__main__':
    globals()[sys.argv[1]]()
