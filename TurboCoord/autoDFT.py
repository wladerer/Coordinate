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


def autoSubmit(dir: str="."):
    """
    Checks each directory to see if a calculation has been completeed, and if so, creates a new directory run the next job
    """
    #iterate over each directory in the current directory
    for subdir in os.listdir(dir):
        #if the directory is a directory
        if os.path.isdir(subdir):
            #change to the directory
            os.chdir(subdir)
            #check if the calculation has converged
            if isConverged():
                #summarize the results
                summarize()
                #archive the files
                archive()
                #if so, create a new directory and run the next job
                nextCycle()
            #change back to the original directory
            os.chdir("..")

def isConverged():
    """
    Checks current folder to see if the file "GEO_OPT_CONVERGED" exists
    """
    if os.path.exists("GEO_OPT_FAILED"):
        print("This geometry optimization failed, please consider the parameters and orginal geometry given at the beginning of the optimization")
    
    if os.path.exists("GEO_OPT_RUNNING"):
        print("Your calculation timed out, consider restarting your calculation")

    return os.path.exists("GEO_OPT_CONVERGED")

def nextCycle():
    """
    Creates a new directory and runs the next job
    """
    os.system("mkdir next_step ; cd next_step ; cp ../coord . ; cp ../*.pbs ; cp ../parameters.yaml . ; python3 ~/.execs/autoDFT.py define")


def archive():
    """
    Archives coord, energy, control, and parameters.yaml files in the current directory
    """
    os.system("mkdir archive ; cp coord energy control parameters.yaml archive")


def summarize():
    """
    Summarizes the results of the optimization
    """
    
    #get the energy of the final geometry in the energy file
    with open("energy", "r") as fh:
        energy = fh.readlines()[-1].split()[0]
    
    #get the name of the current directory
    name = os.getcwd().split("/")[-1]

    #convert the coord file to xyz format and name it after the current directory
    os.system(f"t2x coord > {name}.xyz")

    #get the total scf cycles and time of the calculation from the control file
    with open("control", "r") as fh:
        time = fh.readlines()[0].split()[1]
        cycles = fh.readlines()[1].split()[1]

    #print the results
    print(f"Your final energy is {energy} Hartree")
    print(f"Your final geometry is in the file final_geometry.xyz")
    print(f"Your calculation took {time} seconds and {cycles} cycles")
        

#Allows you to call function in terminal 
if __name__ == '__main__':
    globals()[sys.argv[1]]()
