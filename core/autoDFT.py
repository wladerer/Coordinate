from pexpect import popen_spawn
import yaml 
from turbomoleio.input.define import DefineRunner
import os

def define(file, parameter_file):
    
    os.system(f"x2t {file} > coord")

    with open(parameter_file, "r") as fh:
        dp = yaml.load(fh, Loader=yaml.SafeLoader)
        dr = DefineRunner(parameters=dp)
        dr.run_full()
    
    # os.system(f'sed -i "s/scforbitalshift  closedshell=.05/scforbitalshift  closedshell=.3 /" control')

