from pexpect import popen_spawn
import yaml 
from turbomoleio.input.define import DefineRunner
import os

def prepareTurbo():
    
    with open("parameters.yaml", "r") as fh:
        dp = yaml.load(fh, Loader=yaml.SafeLoader)
        dr = DefineRunner(parameters=dp)
        dr.run_full()
    
    # os.system(f'sed -i "s/scforbitalshift  closedshell=.05/scforbitalshift  closedshell=.3 /" control')
