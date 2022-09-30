from pexpect import popen_spawn
import yaml 
from turbomoleio.input.define import DefineRunner

yaml_file = "parameters.yaml"

with open(yaml_file, "r") as fh:
  dp = yaml.load(fh, Loader=yaml.SafeLoader)

dr = DefineRunner(parameters=dp)
dr.run_full()