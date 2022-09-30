from pexpect import popen_spawn
import yaml 
from autoDFT import define

file = "/home/wladerer/github/Coordinate/core/less-terts-noferr-thf-def2-TZVPP.xyz"
yaml_file = "/home/wladerer/github/Coordinate/utils/parameters.yaml"
define(file, yaml_file)