from ..Coordinate.autoDFT import define

file3 = "/home/wladerer/github/Coordinate/core/damp5_converged/file3_damp5_coord"
file15 = "/home/wladerer/github/Coordinate/core/damp5_converged/file15_damp5_coord"
file18 = "/home/wladerer/github/Coordinate/core/damp5_converged/file18_damp5_coord"
file26 = "/home/wladerer/github/Coordinate/core/damp5_converged/file26_damp5_coord"
files = [file3, file15, file18, file26]

names= ["file3", "file15", "file18", "file26"]

for name, file in zip(names, files):
    os.system(f"mkdir {name}")
    os.system(f"cd {name}")
    os.system(f"cp {file} coord")
    define("/home/wladerer/github/Coordinate/utils/parameters.yaml")
    os.system("cd ..")