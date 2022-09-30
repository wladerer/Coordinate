# Coordinate

Coordinate is a suite of tools used to find empty spaces within coordination complexes. It is a very crude approach -- essentially a heuristic based algorithm to find the least sterically hindered spaces within the first coordination sphere.

![Sampled Sphere](single_sphere.png)
![Grid of generated structures](grid.png)

## Features 

- Find sterically unhindered spaces within a coordination complex
- Create a file to view all potential spaces with dummy atoms (using Avogadro, Molden, etc.)
- Create all acceptable coordination complexes of a desired molecule - ligand pair and write to an .xyz file
- Automate TURBOMOLE dft jobs using yaml files (they're like easier to read JSON files)
  - This utility comes mostly from the [turbomoleio](https://github.com/Matgenix/turbomoleio) package maintained by Matgenix SRL
  
  
