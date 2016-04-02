- [GeauxDock: GPU Accelerated Molecular Docking](#sec-1)
- [Introduction](#sec-2)

# GeauxDock: GPU Accelerated Molecular Docking<a id="orgheadline1"></a>

A LA-SiGMA Software Distribution

Copyright 2016 Louisiana State University

# Introduction<a id="orgheadline2"></a>

GeauxDock is an ultra-fast automated docking program, designed to predict
how small ligands bind to pharmacologically relevant macromolecules.

GeauxDock employs a novel hybrid force field and a Monte Carlo protocol for 
the efficient sampling of conformational space.


A full featured version of GeauxDock runs on the following device:
1.  Multi-core CPU
2.  Xeon Phi
3.  Nvidia GPU (require CC > 3.0), tunned for Kepler GK110 and Maxwell 2.0

This branch of code is significnatly simplified for quickly developing
new features. The GPU computation is not available in this branch.



GeauxDock has applications in:

1.  ligand virtual screening
2.  structure-based drug design
3.  drug repositioning and polypharmacology

Step-by-step setup and operating instructions can be found in 
doc/instructions.txt.

The Official site of GeauxDock project is here
<http://lasigma.loni.org/package/dock/>

For the latest version and other resources visit
<https://github.com/Yeahhhh/geauxdock_cpumic>

LA-SiGMA, the Louisiana Alliance for Simulation-Guided Materials
Applications, is a statewide interdisciplinary collaboration of
material and computer scientists developing computational resources
tuned for the latest generation of accelerator-equipped systems. The
Alliance also develops graduate curricula and is engaged in multiple
outreach activities. Visit us at <http://lasigma.loni.org> .

The accelerator ports of this code were developed by Ye Fang, Yun Ding,
with assistance from Wei Feinstein, David Koppelman, Juana Moreno,
Daniel Case, J. Ramanujam, Michal Brylinski and Mark Jarrell.

This work was supported in part by the National Science Foundation
under the NSF EPSCoR Cooperative Agreement No. EPS-1003897 with
additional support from the Louisiana Board of Regents.
