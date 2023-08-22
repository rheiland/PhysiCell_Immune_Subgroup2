This `upgrade` branch is an attempt to upgrade the main repo (which is itself a fork and some sort of update to https://github.com/adriannejnner/PhysiCell_Immune_Subgroup2). We try to upgrade to a modern version (1.13.1) of PhysiCell.

The steps taken were (approximately):
```
(base) M1P~/git/adrianne_upgrade$ cp Makefile Make-old
(base) M1P~/git/adrianne_upgrade$ cp main.cpp main-old.cpp
(base) M1P~/git/adrianne_upgrade$ cd BioFVM/
(base) M1P~/git/adrianne_upgrade/BioFVM$ cp ~/dev/PhysiCell_V1.13.1/BioFVM/* .
(base) M1P~/git/adrianne_upgrade/BioFVM$ cd ../core   
(base) M1P~/git/adrianne_upgrade/core$ cp ~/dev/PhysiCell_V1.13.1/core/* .
(base) M1P~/git/adrianne_upgrade/core$ cd ../modules/
(base) M1P~/git/adrianne_upgrade/modules$ cp ~/dev/PhysiCell_V1.13.1/modules/* .
(base) M1P~/git/adrianne_upgrade/modules$ cd ..
(base) M1P~/git/adrianne_upgrade$ cp ~/dev/PhysiCell_V1.13.1/main.cpp main-new.cpp
(base) M1P~/git/adrianne_upgrade$ cp ~/dev/PhysiCell_V1.13.1/Makefile .
(base) M1P~/git/adrianne_upgrade$ cp ~/dev/PhysiCell_V1.13.1/VERSION.txt .
```

Edit main.cpp and Makefile. Shockingly, running `make` actually ran without any errors and created the `project` executable.