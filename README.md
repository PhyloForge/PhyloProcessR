# PHYLOCAP

R package For processing high-throughput sequencing data from targeted sequence capture

Main functionality coming soon! Some implemented functions and their purpose are as follows: 

Coming Soon!

# Installation

PHYLOCAP is an R package that has been tested on R version 3.5 or greater. 

PHYLOCAP depends on one R package currently:
  1) data.table (>=1.12)

PHYLOCAP has two additional R package imports (i.e. you need them installed, but not loaded  with library(): 
  ape (>= 5.0)
  stringr (>= 1.4)
  
And to install PHYLOCAP, you can use the R package devtools. Here are step-by-step instructions for installation:

1) Install devtools by typing "install.packages(devtools)" in your R console. 

2) Install PHYLOCAP by typing in your R console: "devtools::install_github("chutter/PHYLOCAP")"

3) Devtools will ask you to install the package dependecies, select "Yes". If devtools asks you to update packages, you may choose to do so. I would recommend not to install packages from source if devtools asks you. Ape is problemic from source and I could not get it to install on my machine. If devtools is still giving you trouble, you can install the dependencies with "install.packages(c("ape", "stringr"))". Then rerun Step 2 and skip package updates. 

4) Devtools should finish and say the package loaded properly. Load the package with library(PHYLOCAP). 

And installation should be done! 


# Examples

1) first install and load the R package. Its a good idea to install new every time as this package is being updated frequently. The warning that the package has no updates can be ignored. 

```r
devtools::install_github("chutter/PHYLOCAP")
library(PHYLOCAP)

```

2) Alignment format conversion 




3) Alignment statistics




4) Alignment conversion within R



5) Alignment trimming 



6) Concatenate alignments 


7) Remove columns consisting of only gaps
