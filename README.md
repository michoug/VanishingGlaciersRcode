# VanishingGlaciersRcode

[![DOI](https://zenodo.org/badge/687518253.svg)](https://zenodo.org/doi/10.5281/zenodo.10648971)

R code for the VanishingGlacier MAGs project with the aim to generate figures and panels for the manuscript.

## Installation

To download all the data from this repository, you first need to have Git-LFS installed. You can download it from [here](https://git-lfs.github.com/).

Then, you can clone the repository using the following command:

```bash

git clone https://github.com/michoug/VanishingGlaciersRcode.git

```

## Usage

The R code is found in the **R** folder, each script can be used independently to generate a figure/panel.
All the figures will be saved in the **Figures** folder that need to be created before running the scripts, it's not included in the repository.
The *renv* package is used to manage the R environment, you can install it using the following command:

```R

install.packages("renv")

```

Then, you can restore the environment using the following command:

```R

renv::restore()

```

If you wish, you can directly open the file **VanishingGlaciersRcode.Rproj** in *RStudio* and the environment will be automatically restored.