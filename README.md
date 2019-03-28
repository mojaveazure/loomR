![Travis](https://img.shields.io/travis/mojaveazure/loomR.svg?style=for-the-badge)
<!-- ![AppVeyor](https://img.shields.io/appveyor/ci/mojaveazure/loomR.svg?style=for-the-badge) -->
<!-- ![CRAN](https://img.shields.io/cran/v/loomR.svg?style=for-the-badge) -->

# loomR

### An R interface for loom files

For more information on loom files, please see the documentation for [loompy](https://github.com/linnarsson-lab/loompy)

## Tutorial

A tutorial for loomR can be found [here](http://satijalab.org/loomR/loomR_tutorial.html). A full function and method reference can be found [here](http://satijalab.org/loomR/loomR.pdf).

## Compatability with loompy

loomR aims to be completely compatible with loompy. Currently, loomR implements the following methods of the loompy API:
 - map/apply
 - create
 - connect
 - combine
 - subset
 - add layer
 - add attriute
 - add graph
 - add cells
 - add loom

## Dependencies

loomR depends on:
 - [R](https://cran.r-project.org/) v3.4.x
 - The [R6](https://cran.r-project.org/package=R6) package
 - The [hdf5r](https://cran.r-project.org/package=hdf5r) package
 - The HDF5 [C++ API](https://support.hdfgroup.org/HDF5/release/obtainsrc.html)

To get the HDF5 C++ API, please see the table below:

| Operating system | Command |
| ---------------- | ------- |
| macOS | Using [Homebrew](https://brew.sh/), `brew install hdf5@1.8` |
| Debian and Debian-based OSes | `sudo apt install libhdf5-dev` |
| Red Hat-based OSes | `sudo dnf install hdf5-devel` or `sudo yum install hdf5-devel` |
| Windows | Download precombiled binaries from Mario Annau [here](https://github.com/mannau/h5-libwin) |

Please note, loomR and hdf5r require HDF5 between version 1.8.13 and 1.10.1
