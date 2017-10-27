# loomR

### An R interface for loom files

For more information on loom files, please see the documentation for [loompy](https://github.com/linnarsson-lab/loompy)

## Dependencies

loomR depends on:
 - [R](https://cran.r-project.org/) v3.4.x
 - The [hdf5r](https://github.com/hhoeflin/hdf5r) package
 - The HDF5 [C++ API](https://support.hdfgroup.org/HDF5/release/obtainsrc.html)

To get the HDF5 C++ API, please see the table below:

| Operating system | Command |
| ---------------- | ------- |
| macOS | Using [Homebrew](https://brew.sh/), `brew install hdf5` |
| Debian and Debian-based OSes | `sudo apt install libhdf5-dev` |
| Red Hat-based OSes | `sudo dnf install hdf5-devel` or `sudo yum install hdf5-devel` |
| Windows | Download precombiled binaries [here](https://github.com/mannau/h5-libwin) |
