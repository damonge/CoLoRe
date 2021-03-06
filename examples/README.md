# Examples

This folder contains a set of examples of use cases for CoLoRe. The main features of each of them are described below:

## 1 `simple`

This example contains all of the current existing features in CoLoRe for a very small simulation box, modest source sample and low map resolution (i.e. something that can easily be run on a laptop). The [param file](./simple/param.cfg) contains descriptions of all the parameters CoLoRe can take. Try running CoLoRe from the root directory with this file as input, and read some of the outputs using the [`read_colore.py`](./read_colore.py) script. You can also read the density field generated by CoLoRe as part of the simulation using the [`read_grid.py`](./read_grid.py) script provided.

The script [`read_skewers.py`](./simple/read_skewers.py) shows how to process skewers along source lines of sight for Lyman-alpha-type studies.


## 2 `cl_test`

This example can also be run on a laptop, and can be used as a quantitative test that the outputs of the code are as expected. Try to run CoLoRe from the root directory with its [param file](./cl_test/param.cfg) as input. Running the script [`cl_test.py`](./cl_test/cl_test.py) from the root directory will compare the angular power spectra of the different tracers generated by CoLoRe with their theoretical expectation. This script also shows how to access the CMB lensing and ISW maps output by CoLoRe.


## 3 `LSST`

The [param file](./LSST/param.cfg) associated with this example shows the typical parameters one would like to run CoLoRe with to generate useful simulations for an LSST-like sample. In this case the memory requirements are extensive (~1 TB).
