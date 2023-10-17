# FOWT_framework

To run the code you need to install the following packages in the conda environment:

- meshmagick - pip install https://github.com/LHEEA/meshmagick/archive/master.zip
- gmsh - pip install --upgrade gmsh
- pygmsh -pip install pygmsh
- h5py - conda install h5py
- h5netcdf - conda install -c conda-forge h5netcdf
- capytaine - conda install -c conda-forge capytaine
- future - conda install future

Also, remember that pyFAST is a third-party package and needs to be updated every now and then when working with new versions of OF

To run simFOWT.py or main.py, the working directory MUST be "FOWT_optim_test"
Please provide at least the .fst or the .sim file names in the definition file, to start the simulation software
