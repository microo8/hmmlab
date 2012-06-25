HMMLab
======

HMMLab is a Hidden Markov Model editor oriented on HMMs for speach recognition. It can create, edit, train and visualize HMMs. HMMLab supports loading/saving HMMs from/to HTK files.

Compilation
-----------

Before compilation you must install `SWIG <http://swig.org/>`, `GSL <http://www.gnu.org/software/gsl/>`
You must use python3.

  cd src/
  python setup_hmmlablib.py build_ext
  cp hmmlablib.py ../bin/
  cp build/lib
