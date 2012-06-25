HMMLab
======

HMMLab is a Hidden Markov Model editor oriented on HMMs for speach recognition. It can create, edit, train and visualize HMMs. HMMLab supports loading/saving HMMs from/to HTK files.

Compilation
```````````

Before compilation you must install:
 - `SWIG <http://swig.org/>`_
 - `GSL <http://www.gnu.org/software/gsl/>`_
 - `Graphviz <http://www.graphviz.org/>`_
<<<<<<< HEAD
 - `PyGObject <https://live.gnome.org/PyGObject/>`_ for **Python3**.
=======
 - `PyGObject <https://live.gnome.org/PyGObject/>` for **Python3**.
>>>>>>> 4dfa0703bb599b596ba22a29eaa4ed7bf622608d

You must use **Python3**.
Go to ``hmmlab/`` directory::
  cd src/
  python setup_hmmlablib.py build_ext
  cp hmmlablib.py ../bin/
  cp build/lib.linux-<arch>-3.X/_hmmlablib.cpython-XXmu.so ../bin/

The ``<arch>`` is your computer architecture and ``XX`` is your **Python3** version.
