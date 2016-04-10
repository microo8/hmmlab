HMMLab
======

HMMLab is a Hidden Markov Model editor oriented on HMMs for speach recognition. It can create, edit, train and visualize HMMs. HMMLab supports loading/saving HMMs from/to HTK files.

Compilation
```````````

Before compilation you must install:

 - `SWIG <http://swig.org/>`_
 - `GSL <http://www.gnu.org/software/gsl/>`_
 - `Graphviz <http://www.graphviz.org/>`_
 - `Libxml2 <http://www.xmlsoft.org/>`_
 - `PyGObject <https://live.gnome.org/PyGObject/>`_ for **Python3**.

You must use **Python3**.

Go to ``hmmlab/`` directory:

::

  $ python3 setup.py build
  $ sudo python3 setup.py install

Now you have a script in ``/usr/bin`` named ``hmmlab``.

Or on a rpm distibution you can simply download the rpm package in the `dist <https://github.com/microo8/hmmlab/tree/master/dist>`_ section. It will install all dependencies.

You must create a configuration file in ``~/.config/hmmlab.conf`` and it might look like:

::

  [mainwindow]
  width = 800
  height = 600

  [model]
  width = 80
  height = 40

Running ``$ hmmlab`` will start the **HMMLab** application. Have fun :)


Donate
``````
bitcoin:1MMV4okYsjDqq2uwG9Z3tiWiTeoVbYiCJ

.. image:: bitcoin.png
