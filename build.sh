#!/usr/bin/env bash

python3 setup.py build && cp build/lib.linux-x86_64-3.3/hmmlab/hmmlablib/_libhmm.cpython-33m.so hmmlab/hmmlablib/ && cp build/lib.linux-x86_64-3.3/hmmlab/hmmlablib/libhmm.py hmmlab/hmmlablib/
