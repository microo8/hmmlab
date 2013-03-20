'''
This file is part of HMMLab.

HMMLab is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

HMMLab is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with HMMLab.  If not, see <http://www.gnu.org/licenses/>.
'''
from distutils.core import setup, Extension

module1 = Extension('hmmlab.hmmlablib._libhmm',
                    sources=['hmmlab/hmmlablib/vmlib.cpp',
                             'hmmlab/hmmlablib/libhmm.cpp',
                             'hmmlab/hmmlablib/gnuplot_pipes.cpp',
                             'hmmlab/hmmlablib/libhmm.i'],
                    swig_opts=['-c++', '-py3'],
                    include_dirs=['/usr/include/graphviz', '/usr/include/gsl/'],
                    runtime_library_dirs=['/usr/lib/graphviz/', '/usr/lib/'],
                    libraries=["gvc","graph","cdt","gsl","gslcblas","m"],
                    extra_compile_args=['-Wno-write-strings'])

setup(  name            = 'hmmlab',
        version         = '0.7',
        author          = 'Bc. Vladimir Magyar',
        author_email    = 'magyarvladimir@gmail.com',
        description     = 'HMMLab is a Hidden Markov Model editor oriented on HMMs for speach recognition. It can create, edit, train and visualize HMMs. HMMLab supports loading/saving HMMs from/to HTK files.',
        license         = 'GPLv3',
        url             = 'https://github.com/microo8/hmmlab',
        ext_modules     = [module1],
        packages        = ['hmmlab', 'hmmlab.hmmlablib'],
        package_dir     = {'hmmlab' : 'hmmlab', 'hmmlab.hmmlablib': 'hmmlab/hmmlablib'},
        package_data    = {'hmmlab' : ['licence', 'glade/*.glade', 'glade/*.png']},
        scripts         = ['hmmlab/hmmlab'])

