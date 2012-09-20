from distutils.core import setup, Extension

module1 = Extension('hmmlab.hmmlablib._libhmm',
                    sources=['hmmlab/hmmlablib/vmlib.cpp',
                             'hmmlab/hmmlablib/libhmm.cpp',
                             'hmmlab/hmmlablib/libhmm.i'],
                    swig_opts=['-c++', '-py3'],
                    include_dirs=['/usr/include/graphviz', '/usr/include/gsl/'],
                    runtime_library_dirs=['/usr/lib/graphviz/', '/usr/lib/'],
                    libraries=["gvc","graph","cdt","gsl","gslcblas", "m"],
                    extra_compile_args=['-Wno-write-strings'])

setup(  name            = 'hmmlab',
        version         = '0.6',
        author          = 'Bc. Vladimir Magyar',
        author_email    = 'magyarvladimir@gmail.com',
        description     = 'HMMLab is a Hidden Markov Model editor oriented on HMMs for speach recognition. It can create, edit, train and visualize HMMs. HMMLab supports loading/saving HMMs from/to HTK files.',
        license         = 'GPLv3',
        url             = 'https://github.com/microo8/hmmlab',
        platforms       = ['x86_64'],
        ext_modules     = [module1],
        packages        = ['hmmlab', 'hmmlab.hmmlablib'],
        package_dir     = {'hmmlab' : 'hmmlab', 'hmmlab.hmmlablib': 'hmmlab/hmmlablib'},
        package_data    = {'hmmlab' : ['glade/*.glade']} )

