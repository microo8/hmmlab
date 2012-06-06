from distutils.core import setup, Extension

module1 = Extension('_hmmlablib',
                    sources=['vmlib.cpp', 'hmmlablib.cpp','hmmlablib.i'],
                    swig_opts=['-c++', '-py3', '-builtin'],
                    include_dirs=['/usr/include/graphviz', '/usr/include/gsl/'],
                    library_dirs=[''],
                    runtime_library_dirs=['/usr/lib/graphviz/', '/usr/lib/'],
                    libraries=["gvc","graph","cdt","gsl","gslcblas", "m"],
                    extra_compile_args=['-Wno-write-strings'])

setup(  name            = 'hmmlablib',
        version         = '0.6',
        author          = 'Bc. Vladimir Magyar',
        author_email    = 'magyarvladimir@gmail.com',
        description     = '',
        ext_modules     = [module1],
        packages        = ['hmmlablib'])

