from setuptools import setup, Extension
import numpy

eigen_include = 'eigen-3.4.0'

setup(name='finemapinf',
    version='1.3',
    author='Zhou Fan',
    author_email='zhou.fan@yale.edu',
    description='Python/C++ implementation of FINEMAP-inf',
    packages=['finemapinf'],
    ext_modules=[Extension('_c_funcs',
      sources=['finemapinf/py_extension.cpp',
        'finemapinf/finemapinf.cpp',
        'finemapinf/hashtables.cpp'],
      include_dirs=[numpy.get_include(), eigen_include],
      extra_compile_args=['-std=c++11'])]
    )

