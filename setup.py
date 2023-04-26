from distutils.core import setup
from distutils.extension import Extension
from Cython.Build import cythonize
import os
import os.path
import shutil


extension = Extension(
    name="pyVCT",
    sources=["libcpmfem.pyx"],
    libraries=["cpmfem"],
    library_dirs=["VCT"],
    include_dirs=["VCT"]
)
setup(
    name="pyVCT",
    ext_modules=cythonize([extension])
)

os.remove('libcpmfem.c')
shutil.rmtree('build')
if not os.path.isdir('output'):
    os.mkdir('output')
if not os.path.isdir('imgs'):
    os.mkdir('imgs')
