from setuptools import setup, find_packages, Extension
from Cython.Build import cythonize

ext_modules = [
    Extension("superintervals.intervalset",
              ["superintervals/intervalset.pyx"],
              language="c++",
              extra_compile_args=["-std=c++17", "-march=native"])
]

setup(
    version='0.1.1',
    name='superintervals',
    description="Rapid interval intersections",
    author="Kez Cleal",
    author_email="clealk@cardiff.ac.uk",
    packages=find_packages(),
    install_requires=['Cython'],
    ext_modules=cythonize(ext_modules),
)