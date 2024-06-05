from setuptools import setup, find_packages, Extension

ext_modules = [
    Extension("matrylist.matrylist",
              ["matrylist/matrylist.pyx"],
              language="c++",
              extra_compile_args=["-std=c++17"])
]

setup(
    version='0.1.0',
    name='matrylist',
    description="Rapid interval intersections",
    author="Kez Cleal",
    author_email="clealk@cardiff.ac.uk",
    packages=find_packages(),
    install_requires=['Cython'],
    ext_modules=ext_modules,
)