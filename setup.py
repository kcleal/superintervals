from setuptools import setup, find_packages, Extension

ext_modules = [
    Extension("matrylist.matrylist",
              ["matrylist/matrylist.pyx"],
              language="c++",
              extra_compile_args=["-std=c++17"])
]

setup(
    version='0.4.2',
    name='sortedintersect',
    description="Interval intersection for sorted query and target intervals",
    author="Kez Cleal",
    author_email="clealk@cardiff.ac.uk",
    packages=find_packages(),
    install_requires=['Cython'],  # Ensure runtime dependencies are listed if any
    ext_modules=ext_modules,
)