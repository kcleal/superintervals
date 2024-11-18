from setuptools import setup, find_packages, Extension
from Cython.Build import cythonize

ext_modules = [
    Extension("superintervals.intervalset",
              ["src/superintervals/intervalset.pyx"],
              include_dirs=["src"],
              language="c++",
              extra_compile_args=["-std=c++17", "-march=native"])
]

print('PAKCAGES', find_packages(where='src'))  # Add this line for debugging

setup(
    name='superintervals',
    description="Rapid interval intersections",
    author="Kez Cleal",
    author_email="clealk@cardiff.ac.uk",
    packages=find_packages(where='src'),
    package_dir={"": "src"},
    install_requires=['Cython'],
    ext_modules=cythonize(ext_modules),
)