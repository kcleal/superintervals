import sys

from setuptools import setup, find_packages, Extension
from Cython.Build import cythonize

MSVC = sys.platform == "win32"

ext_modules = [
    Extension("superintervals.intervalmap",
              ["src/superintervals/intervalmap.pyx"],
              include_dirs=["src"],
              language="c++",
              extra_compile_args=["/std:c++17"] if MSVC else ["-std=c++17"],
              extra_link_args=[] if MSVC else ["-lstdc++"])
]

print('PAKCAGES', find_packages(where='src'))  # Add this line for debugging

setup(
    name='superintervals',
    description="Rapid interval intersections",
    author="Kez Cleal",
    author_email="clealk@cardiff.ac.uk",
    packages=find_packages(where='src'),
    package_dir={"": "src"},
    package_data={"superintervals": ["py.typed", "*.pyi"]},
    ext_modules=cythonize(ext_modules),
)