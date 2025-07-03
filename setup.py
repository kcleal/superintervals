from setuptools import setup, find_packages, Extension
from Cython.Build import cythonize

ext_modules = [
    Extension("superintervals.intervalmap",
              ["src/superintervals/intervalmap.pyx"],
              include_dirs=["src"],
              language="c++",
              extra_compile_args=["-std=c++17"])
]

print('PAKCAGES', find_packages(where='src'))  # Add this line for debugging

setup(
    name='superintervals',
    version='0.3.0',
    description="Rapid interval intersections",
    author="Kez Cleal",
    author_email="clealk@cardiff.ac.uk",
    packages=find_packages(where='src'),
    package_dir={"": "src"},
    ext_modules=cythonize(ext_modules),
)