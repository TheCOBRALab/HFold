import os
import sys
from setuptools import setup, Extension, find_packages
import codecs
import glob
import pybind11
from pybind11.setup_helpers import build_ext


VERSION = "1.0.0"
DESCRIPTION = "HFold: A Python wrapper for the HFold C++ library"

here = os.path.abspath(os.path.dirname(__file__))

# Read the contents of your README file to use as the long description
with codecs.open(os.path.join(here, "README.md"), encoding="utf-8") as fh:
    long_description = "\n" + fh.read()


# Recursively collect all .cpp files from src/ and bindings/, EXCEPT main.cpp
def get_cpp_sources():
    cpp_files = glob.glob("src/**/*.cpp", recursive=True)
    cpp_files += glob.glob("bindings/**/*.cpp", recursive=True)
    cpp_files.remove("src/main.cpp")
    # Exclude main.cpp or anything else
    # cpp_files = [f for f in cpp_files if not os.path.basename(f) == "main.cpp"]
    return cpp_files

def get_c_sources():
    c_files = glob.glob("src/ViennaRNA/**/*.c", recursive=True)
    c_files = [f for f in c_files if "avx512" not in f.lower()]
    # Exclude main.cpp or anything else
    # c_files = [f for f in c_files if not os.path.basename(f) == "main.cpp"]
    return c_files


cpp_sources = get_cpp_sources()
c_sources = get_c_sources()

ext_modules = [
    Extension(
        "hfold",  # This is your pybind11 module, must be C++
        cpp_sources,
        include_dirs=[
            os.path.join(here, "src"),
            os.path.join(here, "bindings"),
            pybind11.get_include(),
            pybind11.get_include(user=True),
        ],
        language="c++",
        extra_compile_args=["-O3", "-std=c++17", "-DHAVE_STRDUP=1"],
    ),
    Extension(
        "hfold_c",  # Dummy C extension just to compile C files (no pybind11)
        c_sources,
        include_dirs=[
            os.path.join(here, "src"),
            os.path.join(here, "bindings"),
        ],
        language="c",
        extra_compile_args=["-O3", "-DHAVE_STRDUP=1"],
    ),
]

# Check if the user is using Python 3.6 or higher
if sys.version_info < (3, 6):
    raise RuntimeError("Python 3.6 or higher is required.")

# Check if the user is using a supported platform
if sys.platform not in ["linux", "darwin"]:
    raise RuntimeError("This package only supports Linux and macOS platforms.")

# Check if the user is using a supported architecture
if sys.maxsize < 2**32:
    raise RuntimeError("This package only supports 64-bit architectures.")


setup(
    name="hfold",
    version=VERSION,
    author="COBRALab",
    author_email="",
    description=DESCRIPTION,
    long_description=long_description,
    long_description_content_type="text/markdown",
    packages=find_packages(exclude=["tests", "docs"]),
    ext_modules=ext_modules,
    install_requires=[
        "numpy",
        "pandas",
    ],
    python_requires=">=3.6",
    classifiers=[
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.6",
        "Programming Language :: Python :: 3.7",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
        "Programming Language :: Python :: 3.11",
        "Programming Language :: Python :: 3.12",
        "Programming Language :: Python :: 3.13"
        "Operating System :: MacOS",
        "Operating System :: Unix",
    ],
    cmdclass={"build_ext": build_ext},
    zip_safe=False,

)