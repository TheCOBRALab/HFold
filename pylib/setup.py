import os
import pybind11
from pybind11.setup_helpers import build_ext
from setuptools import setup, Extension, find_packages
import sysconfig
import sys
import codecs

VERSION = "1.0.0"
DESCRIPTION = "HFold: A Python wrapper for the HFold C++ library"

here         = os.path.abspath(os.path.dirname(__file__))
project_root = os.path.abspath(os.path.join(here, ".."))
src_dir      = os.path.join(project_root, "src")
bindings_dir = os.path.join(here, "bindings")


# Read the contents of your README file to use as the long description
with codecs.open(os.path.join(here, "../README.md"), encoding="utf-8") as fh:
    long_description = "\n" + fh.read()


hfold_sources = [os.path.join(src_dir, f) for f in [
    "W_final.cpp",
    "pseudo_loop.cpp",
    "HFold.cpp",
    "cmdline.cpp",
    "Result.cpp",
    "s_energy_matrix.cpp",
    "Hotspot.cpp",
    "sparse_tree.cpp"
    ]
]

binding_sources = [
    os.path.join(bindings_dir, f) for f in [
    "pybind_module.cpp"
    ]
]

sources = hfold_sources + binding_sources

vienna_lib = os.environ.get("VIENNARNA_LIB")

if not vienna_lib:
    for lib in [
        "/usr/local/lib",
        "/usr/lib",
        "/opt/homebrew/lib",
        "/opt/vienna/lib",
        "/usr/local/lib64",
        "/usr/lib64",
        ]:
         if os.path.exists(os.path.join(lib, "libRNA.so")) or \
            os.path.exists(os.path.join(lib, "libRNA.dylib")) or \
            os.path.exists(os.path.join(lib, "libRNA.a")):
            vienna_lib = lib
            print(f"ViennaRNA library found at {vienna_lib}")
            break
    if not vienna_lib:
        # If ViennaRNA is not found in common locations, raise an error
        raise RuntimeError("Please set VIENNARNA_LIB environment variables")

ext_modules = [
    Extension(
        "hfold",
        sources,
        language="c++",
        include_dirs=[
            os.path.join(here, "../src"),
            os.path.join(here, "bindings"),
            pybind11.get_include(),
            pybind11.get_include(user=True),
        ],
        library_dirs=[vienna_lib],
        libraries=["RNA"],
        extra_compile_args=["-O3", "-std=c++17", "-DHAVE_STRDUP=1"],

    )
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
    install_requires=[],
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
        "Programming Language :: Python :: 3.13",
        "Operating System :: MacOS",
        "Operating System :: Unix",
    ],
    zip_safe=False,
    cmdclass={"build_ext": build_ext},
)
