from setuptools import setup, Extension
import os
import pybind11
from pybind11.setup_helpers import build_ext
from setuptools import setup, Extension, find_packages
import sysconfig
import sys
import codecs

VERSION = "1.0.0"
DESCRIPTION = "HFold: A Python wrapper for the HFold C++ library"

here = os.path.abspath(os.path.dirname(__file__))

# Read the contents of your README file to use as the long description
with codecs.open(os.path.join(here, "README.md"), encoding="utf-8") as fh:
    long_description = "\n" + fh.read()

hfold_sources = [
    "src/W_final.cpp",
    "src/pseudo_loop.cpp",
    "src/HFold.cpp",
    "src/cmdline.cpp",
    "src/Result.cpp",
    "src/s_energy_matrix.cpp",
    "src/Hotspot.cpp",
    "src/sparse_tree.cpp",
    # Add any additional Python bindings here if needed
    "bindings/pybind_module.cpp"
]

vienna_sources = [
    # Constraints
    "src/ViennaRNA/constraints/constraints.c",
    "src/ViennaRNA/constraints/hard.c",
    "src/ViennaRNA/constraints/SHAPE.c",
    "src/ViennaRNA/constraints/soft.c",

    # Datastructures
    "src/ViennaRNA/datastructures/basic_datastructures.c",
    "src/ViennaRNA/datastructures/char_stream.c",
    "src/ViennaRNA/datastructures/lists.c",

    # I/O
    "src/ViennaRNA/io/file_formats.c",
    "src/ViennaRNA/io/io_utils.c",

    # Landscape
    "src/ViennaRNA/landscape/move.c",

    # Loops
    "src/ViennaRNA/loops/external_pf.c",
    "src/ViennaRNA/loops/external.c",
    "src/ViennaRNA/loops/hairpin.c",
    "src/ViennaRNA/loops/internal.c",
    "src/ViennaRNA/loops/multibranch.c",

    # Parameters
    "src/ViennaRNA/params/default.c",
    "src/ViennaRNA/params/params.c",
    "src/ViennaRNA/params/io.c",

    # Utilities
    "src/ViennaRNA/utils/cpu.c",
    "src/ViennaRNA/utils/higher_order_functions.c",
    "src/ViennaRNA/utils/string_utils.c",
    "src/ViennaRNA/utils/structure_utils.c",
    "src/ViennaRNA/utils/utils.c",

    # Core ViennaRNA
    "src/ViennaRNA/alphabet.c",
    "src/ViennaRNA/boltzmann_sampling.c",
    "src/ViennaRNA/centroid.c",
    "src/ViennaRNA/cofold.c",
    "src/ViennaRNA/commands.c",
    "src/ViennaRNA/dp_matrices.c",
    "src/ViennaRNA/equilibrium_probs.c",
    "src/ViennaRNA/eval.c",
    "src/ViennaRNA/fold_compound.c",
    "src/ViennaRNA/fold.c",
    "src/ViennaRNA/gquad.c",
    "src/ViennaRNA/grammar.c",
    "src/ViennaRNA/MEA.c",
    "src/ViennaRNA/mfe.c",
    "src/ViennaRNA/mm.c",
    "src/ViennaRNA/model.c",
    "src/ViennaRNA/part_func.c",
    "src/ViennaRNA/ribo.c",
    "src/ViennaRNA/sequence.c",
    "src/ViennaRNA/subopt.c",
    "src/ViennaRNA/unstructured_domains.c",
]

all_sources = hfold_sources + vienna_sources

ext_modules = [
    Extension(
        "hfold",
        all_sources,
        language="c++",
        include_dirs=[
            os.path.join(here, "src"),
            os.path.join(here, "bindings"),
            pybind11.get_include(),
            pybind11.get_include(user=True),
        ],
        extra_compile_args=["-O3", "-std=c++17", "-DHAVE_STRDUP=1"],
        extra_link_args=[
        "-Wl,--unresolved-symbols=ignore-in-shared-libs",  # optional
        f"-L{sysconfig.get_config_var('LIBDIR')}",
        f"-lpython{sysconfig.get_config_var('VERSION')}",
    ]
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
    zip_safe=False,
    cmdclass={"build_ext": build_ext},
)
