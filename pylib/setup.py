# setup.py  (only needed if you want to keep setup.py around)
from skbuild import setup  # scikit-build or scikit-build-core exports this

setup(
    name="hfold",
    version="0.0.3",
    description="HFold: Python bindings for the HFold C++ library",
    packages=["hfold"],           # any pure-Python helper packages you ship
    python_requires=">=3.8",
    cmake_source_dir="..",
    cmake_args=[
     "-DBUILD_PYTHON=ON"
    ],
)
