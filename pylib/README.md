# How to test the hfold python package

## 1. Make sure you have ViennaRNA installed on your system
### ViennaRNA package installation
The `hfold` package requires ViennaRNA to be installed before compilation.

1. Download the ViennaRNA package: [ViennaRNA-2.7.0.tar.gz](https://github.com/ViennaRNA/ViennaRNA/releases/download/v2.7.0/ViennaRNA-2.7.0.tar.gz)


2. Install the Package:
```bash
tar -zxvf ViennaRNA-2.7.0.tar.gz
cd ViennaRNA-2.7.0
./configure --without-perl
sudo make -j$(nproc)              # Linux
sudo make -j$(sysctl -n hw.ncpu)  # macOS
sudo make install
```
More info: [ViennaRNA GitHub](https://github.com/ViennaRNA/ViennaRNA) 

## 2. Setup a Virtual Environment
Make sure you're in a virtual environment
```bash
python3 -m venv venv
source venv/bin/activate
```

## 3. Install Python Dependencies
Make sure you have pybind11 and setuptools installed
```bash
sudo apt install pybind11-dev # linux
brew install pybind11 # MacOS
python3 -m pip install pybind11 setuptools scikit-build
```

## 4. Build the Python Extension
From the project root:
```bash
pip install -v -e ./pylib/
```

## 5. Start Python
```bash
python3
```

## 6. Run Test Commands
Once in the Python shell, you can test `hfold`:
```python
import hfold

# Basic sequence folding
hfold.hfold("AGGGCUAUCCUU", param_file="../params/rna_DirksPierce09.par")[0][0]

# Folding with suboptimal structures
hfold.hfold("AGGGCUAUCCUU", subopt=5, param_file="../params/rna_DirksPierce09.par")[0]

# Folding sequences from file
hfold.hfold(subopt=2, fileI="./sample_inputs/sample.txt", param_file="../params/rna_DirksPierce09.par")
```

## Full List of Parameters
| Parameter        | Type   | Description                     |
|------------------|--------|---------------------------------|
| `sequence`       | str    | RNA sequence                    |
| `restriction`    | str    | Folding constraint              |
| `pk_free`        | bool   | Exclude pseudoknots             |
| `pk_only`        | bool   | Only allow pseudoknots          |
| `dangles`        | int    | Dangle configuration            |
| `subopt`         | int    | Number of suboptimal structures |
| `fileI`          | str    | Input file path                 |
| `fileO`          | str    | Output file path                |
| `param_file`     | str    | Thermodynamic parameters file   |
| `noConv_given`   | bool   | Do not convert DNA into RNA     |
| `suppress_output`| bool   | Suppress console output         |


Build wheels 
```bash
python -m pip install --upgrade cibuildwheel scikit-build-core twine pybind11 setuptools # install pre-req
cibuildwheel pylib --output-dir dist # build
python -m twine upload dist/* # upload
```