# OFDM for Matlab
[![DOI](https://zenodo.org/badge/141947502.svg)](https://zenodo.org/badge/latestdoi/141947502)

This is an OFDM class for MATLAB. This class allows the user to generate OFDM signals with various subcarrier spacings, cyclic prefix lengths, and number of symbols.

## How to Cite
Please cite the repo if you use it in your projects. An example bibtex entry is below.

```
@misc{TarverOFDMMatlab,
  author       = {Tarver, Chance},
  title        = {OFDM for Matlab v1.0},
  month        = jan,
  year         = 2020,
  doi          = {put apprpriate doi here from current doi above},
}
```

## How to install this: 
### Option 1: Add OFDM.m to the path:
Download this repo. Put the OFDM.m in the folder of your project or in any folder that is in the Matlab path. Then proceed to use the class to use the provided PA model or train a new model. 

### Option 2: Add as a submodule in your git project:
If you already are using git with your project, you can use this as a submodule. In the main directory of your project, run
```
git submodule add https://github.com/ctarver/OFDM-Matlab.git
git submodule update --init --recursive
```
This repo should show up in your project. Then just commit the new submodule to your project like you would commit anything. 
To use the class, you still need to add the OFDM.m to your Matlab path.
```addpath(genpath(OFDM-Matlab))```

## How to use the OFDM class:
Initialize with:
`modulator = OFDM(params);`

The inputs params is a struct with the following properties:
  - n_subcarriers
  - subcarrier_spacing
  - constellation
  - n_symbols
  - use_windowing
  
To modulate data according to the parameters, simply use the `use` method, `[tx_data, tx_constellations] = modulator.use;`. 

To demodulate data according to the parameters, simply use the `demod` method, `rx_constellations = modulator.demod(tx_data);`.

Checkout the example.m to see more.   

 ### Todo:
 - Add option to provide data instead of doing random data. 
 - Have the demodulator also demod the constellation on each subcarrier. 

## History

### v1.0 

January 29, 2020

* Added a DOI
* 1st version with a tag     
