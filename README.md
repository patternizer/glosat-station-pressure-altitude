![image](https://github.com/patternizer/glosat-station-pressure-altitude/blob/master/20crv3-pressure-curve-fit.png)
![image](https://github.com/patternizer/glosat-station-pressure-altitude/blob/master/station-20crv3-pressure-histogram.png)

# glosat-station-pressure-altitude

Python code to look up closest 20CRv3 pressure levels to station pressure altitude. Part of ongoing work for the [GloSAT](https://www.glosat.org) project: www.glosat.org. Pressure altitude is calculated with the international standard atmosphere (ISA) formula and the NIST conversion for comparison.

## Contents

* `glosat-hpa-altitude-converter.py` - python code to locate nearest 20CRv3 pressure level corresponding to station pressure altitute

## Instructions for use

The first step is to clone the latest glosat-staion-pressure-altitude code and step into the installed Github directory: 

    $ git clone https://github.com/patternizer/glosat-station-oressure-altitude.git
    $ cd glosat-station-pressure-altitude

Then create a DATA/ directory and copy to it the required inventories listed in glosat-hpa-altitude-converter.py.

### Using Standard Python

The code is designed to run in an environment using Miniconda3-latest-Linux-x86_64.

    $ python glosat-hpa-altitude-converter.py

## License

The code is distributed under terms and conditions of the [Open Government License](http://www.nationalarchives.gov.uk/doc/open-government-licence/version/3/).

## Contact information

* [Michael Taylor](michael.a.taylor@uea.ac.uk)

