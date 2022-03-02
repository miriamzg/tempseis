# tempseis

Perform higher-order moment tensor inversions

## Requirements

The [Neighbourhood Algorithm](http://www.iearth.org.au/codes/NA/#:~:text=Overview-,The%20neighbourhood%20algorithm%20is%20a%20two%2Dstage%20numerical%20procedure%20for,search%20technique%20for%20global%20optimization.&text=The%20objective%20is%20to%20find,a%20user%20supplied%20objective%20function.) must be obtained from the original authors.

[Python 3](https://www.python.org/) is used for preprocessing and plotting.

## Inversion

This directory contains the problem-specific source code for higher-order moment tensor inversions.  

### Building

In `Inversion/src` there is the main routine `homti_na.f`, some subroutines and a `Makefile`.  To build the code, run

```bash
cd Inversion/src
make NA_DIR=</path/to/NA/package> all
```

The `all` target will build the necessary parts of the NA and link with the routines here.  The default value of `NA_DIR` is `../NA`.

## pytempseis

This python package contains helpful scripts for data processing and plotting.  Install in a `python3` environment using `pip`

```bash
cd tempseis
python3 -m venv .venv
source .venv/bin/activate
pip install .
```

This will make sure all the correct packages and functions are available to the various scripts.  Feel free to use a different environment name.

### Data processing

In `pytempseis/dataprocessing`, edit the `parameters.yml` as required then run

```bash
python process_data.py
```

to run through all the processing steps.  Alternatively you can run the individual `.py` scripts if you wish to process step-by-step.

### Plotting

In `pytempseis/plot` you'll find helpful plotting scripts, most of which will take a result direction as input.  TODO: use `argparse`

## Run Everything

Once the data has been preprocessed, the easiest thing to do is

```bash
./total_proc_script
```

from the top level directory.  This will copy data to the right place so that it can be found by the NA code, run the inversion, and plot the results.  Everything should then be in the `Results` directory.