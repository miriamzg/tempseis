# PYTEMPSEIS - Data Processing

Package for pre-processing of the observed and synthetic seismograms.  Pre-processing includes removal of instrument response, rotating seismograms to radial and transverse directions, automatic phase arrival picking, derivative calculation and windowing.

***

## Usage

Data to be processed must be saved in the following directory structure

```bash
<database>
|   CMTSOLUTION_<event_code>_GCMT
|   |   raw_data
|   |   synthetics
|   |   |   point_source
|   |   |   xm
|   |   |   ym
|   |   |   ...
|   |   |   kernels_info.txt
|   |   CMTSOLUTION_<event_code>_GCMT
|   STATIONS

```

The `<database>` can be anywhere.  Specify its path in `parameters.yml`, along with the event code, filtering periods and wave phases you wish to pick.  Then it is simply a case of running

```bash
python process_data.py
```

This program does require some user input.  In the end, the directory `<database>/processed_data` will be created, containing the filtered seismograms and some comparison figures.  It is recommended that you then copy `parameters.yml` into `<database>/processed_data` for later reference.

If you prefer, steps can be run individually by running, for example

```bash
python automatic_time_windowing.py
```

This will again read from `parameters.yml`

***

## Known Issues

* Some paths are coded using `/` so there may be problems when running on Windows systems
