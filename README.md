# tempseis

Perform higher-order moment tensor inversions

## pytempseis

This python package contains helpful scripts for data processing and plotting.  Install in a `python3` environment using `pip`

```bash
cd tempseis
python3 -m venv .venv
source .venv/bin/activate
pip install .
```

### Data processing

In `pytempseis/dataprocessing`, edit the `parameters.yml` as required then run
```bash
python process_data.py
```
to run through all the processing steps.  Alternatively you can run the individual `.py` scripts if you wish to process step-by-step.