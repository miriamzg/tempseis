from setuptools import setup, find_packages

setup(
    name="pytempseis",
    version="1.0",
    author=["Andrea Berbellini", "Alice Turner", "Auggie Marignier"],
    maintainer_email="augustin.marignier.14@ucl.ac.uk",
    description="Helpful scripts for processing and plotting higher-order moment tensor inversions",
    packages=find_packages(),
    install_requires=["numpy", "obspy", "matplotlib", "pyaml", "cartopy"],
    extras_require={"dev": ["black", "flake8"]},
)
