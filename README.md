# Euclid-Q1

This repository contains tutorial Jupyter Notebooks for accessing and using Euclid-Q1 data on ESA Datalabs. Description of the Euclid Q1 data release is available at: https://www.cosmos.esa.int/web/euclid/euclid-q1-data-release 

| Notebook | Description |
|----------|----------|
| Astroquery | How to access Euclid data using the [**Euclid Astroquery package**](https://astroquery.readthedocs.io/en/latest/esa/euclid/euclid.html#module-astroquery.esa.euclid) pre-installed Python library that interacts with the [Euclid Science Archive System](https://eas.esac.esa.int/sas/).|
| ADQL Queries | Provides examples of querying Euclid Q1 data using [**ADQL**](https://www.ivoa.net/documents/ADQL/). Users can access this data directly on ESA Datalabs (using predefined **datalabs_path**). |
| Cutouts| Generate single **image cutouts** or in bulk from Euclid Q1 data, directly on ESA Datalabs. |
| Source Extraction| Demonstrates a simple source extraction task from Euclid images using the **SEP** Python package. |
| Image Visualisation | Visualise Euclid images using different methods: **static** visualisation using Matplotlib, **interactive** visualisation using Imviz tool from Jdaviz, **sky-based** visualisation using PyESASky |
| Spectra Visualisation | Retrieve and visualise spectra from the Euclid Q1 volume using: Matplotlib (for **static** plots) and Jdaviz (for **interactive** exploration) |

**ERO Notebooks**: 

Additionally, there two notebooks available on the Euclid Early Release Observations (EROs):

| Notebook | Description |
|----------|----------|
| Accessing ERO Data | How to access the data in the Euclid ERO volume.             |
| Image Colourisation | Basic concepts behind astronomical image colourisation using the Euclid ERO image of the Horsehead Nebula.          |

# Disclaimer 

Please note that the tutorials in this repository are subject to change. The current notebooks were developed with the intent to use them inside ESA Datalabs, so all the filepaths point to the Euclid Q1 data volume mounted in ESA Datalabs. It is possible to modify the notebooks to run outside of the platform if needed, but there are a few things to consider:

1. You need to install all the packages used in the notebooks
2. When trying to access any of the data, the products first need to be downloaded using the [**Euclid Astroquery package**](https://astroquery.readthedocs.io/en/latest/esa/euclid/euclid.html#module-astroquery.esa.euclid) get_product method.

After downloading the data you can run the rest of the notebook by substituting the ESA Datalabs file path(s) set in the beginning of the notebook with the one(s) set on the local machine.

The ESA Datalabs platform is still in development, and user access is experimental. Users may encounter downtime and system instability, which could impact user experience. 

