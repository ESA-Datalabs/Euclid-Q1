# Euclid-Q1

This repository contains tutorial Jupyter Notebooks for accessing and using Euclid-Q1 data on ESA Datalabs.

| Notebook | Description |
|----------|----------|
| [Astroquery](../data/euclid_q1/Example_Notebooks/Astroquery.ipynb) | How to access Euclid data using the [**Euclid Astroquery package**](https://astroquery.readthedocs.io/en/latest/esa/euclid/euclid.html#module-astroquery.esa.euclid) pre-installed Python library that interacts with the [Euclid Science Archive System](https://eas.esac.esa.int/sas/).|
| [ADQL Queries](../data/euclid_q1/Example_Notebooks/ADQL_examples.ipynb)  | Provides examples of querying Euclid Q1 data using [**ADQL**](https://www.ivoa.net/documents/ADQL/). Users can access this data directly on ESA Datalabs (using predefined **datalabs_path**). |
| [Cutouts](../data/euclid_q1/Example_Notebooks/Cutouts.ipynb) | Generate single **image cutouts** or in bulk from Euclid Q1 data, directly on ESA Datalabs. |
| [Source Extraction](../data/euclid_q1/Example_Notebooks/Source_extraction.ipynb) | Demonstrates a simple source extraction task from Euclid images using the **SEP** Python package. |
| [Image Visualisation](../data/euclid_q1/Example_Notebooks/Image_visualisation.ipynb) | Visualise Euclid images using different methods: **static** visualisation using Matplotlib, **interactive** visualisation using Imviz tool from Jdaviz, **sky-based** visualisation using PyESASky |
| [Spectra Visualisation](../data/euclid_q1/Example_Notebooks/Spectra_visualisation.ipynb) | Retrieve and visualise spectra from the Euclid Q1 volume using: Matplotlib (for **static** plots) and Jdaviz (for **interactive** exploration) |

**ERO Notebooks**: 

Additionally, there two notebooks available on the Euclid Early Release Observations (EROs):

| Notebook | Description |
|----------|----------|
| [Accessing ERO Data](../data/euclid_ero/Example_Notebooks/ERO_data_access.ipynb) | How to access the data in the Euclid ERO volume.             |
| [Image Colourisation](../data/euclid_ero/Example_Notebooks/ERO_image_colourisation.ipynb) | Basic concepts behind astronomical image colourisation using the Euclid ERO image of the Horsehead Nebula.          |
