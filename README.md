[//]: # (Copyright &#40;c&#41; European Space Agency, 2025.)
[//]: # ()
[//]: # (This file is subject to the terms and conditions defined in file 'LICENCE.txt', which)
[//]: # (is part of this source code package. No part of the package, including)
[//]: # (this file, may be copied, modified, propagated, or distributed except according to)
[//]: # (the terms contained in the file ‘LICENCE.txt’.)
[![License: ESA permissive](https://img.shields.io/badge/ESA%20Public%20License-Permissive-blue.svg)](https://github.com/ESA-Datalabs/Euclid-Q1/blob/main/LICENSE.txt)

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

