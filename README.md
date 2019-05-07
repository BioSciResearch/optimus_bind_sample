<h1 align="center"> Optimus Bind 

[![Build Status](https://travis-ci.org/tcardlab/optimus_bind_sample.png?branch=master)](https://travis-ci.org/tcardlab/optimus_bind_sample) 
[![Website shields.io](https://img.shields.io/website-up-down-green-red/http/domainName.io.svg)](http://shields.io/)
[![PyPI license](https://img.shields.io/pypi/l/ansicolortags.svg)](https://pypi.python.org/pypi/ansicolortags/)
</h1> 

A short description of the project.

<details>
<summary><b>Table of contents </b></summary>

## Table of content
 - [Features](#Features)
	 - [Screenshots](#Screenshots)
 - [Setup](#Setup)
	 - [Prerequisites](#Prerequisites)
	 - [Installing](#Installing)
	 - [Run Tests](#Run-Tests)
 - [For Developers:](#Devs)
	 - <button onclick="myFunction()" href="#Project-Organization" >Project Organization</button>
	 - <a href="#Project-Organization" >Project Organization</a>
	 - [Built With](#Built-With)
	 - [Contributing](#Contributing)
 - [Citations](#Citations)
</details>

## Features

### Screenshots

## Setup
These instructions will get you a copy of the project up and running on your local machine for development and testing purposes.

### Prerequisites 

    enter code here

### Installing

    enter code here

### Run Tests

    enter code here

<details id="Devs">
<summary><h1>For Developers</h1></summary>

Project Organization
------------

    ├── LICENSE
    ├── Makefile           <- Makefile with commands like `make data` or `make train`
    ├── README.md          <- The top-level README for developers using this project.
    ├── data
    │   ├── external       <- Data from third party sources.
    │   ├── interim        <- Intermediate data that has been transformed.   (post pdbfixer)
    │   ├── processed      <- The final, canonical data sets for modeling.   (idk yet)
    │   └── raw            <- The original, immutable data dump.    (OG .pdb's etc.)
    │
    ├── docs               <- A default Sphinx project; see sphinx-doc.org for details
    │
    ├── models             <- Trained and serialized models, model predictions, or model summaries
    │
    ├── notebooks          <- Jupyter notebooks. Naming convention is a number (for ordering),
    │                         the creator's initials, and a short `-` delimited description, e.g.
    │                         `1.0-jqp-initial-data-exploration`.     (contains tests + working toward colab compatability)
    │
    ├── references         <- Data dictionaries, manuals, and all other explanatory materials.
    │
    ├── reports            <- Generated analysis as HTML, PDF, LaTeX, etc.
    │   └── figures        <- Generated graphics and figures to be used in reporting
    │
    ├── requirements.txt   <- The requirements file for reproducing the analysis environment, e.g.
    │                         generated with `pip freeze > requirements.txt`
    │
    ├── setup.py           <- makes project pip installable (pip install -e .) so src can be imported
    ├── src                <- Source code for use in this project.
    │   ├── __init__.py    <- Makes src a Python module
    │   │
    │   ├── data           <- Scripts to download or generate data
    │   │   └── make_dataset.py     (SKEMPI 2 & ZEMu downloads and management)
    │   │
    │   ├── features       <- Scripts to turn raw data into features for modeling
    │   │   └── build_features.py    (need to see example for clarification)
    │   │
    │   ├── models         <- Scripts to train models and then use trained models to make
    │   │   │                 predictions
    │   │   ├── predict_model.py
    │   │   └── train_model.py
    │   │
    │   └── visualization  <- Scripts to create exploratory and results oriented visualizations
    │       └── visualize.py
    │
    └── tox.ini            <- tox file with settings for running tox; see tox.testrun.org

## Built With
too lazy to link them atm...
 - Python Network X library
 - OpenMM/

## Contributing
Please read [CONTRIBUTING.md](https://gist.github.com/PurpleBooth/b24679402957c63ec426)(example) for details on our code of conduct, and the process for submitting pull requests to us.

**Collaboration kit:** 
 - Chat room:  [Slack](https://bioscienceclub.slack.com/messages/CHK7D10MN/details/) by request
 - Ref. library:  [F1000Workspace](https://f1000.com/work/#/items/6730972/detail?collection=321381) by request
 - Repository:  [GitHub](https://github.com/tcardlab/optimus_bind_sample)
 - To-do list:     [Trello, GitHub, other?] Scrum?

**Seeking people who are skilled in:**
 - Bioinformatics 
 - Molecular dynamics and force field development
 - Protein engineering and protein design 
 - Protein Folding and biophysics (computational or experimental) for insight into force field design
 - Basic graph theory and the Python Network X library 
 - Machine learning, particularly on small datasets 
 - Server and website design

</details>

## Citations

 1. "SKEMPI 2.0: An updated benchmark of changes in protein-protein binding energy, kinetics and thermodynamics upon mutation".  Justina Jankauskaitė, Brian Jiménez-García, Justas Dapkūnas, Juan Fernández-Recio, Iain H Moal  _**Bioinformatics**_ (2018), bty635, [https://doi.org/10.1093/bioinformatics/bty635](https://doi.org/10.1093/bioinformatics/bty635)
 2. "etc" et al.

--------

<p><small>Project based on the <a target="_blank" href="https://drivendata.github.io/cookiecutter-data-science/">cookiecutter data science project template</a>. #cookiecutterdatascience</small></p>
<!--stackedit_data:
eyJoaXN0b3J5IjpbNzE5NzA4NTUsLTExMDMyNTE3NTksMTEwND
I5NDI2NiwyMTE2NzM3NzkyLDE3ODA0MzY4OTAsOTQ1Mzg4NzUx
LC0xNjczNDE0MzUzLDIxMjQyNzg1NjUsLTM2OTMwMzU0LC0xNj
cxMTk1OTE5LC0xNDkyMTE5Njk1LDE2NDM0ODgzLDgxNjg4ODkw
OSwtMTg0NDgzOTIxOSw4MTQxMzk2MDRdfQ==
-->