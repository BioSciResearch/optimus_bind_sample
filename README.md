optimus_bind
========================
[![Build Status](https://travis-ci.org/tcardlab/optimus_bind_sample.png?branch=master)
](https://travis-ci.org/tcardlab/optimus_bind_sample)
[![Website shields.io](https://img.shields.io/website-up-down-green-red/http/domainName.io.svg)](http://shields.io/)
A short description of the project.

<details>
<summary>Table of content</summary>

## Table of content
 - [Project Organization](#Project-Organization)
 - 
</details>

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


--------

<p><small>Project based on the <a target="_blank" href="https://drivendata.github.io/cookiecutter-data-science/">cookiecutter data science project template</a>. #cookiecutterdatascience</small></p>
<!--stackedit_data:
eyJoaXN0b3J5IjpbMTAwNTEyMTcxMiwtMTg0NDgzOTIxOSw4MT
QxMzk2MDRdfQ==
-->