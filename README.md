<h1 align="center">
<span style="color:red">
    <a href="https://www.quora.com/q/hxbiokqurmxybuec">
		<img src="https://images-wixmp-ed30a86b8c4ca887773594c2.wixmp.com/f/a00a8432-c08c-46be-9d13-99373ee82e3b/d2e31wf-5e388522-269a-4ae8-9f69-ac56aa48a802.png?token=eyJ0eXAiOiJKV1QiLCJhbGciOiJIUzI1NiJ9.eyJzdWIiOiJ1cm46YXBwOjdlMGQxODg5ODIyNjQzNzNhNWYwZDQxNWVhMGQyNmUwIiwiaXNzIjoidXJuOmFwcDo3ZTBkMTg4OTgyMjY0MzczYTVmMGQ0MTVlYTBkMjZlMCIsIm9iaiI6W1t7InBhdGgiOiJcL2ZcL2EwMGE4NDMyLWMwOGMtNDZiZS05ZDEzLTk5MzczZWU4MmUzYlwvZDJlMzF3Zi01ZTM4ODUyMi0yNjlhLTRhZTgtOWY2OS1hYzU2YWE0OGE4MDIucG5nIn1dXSwiYXVkIjpbInVybjpzZXJ2aWNlOmZpbGUuZG93bmxvYWQiXX0.RnLsTImIZ3RwxiFYUMVhhIjr_V2qg_Shld0T3ZSyWfM" width="100" height="100"></span>
<img src="https://lh3.googleusercontent.com/976GNleJU-C0b-Gu67qotDub8poiFSOrG2IXvDT6WuY2bOy48gC-YmFuP3ZWPG27mmMdsPgF4zzK" width="200" height="100">
</a>

 <!--| Optimus Bind-->

[![Build Status](https://travis-ci.org/tcardlab/optimus_bind_sample.png?branch=master)](https://travis-ci.org/tcardlab/optimus_bind_sample) 
[![Inline docs](http://inch-ci.org/github/tcardlab/optimus_bind_sample.svg?branch=master)](http://inch-ci.org/github/tcardlab/optimus_bind_sample)
[![Website shields.io](https://img.shields.io/website-up-down-green-red/http/domainName.io.svg)](http://shields.io/)
[![PyPI license](https://img.shields.io/pypi/l/ansicolortags.svg)](https://pypi.python.org/pypi/ansicolortags/)
</h1> 

Optimus Bind is a collaborative computational biology project to predict the effects of mutation directly from protein sequences. The Human Genome Project has yielded a wealth of data concerning natural human genetic variation that remains to be fully utilized While genetic sequencing has provided a method of identifying potential subpopulations, the impact of specific mutations is often unknown. The first step to predict a mutation’s effect is to understand how it affects its binding partners in the protein interaction network. This program is intended to scan protein surfaces to evaluate mutations that may affect protein-protein binding. In knowing how the mutation works at the molecular level, you have made the first step to understanding how it work at the cellular and organismal level.

scientists have increasingly turned to computational methods to predict ΔΔG values (changes in the free energy ΔG upon mutation). These methods are computationally expensive for large datasets to the extent that it becomes prohibitive for genome-wide studies or even scanning mutations on a single protein. There is therefore a clear need for new methods that are both fast and accurate.

[Full Summery](https://www.quora.com/Quora-Bioscience-Club-is-considering-collaborative-computational-biology-research-projects-What-topics-are-you-interested-in-and-are-able-to-work-on/answer/Jeffrey-Brender?ch=10&share=fdebe6d2&srid=E3wB) 

### Goals <sup>[1](https://www.quora.com/q/hxbiokqurmxybuec/What-are-the-major-requirements-for-Optimus-Bind-the-collaborative-Quora-project-to-predict-the-impact-of-mutations-on)</sup>
 - Fast – Upper limit of 30 minutes per mutation.  
	 - *Problem:* The most direct approach, accurately simulating the physics of the system to guess the mutation’s effect on the binding free energy (∆∆G), can take as long as 180 hours of computer time for a single mutation.
	 - Impact: For each type of cancer, there can be hundreds of disease associated “driver” mutations.[5] Protein engineering is another case where speed is critical as each amino acid is evaluated at each position in one variant of the procedure. [6] This generates hundreds of mutations if the sites are independent and many, many more if they are not.
 - Accurate – I would like to get this to r>0.9 and an average error of <1 kcal/mol.
 - It should be open source, downloadable, and free
	 - Many computational projects are locked up in web servers. I would like a program anyone can use and, if they wish, build off of.
 - Machine Learning
	 - original program used a random forest model tried to minimize the number of features to avoid overfitting.[3] Later versions[4][5] got rid of machine learning altogether and used a linear sum of two terms.
	 - handle small molecule binding
	 - improve scoring system

### Challenges<sup>[1](https://www.quora.com/q/hxbiokqurmxybuec/Which-is-preferred-genetic-algorithms-neural-networks-or-a-combination-such-as-NEAT), [2](https://www.quora.com/q/hxbiokqurmxybuec?utm_source=quora&utm_medium=referral)</sup>
 - There isn’t a lot of data.
	 - Total number of mutations: 7085 skempi v2.0
	 - You can’t get more of it easily. 
	 - The data is not evenly distributed. 
		 - The SKEMPI database doesn’t evenly sample that mutations we want to consider. While for a few proteins there are multiple entries, for others there is nothing at all. Overall the coverage is pretty sparse. As described in more detail here, this sort of imbalanced dataset can skew the machine learning process. The model bases its predictions on the data available. When coverage is heavy in some areas and sparse in others, the accuracy of the model ends up skewed towards the subpopulations where the coverage is heaviest.
	 - While machine learning based methods that do not use physics approaches are fast and appear to have good accuracy, they are often overtrained and fall apart when confronted with new data. 
 - Less than 10% of protein complexes have structures.
	 - While it is possible in most cases to make a model of the protein complex,[13] the accuracy of the model is not perfect. 
 - Large mutation space to explore 
	 - (20 amino acids)^(the proteins length)
	 - Even if we restrict the search to single mutations, as we would if we are looking at the possible effect of SNPs, this still comes out to hundreds of mutations that need to be evaluated for each protein complex. 
 - molecular dynamics is slow and dependent upon the structure’s resolution

### Combative Design

 - The solution is to use stratified sampling so that we sample all possible cases evenly. To do this we need a feature list that defines the effectively describes the different types of proteins we might encounter.
 - Accept some errors are going to exist and at least have the option of low resolution scoring, using both residue level scoring and local estimates of the accuracy and precision of the structure (either real or predicted) as a feature in machine learning.
 - One solution is to infer the dynamics of the protein by scaling interactions directly from the sequence (similar to DynaMine and antecedent to how FoldX’s operates).

<details>
<summary><b>Table of contents </b></summary>

## Table of content
 - [Features](#Features)
	 - [Screenshots](#Screenshots)
 - [Setup](#Setup)
	 - [Prerequisites](#Prerequisites)
	 - [Installing](#Installing)
	 - [Run Tests](#Run-Tests)
 - [For Developers](#For-Developers)
	 - [Project Organization](#Project Organization)
	 - [Built With](#Built With)
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

<br>
<br>

# For Developers
<details id="Devs" open>
<summary> <strong>Toggle</strong> </summary>

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
 - OpenMM
 - pdbfixer
 - Python Network X library
 - etc.

## Contributing
Please read [CONTRIBUTING.md](https://gist.github.com/PurpleBooth/b24679402957c63ec426)(example) for details on our code of conduct, and the process for submitting pull requests to us.

**Collaboration kit:** 
 - Dev Space:  [Quora](https://www.quora.com/q/hxbiokqurmxybuec) open for Q&A and Announcements
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
eyJoaXN0b3J5IjpbMzExNTA2Mzg2LC05NDkzMDgwMDgsMTQyND
cxNDkyMCwxOTI0NDMyMTM2LC03OTM2MzQ1NjUsMjA2ODk5MjY3
MiwtNTc4NDQ2ODcsMjk4MzY1NzIwLC0xODk5MDEwMjMwLDk1Nz
g2OTQyNyw2NjM2MDI5MTYsMTA3NzQxMDIzOCwtMTk2Mjg3OTI1
LDE5NTk0NTIyMzgsLTE4MzA2ODM0MCwxNzE0ODUzNTMzLC0xMz
QzODA1MDY5LC0xMTAzMjUxNzU5LDExMDQyOTQyNjYsMjExNjcz
Nzc5Ml19
-->