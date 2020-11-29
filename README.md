<!-- README TITLE -->
# BIOF501A Final Project: Automating lineage recognition from COVID19 reads 

<!-- TABLE OF CONTENTS -->
## Table of Contents

* [Project Information](#project-information)
  + [Dependencies](#dependencies)
* [Getting Started](#getting-started)
  + [Installation](#installation)
* [Usage](#usage)
* [Expected Output](#expected-output)
* [Team Members](#Team-members)

***
<!-- PROJECT INFORMATION -->
## Project Information

This is the repository for the *BIOF501A Final Project* . 

### Background Information and Hypothesis

idea of this project is to offer a dynamic workflow to call lineages
clinically relevant because physicians might notice everybody has a lineage from this one area of the world - perhaps in the future informing of local outbreaks 
currently, only 1 sample is processed, but for the imaginable future, more can be used. 

### Dependencies

Currently, [https://ubc-mds.github.io/resources_pages/install_ds_stack_mac/#git](git), [https://ubc-mds.github.io/resources_pages/install_ds_stack_mac/#r-xquartz-irkernel-and-rstudio](R), [https://ubc-mds.github.io/resources_pages/install_ds_stack_mac/#latex](LaTeX), and [https://docs.conda.io/en/latest/miniconda.html](Miniconda) are required.

***
<!-- GETING STARTED -->

## Getting Started

All milestones will be titled *milestones1*, *milestones2*, etc. for milestone projects. Worksheets will be titled *worksheet1*, *worksheet2*, etc. for worksheet projects. 

### Installation

To clone the repository, run the following shell command: 
```sh
git clone https://github.com/github_username/repo_name.git
```

***
<!-- USAGE -->

## Usage 

Once the workflow has been acquired, more effort needs to be done to get the full installation. 

Basically just follow the normal pangolin installation:

```sh
git clone https://github.com/cov-lineages/pangolin.git 
cd pangolin
conda env create -f environment.yml
conda activate pangolin
python setup.py install
```

conda activate /Users/owner/Documents/UBC/BIOF501A/conda-env

***
<!-- EXPECTED OUTPUT -->
## Expected Output 

![Expected histogram](expected_results/covid_histogram.png "Histogram illustrating lineages among all samples")

***
<!-- TEAM -->
## Team-members

**Team Member** | **Degree** | **PI** | **Hobbies** 
------ | ---------- | -------- | ------
Jeremy Fan | Bioinformatics | Steven Jones | Weightlifting, annoying my roommate by cooking instant noodles at 3 AM 
***


