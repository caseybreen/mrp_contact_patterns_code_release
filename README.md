## Estimating Subnational Contact Patterns using Multilevel Regression with Poststratification

This repository contains code and materials to replicate ["Estimating Subnational Contact Patterns using Multilevel Regression with Poststratification."](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1010742)

### Replication Package

The this repository includes code to replicate all figures and tables in the paper. There are three steps to running the replication code: 

1. Clone this repository
2. Download the `data.zip` file from the [accompanying OSF project](https://osf.io/aecwn/), unzip the data file, and move it to the top level of the repository. 
3. Run the `00_run_all.Rmd` script, which will run all code (or run scripts `01` to `05` individually)


#### Data 

Please download all data for replication from the project's [OSF project](https://osf.io/aecwn/). The data were originally obtained from: 

- ACS 2014-2018 5-year, csv and .xml [[link](https://usa.ipums.org/usa/)]
- Google Mobility Reports [[link](https://www.google.com/covid19/mobility/)]

#### Code 

After downloading the required data and moving into the top level of the replication repository, researchers can run the following script to replicate all figures and tables: 

- `00_run_all.Rmd` - this file runs all scripts. 

Alternatively, researchers can run the following files individually in order: 

- `01_prep_google_pca.Rmd` - this file prepares the mobility data from Google Mobility Reports
- `02_prep_acs_targets.Rmd` - this file calculates post-stratification targets from the 5-year American Community Survey 
- `03_run_models.Rmd` - this file runs the two Bayesian hierarchical models 
- `04_predict_and_poststratify.Rmd` - this file performs prediction from the posterior predictive distribution and then poststratifies to produce daily estimates of contact patterns 
- `05_create_figures.Rmd` - this script constructs all figures and tables for paper 
- `06_model_selection_loo.Rmd` - this script runs a few models with different specifications and performs model selection  

### Authors

- [Casey Breen](caseybreen.com)
- [Ayesha Mahmud](https://ayeshamahmud.github.io/)
- [Dennis Feehan](https://dennisfeehan.org/)

