# LSLOD Query
This repository contains the set of scripts to analyze the LSLOD Cloud -- understand the characteristics of data and knowledge stored in different LSLOD sources, as well as to understand and quantify the semantic heterogeneity in the LSLOD Sources (excluding the Biomedical Ontologies from BioPortal -- please check the OntoReuse repository for the corresponding scripts related to BioPortal). 

Scripts for the following purposes are present here:
 - Automate querying of the LSLOD cloud, generate schemas, and instance characteristics in a data-driven manner.
 - Increase in assertion count with respect to the addition of subsequent sources to querying scheme.
 - Analyze the semantic heterogeneity of the LSLOD Cloud
 - Generate a Roadmap visualization of the LSLOD Cloud (similar to Google Maps view)

This repository is a mix of Python Jupyter notebooks, Python scripts, R scripts, and JavaScript scripts.

---
### To setup the virtual environment

The scripts require Python 3.6 and above. Execute the following commands:
```
virtualenv venv
source venv/bin/activate
pip install -r requirements.txt
```

---
### Generate the LSLOD network

Execute the scripts for one SPARQL endpoint
```
cd Extractors/
python schema_extractor.py --qe http://sparql.wikipathways.org/sparql --name WikiPathways
('Executing for', 'http://sparql.wikipathways.org/sparql', 'WikiPathways')
python instance_ret_main.py
```



