# Mapping between ENA specimen information and GBIF records

## Introduction
In working group 23 of the Elixir BioHackathon, the BioSyncLink project aims at developing tools that facilitate the linking between sequence data and specimens/sample records. 
Using, the specimen source information in sequence records, it may be possible to search and link to the public collection records describing the specimens or samples from which the sequences were derived.
During the hackathon, a script was implemented that maps GGBN records to ENA accessions and vice versa (bidirectional linking). 
Further, another Python script was implemented to demonstrate the mapping of ENA accessions to GBIF occurrence records by using the specimen search API SpASe which was developed within the BiCIKL project. This script is described in this document.



## What does this script do?
This section describes the functionality of the script. Starting with the input data, the output data and the SpASe API, the workflow of the script is described.
### Input
Since this work revolves around the mapping of ENA accessions to GBIF occurrence records, the input data is a TSV file containing the ENA accessions and the corresponding specimen information. The file was obtained by [querying](https://www.ebi.ac.uk/ena/portal/api/search?result=sequence&query=specimen_voucher=%22*%22&fields=accession,scientific_name,specimen_voucher,country,collected_by,collection_date&format=tsv&limit=1000000) the ENA API for the fields accession, scientific_name, specimen_voucher, country, collected_by and collection_date. Further, as precondition for the mapping of sequences derived from specimens, only records with a specimen_voucher were considered. As this script is a proof of concept, the input data is limited to 1 million records.

### Output
Two CSV files are generated as output. The file **xref_file.csv** contains the mapping between ENA accessions and GBIF occurrence records. The file **found_records.csv** contains the specimen information found in the occurrence records for those accessions that could be mapped. Hence, each row contains the corresponding accession number.

### SpASe
[SpASe](https://services.bgbm.org/spase/api/specimens/file) is a web service developed in the course of the BiCIKL project and consists of an API and a GUI. It was designed to search specimen information systems for a set of parameter values. Further, it evaluates and assesses the similarity of found records compared to the search parameter values and provides it in the form of a similarity score. Besides DiSSCo and iDigBio, it queries GBIF and can thus be used to map ENA accessions to GBIF occurrence records. The web service is described in detail [here](https://docs.google.com/document/d/1vb9JlHm4DK-Z9V6BxnYC2mnsQaqWJLNR/edit?usp=sharing&ouid=112373834655561411383&rtpof=true&sd=true), the source code can be found [here](https://git.bgbm.org/bicikl/wp7-web-service). The allowed input parameters follow the source feature qualifiers defined in the [INSDC Feature Table Specification](https://www.insdc.org/documents/feature-table), which makes SpASE a suitable tool for this task. To perform the mapping, SpASe searches only GBIF, excludes INSDC and BOLD dataset records and returns only records with a similarity score of 1.

## How to run this script?
**When using PyCharm:**
1. Open the project in PyCharm.
2. Add a new interpreter (File > Settings > Project: bio-sync-link > Project Interpreter > Add) and select the Python 3.11 interpreter.
3. Add a new run configuration (Run > Edit Configurations > Add New Configuration > Python) and select the previously added interpreter.
4. Set the script path to the run.py file.
5. Set the working directory to ENA_GBIF_Mapper.
6. Open requirements.txt and install the required packages.
7. Run the script. 

**When using the command line:**
1. Open the command line.
2. Navigate to the ENA_GBIF_Mapper folder.
3. Create a virtual environment: `python -m venv venv`
4. Activate the virtual environment: `venv\Scripts\activate.bat`
5. Install the required packages: `pip install -r requirements.txt`
6. Run the script: `python run.py`.

In either case, the results will be written into the corresponding CSV files in the outputs folder.
