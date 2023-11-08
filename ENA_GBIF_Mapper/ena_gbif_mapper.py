# std
from dateutil.parser import parse
import re
import requests as requests

# 3rd party
import pandas as pd

PARAMS = [
        'scientific_name',
        'specimen_voucher',
        'country',
        'collected_by',
        'collection_date'
]


def read_ena_dump_pd(ena_dump_file):
    ena_dump = pd.read_csv(ena_dump_file, sep='\t')
    return ena_dump


def process(data):
    xref_file = load_xref_file()
    found_records_list = list()
    sample_size: int = int(input("Enter sample size: "))
    if sample_size > len(data):
        print("Sample size is larger than data size. Using data size instead.")
        sample_size = len(data)
    else:
        print("Using sample size of ", sample_size)
    for index, row in data.sample(sample_size).iterrows():
        parameters = extract_parameters(row)
        try:
            found_records = request(parameters)
            for found_record in found_records['result']:
                if found_record['score'] == 1:
                    found_records_list.append({'accession': row['accession'], **found_record})
                    xref_entry = {'SOURCE_PRIMARY_ID': found_record['sourceURL'].replace(
                        'https://api.gbif.org/v1/occurrence/', ''), 'TARGET_PRIMARY_ACC': row['accession']}
                    xref_file.loc[len(xref_file)] = xref_entry
        except Exception as e:
            print(e)
            print(parameters)
    found_records_df = pd.DataFrame(found_records_list)
    found_records_df.to_csv('output/found_records.csv', sep='\t', index=False)
    xref_file.to_csv('output/xref_file.csv', sep='\t', index=False)


def request(parameters):
    url = 'https://services.bgbm.org/spase/api/specimen/guid?'
    try:
        # converting dates avoids format issues with SpASe requests
        if 'collection_date' in parameters.keys():
            if re.match(r'^[A-Z][a-z|A-Z]{2}-\d{4}$', parameters['collection_date']):
                parameters['collection_date'] = parse(parameters['collection_date']).strftime('%Y-%m')
            elif re.match(r'^\d{2}-[A-Z][a-z|A-Z]{2}-\d{4}$', parameters['collection_date']):
                parameters['collection_date'] = parse(parameters['collection_date']).strftime('%Y-%m-%d')
        print(parameters)
        r = requests.get(url, params={**{'exactMatch': 'true', 'quickSearch':'true', 'gbifOnly': 'true'},
                                      **parameters})
        return r.json()
    except Exception as e:
        print(e)
        print(parameters)


def load_xref_file():
    xref_file = pd.read_csv('input/Ggbn_ena_template.csv', sep='\t')
    return xref_file


def extract_parameters(row):
    parameters = dict()
    for col in row.index:
        if col in PARAMS and not pd.isnull(row[col]):
            if col == 'scientific_name':
                parameters['organism'] = row[col]
            else:
                parameters[col] = row[col]
    return parameters
