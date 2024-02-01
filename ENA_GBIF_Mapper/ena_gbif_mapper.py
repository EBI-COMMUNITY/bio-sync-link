# std
from dateutil.parser import parse
import re
import requests as requests
import csv
import random

# 3rd party
import pandas as pd

PARAMS = [
        'scientific_name',
        'specimen_voucher',
        'country',
        'collected_by',
        'collection_date'
]


def count_lines(file: str) -> int:
    with open(file, 'rb') as f:
        num_lines = sum(1 for _ in f) - 1
    return num_lines


def process(file_path: str, sample_size: int) -> None:
    found_records_list = list()
    number_of_lines = count_lines(file_path)
    clear_file('output/xref_file.csv', 1)
    clear_file('output/found_records.csv', 1)
    with open(file_path, 'r') as f:
        try:
            random_numbers = set(random.sample(range(2, number_of_lines + 1), sample_size))
            reader = csv.reader(f, delimiter='\t')
            header = next(reader)
            for row in reader:
                if reader.line_num in random_numbers:
                    parameters = dict()
                    parameters = extract_parameters(pd.Series(row, index=header))
                    found_records = request(parameters)
                    for found_record in found_records['result']:
                        if found_record['score'] == 1:
                            found_record = {'accession': row[0], **found_record}
                            found_records_list.append({'accession': row[0], **found_record})
                            xref_entry = {'SOURCE_PRIMARY_ID': found_record['sourceURL'].replace(
                                'https://api.gbif.org/v1/occurrence/', ''), 'SOURCE_SECONDARY_ID': '',
                                'TARGET_PRIMARY_ACC': row[0], 'TARGET_SECONDARY_ACC': ''}
                            with open('output/xref_file.csv', 'a', newline='') as xref_file:
                                writer = csv.DictWriter(xref_file, fieldnames=xref_entry.keys(), delimiter='\t')
                                writer.writerow(xref_entry)
                            with open('output/found_records.csv', 'a', newline='') as found_records_file:
                                writer = csv.DictWriter(found_records_file, fieldnames=found_record.keys(),
                                                        delimiter='\t')
                                writer.writerow(found_record)
        except Exception as e:
            print(e)
            print(parameters) if 'parameters' in locals() else print('parameters not defined')


def request(parameters: dict) -> dict:
    url = 'https://services.bgbm.org/spase/api/specimen/guid?'
    try:
        # converting dates avoids format issues with SpASe requests
        if 'collection_date' in parameters.keys():
            if re.match(r'^[A-Z][a-z|A-Z]{2}-\d{4}$', parameters['collection_date']):
                parameters['collection_date'] = parse(parameters['collection_date']).strftime('%Y-%m')
            elif re.match(r'^\d{2}-[A-Z][a-z|A-Z]{2}-\d{4}$', parameters['collection_date']):
                parameters['collection_date'] = parse(parameters['collection_date']).strftime('%Y-%m-%d')
        print(parameters)
        r = requests.get(url, params={**{'exactMatch': 'true', 'quickSearch': 'true', 'gbifOnly': 'true'},
                                      **parameters})
        return r.json()
    except Exception as e:
        print(e)
        print(parameters)


def clear_file(file_path: str, start_line: int) -> None:
    with open(file_path, 'r') as f:
        lines = f.readlines()
    with open(file_path, 'w') as f:
        f.writelines(lines[:start_line])


def extract_parameters(row: pd.Series) -> dict:
    parameters = dict()
    for col in row.index:
        if col in PARAMS and row[col] != "":
            if col == 'scientific_name':
                parameters['organism'] = row[col]
            else:
                parameters[col] = row[col]
    return parameters
