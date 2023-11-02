# This is a sample Python script.
import csv
import json

import requests
import tripletSearcher


def call_ena_api(row):
    uuid = row[22]
    cleaned_uuid = (uuid.replace('DNA-', '')
                    .replace('TIS-', '')
                    .replace('DNA Prep.', '')
                    .replace('Tissue', '')
                    .replace('DNA_Moll_', '')
                    .replace('DSM ', '').strip())
    try:
        response = requests.get(
            f'https://www.ebi.ac.uk/ena/portal/api/search?result=sequence&fields=all&limit=10&format=json&query=specimen_voucher="${cleaned_uuid}"')
        response_json = json.loads(response.content)
        if len(response_json) != 0:
            print(json.dumps(response_json))
            return response_json
        else:
            print(f'No record found for {row[22]}')
    except Exception as e:
        print(f'Error while calling ENA API for {row[22]}')


def process_row(row):
    print(f'This record has a UnitID: {row[22]}')
    return call_ena_api(row)


def write_result_to_file(writer, row, result):
    for result_row in result:
        if row[23] == result_row['scientific_name']:
            result = {
                'ggbn_unitid': row[22],
                'ena_specimen_voucher': result_row['specimen_voucher'],
                'ggbn_scietific_name': row[23],
                'ena_scientific_name': result_row['scientific_name'],
                'ggbn_country': row[8],
                'ena_country': result_row['country'],
                'ggbn_collection_date': row[4],
                'ena_collection_date': result_row['collection_date'],
                'ggbn_collector': row[6],
                'ena_collector': result_row['collected_by']
            }
            writer.writerow(result)
            return result


def main_method():
    with open('dump_100k.csv', 'r', encoding='ISO-8859-1') as csvfile:
        datareader = csv.reader(csvfile, delimiter=';', quotechar='"')
        with open('0.4-switch-to-comma.csv', 'w', newline='') as outputfile:
            writer = csv.DictWriter(outputfile, delimiter=',', quotechar='"',
                                    fieldnames=['ggbn_unitid', 'ena_specimen_voucher', 'ggbn_scietific_name',
                                                'ena_scientific_name', 'ggbn_country', 'ena_country',
                                                'ggbn_collection_date', 'ena_collection_date', 'ggbn_collector',
                                                'ena_collector'])
            writer.writeheader()
            for i, row in enumerate(datareader):
                if i == 0:
                    continue
                elif i < 10000:
                    result = process_row(row)
                    if result:
                        r = write_result_to_file(writer, row, result)
                    else:
                        r = tripletSearcher.search_triplets(row, writer)
                    if r:
                        tripletSearcher.annotate_triplet(r, row)
                else:
                    break


if __name__ == '__main__':
    main_method()
