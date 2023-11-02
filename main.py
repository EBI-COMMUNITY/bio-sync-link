# This is a sample Python script.
import concurrent.futures
import csv
import json
import logging
import os
import uuid

import boto3

import requests

logging.basicConfig(format='%(asctime)s - %(threadName)s - %(message)s', level=logging.INFO)


def call_ena_api(row):
    uuid = row[22]
    cleaned_uuid = (uuid.replace('DNA-', '')
                    .replace('TIS-', '')
                    .replace('DNA Prep.', '')
                    .replace('Tissue', '')
                    .replace('DNA_Moll_', '')
                    .replace('DSM ', '').strip())
    try:
        response_json = cal_ena_api_with_uuid('sequence', cleaned_uuid)
        if len(response_json) != 0:
            logging.info(json.dumps(response_json))
            for result_row in response_json:
                result_row['source'] = 'uuid'
                result_row['api'] = 'sequence'
            return response_json
        else:
            response_json = cal_ena_api_with_uuid('sample', cleaned_uuid)
            if len(response_json) != 0:
                logging.info(json.dumps(response_json))
                for result_row in response_json:
                    result_row['source'] = 'uuid'
                    result_row['api'] = 'sample'
                return response_json
    except Exception as e:
        logging.error(f'Error while calling ENA API for {row[22]}', e)


def cal_ena_api_with_uuid(result, cleaned_uuid):
    response = requests.get(
        f'https://www.ebi.ac.uk/ena/portal/api/search?result={result}&fields=all&limit=10&format=json&'
        f'query=specimen_voucher="{cleaned_uuid}" '
        f'OR bio_material="{cleaned_uuid}" '
        f'OR culture_collection="{cleaned_uuid}" '
        f'OR isolation_source="{cleaned_uuid}"')
    response_json = json.loads(response.content)
    return response_json


def call_ena_with_accession(accession_numbers):
    try:
        response_json = call_ena_api_with_accession('sequence', accession_numbers)
        if len(response_json) != 0:
            logging.info(json.dumps(response_json))
            for result_row in response_json:
                result_row['source'] = 'accession'
                result_row['api'] = 'sequence'
            return response_json
        else:
            response_json = call_ena_api_with_accession('sample', accession_numbers)
            if len(response_json) != 0:
                logging.info(json.dumps(response_json))
                for result_row in response_json:
                    result_row['source'] = 'accession'
                    result_row['api'] = 'sample'
                return response_json
        return None
    except Exception as e:
        logging.error(f'Error while calling ENA API for {accession_numbers}', e)


def call_ena_api_with_accession(result, accession_numbers):
    query_string = f'https://www.ebi.ac.uk/ena/portal/api/search?result={result}&fields=all&limit=10&format=json&query='
    for i, accession_number in enumerate(accession_numbers):
        if i == 0:
            query_string = query_string + f'accession="{accession_number.strip()}"'
        else:
            query_string = query_string + f' OR accession="{accession_number.strip()}"'
    response = requests.get(query_string)
    response_json = json.loads(response.content)
    return response_json


def write_unmatched_accession(row, writer_invalid_accession):
    result = {
        'ggbn_guid': row[19],
        'ggbn_unitid': row[22],
        'ggbn_scietific_name': row[23],
        'ggbn_country': row[8],
        'ggbn_collection_date': row[4],
        'ggbn_collector': row[6]
    }
    writer_invalid_accession.writerow(result)


def process_row(row, writer_invalid_accession):
    accession_numbers = row[18].split('|')
    result = None
    if len(accession_numbers) > 0 and accession_numbers[0] != '\\N':
        logging.info(f'This record has accession_numbers {accession_numbers} checking accession numbers')
        result = call_ena_with_accession(accession_numbers)
        if not result and accession_numbers[0] != '\\N':
            write_unmatched_accession(row, writer_invalid_accession)
        elif len(accession_numbers) != len(result) and accession_numbers[0] != '\\N':
            for number in accession_numbers:
                if number not in result:
                    logging.info(f'Accession number {number} not found in ENA')
                    write_unmatched_accession(row, writer_invalid_accession)
    if not result:
        logging.info(f'No accession_numbers found, checking uuid {row[22]}')
        result = call_ena_api(row)
    return result


def write_result_to_file(writer, row, result):
    for result_row in result:
        result = {
            'source': result_row['source'],
            'ggbn_unitid': row[22],
            'ena_hit_on': check_hit(result_row),
            'ggbn_scietific_name': row[23],
            'ena_scientific_name': result_row['scientific_name'],
            'ggbn_country': row[8],
            'ena_country': result_row['country'],
            'ggbn_collection_date': row[4],
            'ena_collection_date': result_row['collection_date'],
            'ggbn_collector': row[6],
            'ena_collector': result_row['collected_by'],
            'ena_id': pick_correct_accession_field(result_row),
            'ena_api': result_row['api'],
            'gbbn_type': row[2],
            'ggbn_guid': row[19],
            'ggbn_id': row[0]
        }
        split_scientific_name_ena = row[23].split(' ')
        for name_part in split_scientific_name_ena:
            if name_part in result_row['scientific_name']:
                result['tax_match'] = 'True'
        if 'tax_match' not in result:
            result['tax_match'] = 'False'
        if result['source'] == 'uuid' and result['tax_match'] == 'False':
            logging.info(f'Skipping row as the tax match is false {result}')
        elif result['gbbn_type'].lower() not in ['dna', 'tissue',
                                                 'specimen',
                                                 'culture',
                                                 'environmental sample']:
            logging.info(f'Skipping row as the type is not accepted {result}')
        else:
            writer.writerow(result)


def pick_correct_accession_field(result_row):
    value = result_row.get('sequence_accession')
    if value:
        return value
    else:
        return result_row.get('sample_accession')


def check_hit(result_row):
    if result_row['specimen_voucher']:
        return result_row['specimen_voucher']
    if result_row['bio_material']:
        return result_row['bio_material']
    if result_row['culture_collection']:
        return result_row['culture_collection']
    if result_row['isolation_source']:
        return result_row['isolation_source']


def main_method():
    with open('new-full-dump.csv', 'r', encoding='ISO-8859-1') as csvfile:
        datareader = csv.reader(csvfile, delimiter='\t', quotechar='"')
        with open('output.csv', 'w', newline='') as outputfile:
            writer = csv.DictWriter(outputfile, delimiter=',', quotechar='"',
                                    fieldnames=['source', 'tax_match', 'ggbn_unitid', 'ena_hit_on',
                                                'ggbn_scietific_name',
                                                'ena_scientific_name', 'ggbn_country', 'ena_country',
                                                'ggbn_collection_date', 'ena_collection_date', 'ggbn_collector',
                                                'ena_collector', 'ena_id', 'ggbn_guid', 'ena_api', 'gbbn_type',
                                                'ggbn_guid', 'ggbn_id'])
            with open('unmatched-accession.csv', 'w', newline='') as unmatched_accession_outputfile:
                writer_invalid_accession = csv.DictWriter(unmatched_accession_outputfile, delimiter=',', quotechar='"',
                                                          fieldnames=['ggbn_guid', 'ggbn_unitid', 'ggbn_scietific_name',
                                                                      'ggbn_country', 'ggbn_collection_date',
                                                                      'ggbn_collector'])
                # with concurrent.futures.ThreadPoolExecutor(max_workers=10) as executor:
                writer.writeheader()
                for i, row in enumerate(datareader):
                    if i == 0:
                        continue
                    elif i < int(os.environ.get('RECORD_LIMIT')):
                        # executor.submit(process_row_in_ggbn, row, writer, writer_invalid_accession)
                        process_row_in_ggbn(row, writer, writer_invalid_accession)
                    else:
                        break
    s3_client = boto3.client('s3',
                             aws_access_key_id=os.environ.get('ACCESS_KEY'),
                             aws_secret_access_key=os.environ.get('SECRET_KEY'))
    s3_client.upload_file('output.csv', 'ggbn-ena-mapping', f'output-{str(uuid.uuid4())}.csv')
    s3_client.upload_file('unmatched-accession.csv', 'ggbn-ena-mapping', f'unmatched-accession-{str(uuid.uuid4())}.csv')


def process_row_in_ggbn(row, writer, writer_invalid_accession):
    result = process_row(row, writer_invalid_accession)
    if result:
        write_result_to_file(writer, row, result)


# Press the green button in the gutter to run the script.
if __name__ == '__main__':
    main_method()
