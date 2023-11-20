import concurrent.futures
import csv
import json
import logging
import os
import uuid
from typing import List, Dict, TextIO, Union

import boto3

import requests
from requests import RequestException
import tripletManager

logging.basicConfig(format='%(asctime)s - %(threadName)s - %(message)s', level=logging.INFO)


def call_ena_api_with_unitid(row: List[str]) -> Union[List[Dict[str, str]], None]:
    """
    Try to retrieve a match from ENA based on the unitid.
    First try the sequence endpoint, if no match is found try the sample endpoint.
    :param row: The GGBN record
    :return: A list of results from ENA or None if no match was found
    """
    unitid = row[22]
    cleaned_unitid = clean_uuid(unitid)
    try:
        response_json = cal_ena_api_with_unitid('sequence', cleaned_unitid)
        if len(response_json) != 0:
            logging.info(json.dumps(response_json))
            for result_row in response_json:
                result_row['source'] = 'unitid'
                result_row['api'] = 'sequence'
            return response_json
        else:
            response_json = cal_ena_api_with_unitid('sample', cleaned_unitid)
            if len(response_json) != 0:
                logging.info(json.dumps(response_json))
                for result_row in response_json:
                    result_row['source'] = 'unitid'
                    result_row['api'] = 'sample'
                return response_json
    except RequestException:
        logging.error(f'Error while calling ENA API for {unitid}')
    return None


def clean_uuid(unitid: str) -> str:
    """
    Removes any prefixes from the unitid.
    The cleaned unitid will be used to call the ENA API.
    :param unitid: unitid from the GGBN record
    :return: Cleaned unitid
    """
    return (unitid.replace('DNA-', '')
            .replace('TIS-', '')
            .replace('DNA Prep.', '')
            .replace('Tissue', '')
            .replace('DNA_Moll_', '')
            .replace('DSM ', '').strip())


def cal_ena_api_with_unitid(result: str, cleaned_unitid: str) -> List[Dict[str, str]]:
    """
    Create the query string based on the cleaned unitid and make the API call.
    :param result: Which of the ENA endpoints need to be called (sequence or sample)
    :param cleaned_unitid: The cleaned unit id from the GGBN record
    :return: A list of results from the ENA API
    """
    response = requests.get(
        f'https://www.ebi.ac.uk/ena/portal/api/search?result={result}&fields=all&limit=10&format=json&'
        f'query=specimen_voucher="{cleaned_unitid}" '
        f'OR bio_material="{cleaned_unitid}" '
        f'OR culture_collection="{cleaned_unitid}" '
        f'OR isolation_source="{cleaned_unitid}"')
    response_json = json.loads(response.content)
    return response_json


def call_ena_with_accession(accession_numbers: List[str]) -> Union[List[Dict[str, str]], None]:
    """
    Try to retrieve a match from ENA based on accession numbers.
    First check the sequence endpoint, if no match is found check the sample endpoint.
    On exception the exception will be logged and a None will be thrown
    :param accession_numbers: List of the accession numbers
    :return: A list of results from ENA or None if no match was found
    """
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
    except RequestException:
        logging.error(f'Error while calling ENA API for {accession_numbers}')
        return None


def call_ena_api_with_accession(result: str, accession_numbers: List[str]) -> List[Dict[str, str]]:
    """
    Create the query string based on the accession numbers and make the API call.
    :param result: Which of the ENA endpoints need to be called (sequence or sample)
    :param accession_numbers: List of accession numbers for the GGBN record
    :return: Returns the result from the ENA API
    """
    query_string = f'https://www.ebi.ac.uk/ena/portal/api/search?result={result}&fields=all&limit=10&format=json&query='
    for i, accession_number in enumerate(accession_numbers):
        if i == 0:
            query_string = query_string + f'accession="{accession_number.strip()}"'
        else:
            query_string = query_string + f' OR accession="{accession_number.strip()}"'
    response = requests.get(query_string)
    response_json = json.loads(response.content)
    return response_json


def write_unmatched_accession(row: List[str], writer_invalid_accession: csv.DictWriter, accession_number: str) -> None:
    """
    Write the non-matching accession numbers to file.
    :param row: The GGBN record
    :param writer_invalid_accession: Writer for the non-matching accession numbers
    :param accession_number: The accession number that did not match
    :return:
    """
    result = {
        'accession_number': accession_number,
        'ggbn_guid': row[19],
        'ggbn_unitid': row[22],
        'ggbn_scietific_name': row[24],
        'ggbn_country': row[8],
        'ggbn_collection_date': row[5],
        'ggbn_collector': row[6],
        'ggbn_full_id': row[32]
    }
    writer_invalid_accession.writerow(result)


def process_row(row: List[str], writer_invalid_accession: csv.DictWriter) -> Union[List[Dict[str, str]], None]:
    """
    Try to find a match for a given GGBN record.
    We first check if there are accession_numbers present in the record.
    If so we will use this numbers to call the ENA API with the accession numbers.
    If there is a match this will be returned, otherwise the non-matching accession numbers are writer to file.
    If no accession_numbers are present or if they did not match we will try to match on uuid.
    :param row: GGBN record row
    :param writer_invalid_accession: Writer for the non-matching accession numbers
    :return: A List of Dicts with the match from ENA. There could be a match to multiple results therefore it is a list.
    """
    accession_numbers = row[18].split('|')
    result = None
    if len(accession_numbers) > 0 and accession_numbers[0] != '\\N':
        logging.info(f'This record has accession_numbers {accession_numbers} checking accession numbers')
        result = call_ena_with_accession(accession_numbers)
        if not result:
            for accession_number in accession_numbers:
                write_unmatched_accession(row, writer_invalid_accession, accession_number)
        if result:
            found_accessions = [pick_correct_accession_field(item) for item in result]
            not_found_accessions = list(
                filter(lambda accession: accession.strip() not in found_accessions, accession_numbers))
            for accession_number in not_found_accessions:
                write_unmatched_accession(row, writer_invalid_accession, accession_number)
    if not result:
        logging.info(f'No accession_numbers found, checking unitid {row[22]}')
        result = call_ena_api_with_unitid(row)
    return result


def write_result_to_file(writer: csv.DictWriter, row: List[str], result: List[Dict[str, str]]) -> bool:
    """
    Map the result from ENA to the GGBN record to a result record.
    Check if the taxonomy matches.
    If the match is based on an unitid and the taxonomy does not match ignore the result.
    If the GGBN type for the match is not accepted ignore the result.
    :param writer: The writer for the matches file
    :param row: GGBN record
    :param result: ENA API result
    :return: Will not return anything but write the result to file
    """
    for result_row in result:
        mapped_result = map_match_to_result(result_row, row)
        split_scientific_name_ena = row[24].split(' ')
        for name_part in split_scientific_name_ena:
            if name_part in result_row['scientific_name']:
                mapped_result['tax_match'] = 'True'
        if 'tax_match' not in mapped_result:
            mapped_result['tax_match'] = 'False'
        if mapped_result['source'] == 'unitid' and mapped_result['tax_match'] == 'False':
            logging.info(f'Skipping row as the tax match is false {mapped_result}')
            return False
        elif mapped_result['gbbn_type'].lower() not in ['dna', 'tissue',
                                                        'specimen',
                                                        'culture',
                                                        'environmental sample']:
            logging.info(f'Skipping row as the type is not accepted {mapped_result}')
            return False
        else:
            writer.writerow(mapped_result)
            return True


def map_match_to_result(result_row: Dict[str, str], row: List[str]) -> Dict[str, str]:
    """
    Map the information from GGBN and ENA to a result record.
    :param result_row: ENA result
    :param row: GGBN record
    :return: Formatted Dict that contains all necessary information
    """
    result = {
        'source': result_row['source'],
        'ggbn_unitid': row[22],
        'ena_hit_on': check_hit(result_row),
        'ggbn_scietific_name': row[24],
        'ena_scientific_name': result_row['scientific_name'],
        'ggbn_country': row[8],
        'ena_country': result_row['country'],
        'ggbn_collection_date': row[5],
        'ena_collection_date': result_row['collection_date'],
        'ggbn_collector': row[6],
        'ena_collector': result_row['collected_by'],
        'ena_id': pick_correct_accession_field(result_row),
        'ena_api': result_row['api'],
        'gbbn_type': row[2],
        'ggbn_guid': row[19],
        'ggbn_short_id': row[23],
        'ggbn_full_id': row[32]
    }
    return result


def pick_correct_accession_field(result_row: Dict[str, str]):
    """
    Retrieve the accession field from the ENA API result.
    :param result_row: ENA API result
    :return: accession number
    """
    return result_row.get('accession')


def check_hit(result_row: Dict[str, str]) -> Union[str, None]:
    """
    Collect the identifier on which the ENA API produced a hit.
    :param result_row: ENA API result
    :return: The identifier on which the ENA API produced a hit.
    """
    if result_row['specimen_voucher']:
        return result_row['specimen_voucher']
    if result_row['bio_material']:
        return result_row['bio_material']
    if result_row['culture_collection']:
        return result_row['culture_collection']
    if result_row['isolation_source']:
        return result_row['isolation_source']
    return None


def main_method() -> None:
    """
    Main method that will read the CSV file and process each row.
    Each matching result will be writing to the output.csv file.
    Any accession which wasn't found in ENA is writen to unmatched-accession.csv.
    When finished the output.csv and unmatched-accession.csv will be uploaded to S3.
    :return: Will not produce a result, results can be found in the S3 bucket
    """
    with (open('new-full-dump.csv', 'r', encoding='ISO-8859-1') as csvfile):
        datareader = csv.reader(csvfile, delimiter='\t', quotechar='"')
        with open('output.csv', 'w', newline='') as outputfile, open('annotationOutput.csv', 'w',
                                                                     newline='') as annotationoutfile:
            writer = init_match_writer(outputfile)
            annotation_writer = csv.writer(annotationoutfile)
            with open('unmatched-accession.csv', 'w', newline='') as unmatched_accession_outputfile:
                writer_invalid_accession = init_invalid_accession_writer(unmatched_accession_outputfile)
                with concurrent.futures.ThreadPoolExecutor(max_workers=10) as executor:
                    writer.writeheader()
                    writer_invalid_accession.writeheader()
                    annotation_writer.writerow(['ena_id', 'triplet', 'annotation'])
                    for i, row in enumerate(datareader):
                        if i == 0:
                            continue
                        elif i < int(os.environ.get('RECORD_LIMIT')):
                            executor.submit(process_row_in_ggbn, row, writer, writer_invalid_accession,
                                            annotation_writer)
                            # process_row_in_ggbn(row, writer, writer_invalid_accession, annotation_writer)
                        else:
                            break
    write_result_to_s3()


def write_result_to_s3() -> None:
    """
    Will write both result files to S3 buckets.
    Secrets need to be provided through environmental variables
    :return:
    """
    s3_client = boto3.client('s3',
                             aws_access_key_id=os.environ.get('ACCESS_KEY'),
                             aws_secret_access_key=os.environ.get('SECRET_KEY'))
    s3_client.upload_file('output.csv', 'ggbn-ena-mapping', f'output-{str(uuid.uuid4())}.csv')
    s3_client.upload_file('unmatched-accession.csv', 'ggbn-ena-mapping', f'unmatched-accession-{str(uuid.uuid4())}.csv')


def init_invalid_accession_writer(unmatched_accession_outputfile: TextIO) -> csv.DictWriter:
    """
    Initialize the writer for the unmatched accessions output file
    :param unmatched_accession_outputfile: File for the unmatched accessions
    :return: Returns the writer
    """
    return csv.DictWriter(unmatched_accession_outputfile, delimiter=',', quotechar='"',
                          fieldnames=['accession_number',
                                      'ggbn_guid',
                                      'ggbn_unitid',
                                      'ggbn_scietific_name',
                                      'ggbn_country',
                                      'ggbn_collection_date',
                                      'ggbn_collector',
                                      'ggbn_full_id'])


def init_match_writer(outputfile: TextIO) -> csv.DictWriter:
    """
    Initialize the writer for the matches output file
    :param outputfile: File for the matches file
    :return: Returns the writer
    """
    return csv.DictWriter(outputfile, delimiter=',', quotechar='"',
                          fieldnames=['source', 'tax_match', 'ggbn_unitid', 'ena_hit_on',
                                      'ggbn_scietific_name',
                                      'ena_scientific_name', 'ggbn_country', 'ena_country',
                                      'ggbn_collection_date', 'ena_collection_date', 'ggbn_collector',
                                      'ena_collector', 'ena_id', 'ggbn_guid', 'ena_api', 'gbbn_type',
                                      'ggbn_guid', 'ggbn_short_id', 'ggbn_full_id'])


def process_row_in_ggbn(row: List[str], writer: csv.DictWriter, writer_invalid_accession: csv.DictWriter,
                        annotation_writer: csv.writer) -> None:
    """
    Starts processing a single row from the GGBN dump
    If there is a result write it to file, otherwise try to make a triplet and call the ENA API with the triplet.
    If the triplet call produces a result add it as annotation to the triplet file.
    :param row: Single row from the GGBN dump
    :param writer: Writer to the matches file
    :param writer_invalid_accession: Writer to the unmatched accessions file
    :param annotation_writer: Writer to annotation file
    :return: No return value
    """
    result = process_row(row, writer_invalid_accession)
    r = False
    if result:
        r = write_result_to_file(writer, row, result)
    # if not r:
    #     r = tripletManager.call_ena_api_triplets(row, writer)
    # if r:
    #     tripletManager.annotate_triplet(result, row, annotation_writer)


if __name__ == '__main__':
    main_method()
