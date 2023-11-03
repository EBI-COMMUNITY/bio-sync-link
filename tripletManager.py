import csv

import requests
import re
import logging

sah_endpoint = 'https://www.ebi.ac.uk/ena/sah/api/'
ena_endpoint = 'https://www.ebi.ac.uk/ena/portal/api/search?result=sequence&fields=all&limit=10&format=json&query=specimen_voucher='
annotation_endpoint = 'https://www.ebi.ac.uk/ena/clearinghouse/api/curations'
institution_dict = {}

logging.basicConfig(format='%(asctime)s - %(threadName)s - %(message)s', level=logging.INFO)
send_to_api = False

def call_ena_api_triplets(row, writer):
    endpoint = 'sequence'
    results = search_triplets(row, endpoint)
    if not results:
        endpoint = 'sample'
        results = search_triplets(row, endpoint)
    if results:
        return write_positive_match(results, row, writer, endpoint)


def search_triplets(row, endpoint):
    col_triplet, space_triplet = build_searchable_triplet(row, True)
    if col_triplet and validate_triplet(col_triplet):
        logging.info("No match with accession/uuid found. checking triplet " + col_triplet)
        results = None
        try:
            results = requests.get(build_ena_request(col_triplet, space_triplet, endpoint)).json()
        except Exception:
            logging.error(f'Error while calling ENA API for triplet {col_triplet}')
            return
        if results and len(results) > 0:
            return results
    else:
        col_doublet, space_doublet = build_searchable_triplet(row, False)
        if col_doublet and validate_triplet(col_doublet):
            logging.info("No match with accession/uuid found. checking triplet " + col_doublet)
            results = None
            try:
                results = requests.get(build_ena_request(col_doublet, space_doublet, endpoint)).json()
            except Exception:
                logging.error(f'Error while calling ENA API for {col_doublet}')
                return
            if results and len(results) > 0:
                return results
            else:
                logging.info("No match found for dwc triplet " + col_triplet + " or " + col_doublet)


def build_ena_request(col_triplet, space_triplet, endpoint):
    return (f'https://www.ebi.ac.uk/ena/portal/api/search?result={endpoint}&fields=all&limit=10&format=json&'
            f'query=specimen_voucher="{col_triplet}" '
            f'OR bio_material="{col_triplet}" '
            f'OR culture_collection="{col_triplet}" '
            f'OR isolation_source="{col_triplet}" '
            f'OR specimen_voucher="{space_triplet}" '
            f'OR bio_material="{space_triplet}" '
            f'OR culture_collection="{space_triplet}" '
            f'OR isolation_source="{space_triplet}"')


def write_positive_match(result, row, writer, endpoint):
    for result_row in result:
        result = {
            'source': 'triplet',
            'ggbn_unitid': row[22],
            'ena_hit_on': check_hit(result_row),
            'ena_specimen_voucher': result_row['specimen_voucher'],
            'ggbn_scientific_name': row[23],
            'ena_scientific_name': result_row['scientific_name'],
            'ggbn_country': row[8],
            'ena_country': result_row['country'],
            'ggbn_collection_date': row[4],
            'ena_collection_date': result_row['collection_date'],
            'ggbn_collector': row[6],
            'ena_collector': result_row['collected_by'],
            'ena_id': pick_correct_accession_field(result_row),
            'ena_api': endpoint,
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
        if result['source'] == 'uuid' and result['tax_match'] == 'False' or result['gbbn_type'] not in ['DNA', 'Tissue',
                                                                                                        'specimen',
                                                                                                        'culture',
                                                                                                        'environmental sample']:
            logging.info(f'Skipping row as the tax match is false {result}')
        else:
            logging.info("A match has been found with DWC triplet" + str(result))
            writer.writerow(result)
            return result


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


def validate_triplet(triplet):
    params = {
        "value": triplet,
        "qualifier_type": "specimen_voucher"
    }
    try:
        r = requests.get(sah_endpoint + 'validate', params).json()
    except Exception:
        logging.info(f'Error while calling ENA SAH API for {triplet}')
        return
    return r["success"]


def build_searchable_triplet(row, include_collection):
    institution = translate_institution(row[20].replace('"', '').replace('\\N', ''))
    if not institution:
        return '', ''
    collection = False
    if include_collection:
        collections = get_collection_codes(institution)
        if len(collections) == 1:
            collection = collections[0]
    if not collection:
        collection = row[21].replace('"', '').replace('\\N', '')
    unit = remove_institution_from_voucher_id(row[22].replace('"', '').replace('\\N', ''), institution)
    return (assemble_triplet(institution, collection, unit, ":", include_collection),
            assemble_triplet(institution, collection, unit, " ", include_collection))


def assemble_triplet(institution, collection, unit, delimiter, include_collection):
    institution = translate_institution(institution.replace('\\N', ''))
    if not institution or not unit:
        return ''

    if collection and include_collection:
        return institution + delimiter + collection + delimiter + unit
    return institution + delimiter + unit


def set_institution_dict():
    global institution_dict
    if len(institution_dict) == 0:
        with open("institutions.csv", 'r') as file:
            reader = csv.reader(file)
            first = True
            for row in reader:
                if not first:
                    k = row[0]
                    v = row[6]
                    if k not in institution_dict or (k in institution_dict and institution_dict[k] == '#N/A'):
                        institution_dict[k] = v
                else:
                    first = False


def translate_institution(institution):
    set_institution_dict()
    global institution_dict
    translated = institution_dict.get(institution)
    if institution not in institution_dict or translated == '#N/A':
        return ''
    return institution_dict.get(institution)


def read_csv():
    with open('dump_100k.csv', 'r', encoding='ISO-8859-1') as f:
        reader = csv.reader(f, delimiter=';')
        for row in reader:
            yield row


def get_collection_codes(institution_ena):
    try:
        r = requests.get(sah_endpoint + 'institution/' + institution_ena + '/collection?').json()
    except Exception:
        print(f'Error while calling ENA SAH API for {institution_ena}')
        return []
    if not r or not r.get("institutions"):
        print(r, " ", institution_ena)
    collections = r.get("institutions")[0].get("collections")  # this endpoint will return unique institutions
    col_list = []
    for collection in collections:
        col_list.append(collection.get("collectionCode"))
    logging.debug(institution_ena, " has collections", col_list)
    return col_list


def annotate_triplet(result, ggbn_row, annotation_writer):
    global institution_dict
    set_institution_dict()
    voucher_id = result.get('ena_specimen_voucher')
    if is_triplet(voucher_id):
        return
    ena_institution = institution_dict.get(ggbn_row[20].replace('"', '').replace('\\N', ''))
    voucher_id = remove_institution_from_voucher_id(voucher_id, ena_institution)
    col_list = get_collection_codes(ena_institution)
    if not col_list:
        triplet = assemble_triplet(ena_institution, '', voucher_id, ":", False)
        logging.info("Voucher ID for this specimen should be annotated as triplet: " + triplet)
        request = write_annotation_to_file(triplet, voucher_id, annotation_writer)
        if send_to_api:
            requests.post(annotation_endpoint, json=request)
    elif len(col_list) == 1:
        triplet = assemble_triplet(ena_institution, col_list[0], voucher_id, ":", True)
        logging.info("Voucher ID for this specimen should be annotated as triplet: ", triplet)
        request = write_annotation_to_file(triplet, voucher_id, annotation_writer)
        if send_to_api:
            requests.post(annotation_endpoint, json=request)

    else:
        logging.warning("Too many collection codes to construct triplet" + str(col_list))
        annotation_writer.writerow([voucher_id, "",
                                    f'"Unable to construct triplet. Institution {ena_institution} has too many collections {col_list}'])


def write_annotation_to_file(triplet, voucher_id, annotation_writer):
    request = get_annotation_request(triplet, voucher_id)
    annotation_writer.writerow([voucher_id, triplet, request])
    return request


def get_annotation_request(triplet, voucher_id):
    params = {
        "recordType": "sequence",
        "record_id": voucher_id,
        "attributePost": voucher_id,
        "valuePost": triplet,
        "assertionMethod": "automaticAssertion",
        "assertionEvidence": [{
            "label": "inference based on the construction of a DWC triplet"
        }],
        "providerSource": "https://github.com/EBI-COMMUNITY/bio-sync-link",
        "providerName": "BioSyncLink",
        "assertionAdditionalInfo": "This assertion was made by the Bio-Sync-Link project, Elixir Biohackathon 2023",
    }
    return str(params)


def remove_institution_from_voucher_id(voucher_id, ena_institution):
    pattern = re.escape(
        ena_institution) + r'[\s\-_]*'  # Match the substring with trailing spaces, hyphens, or underscores
    return re.sub(pattern, "", voucher_id)


def is_triplet(ena_voucher):
    pattern = r'^(\w+):(\w+)(?::(\w+))?$'  # Find items either {institution}:{id} or {institution}:{collection}:{id}
    if not ena_voucher:
        print("")
    match = re.match(pattern, ena_voucher)
    if match:
        return True
    else:
        return False


def triplet_workflow():
    with open('new-full-dump.csv', 'r', encoding='ISO-8859-1') as csvfile:
        datareader = csv.reader(csvfile, delimiter='\t', quotechar='"')
        with open("tripletOutput.csv", 'w', newline='') as file, open('annotationOutput.csv', 'w',
                                                                      newline='') as annotationoutfile:
            annotation_writer = csv.writer(annotationoutfile)
            writer = csv.DictWriter(file, delimiter=',', quotechar='"',
                                    fieldnames=['source', 'tax_match', 'ggbn_unitid', 'ena_hit_on',
                                                'ena_specimen_voucher',
                                                'ggbn_scientific_name',
                                                'ena_scientific_name', 'ggbn_country', 'ena_country',
                                                'ggbn_collection_date', 'ena_collection_date', 'ggbn_collector',
                                                'ena_collector', 'ena_id', 'ggbn_guid', 'ena_api', 'gbbn_type',
                                                'ggbn_guid', 'ggbn_id'])
            writer.writeheader()
            for i, row in enumerate(datareader):
                if i == 0:
                    continue
                elif i < 1000:
                    triplet_match = call_ena_api_triplets(row, writer)
                    if triplet_match:
                        annotate_triplet(triplet_match, row, annotation_writer)
                else:
                    break


if __name__ == '__main__':
    triplet_workflow()
