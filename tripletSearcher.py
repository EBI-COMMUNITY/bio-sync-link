import csv
import requests
import re
import logging

sah_endpoint = 'https://www.ebi.ac.uk/ena/sah/api/'
ena_endpoint = 'https://www.ebi.ac.uk/ena/portal/api/search?result=sequence&fields=all&limit=10&format=json&query=specimen_voucher='
annotation_endpoint = 'https://www.ebi.ac.uk/ena/clearinghouse/api/curations'
institution_dict = {}

logging.basicConfig(format='%(asctime)s - %(threadName)s - %(message)s', level=logging.INFO)


def search_triplets(row, writer):
    col_triplet, space_triplet = build_searchable_triplet_from_unit_id(row)
    if col_triplet:
        if validate_triplet(col_triplet):
            print("This record has a triplet: ", col_triplet)
            try:
                results = requests.get(ena_endpoint + '"' + col_triplet + '"').json()
            except Exception:
                print(f'Error while calling ENA API for {col_triplet}')
                return
            if len(results) == 1:
                return write_positive_match(results[0], row, writer)
            for result in results:
                if result.get('specimen_voucher') == col_triplet or result.get('specimen_voucher') == space_triplet:
                    return write_positive_match(result, row, writer)
        else:
            col_doublet, space_doublet = build_searchable_triplet_from_unit_id(row)
            if validate_triplet(col_doublet):
                print("This record has a triplet: ", col_doublet)
                try:
                    results = requests.get(ena_endpoint + '"' + col_doublet + '"').json()
                except Exception:
                    print(f'Error while calling ENA API for {col_doublet}')
                    return
                if len(results) == 1:
                    return write_positive_match(results[0], row, writer)
                for result in results:
                    if result.get('specimen_voucher') == col_doublet or result.get('specimen_voucher') == space_doublet:
                        return write_positive_match(result, row, writer)
        print("\tNo match found for triplet")


def write_positive_match(result, row, writer):
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
            'ena_id': result_row['sequence_accession'],
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
        if result['source'] == 'uuid' and result['tax_match'] == 'False' or result['gbbn_type'] not in ['DNA', 'Tissue',
                                                                                                        'specimen',
                                                                                                        'culture',
                                                                                                        'environmental sample']:
            logging.info(f'Skipping row as the tax match is false {result}')
        else:
            writer.writerow(result)
            return result


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
        print(f'Error while calling ENA SAH API for {triplet}')
        return
    return r["success"]


def build_searchable_triplet_from_unit_id(row):
    institution = translate_institution(row[20].replace('"', '').replace('\\N', ''))
    collection = row[21].replace('"', '').replace('\\N', '')
    unit = remove_institution_from_voucher_id(row[22].replace('"', '').replace('\\N', ''), institution)
    return (assemble_triplet(institution, collection, unit, ":"),
            assemble_triplet(institution, collection, unit, " "))


def assemble_triplet(institution, collection, unit, delimiter):
    institution = translate_institution(institution.replace('\\N', ''))
    if not institution or not unit:
        return ''

    if collection:
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
                    v = row[4]
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
    if institution != translated:
        print(institution, " -> ", translated)
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
    collections = r.get("institutions")[0].get("collections")  # this endpoint will return unique institutions
    col_list = []
    for collection in collections:
        col_list.append(collection.get("collectionCode"))
    print(institution_ena, " has collections", col_list)
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
        triplet = assemble_triplet(ena_institution, '', voucher_id, ":")
        print("triplet for this specimen should be annotated as: ", triplet)
        write_annotation_to_file(triplet, voucher_id, annotation_writer)
    elif len(col_list) == 1:
        triplet = assemble_triplet(ena_institution, col_list[0], voucher_id, ":")
        print("triplet for this specimen should be annotated as: ", triplet)
        write_annotation_to_file(triplet, voucher_id, annotation_writer)
    else:
        print("** too many collection codes!", col_list)


def write_annotation_to_file(triplet, voucher_id, annotation_writer):
    annotation_writer.writeRow(voucher_id, triplet, get_annotation_request(triplet, voucher_id))


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
    return params
    # requests.post(annotation_endpoint, json=params)


def remove_institution_from_voucher_id(voucher_id, ena_institution):
    pattern = re.escape(
        ena_institution) + r'[\s\-_]*'  # Match the substring with trailing spaces, hyphens, or underscores
    return re.sub(pattern, "", voucher_id)


def is_triplet(ena_voucher):
    pattern = r'^(\w+):(\w+)(?::(\w+))?$'  # Find items either {institution}:{id} or {institution}:{collection}:{id}
    match = re.match(pattern, ena_voucher)
    if match:
        return True
    else:
        return False
