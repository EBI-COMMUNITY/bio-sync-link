import csv
import requests
import re

sah_endpoint = 'https://www.ebi.ac.uk/ena/sah/api/'
ena_endpoint = 'https://www.ebi.ac.uk/ena/portal/api/search?result=sequence&fields=all&limit=10&format=json&query=specimen_voucher='
pattern = r'^(\w+):(\w+)(?::(\w+))?$'
institution_dict = {}


def search_triplets(row, writer):
    col_triplet, space_triplet = build_searchable_triplet(row, True)
    if col_triplet:
        if validate_triplet(col_triplet, row[22]):
            print("This record has a triplet: ", col_triplet)
            try:
                results = requests.get(ena_endpoint + '"' + col_triplet + '"').json()
            except Exception:
                print("Unable to send request")
                return
            if len(results) == 1:
                write_positive_match(results[0], row, writer)
                return
            for result in results:
                if result.get('specimen_voucher') == col_triplet or result.get('specimen_voucher') == space_triplet:
                    write_positive_match(result, row, writer)
                    return
        else:
            col_doublet, space_doublet = build_searchable_triplet(row, False)
            if validate_triplet(col_doublet, row[22]):
                print("This record has a triplet: ", col_doublet)
                try:
                    results = requests.get(ena_endpoint + '"' + col_doublet + '"').json()
                except Exception:
                    print("Unable to send request")
                    return
                if len(results) == 1:
                    write_positive_match(results[0], row, writer)
                    return
                for result in results:
                    if result.get('specimen_voucher') == col_doublet or result.get('specimen_voucher') == space_doublet:
                        write_positive_match(result, row, writer)
                        return
        print("\tNo match found for triplet")


def write_positive_match(result, row, writer):
    if row[23] == result.get('scientific_name'):
        write_row = {
            'ggbn_unitid': row[22],
            'ena_specimen_voucher': result.get('specimen_voucher'),
            'ggbn_scietific_name': row[23],
            'ena_scientific_name': result.get('scientific_name'),
            'ggbn_country': row[8],
            'ena_country': result.get('country'),
            'ggbn_collection_date': row[4],
            'ena_collection_date': result.get('collection_date'),
            'ggbn_collector': row[6],
            'ena_collector': result.get('collected_by')
        }
        print("result found: ", write_row)
        writer.writerow(write_row)


def validate_triplet(triplet, ggbn_id):
    params = {
        "value": triplet,
        "qualifier_type": "specimen_voucher"
    }
    try:
        r = requests.get(sah_endpoint + 'validate', params).json()
    except Exception:
        print("Unable to send request")
        return
    return r["success"]


def build_searchable_triplet(row, include_collection):
    institution = translate_institution(row[20].replace('"', '').replace('\\N', ''))
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
        with open("../institutions.csv", 'r') as file:
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
    with open('../ena-caller/dump_100k.csv', 'r', encoding='ISO-8859-1') as f:
        reader = csv.reader(f, delimiter=';')
        for row in reader:
            yield row


def get_collection_codes(institution_ena):
    try:
        r = requests.get(sah_endpoint + 'institution/' + institution_ena + '/collection?').json()
    except Exception:
        print("Unable to send request")
        return []
    collections = r.get("institutions")[0].get("collections")  # this endpoint will return unique institutions
    col_list = []
    for collection in collections:
        col_list.append(collection.get("collectionCode"))
    print(institution_ena, " has collections", col_list)
    return col_list


def annotate_triplet(result, ggbn_row):
    set_institution_dict()
    global institution_dict
    voucher_id = result.get('ena_specimen_voucher')
    if is_triplet(voucher_id):
        return
    ena_institution = institution_dict.get(ggbn_row[20].replace('"', '').replace('\\N', ''))
    voucher_id = remove_institution_from_voucher_id(voucher_id, ena_institution)
    col_list = get_collection_codes(ena_institution)
    if not col_list:
        triplet = assemble_triplet(ena_institution, '', voucher_id, ":", False)
        print("triplet for this specimen should be annotated as: ", triplet)
    elif len(col_list) == 1:
        triplet = assemble_triplet(ena_institution, col_list[0], voucher_id, ":", True)
        print("triplet for this specimen should be annotated as: ", triplet)
    else:
        print("** too many collection codes!")


def send_annotation_request(triplet, voucher_id):
    params = {
        "recordType": "sequence",
        "record_id": voucher_id,
        "attributePost": voucher_id,
        "valuePost": triplet,
        "assertionMethod": "automaticAssertion",
        "assertionEvidence": "identifier",  # ?
        "assertionSource": "https://github.com/EBI-COMMUNITY/bio-sync-link",
        "assertionAdditionalInfo": "This assertion was made by the Bio-Sync-Link project",

    }


def remove_institution_from_voucher_id(voucher_id, ena_institution):
    pattern = re.escape(
        ena_institution) + r'[\s\-_]*'  # Match the substring with trailing spaces, hyphens, or underscores
    return re.sub(pattern, "", voucher_id)


def is_triplet(ena_voucher):
    match = re.match(pattern, ena_voucher)
    if match:
        return True
    else:
        return False


def main():
    with open('triplet_matches.csv', 'w', newline='') as outputfile:
        writer = csv.DictWriter(outputfile, delimiter=',', quotechar='"',
                                fieldnames=['ggbn_unitid', 'ena_specimen_voucher', 'ggbn_scietific_name',
                                            'ena_scientific_name', 'ggbn_country', 'ena_country',
                                            'ggbn_collection_date', 'ena_collection_date', 'ggbn_collector',
                                            'ena_collector'])

        i = 1
        for row in read_csv():
            if i > 1:
                search_triplets(row, writer)
            i = i + 1
            if i > 100000:
                break


if __name__ == '__main__':
    ena_institution = "NHMUK"
    voucher_id = "NHMUK_2018"
    pattern = re.escape(ena_institution) + r'[\s\-_]*'
    newStr = re.sub(pattern, '', voucher_id)
    print(newStr)
