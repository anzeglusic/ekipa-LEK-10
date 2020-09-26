import requests
import base64, zlib
import pandas
import csv


HEADERS = {'apikey': 'tjnOWFH3TwnZwbNgVCvFGGReW9DvO0qj'}
BASE_URL = "https://api.rsc.org/compounds/v1"


def get_query_id(smiles):
    url = f"{BASE_URL}/filter/smiles"
    data = {"smiles": smiles}
    resp = requests.post(url, json=data, headers=HEADERS)
    if resp.ok:
        return resp.json()["queryId"]
    else:
        return f"ERROR {resp.status_code}"


def get_results_sdf(query_id):
    url = f"{BASE_URL}/filter/{query_id}/results/sdf"
    resp = requests.get(url, headers=HEADERS)
    return resp.json()["results"]


def save_sdf_file(smiles, filename):
    query_id = get_query_id(smiles)
    sdf_zipped_contents = get_results_sdf(query_id)
    f = open(filename, "w")
    f.write(zlib.decompress(base64.b64decode(sdf_zipped_contents), 16 + zlib.MAX_WBITS).decode('utf-8'))
    f.close()


def get_results(query_id):
    url = f"{BASE_URL}/filter/{query_id}/results"
    resp = requests.get(url, headers=HEADERS, params={"start":0, "count":100})
    return resp.json()["results"]




def save_all_query_IDs():
    excel = pandas.ExcelFile("BA lista za Hackatlon.xlsx", )
    data_frame = excel.parse("Sheet1",0)
    skip = True
    for smiles in data_frame["Updated SMILES"]:
        if smiles == "O=C(NC1CC2C(N(C1)C)Cc1c3c2cccc3[nH]c1)N(CC)CC":
            skip = False
        if skip:
            continue
        query_id = get_query_id(smiles)
        f = open("query_ids.txt", "a")
        f.write(f"{smiles},{query_id}\n")
        f.close()



def save_all_results():
    with open('employee_birthday.txt') as csv_file:
        csv_reader = csv.reader(csv_file, delimiter=',')
        line_count = 0
        for row in csv_reader:
            if line_count == 0:
                print(f'Column names are {", ".join(row)}')
                line_count += 1
            else:
                print(f'\t{row[0]} works in the {row[1]} department, and was born in {row[2]}.')
                line_count += 1
        print(f'Processed {line_count} lines.')


save_sdf_file("OCC1CC(n2c3nc(nc(NC4CC4)c3nc2)N)C=C1", "Abacavir")