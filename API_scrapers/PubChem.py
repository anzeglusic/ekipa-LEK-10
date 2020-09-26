import requests
import urllib.parse
import os
from pandas import ExcelFile

BASE_URL = "https://pubchem.ncbi.nlm.nih.gov/rest/pug"


def get_json(smiles):
    url = f"{BASE_URL}/compound/smiles/{urllib.parse.quote(smiles)}/JSON"
    request = requests.get(url)
    return request.json()


def get_SDF(smiles, filename):
    url = f"{BASE_URL}/compound/smiles/{urllib.parse.quote(smiles)}/SDF?record_type=3d"
    request = requests.get(url)
    if request.ok:
        with open(os.path.join("SDF", filename), 'wb') as file:
            file.write(request.content)
    else:
        print(f"ERROR {request.status_code} at {smiles}")


def get_all_SDFs(excel_path):
    excel = ExcelFile(excel_path)
    excel_df = excel.parse("Sheet1", 0)
    for name, smiles in zip(excel_df['Name'], excel_df['Updated SMILES']):
        get_SDF(smiles, f"{name}.sdf")


def parse_json(json):
    dict = {}
    props = json["PC_Compounds"][0]["props"]
    dict["Compound Complexity"] = find_by_x(props, "label", "Compound Complexity")
    dict["Hydrogen Bond Acceptor"] = find_by_x(props, "name", "Hydrogen Bond Acceptor")
    dict["Hydrogen Bond Donor"] = find_by_x(props, "name", "Hydrogen Bond Donor")
    dict["Rotatable Bond"] = find_by_x(props, "name", "Rotatable Bond")
    dict["Log P"] = find_by_x(props, "label", "Log P")
    dict["Mass"] = find_by_x(props, "label", "Mass")
    dict["Molecular Weight"] = find_by_x(props, "label", "Molecular Weight")
    dict["Polar Surface Area"] = find_by_x(props, "name", "Polar Surface Area")
    dict["MonoIsotopic Weight"] = find_by_x(props, "name", "MonoIsotopic")
    count = json["PC_Compounds"][0]["count"]
    dict["heavy_atom"] = count["heavy_atom"]
    dict["atom_chiral"] = count["atom_chiral"]
    dict["atom_chiral_def"] = count["atom_chiral_def"]
    dict["atom_chiral_undef"] = count["atom_chiral_undef"]
    dict["bond_chiral"] = count["bond_chiral"]
    dict["bond_chiral_def"] = count["bond_chiral_def"]
    dict["bond_chiral_undef"] = count["bond_chiral_undef"]
    dict["isotope_atom"] = count["isotope_atom"]
    dict["covalent_unit"] = count["covalent_unit"]
    dict["tautomers"] = count["tautomers"]
    return dict

def find_by_x(array, x, label):
    found = list(filter(lambda el : (x in el["urn"] and el["urn"][x] == label), array))
    return list(found[0]["value"].values())[0]


def get_all_descriptors(excel_path):
    excel = ExcelFile(excel_path)
    excel_df = excel.parse("Sheet1", 0)
    descriptors = []
    for name, smiles in zip(excel_df['Name'], excel_df['Updated SMILES']):
        try:
            desc = parse_json(get_json(smiles))
            desc["name"] = name
            desc["smiles"] = smiles
            descriptors.append(desc)
        except:
            print(f"Cannot parse descriptors for {smiles}")
    #TODO return DataFrame(descriptors)
    return descriptors


if __name__ == '__main__':
    #get_all_SDFs("BA lista za Hackatlon.xlsx")
    d = get_all_descriptors("BA lista za Hackatlon.xlsx")
