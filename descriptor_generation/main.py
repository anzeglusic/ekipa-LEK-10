import os
from rdkit import Chem
from mordred import Calculator, descriptors, get_descriptors_in_module
from pandas import ExcelFile
from importlib import import_module

USEFUL_DESCRIPTORS = ['nAtom',
 'nHeavyAtom',
 'nH',
 'nO',
 'nBonds',
 'nBondsO',
 'nBondsS',
 'nBondsKS',
 'C2SP3',
 'NsOH',
 'nHBAcc',
 'MPC2',
 'ATS0se',
 'ATS1se',
 'ATS1pe',
 'ATS1are',
 'ATSC0p',
 'ETA_eta_L',
 'ZMIC0',
 'PEOE_VSA1',
 'SlogP_VSA2',
 'VSA_EState9',
 'bpol',
 'MW',
 'nRot',
 'FilterItLogS']

def sdf_parser(path):
    mols = list(Chem.SDMolSupplier(path, removeHs=False))
    if len(mols) != 1 or mols[0] is None:
        raise Exception(f"ERROR reading file {path}")
    mol = mols[0]
    mol.SetProp("_Name", os.path.basename(path))
    return mol

def calculate_descriptors_from_molfiles(molfiles_dir):
    calc = Calculator(descriptors, ignore_3D=False)
    mols = list(map(lambda filename: sdf_parser(os.path.join(molfiles_dir, filename)), os.listdir(molfiles_dir)))
    df = calc.pandas(mols)
    df["name"] = [mol.GetProp("_Name") for mol in mols]
    df.to_excel("mordred_descriptors.xlsx")

def calculate_descriptors_from_smiles(excel_filename):
    calc = Calculator(descriptors, ignore_3D=True)

    excel = ExcelFile(excel_filename)
    excel_data_frame = excel.parse("Sheet1", 0)

    mols = [Chem.MolFromSmiles(smi) for smi in excel_data_frame["Updated SMILES"]]

    df = calc.pandas(mols)
    df = df[USEFUL_DESCRIPTORS]
    df["name"] = excel_data_frame["Name"]
    df["Updated SMILES"] = excel_data_frame["Updated SMILES"]
    #df["BA"] = excel_data_frame["Biološka uporabnost (%)"]
    df.to_excel("mordred_descriptors.xlsx")


def calculate_specific_descriptors_from_smiles(excel_filename, descriptors_names):
    descs = (
                d
                for m in descriptors_names
                for d in get_descriptors_in_module(import_module("mordred."+m))
            )
    calc = Calculator(descs, ignore_3D=True)

    excel = ExcelFile(excel_filename)
    excel_data_frame = excel.parse("Sheet1", 0)

    mols = [Chem.MolFromSmiles(smi) for smi in excel_data_frame["Updated SMILES"]]

    df = calc.pandas(mols)
    df["name"] = excel_data_frame["Name"]
    df["BA"] = excel_data_frame["Biološka uporabnost (%)"]
    df.to_excel("mordred_descriptors.xlsx")


def fix_molfiles(molfiles_dir):
    for filename in os.listdir(molfiles_dir):
        filepath = os.path.join(molfiles_dir,filename)
        with open(filepath, 'r', errors='ignore') as file:
            if "V2000" not in file.readlines()[3]:
                with open(filepath, 'r+', errors='ignore') as f:
                    content = f.read()
                    f.seek(0)
                    f.write('\n' + content)


if __name__ == '__main__':
    #molfiles_dir = "SDF"
    #fix_molfiles(molfiles_dir)
    #calculate_descriptors_from_molfiles(molfiles_dir)
    calculate_descriptors_from_smiles("BA lista za Hackatlon - test set1.xlsx")
    # calculate_specific_descriptors_from_smiles("BA lista za Hackatlon.xlsx", [
    #    "ABCIndex"
    # ])



