from fastapi import FastAPI
from pydantic import BaseModel
import pandas as pd
from fastapi.middleware.cors import CORSMiddleware

from app.modules.hackathon_make_prediction import ML_model

app = FastAPI()

origins = ["*"]

app.add_middleware(
    CORSMiddleware,
    allow_origins=origins,
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)

# STRENIRAMO MODEL
model = ML_model("app/data/training.csv")

# NALOŽIMO PODATKE
excel = pd.ExcelFile('app/data/BA lista za Hackatlon COMBINED.xlsx')
metapodatki = excel.parse("Sheet1",0)
metapodatki.set_index('Updated SMILES', inplace=True)

class PredictResponse(BaseModel):
    smiles: str
    drug_name: str
    bioavailability: float

@app.get("/smiles/{smiles_id}", response_model=PredictResponse)
def predict_by_smiles_formula(smiles_id: str):
    smiles_id = smiles_id.strip()
    drug_name = get_drug_name(smiles_id, test_mode=True)
    ba_value = get_ba(smiles_id)
    #ba_value = 45.55

    result = {'smiles': smiles_id, 'drug_name': drug_name, 'bioavailability': ba_value}
    return result


def get_drug_name(smiles, test_mode=False):
    try:
        podatki = metapodatki.loc[smiles]
    except:
        return 'Drug not existing'

    if test_mode:
        print(f"PRAVILNA BA: {podatki['Biološka uporabnost (%)']}")
        return podatki['Name']
    else:
        return podatki['Name']
    
def get_ba(smiles_id):
    #NALOŽIMO ZNAČILKE ZA PRODUKCIJO
    znacilke = pd.read_csv('app/data/predict.csv')
    znacilke.drop(columns= ['name'], inplace=True)
    znacilke.set_index('Updated SMILES', inplace=True)
    try:
        znacilke = znacilke.loc[smiles_id]
        BA = model.make_prediction(znacilke)
        return BA
    except:
        print('ERROR!!!!! Value not found')
        return 0.0



