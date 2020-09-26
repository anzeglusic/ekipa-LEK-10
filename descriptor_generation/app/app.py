from rdkit import Chem
from mordred import Calculator, descriptors
from flask import Flask, Response

app = Flask(__name__)

@app.route('/')
def mordred():
    calc = Calculator(descriptors, ignore_3D=True)
    mol = Chem.MolFromSmiles('c1ccccc1')
    df = calc.pandas([mol])

    return Response(df.to_csv(), mimetype='text/csv')

if __name__ == '__main__':
    app.run()