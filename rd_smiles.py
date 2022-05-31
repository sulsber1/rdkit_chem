from xmlrpc.client import Boolean
from google.cloud import datastore
from flask import Flask, render_template, request, redirect, session, make_response
from rdkit import Chem
import classy_kit


app = Flask(__name__)

# smile = 'CC2(C)CCCC(\C)=C2\C=C\C(\C)=C\C=C\C(\C)=C\C=C\C=C(/C)\C=C\C=C(/C)\C=C\C1=C(/C)CCCC1(C)C'

@app.route('/')
def index():
    return "Please navigate to /smiles to use this API"

@app.route('/<smiles>', methods=['GET', 'POST'])
def get_properties(smiles):

    if request.method == 'GET':
        res = make_response(classy_kit.Molecule_obj(smiles).to_json())
        res.headers.set('Content-Type', 'application/json')
        if 'SMILES' in res.json['errors']:
            err = "Error Parsing SMILES, try POST or review SMILES Submission"
            return (err, 400)
        else:
            return (res, 200)

    if request.method == 'POST':
        content = request.get_json()
        res = make_response(classy_kit.Molecule_helper(content["smiles"]).to_json())
        res.headers.set('Content-Type', 'application/json')
        return (res, 200)

if __name__ == '__main__':
    app.run(host='127.0.0.1', port=8080, debug=True)