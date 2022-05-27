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

@app.route('/<smiles>', methods=['GET'])
def get_properties(smiles):
    if request.method == 'GET':
        res = make_response(classy_kit.Molecule_obj(smiles).to_json())
        res.headers.set('Content-Type', 'application/json')
        res.status_code = 200
        return res


if __name__ == '__main__':
    app.run(host='127.0.0.1', port=8080, debug=True)