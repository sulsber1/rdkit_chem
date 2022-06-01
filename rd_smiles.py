from xmlrpc.client import Boolean
from google.cloud import datastore
from flask import Flask, render_template, request, redirect, session, make_response
from rdkit import Chem
import classy_kit, json


app = Flask(__name__)

# smile = 'CC2(C)CCCC(\C)=C2\C=C\C(\C)=C\C=C\C(\C)=C\C=C\C=C(/C)\C=C\C=C(/C)\C=C\C1=C(/C)CCCC1(C)C'

@app.route('/')
def index():
    return "Please navigate to /smiles to use this API"

@app.route('/<smiles>', methods=['GET', 'POST'])
def get_properties(smiles):

    if request.method == 'GET':
        try:
            res = make_response(classy_kit.Molecule_helper([smiles]).to_json())
            res.headers.set('Content-Type', 'application/json')
            return (res, 200)
        except:
            return ('', 404)
    if request.method == 'POST':
        try:
            content = request.get_json()
            res = make_response(classy_kit.Molecule_helper(content["smiles"]).to_json())
            res.headers.set('Content-Type', 'application/json')
            return (res, 200)
        except:
            return ('', 404)

if __name__ == '__main__':
    app.run(host='127.0.0.1', port=8080, debug=True)