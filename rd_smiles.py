from flask import Flask, request, make_response
import classy_kit_concurrent
import json

app = Flask(__name__)


@app.route('/smiles', methods=['GET', 'POST'])
def get_properties():
    if request.method == 'POST':
        try:
            content = request.get_json()
            res = make_response(json.dumps(classy_kit_concurrent.get_descriptors(content)))
            res.headers.set('Content-Type', 'application/json')
            return (res, 200)
        except:
            return ('', 404)

if __name__ == '__main__':
    app.run(host='127.0.0.1', port=8080, debug=True)