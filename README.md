<div id="top"></div>

<h3 align="center">RD_SMILES and Classy_Kit</h3>

  <p align="center">
    project_description (WIP)
    <br />
    <a href="https://github.com/github_username/repo_name"><strong>Explore the docs »</strong></a>
    <br />
    <br />
    <a href="https://github.com/github_username/repo_name">View Demo</a>
    ·
    <a href="https://github.com/github_username/repo_name/issues">Report Bug</a>
    ·
    <a href="https://github.com/github_username/repo_name/issues">Request Feature</a>
  </p>
</div>


<!-- ABOUT THE PROJECT -->
## About The Project

<p align="right">(<a href="#top">back to top</a>)</p>



### Built With

* [python](https://python.org/)
* [rdkit](https://rdkit.org/)

<p align="right">(<a href="#top">back to top</a>)</p>



<!-- GETTING STARTED -->
### Installation and starting the server

1. Clone the repo:
   ``` sh
   git clone https://github.com/sulsberd/rdkit_chem.git
   ```
   
2. Initialize a virtual environment:
   ``` sh
   python -m venv rd_smiles .
   ```
   
3. Activate the virtual environment:
   ``` sh
   Scripts\activate - Linux
   Scripts\activate - Windows
   ```
   
4. Install dependencies within the virtual environment:
   ``` sh
   pip install -r requirements.txt
   ```
   
5. Activate a local instance of the server, the server will default to localhost:8080:
   ``` sh
   python rd_smiles.py
   ```

<p align="right">(<a href="#top">back to top</a>)</p>



<!-- USAGE EXAMPLES -->
## Usage

Both GET and POST requests are handled at localhost:8080/<smiles> where <smiles> is a simplified molecular-input line-entry system (SMILES) of a molecule.  

GET requests can be riskier and are not recommended, the parsing of the parameters from the URL isn't robust.
POST requests can handle 0-n requests
  
POST body must contain a "smiles" key, followed by an array of comma seperated SMILES strings.
  
  ``` sh
  {
    "smiles": ["CC(=O)NCCc1c[nH]c2ccc(OC)cc12", "CC1(CCC(=C(C1)C2=CC=C(C=C2)Cl)CN3CCN(CC3)C4=CC(=C(C=C4)C(=O)NS(=O)(=O)C5=CC(=C(C=C5)NCC6CCOCC6)[N+](=O)[O-])OC7=CN=C8C(=C7)C=CN8)C",
    "C=CC(=O)N1CCC[C@H](C1)N2C3=C(C(=N2)C4=CC=C(C=C4)OC5=CC=CC=C5)C(=NC=N3)N"]
}
```

_For more examples, please refer to the [Documentation](https://example.com)_ (WIP)

<p align="right">(<a href="#top">back to top</a>)</p>





