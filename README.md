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



<!-- TABLE OF CONTENTS -->
<details>
  <summary>Table of Contents</summary>
  <ol>
    <li>
      <a href="#about-the-project">About The Project</a>
      <ul>
        <li><a href="#built-with">Built With</a></li>
      </ul>
    </li>
    <li>
      <a href="#getting-started">Getting Started</a>
      <ul>
        <li><a href="#prerequisites">Prerequisites</a></li>
        <li><a href="#installation">Installation</a></li>
      </ul>
    </li>
    <li><a href="#usage">Usage</a></li>
    <li><a href="#roadmap">Roadmap</a></li>
    <li><a href="#contributing">Contributing</a></li>
    <li><a href="#license">License</a></li>
    <li><a href="#contact">Contact</a></li>
    <li><a href="#acknowledgments">Acknowledgments</a></li>
  </ol>
</details>



<!-- ABOUT THE PROJECT -->
## About The Project

[![Product Name Screen Shot][product-screenshot]](https://example.com)


<p align="right">(<a href="#top">back to top</a>)</p>



### Built With

* [python](https://python.org/)
* [rdkit](https://rdkit.org/)

<p align="right">(<a href="#top">back to top</a>)</p>



<!-- GETTING STARTED -->
### Installation and starting the server

1. Clone the repo:
   ''' sh
   git clone https://github.com/sulsberd/rdkit_chem.git
   '''
   
2. Initialize a virtual environment:
   ''' sh
   python -m venv rd_smiles .
   '''
   
3. Activate the virtual environment:
   ''' sh
   Scripts\activate - Linux
   Scripts\activate - Windows
   '''
   
4. Install dependencies within the virtual environment:
   ''' sh
   pip install -r requirements.txt
   '''
   
5. Activate a local instance of the server, the server will default to localhost:8080:
   ''' sh
   python rd_smiles.py
   '''

<p align="right">(<a href="#top">back to top</a>)</p>



<!-- USAGE EXAMPLES -->
## Usage

Both GET and POST requests are handled at localhost:8080/<smiles> where <smiles> is a simplified molecular-input line-entry system (SMILES) of a molecule.  

GET requests can be riskier and are not recommended, the parsing of the parameters from the URL isn't robust.
POST requests can handle 0-n requests

_For more examples, please refer to the [Documentation](https://example.com)_ (WIP)

<p align="right">(<a href="#top">back to top</a>)</p>





