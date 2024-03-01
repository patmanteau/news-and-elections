MBB 2022 K -- Big Data Analyseprojekt
=====================================

Getting started
---------------

Ensure that `grep`, `jq` and `gzip` are installed and on your PATH.

Clone this repo:
~~~bash
git clone https://github.com/patmanteau/news-and-elections.git
~~~

Create and activate the conda environment (if required, rename before creating 
by editing the first line of `environment.yml`):
~~~bash
cd news-and-elections
conda env create -f environment.yml
conda activate news-and-elections
~~~

Prepare spaCy:
~~~bash
python -m spacy download de_core_news_lg
~~~

To use our scrapes, you need to request access to our [Zenodo repository](https://doi.org/10.5281/zenodo.10729326) and we'll provide you with a download token. Then, download and unpack the scrapes: 
~~~bash
curl -o files-archive.zip 'https://zenodo.org/api/records/10729326/files-archive?token=<YOUR_ZENODO_TOKEN>'
unzip files-archive.zip
~~~

Now set up *entity-fishing* (refer to its documentation at <https://nerd.readthedocs.io/en/latest/docker.html#running-entity-fishing-with-docker>). First prepare knowledge data:
~~~bash
cd models/entity-fishing
mkdir -p data/db
cd data/db
curl -O https://science-miner.s3.amazonaws.com/entity-fishing/0.0.6/db-kb.zip
curl -O https://science-miner.s3.amazonaws.com/entity-fishing/0.0.6/db-en.zip
curl -O https://science-miner.s3.amazonaws.com/entity-fishing/0.0.6/db-de.zip
unzip db-kb.zip
unzip db-en.zip
unzip db-de.zip
cd ../..
~~~

Launch *entity-fishing* like so:
~~~bash
cd models/entity-fishing
sudo docker compose pull
sudo docker compose up -d
cd ../..
~~~

If everything worked, you should now be able to access ``entity-fishing`` at
http://localhost:8090


Then ask Snakemake for output and grab a coffee or two:
~~~bash
snakemake -c1 all_entities
snakemake -c1 all
~~~

Project Organization
------------

    ├── LICENSE
    ├── Makefile           <- Makefile with commands like `make data` or `make train`
    ├── README.md          <- The top-level README for developers using this project.
    ├── data
    │   ├── external       <- Data from third party sources.
    │   ├── interim        <- Intermediate data that has been transformed.
    │   ├── processed      <- The final, canonical data sets for modeling.
    │   └── raw            <- The original, immutable data dump.
    │
    ├── docs               <- A default Sphinx project; see sphinx-doc.org for details
    │
    ├── .env               <- Credentials, tokens and the like go here. Don't ever add this to git!
    │
    ├── models             <- Trained and serialized models, model predictions, or model summaries
    │
    ├── notebooks          <- Jupyter notebooks. Naming convention is a number (for ordering),
    │                         the creator's initials, and a short `-` delimited description, e.g.
    │                         `1.0-jqp-initial-data-exploration`.
    │
    ├── references         <- Data dictionaries, manuals, and all other explanatory materials.
    │
    ├── reports            <- Generated analysis as HTML, PDF, LaTeX, etc.
    │   └── figures        <- Generated graphics and figures to be used in reporting
    │
    ├── requirements.txt   <- The requirements file for reproducing the analysis environment, e.g.
    │                         generated with `pip freeze > requirements.txt`
    │
    ├── setup.py           <- makes project pip installable (pip install -e .) so src can be imported
    ├── Snakefile          <- Snakemake input file with data workflow definitions
    ├── snakemake.yaml     <- Snakemake configuration file
    ├── src                <- Source code for use in this project.
    │   ├── __init__.py    <- Makes src a Python module
    │   │
    │   ├── data           <- Scripts to download or generate data
    │   │   └── make_dataset.py
    │   │
    │   ├── features       <- Scripts to turn raw data into features for modeling
    │   │   └── build_features.py
    │   │
    │   ├── models         <- Scripts to train models and then use trained models to make
    │   │   │                 predictions
    │   │   ├── predict_model.py
    │   │   └── train_model.py
    │   │
    │   └── visualization  <- Scripts to create exploratory and results oriented visualizations
    │       └── visualize.py
    │
    └── tox.ini            <- tox file with settings for running tox; see tox.readthedocs.io

Credentials, tokens and the like go into a top-level `.env` file:

~~~sh
# Environment variables go here, can be read by `python-dotenv` package:
#
#   `src/script.py`
#   ----------------------------------------------------------------
#    import dotenv
#
#    project_dir = os.path.join(os.path.dirname(__file__), os.pardir)
#    dotenv_path = os.path.join(project_dir, '.env')
#    dotenv.load_dotenv(dotenv_path)
#   ----------------------------------------------------------------
#
# DO NOT ADD THIS FILE TO VERSION CONTROL!
~~~




--------

<p><small>Project based on the <a target="_blank" href="https://drivendata.github.io/cookiecutter-data-science/">cookiecutter data science project template</a>. #cookiecutterdatascience</small></p>
