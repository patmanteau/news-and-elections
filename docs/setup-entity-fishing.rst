Setting up ``entity-fishing``
=============================

(For ``entity-fishing`` documentation, see <https://nerd.readthedocs.io/en/latest/docker.html#running-entity-fishing-with-docker>

Prepare knowledge data::

  $ mkdir -p data/db
  $ cd data/db
  $ curl -O https://science-miner.s3.amazonaws.com/entity-fishing/0.0.6/db-kb.zip
  $ curl -O https://science-miner.s3.amazonaws.com/entity-fishing/0.0.6/db-en.zip
  $ curl -O https://science-miner.s3.amazonaws.com/entity-fishing/0.0.6/db-de.zip
  $ unzip db-kb.zip
  $ unzip db-en.zip
  $ unzip db-de.zip
  $ cd ../..


Link compose file::

  $ ln -s ../news_and_elections/docs/docker-compose.yml .

Launch ``entity-fishing``::

  $ sudo docker-compose pull
  $ sudo docker-compose up -d

If everything worked, you should now be able to access ``entity-fishing`` at
http://localhost:8090
