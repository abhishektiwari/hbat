CCD Package
===========

.. toctree::
   :maxdepth: 2
   
   ccd_analyzer
   constants_generator
   generate_ccd_constants

Data Storage
------------

CCD data files are automatically downloaded and stored in:

- Default location: ``~/.hbat/ccd-data/``
- Files: ``cca.bcif`` (atoms) and ``ccb.bcif`` (bonds)
- Source: https://models.rcsb.org/

The data is cached locally for offline use and only downloaded when not present.