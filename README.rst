MetFrag for Galaxy
==================
|Build Status (Travis)| |Git| |Bioconda|

Galaxy tool wrapper for MetFrag.

Website: http://c-ruttkies.github.io/MetFrag

Source code: https://github.com/c-ruttkies/MetFrag


Version
------

0.1.4


Galaxy
------
`Galaxy <https://galaxyproject.org>`_ is an open, web-based platform for data intensive biomedical research. Whether on the free public server or your own instance, you can perform, reproduce, and share complete analyses. 


TODO
----
- Different data types for input
- Replace python file with `configfile option <https://docs.galaxyproject.org/en/latest/dev/schema.html#tool-configfiles-configfile>`_
- Read Adduct annotations, calculate MZs and recursorIonMode
- Run in parallel


Suspect list
------------

The list of suspects is an aggregated list of in silico predicted MS/MS spectra of natural products from the Universal Natural Products Database (http://pkuxxj.pku.edu.cn/UNPD/index.php). The list is an aggregated version of the github repository https://github.com/oolonek/ISDB/tree/master/Data/dbs.


Developers & Contributors
-------------------------
 - Jordi Capellades (j.capellades.to@gmail.com) - Universitat Rovira i Virgili (Tarragona, Spain)
 - Julien Saint-Vanne (julien.saint-vanne@sb-roscoff.fr) - ABiMS (France)
 - Tom Lawson (t.n.lawson@bham.ac.uk) - University of Birmingham (UK)
 - Kristian Peters (kpeters@ipb-halle.de) - IPB Halle (Germany)
 - Payam Emami (payam.emami@medsci.uu.se) - Uppsala Universitet (Sweden)
 - Steffen Neumann (sneumann@ipb-halle.de) - IPB Halle (Germany)
 - Christoph Ruttkies (christoph.ruttkies@ipb-halle.de) - IPB Halle (Germany)
 - Ralf J. M. Weber (r.j.weber@bham.ac.uk) - `University of Birmingham (UK) <http://www.birmingham.ac.uk/index.aspx>`_


Changes
-------

Version 0.1.4:
 - added UNPD InCHIkey database to be used for automated testing
 - acknowledge additional contributors


License
-------
Released under the GNU General Public License v3.0 (see LICENSE file)


.. |Build Status (Travis)| image:: https://img.shields.io/travis/computational-metabolomics/metfrag-galaxy.svg?style=flat&maxAge=3600&label=Travis-CI
   :target: https://travis-ci.org/computational-metabolomics/metfrag-galaxy

.. |Git| image:: https://img.shields.io/badge/repository-GitHub-blue.svg?style=flat&maxAge=3600
   :target: https://github.com/c-ruttkies/MetFrag

.. |Bioconda| image:: https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg?style=flat&maxAge=3600
   :target: http://bioconda.github.io/recipes/metfrag/README.html
