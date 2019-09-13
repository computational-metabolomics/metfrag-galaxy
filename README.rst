MetFrag for Galaxy
==================
|Build Status (Travis)| |Git| |Bioconda|

Galaxy tool wrapper for MetFrag.

Website: http://c-ruttkies.github.io/MetFrag

Source code: https://github.com/c-ruttkies/MetFrag


Version
------

v2.4.5+galaxy0.1.13

(Using `MetFrag v2.4.5 <https://anaconda.org/bioconda/metfrag>`_)


Galaxy
------
`Galaxy <https://galaxyproject.org>`_ is an open, web-based platform for data intensive biomedical research. Whether on the free public server or your own instance, you can perform, reproduce, and share complete analyses. 


TODO
----
- Additional unit-tests



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
Version 2.4.5+galaxy0.1.13
 - Improved reporting of adducts
 - Removed print output of settings that might reveal sensitive details regarding the database
 - Added missing inchikey column to metchem output

Version 2.4.5+galaxy0.1.12
 - Allow metchem database to be configured via a config.ini file (allows more flexibility)

Version 2.4.5+galaxy0.1.11
 - Update to MetFrag version 2.4.5
 - Boolean for default suspect list fixed
 - Check for no results added
 - Doc updates

Version 2.4.2+galaxy0.1.10
 - Add default option for suspectlist (so user does not need to upload the GNPS list)
 - Changed "Minimum percentage of explain peaks" argument to float
 - Added "skip invalid adducts" option
 - Updated unit-tests
 - Check for empty input file added

Version 2.4.2+galaxy0.1.9
 - Bug fix for neutral mass calculation (Introduced from 0.1.3)
 - Additional adduct forms for negative added
 - Version system changed to IUC style (TOOL_VERSION+GALAXY_TOOL_VERSION)
 - Memory regex check added to stdio - should now repeat at higher memory if low memory error found
 - Update unit-tests to use RP022611 MassBank data (Glucose)

Version 0.1.8:
 - removed quotes for output columns (previously double quotes " were being added for inchikeys and scores
   causing problems for downstream processing.

Version 0.1.7:
 - Bug Fix for ID tracking
 - Bug Fix for skipping final MSP spectra if two empty lines not present

Version 0.1.6:
 - Bug Fix for when NoExplPeaks is zero

Version 0.1.5:
 - Fix to add MetChem command line param
 - Added auto select for MSP schema
 - Added unit test for suspect list

Version 0.1.4:
 - added UNPD InCHIkey database to be used for automated testing
 - acknowledge additional contributors

Version 0.1.3:
 - merge with the latest PhenoMeNal develop version of the module, based on https://github.com/korseby/container-msnbase
 - merge with changes of Julien Saint-Vanne


License
-------
Released under the GNU General Public License v3.0 (see LICENSE file)


.. |Build Status (Travis)| image:: https://img.shields.io/travis/computational-metabolomics/metfrag-galaxy/master.svg?style=flat&maxAge=3600&label=Travis-CI
   :target: https://travis-ci.org/computational-metabolomics/metfrag-galaxy

.. |Git| image:: https://img.shields.io/badge/repository-GitHub-blue.svg?style=flat&maxAge=3600
   :target: https://github.com/c-ruttkies/MetFrag

.. |Bioconda| image:: https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg?style=flat&maxAge=3600
   :target: http://bioconda.github.io/recipes/metfrag/README.html
