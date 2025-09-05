.. SMiPoly_docs documentation master file, created by
   sphinx-quickstart on Tue May  6 20:54:12 2025.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

SMiPoly documentation
==========================

Current version: |version|, release: |release|


About SMiPoly  
==========================

"SMiPoly (\ **S**\ mall **M**\ olecules **i**\ nto **Poly**\ mers)" is rule-based virtual library generator for discovery of functional polymers. It is consist of two submodules, "monc.py" and "polg.py".  
"monc.py" is a monomer classifier from a list of small molecules, and "polg.py" is a polymer repeating unit generator from the classified monomer list.  

| **How To Cite (publications)** 
| SMiPoly: Generation of a Synthesizable Polymer Virtual Library Using Rule-Based Polymerization Reactions  
| Mitsuru Ohno, Yoshihiro Hayashi, Qi Zhang, Yu Kaneko, and Ryo Yoshida  
| *Journal of Chemical Information and Modeling* **2023** *63* (17), 5539-5548  
| DOI: 10.1021/acs.jcim.3c00329  
| https://pubs.acs.org/doi/full/10.1021/acs.jcim.3c00329
| (version 0.0.1 was used)

**Project Link** https://github.com/PEJpOhno/SMiPoly  

**Demo video**  

.. image:: ./_static/SMiPoly_movie1.jpg
   :target: https://youtu.be/ilzYwNWvTeQ
   :alt: Demo video
   :width: 300px
   
.. toctree::
   :maxdepth: 2
   :caption: Contents:
   :hidden:

   1_Installing_SMiPoly
   2_Getting_start
   3_Modules_Overview
   smipoly.smip.monc
   smipoly.smip.polg
   smipoly.smip.funclib
   4_utilities
   41_MonomerDefiner
   42_Ps_rxnL
   43_Ps_GenL
   5_Postface


