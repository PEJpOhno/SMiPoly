# 5. utilities; Overview

.ipynb files in directiry "utirities" are tools to imprementing the rules of monomer extraction-cllasification , polymerization reaction definition and polymer generation.  

By using the files in './utilities' directory, one can modify or add the definition of monomers, the rules of polymerization reactions and polymer classes.  
To apply the new rule(s), replace the old './smipoly/rules' directory by the new one. The files must be run according to the number assigned the head of the each filename.  

  - 1_MonomerDefiner.ipynb: definitions of monomers  
  - 2_Ps_rxnL.ipynb: rules of polymerization reactions    
  - 3_Ps_GenL.ipynb: definitions of polymer classes with combinations of starting monomer(s) and polymerization reaction  

Each definations were written by [SMARTS (SMiles ARbitrary Target Specification](https://www.daylight.com/dayhtml_tutorials/languages/smarts/)) strings.  



| Polymer Class   | Monomer Class 1           | Monomer Class 2           | Reaction Type                 |
|-----------------|---------------------------|---------------------------|-------------------------------|
| Polyolefin      | Vinyl                    |                           | Addition Chain Polymerization |
| Polyolefin      | Cyclic Olefin            |                           | Addition Chain Polymerization |
| Polyolefin      | Vinyl                    | Vinyl                     | Addition Chain Polymerization |
| Polyolefin      | Vinyl                    | Cyclic Olefin             | Addition Chain Polymerization |
| Polyolefin      | Cyclic Olefin            | Cyclic Olefin             | Addition Chain Polymerization |
| Polyester       | Lactone                  |                           | Ring-Opening Chain Polymerization |
| Polyester       | Hydroxy Carboxylic Acidc |                           | Polycondensation             |
| Polyester       | Hydroxy Carboxylic Acidc | Hydroxy Carboxylic Acidc  | Polycondensation             |
| Polyester       | Di/Polycarboxylic Acid   | Di/Polyol                 | Polycondensation             |
| Polyester       | Di/Polyol                | Carbon Monoxidee          | Polycondensationf            |
| Polyester       | Cyclic Anhydride         | Epoxide                   | Ring-Opening Chain Polymerization |
| Polyether       | Epoxide                  |                           | Ring-Opening Chain Polymerization |
| Polyether       | Hindered Phenol          |                           | Polycondensationg            |
| Polyetherh      | Bis(p-Halogenated Aryl)Sulfone | Di/Polyol (Without Thiol) | Polycondensation             |
| Polyetheri      | Bis(p-Fluoroaryl)Ketone  | Di/Polyol (Without Thiol) | Polycondensation             |
| Polyamide       | Lactam                  |                           | Ring-Opening Chain Polymerization |
| Polyamide       | Amino Acidc             |                           | Polycondensation             |
| Polyamide       | Amino Acidc             | Amino Acidc               | Polycondensation             |
| Polyamide       | Di/Polycarboxylic Acid  | Di/Polyamine              | Polycondensation             |
| Polyimide       | Di/Polycyclic Anhydride | Primary Di/Polyamine      | Polycondensation             |
| Polyurethane    | Di/Polyisocyanate       | Di/Polyol                 | Polyaddition                 |
| Polyoxazolidone | Di/Polyepoxide          | Di/Polyisocyanate         | Polyaddition                 |

