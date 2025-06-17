# 4-3. 3_Ps_GenL  

Polymer Class are same as key of Ps_GenL.  
The alocated number of polymerization reaction are key of Ps_rxnL.  
The numbers within the square brackets in the monomer class correspond to the values of mon_dic and key of mon_dic_inv.  

## Correspondence Table of Numbers and Content  

**Table 4-8. For the function "bipolym"**  
Defined polymer class, assigned polymerization reactions and corresponding monomer classes.  

| Polymer Class  | Polymerization Reaction | Monomer Class 1    | Monomer Class 2    | Reaction Type         | Product          |
|----------------|-------------------------|--------------------|--------------------|-----------------------|------------------|
| Polyolefin     | 1                       | Vinyl \[1\]        |                    | Addition Chain Polymerization | Homopolymer |
| Polyolefin     | 3                       | Cyclic Olefin \[3\] |                   | Addition Chain Polymerization | Homopolymer |
| Polyolefin     | 101                     | Vinyl \[1\]        | Vinyl \[1\]        | Addition Chain Polymerization | Alternating copolymer |
| Polyolefin     | 102                     | Vinyl \[1\]        | Cyclic Olefin \[3\] | Addition Chain Polymerization | Alternating copolymer |
| Polyolefin     | 103                     | Cyclic Olefin \[3\] | Cyclic Olefin \[3\] | Addition Chain Polymerization | Alternating copolymer |
| Polyester      | 4                       | Lactone \[4\]      |                    | Ring-Opening Chain Polymerization | Homopolymer       |
| Polyester      | 6                       | Hydroxy Carboxylic Acid \[6\] |         | Polycondensation | Homopolymer |
| Polyester      | 105                     | Hydroxy Carboxylic Acid \[6\] | Hydroxy Carboxylic Acid \[6\] | Polycondensation | Alternating copolymer |
| Polyester      | 104                     | Di/Polycarboxylic Acid \[52\] | Di/Polyol \[53\]  | Polycondensation | Homopolymer |
| Polyester (polycarbonate) | 106          | Di/Polyol \[53\]   | Carbon Monoxide \[10\]       | Polycondensation | Homopolymer |
| Polyester      | 112                     | Cyclic Anhydride \[9\] | Epoxide \[2\]  | Ring-Opening Chain Polymerization | Homopolymer |
| Polyether      | 2                       | Epoxide \[2\]       |                   | Ring-Opening Chain Polymerization | Homopolymer |
| Polyether      | 8                       | Hindered Phenol \[8\] |                 | Polycondensation               | Homopolymer    |
| Polyetherh     | 114                     | Bis(p-Halogenated Aryl)Sulfone \[12\] | Di/Polyol (Without Thiol) \[58\] | Polycondensation | Homopolymer |
| Polyetheri     | 115                     | Bis(p-Fluoroaryl)Ketone \[13\] | Di/Polyol (Without Thiol) \[58\] | Polycondensation | Homopolymer |
| Polyamide      | 5                       | Lactam \[5\]       |                    | Ring-Opening Chain Polymerization | Homopolyme   |
| Polyamide      | 7                       | Amino Acid \[7\]   |                    | Polycondensation                  | Homopolymer |
| Polyamide      | 109                     | Amino Acid \[7\]   | Amino Acid \[7\]   | Polycondensation
| Polyamide      | 108                     | Di/Polycarboxylic Acid \[52\] | Di/Polyamine \[54\] | Polycondensation | Homopolymer |
| Polyimide      | 110                     | Di/Polycyclic Anhydride \[56\] | Primary Di/Polyamine \[57\] | Polycondensation | Homopolymer |
| Polyurethane   | 111                     | Di/Polyisocyanate \[55\]       | Di/Polyol \[53\]            | Polyaddition | Homopolymer |
| Polyoxazolidone | 113                    | Di/Polyepoxide \[51\] | Di/Polyisocyanate \[55\]    | Polyaddition | Homopolymer |


**Table 4-9. For the unction "ole_copolym"**  
Defined polymer class, assigned polymeriztion reactions and corresponding monomer classes for the detailed olefin polymerization.  
| Polymer Class | Polymerization Reaction | Monomer Class      | 
|---------------|-------------------------|--------------------|
| ROMP          | 1050                    | cycCH \[1050\]        |
| COC           | none                    | cycCH \[1050\] + aliphCH \[1052\] |
| rec:radi      | -                       | acryl \[1001\], bEWole \[1002\], styryl \[1003\], allyl \[1004\], haloCH \[1005\], vinylester \[1006\], malei \[1007\], conjdiene \[1020\] |
| rec:cati *    | -                       | styryl \[1003\], vinylether \[1030\], tertcatCH \[1031\] |
| rec:ani  *    | -                       | acryl \[1001\], styryl \[1003\], conjdiene \[1020\] |
| rec:coord *   | ROMP, ROMPH, COC        | cycCH \[1050\], aliphCH \[1052\] (for COC) |

* When all the components of the polymerization system correspond to the monomer class, it is recommended as a polymerization method.  