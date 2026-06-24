# Kincore-standalone
## Installation and environment setup
First, download the repository:
```
git clone https://github.com/DunbrackLab/Kincore-standalone3
```
Then, create a virtual environment with the necessary packages to run kincore.
### Option 1: use pixi
Install pixi if you haven't already:
```
curl -fsSL https://pixi.sh/install.sh | bash
```
```
cd Kincore-standalone3
```
That's it. You should now be able to run kincore inside the Kincore-standalone3 directory using `pixi run kincore` (see below for more details). 

*Note: the first time you run this command, pixi will create the environment and install all the necessary packages in the Kincore-standalone3 directory.*

### Option 2: use Anaconda
You can also create an anconda enviroment and install the necessary packages.
```
conda create --name 'kincore-standalone3' python=3.8 pandas numpy biopython hmmer --channel conda-forge --channel bioconda
```
Activate the virtual environment
```
conda activate kincore-standalone3
```
## Running Kincore
Go to the Kincore-standalone3 folder if you haven't already:
```
cd Kincore-standalone3
```
Note: any of the following commands can also be run using ```python kinase_state.py``` instead of ```pixi run kincore```, provided you have the necessary packages installed (such as hmmsearch) and/or virtual environment active.

Run the help command to see available options:
```
pixi run kincore -h
```
Kincore can be run on a single mmCIF or PDB file (can contain multiple models in a single file, e.g., MD trajectories or NMR structures).
```
pixi run kincore 1OL5.cif
pixi run kincore 1OL5.pdb
```
Also works on gzipped files:
```
pixi run kincore 1OL5.cif.gz
pixi run kincore 1OL5.pdb.gz
```
To run Kincore on more than one structure, give a list of structure filenames instead (full path names recommended for robustness).
```
pixi run kincore list.txt
```
Recommended: direct the output of Kincore to a text file so you can access it later. Also, attaching "&" at the end allows kincore to run in the background.
```
pixi run kincore list.txt > output.txt &
```
### Other ways of running kincore
First, activate your virtual environment inside the Kincore-standalone3 folder
```
pixi shell
```
Or (if using Anaconda): 
```
conda activate kincore-standalone3
```
Now you can run kincore using ```python kinase_state.py```. For example:
```
python kinase_state.py list.txt > output.txt &
```
## Understanding the output
Each kinase chain occupies at-least 17 lines, with different data listed on each line (18 lines if an autoinhibitory interaction with HRDAsp is detected):
```
1OL5A               0 A    Active_DFGin_BLAminus_SBin_NTin_CTin          Sequence         family CAMK  hmm CAMK     score  273.7   length 405
1OL5A               0 A    Active_DFGin_BLAminus_SBin_NTin_CTin          Labels           Active   DFGin    DFGAsp-rot-in BLAminus Chelix-in   Saltbr-in   ActLoopNT-in   ActLoopCT-in   HRD-in Spine-in  
1OL5A               0 A    Active_DFGin_BLAminus_SBin_NTin_CTin          APEdihe_Labels   APE10-dihe-na APE9-dihe-na APE8-dihe-in APE8-rot-in APE67-dihe-in
1OL5A               0 A    Active_DFGin_BLAminus_SBin_NTin_CTin          APEdist_Labels   APE12-dist-in APE11-dist-in APE10-dist-in APE9-dist-in
1OL5A               0 A    Active_DFGin_BLAminus_SBin_NTin_CTin          NTdom_Residues   Lys.K162     Glu.E181     Glu4.Q185    HPN7.L196 
1OL5A               0 A    Active_DFGin_BLAminus_SBin_NTin_CTin          CTdom_Residues   XHRD.I253    HRD.H254     Arg.R255     HRDAsp.D256  aFasp.D311 
1OL5A               0 A    Active_DFGin_BLAminus_SBin_NTin_CTin          DFG_Residues     XDFG.A273    Asp.D274     Phe.F275     Gly.G276     DFG4.W277   DFG6.V279 
1OL5A               0 A    Active_DFGin_BLAminus_SBin_NTin_CTin          APE_Residues     APE.E299     APE6.D294    APE7.L293    APE8.T292    APE9.G291   APE10.C290 APE11.L289 APE12.T288
1OL5A               0 A    Active_DFGin_BLAminus_SBin_NTin_CTin          HRD_bbDihedrals  HRDHis    -63.28  -60.88 HRDArg   73.60    1.59
1OL5A               0 A    Active_DFGin_BLAminus_SBin_NTin_CTin          XDFG_bbDihedrals XDFG     -139.20 -172.08 DFGAsp   52.22   80.16 DFGPhe  -90.83   24.66 DFGGly  -49.87  -44.18
1OL5A               0 A    Active_DFGin_BLAminus_SBin_NTin_CTin          APE_bbDihedrals  APE6      -62.23  -14.30 APE7    -43.16  -62.50 APE8   -110.77  140.24 APE9    119.45 -153.19 APE10  -153.42  162.50
1OL5A               0 A    Active_DFGin_BLAminus_SBin_NTin_CTin          scDihedrals      DFGAsp    196.46   -3.74 DFGPhe  282.61   72.55 APE8    313.06
1OL5A               0 A    Active_DFGin_BLAminus_SBin_NTin_CTin          APE_Distances    APE9_Arg   3.68 APE10_DFG4   5.27 APE11_DFG4  10.96 APE12_DFG4   8.17
1OL5A               0 A    Active_DFGin_BLAminus_SBin_NTin_CTin          Spine_Distances  Spine1     3.59 Spine2       3.51 Spine3       3.91 MaxSpine     3.91
1OL5A               0 A    Active_DFGin_BLAminus_SBin_NTin_CTin          Other_Distances  Glu4_Phe   5.58 Lys_Phe     15.10 Lys_Glu      9.45 SaltBr       2.91 DFG6_XHRD   2.93
1OL5A               0 A    Active_DFGin_BLAminus_SBin_NTin_CTin          Ligands          ADP:1388,MG:1389,MG:1390,MG:1394    ATPlike,Allosteric,Allosteric,Allosteric
1OL5A               0 A    Active_DFGin_BLAminus_SBin_NTin_CTin          Actloop          Len   26    Ave 17.70   Min 12.06   Max 28.30
```

### Headers
Each line starts with a row header that identifies the structure file, model number, and chain ID, followed by the conformational state. This information is repeated on each line.
```
1OL5      0 A    Active_DFGin_BLAminus_SBin_NTin_CTin
```
| Header  | Meaning |
| ------------- | ------------- |
| 1OL5  | Structure file  |
| 0  | Model number  |
| A | Chain ID |

### Conformational label
The long string in the 4th column **Active_DFGin_BLAminus_SBin_NTin_CTin** is the overall conformational state of the kinase chain.

| Label | Meaning |
| --- | ---|
| Active | The activation loop and αC-helix pass all of the criteria for active kinases |
| DFGin | The DFG Phe (or equivalent residue) is flipped "in" towards the αC-helix |
| BLAminus | The backbone dihedrals of the first three residues of the XDFG motif occupy the Ramchandran regions B, L, and A, and the χ1 dihedral of the DFG Phe sidechain (or equivalent residue) is gauche-minus|
| SBin | A salt bridge is formed between the αC-helix Glu and β3-strand Lys |
| NTin | The N-terminus of the activation loop (near the DFG motif) is hydrogen bonded to the HRD loop |
| CTin | The C-terminus of the activation loop (near the APE motif) adopts the appropriate structure for binding protein substrates |

### Data
The 5th column contains another header that tells you what kind of data follows it. Residues with structural parameters that go into determining the state of the activation loop C-terminus are labeled APE*i* where *i* is the *i*th residue counting backwards from the end of the APE motif.

#### Row 1: ``` Sequence ``` 
| Label  | Meaning |
| ------------------- | ------------------- |
| family CAMK  | 1OL5 chain A belongs to the CAMK family |
| hmm CAMK  | Membership to the CAMK family was determined via the HMM filename "CAMK.hmm" located in the "HMMs" subdirectory. There are extra HMMs for unusual kinases, such as BUB, PEAK, and TP53RK that are members of the OTHER family.|
| score  273.7 | HMM score of the chain A sequence for CAMK.hmm |
| length 405 | Length of the sequence | 
 
#### Row 2: ``` Labels ``` 
| Label  | Meaning |
| ------------------- | ------------------- |
| Active | The activation loop and αC-helix pass all of the criteria for active kinases |
| DFGin | The DFG Phe (or equivalent residue) is flipped "in" towards the αC-helix |
| DFGAsp-rot-in | The χ<sub>1</sub> rotamer of the DFG Asp residue is trans, allowing it to bind ATP and Mg |
| BLAminus | The backbone dihedrals of the first three residues of the XDFG motif occupy the Ramachandran regions B, L, and A, and the χ1 dihedral of the DFG Phe sidechain (or equivalent residue) is gauche-minus|
| Saltbr-in | A salt bridge is formed between the αC-helix Glu and β3-strand Lys |
| ActLoopNT-in | The N-terminus of the activation loop (near the DFG motif) is hydrogen bonded to the HRD loop |
| ActLoopCT-in | The N-terminus of the activation loop (near the DFG motif) is hydrogen bonded to the HRD loop |
| HRD-in | Backbone dihedral states for the HRD His and HRD Arg residues are active-like. |
| Spine-in | The regulatory spine is formed. |

### Row 3: ``` APEdihe_Labels ```
| Label  | Meaning |
| ------------------- | ------------------- |
| APE10dihe-na | Backbone dihedral state of the 10th residue from the end of the APE motif (in the C-terminus of the activation loop). 1OL5 chain A belongs to the CAMK family, not TYR, so dihedral criteria for this residue are not applicable to determining the Active/Inactive state. |
| APE9dihe-na | Backbone dihedral state for the 9th residue from the end of the APE motif. Labeled "NA" for the same reason as APE10dihe. |
| APE8dihe-in | Backbone dihedral state for the 8th residue from the end of the APE motif is active-like, hence the label "APE8dihe-in" (in analogy to how DFG-"in" is an active-like conformation of the DFG motif). |
| APE8rot-in | Sidechain rotamer state for the 8th residue from the end of the APE motif is active-like. |
| APE67dihe-in | Backbone dihedral states for the 6th and 7th residues from the end of the APE motif are active-like. |

### Row 4: ``` APEdist_Labels ```
| Label  | Meaning |
| ------------------- | ------------------- |
| APE12-dist-in | Distance between CB of the APE12 residue and CA of the DFG4 residue is active-like. |
| APE11-dist-in | Distance between CB of the APE11 residue and CA of the DFG4 residue is active-like. |
| APE10-dist-in | Distance between CB of the APE10 residue and CA of the DFG4 residue is active-like. |
| APE9-dist-in | Distance between CA of the APE9 residue and backbone carbonyl oxygen of the HRD Arg residue is active-like. |

### Row 5: ``` NTdom_Residues ```
Residues in the N-terminal domain (N-terminal lobe)
| Nomenclature | Description |
| ------------------- | ------------------- |
| Lys | The conserved Lys in the β3-strand (or equivalent residue). |
| Glu | The conserved Glu in the αC-helix (or equivalent residue). |
| Glu4 | The 4th residue from the Glu in the αC-helix. |
| HPN7 | The 7th residue from the beginning of the HPN motif, located in the β4-strand (part of the regulatory spine). |

### Row 6: ``` CTdom_Residues ```
Residues in the C-terminal domain (C-terminal lobe)
| Nomenclature | Description |
| ------------------- | ------------------- |
| XHRD | The "X" residue immediately before the HRD motif. |
| HRD | The conserved His (or sometimes Tyr) residue at the beginning of the HRD motif. |
| Arg | The 2nd residue from the beginning of the HRD motif / the "R" of the HRD motif (or equivalent residue). |
| HRDAsp | The 3nd residue from the beginning of the HRD motif / the "D" of the HRD motif (or equivalent residue). |
| aFasp | The conserved Asp in the αF-helix.  |

### Row 7: ``` DFG_Residues ```
Residues in the Activation Loop N-terminal segment (ActLoopNT), starting with the XDFG motif (X-Asp-Phe-Gly)
| Nomenclature | Description |
| ------------------- | ------------------- |
| XDFG | The "X" residue immediately before the DFG motif. |
| Asp | The beginning of the DFG motif / the "D" of the DFG (or equivalent residue). |
| Phe | The 2nd residue from the beginning of the DFG motif / the "F" of the DFG (or equivalent residue). |
| Gly | The 3nd residue from the beginning DFG motif / the "G" of the DFG (or equivalent residue). |
| DFG4 | The 4th residue from the beginning DFG motif, counting forwards from the "D" of the DFG. |
| DFG6 | The 6th residue from the beginning DFG motif. |

### Row 8: ``` APE_Residues ```
Residues in the Activation Loop C-terminal segment (ActLoopCT), ending with the APE motif (Ala-Pro-Glu). 
ActLoopCT residues are enumerated backwards (C-term. to N-term.) starting with Glu of APE ("APE1" or just "APE").
| Nomenclature | Description |
| ------------------- | ------------------- |
| APE | The end of the APE motif / the Glu ("E") of the APE (or equivalent residue). |
| APE6 | The 6th residue from the end of the APE motif, counting backwards starting from the "E" of the APE. |
| APE7 | The 7th residue from the end of the APE motif. |
| APE8 | The 8th residue from the end of the APE motif. |
| APE9 | The 9th residue from the end of the APE motif. |
| APE10 | The 10th residue from the end of the APE motif. |
| APE11 | The 11th residue from the end of the APE motif. |
| APE12 | The 12th residue from the end of the APE motif. |

### Row 9: ``` HRD_bbDihedrals ```
Backbone (φ, ψ) dihedral angles of the HRD loop (also called "catalytic loop").
| Nomenclature | Description |
| ------------------- | ------------------- |
| HRD |	φ,ψ	of HRD-His in A region	of Ramachandran	map for	active kinases
| Arg | φ,ψ of HRD-Arg in L region	of Ramachandran	map for	active kinases

### Row 10: ``` XDFG_bbDihedrals ```
Backbone (φ, ψ) dihedral angles of the XDFG (X-Asp-Phe-Gly) motif at the N-terminus of the activation loop.
| Nomenclature | Description |
| ------------------- | ------------------- |
| XDFG | φ,ψ of X of XDFG motif. B region for BLAminus state of active kinases
| DFGAsp |  φ,ψ,χ<sub>1</sub>,χ<sub>2</sub> of Asp	of DFG motif. L region for BLAminus state of active kinases (and trans χ<sub>1</sub> rotamer)
| DFGPhe |  φ,ψ,χ<sub>1</sub>,χ<sub>2</sub> of Phe of DFG	motif. A region for BLAminus state of active kinases
| DFGGly |  φ,ψ of Gly of DFG	motif. A region for BLA(A)minus state of active kinases

### Row 11: ``` APE_bbDihedrals ```
Backbone (φ, ψ) dihedral angles of the Activation Loop C-terminal segment (ActLoopCT).
| Nomenclature | Description |
| ------------------- | ------------------- |
| APE6 |  φ,ψ of APE6 residue (6th residue from end of activation loop: X in XxxAPE). APE7,APE6 == AA or BL for active kinases
| APE7 |  φ,ψ of APE7 residue (6th residue from end of activation loop: X in XxxxAPE). APE7,APE6 == AA or BL for active kinases
| APE8 |  φ,ψ,χ<sub>1</sub> of APE8 residue (6th residue from end of activation loop: X in XxxxxAPE). B region for active	kinases. χ<sub>1</sub> in (-120°,0°).
| APE9 |  φ,ψ of APE9 residue (6th residue from end of activation loop: X in XxxxxxAPE). B region for active TYR kinases
| APE10 |  φ,ψ of APE10 residue (6th residue from end of activation loop: X in XxxxxxxAPE).	B region for active TYR kinases

### Row 12: ``` scDihedrals ```
Sidechain dihedral angles of select residues involved in conformational assignments (χ<sub>1</sub> and, for certain residues, χ<sub>2</sub>).
| Nomenclature | Description |
| ------------------- | ------------------- |
| DFGAsp |  χ<sub>1</sub>,χ<sub>2</sub> of Asp	of DFG motif. trans χ<sub>1</sub> rotamer for active kinases.
| DFGPhe |  χ<sub>1</sub>,χ<sub>2</sub> of Phe of DFG	motif. gauche-minus χ<sub>1</sub> rotamer for active kinases.
| APE8 |  χ<sub>1</sub> of APE8 residue (6th residue from end of activation loop: X in XxxxxAPE). gauche-minus χ<sub>1</sub> rotamer in (-120°,0°) for active non-TYR kinases.

### Row 13: ``` APE_Distances ```
Distances involving ActLoopCT residues used for conformational assignments of non-TYR kinases.
| Nomenclature | Description |
| ------------------- | ------------------- |
| APE9_Arg | APE9-Cα / hRd-Arg O for ActLoopCT-in/ActLoopCT-out calculation for nonTYR kinases
| APE10_DFG4 | APE10-Cβ / DFG4-Cα for ActLoopCT-in/ActLoopCT-out calculation for nonTYR kinases
| APE11_DFG4 | APE10-Cβ / DFG4-Cα for ActLoopCT-in/ActLoopCT-out calculation for nonTYR kinases
| APE12_DFG4 | APE10-Cβ / DFG4-Cα for ActLoopCT-in/ActLoopCT-out calculation for nonTYR kinases

### Row 14: ``` Other_Distances ```
Distances involving the DFG motif, ActLoopNT, and C-helix residues used for conformational assignments.
| Nomenclature | Description |
| ------------------- | ------------------- |
| Glu4_Phe | Glu4-Cα / DFG-Phe Cζ distance for DFGin/DFGout/DFGinter calculation. Glu4 is 4 residues after the salt-bridge Glu
| Lys_Phe | Lys-Cα / DFG-Phe Cζ distance for DFGin/DFGout/DFGinter calculation. Lys is the salt-bridge Lys
| Lys_Glu | Lys-Cβ / Glu-Cβ distance for Chelix-in/Chelix-out calculation.
| SaltBr | Lys-Nζ / Glu OE1,OE2 distance for SaltBr-in/SaltBr-out calculation. Minimum of distance to OE1 and OE2
| DFG6_XHRD | DFG6-N/O / Xhrd O/N distance for ActLoopNT-in/ActLoopNT-out calculation. Minimum of two backbone-backbone hydrogen bonds distances.

### Row 15: ``` Spine_Distances ```
Distances used to measure the state of the regulatory spine (Spine-in/Spine-out). 
 
| Spine1 | Nearest sidechain-atom distance between HRD-His and DFG-Phe
| Spine2 | Nearest sidechain-atom distance between DFG-Phe and Glu4
| Spine3 | Nearest sidechain-atom distance between Glu4 and HPN7
| MaxSpine | Maximum value of the three spine distances to determine if spine is broken (Spine-out when MaxSpine>4.5 Å)

#### Row 16: ``` Ligands ``` 
List of ligands and their types (3 or 5 letter codes from PDB, residue numbers, and types
| Nomenclature | Description |
| ------------------- | ------------------- |
| ATPlike | Any ATP-like PDB ligand (ATP, ACP, ANP, ADP, AGS)
| Type1 | Ligand occupies ATP-binding site
| Type1.5 | Ligand occupies ATP-binding site and part of Chelix site
| Type2 | Ligand occupies both ATP-binding site and Chelix site
| Type3 | Ligand occupies Chelix site
| Allosteric | Ligand is elsewhere

#### Row 17: ``` Actloop ``` 
Minimum, maximum, and average of B-factors of Ca atoms of activation loop. Useful for calculating min(pLDDT) of activation loop of AlphaFold/Boltz models.

#### Row 18: ``` Autoinhibit ``` 
Present if there is any Ser, Thr, or Tyr in activation loop in hydrogen bonding distance of HRD-Asp
