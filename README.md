# Kincore-standalone
## Installation and environment setup
First, download the repository:
```
git clone https://github.com/DunbrackLab/Kincore-standalone3
```
Then, create a virtual environment with the necessary packages to run kincore.
### Option 1: use pixi (recommended)
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
Optionally, you can use Anaconda to create a virtual enviroment and install the necessary packages.
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
*Note: any of the following commands can also be run using ```python kinase_state.py``` instead of ```pixi run kincore```, provided you have the necessary packages installed and/or virtual environment active.*

Run the help command to see available options:
```
pixi run kincore -h
```
Kincore can be run on a single mmCIF or PDB file (can contain multiple models in a single file, e.g., MD trajectories or PDB structures).
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
Now you can run kincore using ```python kinase_state.py``` instead of ```pixi run kincore```. For example:
```
python kinase_state.py list.txt > output.txt &
```
## Understanding the output
Each kinase chain occupies 6 lines, with different data listed on each line:
```
1OL5      0 A    Active_DFGin_BLAminus_SBin_NTin_CTin          Labels    family CAMK  hmm CAMK     score  273.7   Active   DFGin    BLAminus Chelix-in   Saltbr-in   ActLoopNT-in   ActLoopCT-in   Spine-in   APE10-dihe-na APE9-dihe-na APE8-dihe-in APE8-rot-in APE67-dihe-in APE12-dist-in APE11-dist-in APE10-dist-in APE9-dist-in HRD-in
1OL5      0 A    Active_DFGin_BLAminus_SBin_NTin_CTin          Residues  APEp3.E302 APEp2.I301 APE.E299  APE6.D294 APE7.L293 APE8.T292 APE9.G291 APE10.C290 APE11.L289 APE12.T288 XDFG.A273  Asp.D274  Phe.F275  Gly.G276  DFG4.W277 DFG6.V279 Lys.K162  Glu.E181  Glu4.Q185 HPN7.L196  XHRD.I253  HRD.H254  Arg.R255  HRDAsp.D256 aFasp.D311
1OL5      0 A    Active_DFGin_BLAminus_SBin_NTin_CTin          Distances Glu4_Phe   5.58 Lys_Phe  15.10 Lys_Glu   9.45 SaltBr   2.91 DFG6_XHRD   2.93 APE9_Arg   3.68 APE3_aFasp3   6.13 APE10_DFG4   5.27 APE11_DFG4  10.96 APE12_DFG4   8.17 Spine    3.59   3.51   3.91   3.91
1OL5      0 A    Active_DFGin_BLAminus_SBin_NTin_CTin          Dihedrals HRD  -63.28  -60.88 Arg   73.60    1.59 XDFG -139.20 -172.08 DFG   52.22   80.16 -163.54   -3.74 DFG2  -90.83   24.66  282.61   72.55 DFG3  -49.87  -44.18 APE6  -62.23  -14.30 APE7  -43.16  -62.50 APE8 -110.77  140.24  313.06 APE9  119.45 -153.19 APE10 -153.42  162.50
1OL5      0 A    Active_DFGin_BLAminus_SBin_NTin_CTin          Ligands   ADP:1388,MG:1389,MG:1390,MG:1394    ATPlike,Allosteric,Allosteric,Allosteric
1OL5      0 A    Active_DFGin_BLAminus_SBin_NTin_CTin          Actloop     26    17.70   12.06   28.30
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

#### Row 1: ``` Labels ``` 

| Label  | Meaning |
| ------------------- | ------------------- |
| family CAMK  | 1OL5 chain A belongs to the CAMK family |
| hmm CAMK  | Membership to the CAMK family was determined via the HMM filename "CAMK.hmm" located in the "HMMs" subdirectory. There are extra HMMs for unusual kinases, such as BUB, PEAK, and TP53RK that are members of the OTHER family.|
| A | Chain ID |
| score  273.7 | HMM score of the chain A sequence for CAMK.hmm |
| Active | The activation loop and αC-helix pass all of the criteria for active kinases |
| DFGin | The DFG Phe (or equivalent residue) is flipped "in" towards the αC-helix |
| BLAminus | The backbone dihedrals of the first three residues of the XDFG motif occupy the Ramchandran regions B, L, and A, and the χ1 dihedral of the DFG Phe sidechain (or equivalent residue) is gauche-minus|
| Saltbr-in | A salt bridge is formed between the αC-helix Glu and β3-strand Lys |
| ActLoopNT-in | The N-terminus of the activation loop (near the DFG motif) is hydrogen bonded to the HRD loop |
| ActLoopCT-in | The N-terminus of the activation loop (near the DFG motif) is hydrogen bonded to the HRD loop |
| Spine-in | The regulatory spine is formed. |
| APE10dihe-na | Backbone dihedral state of the 10th residue from the end of the APE motif (in the C-terminus of the activation loop). 1OL5 chain A belongs to the CAMK family, not TYR, so dihedral criteria for this residue are not applicable to determining the Active/Inactive state. |
| APE9dihe-na | Backbone dihedral state for the 9th residue from the end of the APE motif. Labeled "NA" for the same reason as APE10dihe. |
| APE8dihe-in | Backbone dihedral state for the 8th residue from the end of the APE motif is active-like, hence the label "APE8dihe-in" (in analogy to how DFG-"in" is an active-like conformation of the DFG motif). |
| APE8rot-in | Sidechain rotamer state for the 8th residue from the end of the APE motif is active-like. |
| APE67dihe-in | Backbone dihedral states for the 6th and 7th residues from the end of the APE motif are active-like. |
| APE12-dist-in | Distance between CB of the APE12 residue and CA of the DFG4 residue is active-like. |
| APE11-dist-in | Distance between CB of the APE11 residue and CA of the DFG4 residue is active-like. |
| APE10-dist-in | Distance between CB of the APE10 residue and CA of the DFG4 residue is active-like. |
| APE9-dist-in | Distance between CA of the APE9 residue and backbone carbonyl oxygen of the HRD Arg residue is active-like. |
| HRD-in | Backbone dihedral states for the HRD His and HRD Arg residues are active-like. |

#### Row 2: ``` Residues ``` 
List of residues used for conformational assignments. Kincore nomenclature followed by the residue type and number in the sequence (separated by ".").
| Nomenclature | Description |
| ------------------- | ------------------- |
| APEp3  | "APE plus 3", i.e., the third residue N-terminal to the APE motif. |
| APEp2  | "APE plus 2", i.e., the 2nd residue N-terminal to the APE motif. |
| APE | The end of the APE motif / the "E" of the APE (or equivalent residue). |
| APE6 | The 6th residue from the end of the APE motif, counting backwards starting from the "E" of the APE. |
| APE7 | The 7th residue from the end of the APE motif. |
| APE8 | The 8th residue from the end of the APE motif. |
| APE9 | The 9th residue from the end of the APE motif. |
| APE10 | The 10th residue from the end of the APE motif. |
| APE11 | The 11th residue from the end of the APE motif. |
| APE12 | The 12th residue from the end of the APE motif. |
| XDFG | The "X" residue immediately before the DFG motif. |
| Asp | The beginning of the DFG motif / the "D" of the DFG (or equivalent residue). |
| Phe | The 2nd residue from the beginning of the DFG motif / the "F" of the DFG (or equivalent residue). |
| Gly | The 3nd residue from the beginning DFG motif / the "G" of the DFG (or equivalent residue). |
| DFG4 | The 4th residue from the beginning DFG motif, counting forwards from the "D" of the DFG. |
| DFG6 | The 6th residue from the beginning DFG motif. |
| Lys | The conserved Lys in the β3-strand (or equivalent residue). |
| Glu | The conserved Glu in the αC-helix (or equivalent residue). |
| Glu4 | The 4th residue from the Glu in the αC-helix. |
| HPN7 | The 7th residue from the beginning of the HPN motif, located in the β4-strand (part of the regulatory spine). |
| XHRD | The "X" residue immediately before the HRD motif. |
| HRD | The conserved His (or sometimes Tyr) residue at the beginning of the HRD motif. |
| Arg | The 2nd residue from the beginning of the HRD motif / the "R" of the HRD motif (or equivalent residue). |
| HRDAsp | The 3nd residue from the beginning of the HRD motif / the "D" of the HRD motif (or equivalent residue). |
| aFasp | The conserved Asp in the αF-helix.  |

#### Row 3: ``` Distances ``` 
