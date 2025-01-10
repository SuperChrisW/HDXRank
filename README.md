# HDXRank
**HDXRank is an open-source pipeline to apply HDX-MS restraints to protein-protein complex prediction ranking.**

## Method overview:
<img src="figures/HDXRank_overview.png" style="width:100%;">
Integrating sparse experimental data into protein complex modeling workflows can significantly improve the accuracy and reliability of predicted models. Despite the valuable insights that hydrogen-deuterium exchange (HDX) data provide about protein binding interfaces, there is currently no standard protocol for incorporating this information into the complex model selection process. Inspired by advances in graph-based deep learning for protein representation, we utilize it as a backbone for a flexible scoring framework in protein-protein complex model ranking based on their alignment with experimental HDX profiles. It offers a robust, HDX-informed selection protocol with improved prediction accuracy.

## Installation:
clone the repository and Use the `HDXRank_minimum.yml` file to create a Conda environment with all necessary dependencies:
```
git clone https://github.com/SuperChrisW/HDXRank.git
cd HDXRank
conda env create -f ./HDXRank_minimum.yml
conda activate HDXRank
```

## Preprocessing
* temporarily refers to AI-HDX document(https://github.com/Environmentalpublichealth/AI-HDX/blob/main/Documentations/MSA_embedding.md)

## Getting Started
HDXRank requires three input files:

1. **Protein structure file** (`.pdb`)  
2. **MSA file** (`.hhm`)  
3. **HDX-MS file** (`.xlsx`)  

Additionally, HDXRank uses a settings file (`.xml`) to control the pipeline.

### Workflow:

1. **Protein embedding**: HDXRank extracts embeddings from `.pdb` and `.hhm` files.  
2. **Protein graph construction**: Constructs a protein graph from the `.pdb` file.
3. **Peptide graph splitting**: Splits the protein graph into peptide graphs based on the provided HDX-MS file.

### Execution:
With all input files prepared, run the following command to start the pipeline:
```bash
python main.py -input ./settings/BatchTable_setting.xml
```

## Citing HDXRank
if you use HDXRank, please cite the following paper: 
[submit for publication]