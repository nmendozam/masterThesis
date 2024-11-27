# Genome-Scale Metabolic Modeling of a Human Astrocyte

This repository contains the data, input models, and code to construct and analyze the genome-scale metabolic model (GEM) of a human astrocyte under three conditions: **healthy**, **inflamed with palmitate**, and **treated with tibolone**.

The final models are available in the [gibbslab/GSMs](https://github.com/gibbslab/GSMs) repository, which also includes a collection of other GEMs.

## Cloning and Running the Repository

To clone this repository and ensure all submodules are downloaded correctly (including the COBRA Toolbox), run the following commands:

```bash
git clone --recurse-submodules https://github.com/nmendozam/masterThesis.git
cd masterThesis
git clone https://github.com/gibbslab/GSMs.git
```

## Description of Scripts

The repository includes the following scripts:

1. **`DataPrep.m`**  
   Prepares the omic data used for contextualizing the GEM model.

2. **`BuildAstrocyteModel.m`**  
   Builds the baseline GEM of the human astrocyte.

3. **`GapFilling.m`**  
   Fills the gaps in the astrocyte model to ensure metabolic network connectivity.

4. **`leakSiphonModes.m`**  
   Semi-automated script to identify and fix leak and siphon reactions in the model.

5. **`FBAanalysis.m`**  
   Performs Flux Balance Analysis (FBA) for the different astrocyte conditions.  
   **Note**: This script requires cloning the [gibbslab/GSMs](https://github.com/gibbslab/GSMs) repository.



## Models

The four astrocyte models generated in this project are listed in this [pull request](https://github.com/gibbslab/GSMs/pull/1/files) and can be downloaded by cloning the [gibbslab/GSMs](https://github.com/gibbslab/GSMs) repository.

1. **[Astrocyte_Mendoza2022.xml.gz](https://github.com/gibbslab/GSMs/blob/master/Astrocyte_Mendoza2022.xml.gz)**  
   A reduced model containing only astrocyte reactions, with all exchange reactions open. This serves as the base model for the three contextualized models.

2. **[Astrocyte_Healthy_Mendoza2022.xml.gz](https://github.com/gibbslab/GSMs/blob/master/Astrocyte_Healty_Mendoza2022.xml.gz)**  
   The model contextualized for a healthy astrocyte.

3. **[Astrocyte_InflamedPalmitate_Mendoza2022.xml.gz](https://github.com/gibbslab/GSMs/blob/master/Astrocyte_InflamedPalmitate_Mendoza2022.xml.gz)**  
   The model contextualized for an astrocyte inflamed with palmitate.

4. **[Astrocyte_TreatedTibolone_Mendoza2022.xml.gz](https://github.com/gibbslab/GSMs/blob/master/Astrocyte_TreatedTibolone_Mendoza2022.xml.gz)**  
   The model contextualized for an astrocyte inflamed with palmitate and treated with tibolone.

---

## Citation

If you use this repository or any of the associated models, please cite the following manuscript (currently under review):

Angarita-Rodríguez, A., Mendoza-Mejía, N., González, J., Papin, J. A., & Pinzón, A. (2024). *Improvement in the prediction power of an astrocyte genome-scale metabolic model using multi-omic data*. Manuscript under review, *Frontiers*.


