# Genome-scale metabolic modeling of astrocytes

This repository contains the code and input models to produce the genome-scale metabolic models of astrocytes described in my master thesis. However, the final models are located in [gibbslab/GSMs](https://github.com/gibbslab/GSMs) rather than here. The 4 models produced for the astrocyte are listed in this [pull request](https://github.com/gibbslab/GSMs/pull/1/files) and can be downloaded by cloning the repository [gibbslab/GSMs](https://github.com/gibbslab/GSMs).

1.  [Astrocyte_Mendoza2022.xml.gz](https://github.com/gibbslab/GSMs/blob/master/Astrocyte_Mendoza2022.xml.gz "Astrocyte_Mendoza2022.xml.gz") It is the reduced model (which has only the reactions of the astrocyte) but with all the exchange reactions open. This is the base model for the other three contextualized models.
2. [Astrocyte_Healty_Mendoza2022.xml.gz](https://github.com/gibbslab/GSMs/blob/master/Astrocyte_Healty_Mendoza2022.xml.gz "Astrocyte_Healty_Mendoza2022.xml.gz") It is the model contextualized to a healthy astrocyte.
3. [Astrocyte_InflamedPalmitate_Mendoza2022.xml.gz](https://github.com/gibbslab/GSMs/blob/master/Astrocyte_InflamedPalmitate_Mendoza2022.xml.gz "Astrocyte_InflamedPalmitate_Mendoza2022.xml.gz") It is the model contextualized to an inflamed astrocyte.
4. [Astrocyte_TreatedTibolone_Mendoza2022.xml.gz](https://github.com/gibbslab/GSMs/blob/master/Astrocyte_TreatedTibolone_Mendoza2022.xml.gz "Astrocyte_TreatedTibolone_Mendoza2022.xml.gz") It is the model contextualized to an inflamed astrocyte and then treated with tibolone.

# Running this code
To clone this repository you should run as it will make sure to download the cobratoolbox submodule.
```bash
git clone --recurse-submodules https://github.com/nmendozam/masterThesis.git
```


