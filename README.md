# Code accompanying Portraits of genetic intra-tumour heterogeneity and subclonal selection across cancer types
This repository contains code corresponding to analysis reported in [Portraits of genetic intra-tumour heterogeneity and subclonal selection across cancer types](https://www.biorxiv.org/content/early/2018/05/07/312041)

The submodules in this repository will link to the version of the code that was used to produce the results and are provided as a permanent archive. We are currently working on organising and documenting all the code and, once ready, more links will be added.

Some of the methods may see further development, therefore, if you're interested in applying a method it is encouraged to go to the linked repositories and check out the latest version there.


## Contents

| Repository | Description |
| --- | --- |
| CICC | Method to construct a consensus subclonal architecture |
| CSR | Method to construct a consensus subclonal architecture |
| SimClone1000 | Pipeline that created the SimClone1000 simulated data set |
| icgc\_consensus\_copynumber | Pipeline that constructs the consensus copy number profiles |
| icgc\_consensus\_purity | Pipeline that establishes a consensus purity estimate |
| icgc\_consensus\_clustering\_assignment | Pipeline that assigns all SNVs, indels and SVs to consensus mutation clusters |

## Dependencies

Each of the above named repositories contains a description of the depencies required for that component.

## Used input methods

### Copy number
* ABSOLUTE - [code](https://software.broadinstitute.org/cancer/cga/absolute) - [paper](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4383288/)
* Aceseq - [code](https://github.com/DKFZ-ODCF/ACEseqWorkflow) - [paper](https://www.biorxiv.org/content/early/2017/10/29/210807)
* Battenberg - [code](https://github.com/Wedge-Oxford/battenberg) - [paper](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3428864/)
* cloneHD - [code](https://github.com/andrej-fischer/cloneHD) - [paper](https://www.ncbi.nlm.nih.gov/pubmed/24882004)
* JaBbA - code - paper
* Sclust - code - [paper](https://www.ncbi.nlm.nih.gov/pubmed/29844525)

### Subclonal architecture
* Bayclone - code - [paper](https://www.ncbi.nlm.nih.gov/pubmed/25592605)
* Ccube - [code](https://github.com/keyuan/ccube) - paper
* CliP - [code](https://github.com/wwylab/CliP) - paper
* cloneHD - [code](https://github.com/andrej-fischer/cloneHD) - [paper](https://www.ncbi.nlm.nih.gov/pubmed/24882004)
* CTPSingle - [code](https://github.com/nlgndnmz/CTPsingle) - [paper](https://www.ncbi.nlm.nih.gov/pubmed/28056180)
* DPClust - code - [paper](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3428864/)
* PhylogicNDT - code - paper
* PhyloWGS - [code](https://github.com/morrislab/phylowgs) - [paper](https://www.ncbi.nlm.nih.gov/pubmed/25786235)
* Pyclone - [code](https://bitbucket.org/aroth85/pyclone/) - [paper](https://www.ncbi.nlm.nih.gov/pubmed/24633410)
* Sclust - code - [paper](https://www.ncbi.nlm.nih.gov/pubmed/29844525)
* SVclone - [code](https://github.com/mcmero/SVclone) - [paper](https://www.biorxiv.org/content/early/2017/08/04/172486)

### Other
* TrackSig - [code](https://github.com/YuliaRubanova/TrackSig) - [paper](https://www.biorxiv.org/content/early/2018/03/29/260471)
* MutationTime.R - [code](https://github.com/gerstung-lab/MutationTimeR) - [paper](https://www.biorxiv.org/content/early/2017/08/30/161562)
