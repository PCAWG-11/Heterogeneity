# Code accompanying Characterizing genetic intra-tumor heterogeneity across 2,658 human cancer genomes
This repository contains code corresponding to analysis reported in [Characterizing genetic intra-tumor heterogeneity across 2,658 human cancer genomes](https://www.biorxiv.org/content/10.1101/312041v5)

The submodules in this repository link to the version of the code that was used to produce the results and are provided as a permanent archive. 

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
| MutationTimeR | Version of MutationTimeR used to assign mutations to clusters |
| weme | Method to construct a consensus subclonal architecture |
| SpoilSport | Method to correct cluster positions and sizes for missed SNVs (i.e. winner's curse) |

## Dependencies

Each of the above named repositories contains a description of the depencies required for that component.

## Used input methods

### Copy number
* ABSOLUTE - [code](https://software.broadinstitute.org/cancer/cga/absolute) - [paper](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4383288/)
* Aceseq - [code](https://github.com/DKFZ-ODCF/ACEseqWorkflow) - [paper](https://www.biorxiv.org/content/early/2017/10/29/210807)
* Battenberg - [code](https://github.com/Wedge-Oxford/battenberg) - [paper](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3428864/)
* cloneHD - [code](https://github.com/andrej-fischer/cloneHD) - [paper](https://www.ncbi.nlm.nih.gov/pubmed/24882004)
* JaBbA - [code](https://github.com/mskilab/JaBbA) - [paper](https://doi.org/10.1016/j.cell.2020.08.006)
* Sclust - [code](http://www.uni-koeln.de/med-fak/sclust/Sclust.tgz) - [paper](https://www.ncbi.nlm.nih.gov/pubmed/29844525)

### Subclonal architecture
* Bayclone - code - [paper](https://www.ncbi.nlm.nih.gov/pubmed/25592605)
* Ccube - [code](https://github.com/keyuan/ccube) - [paper](https://www.biorxiv.org/content/10.1101/484402v1)
* CliP - [code](https://github.com/wwylab/CliP) - paper
* cloneHD - [code](https://github.com/andrej-fischer/cloneHD) - [paper](https://www.ncbi.nlm.nih.gov/pubmed/24882004)
* CTPSingle - [code](https://github.com/nlgndnmz/CTPsingle) - [paper](https://www.ncbi.nlm.nih.gov/pubmed/28056180)
* DPClust - [code](https://github.com/Wedge-Oxford/dpclust) - [paper](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3428864/)
* PhylogicNDT - [code](https://github.com/broadinstitute/PhylogicNDT/) - [paper](https://www.biorxiv.org/content/10.1101/508127v2)
* PhyloWGS - [code](https://github.com/morrislab/phylowgs) - [paper](https://www.ncbi.nlm.nih.gov/pubmed/25786235)
* Pyclone - [code](https://bitbucket.org/aroth85/pyclone/) - [paper](https://www.ncbi.nlm.nih.gov/pubmed/24633410)
* Sclust - code - [paper](https://www.ncbi.nlm.nih.gov/pubmed/29844525)
* SVclone - [code](https://github.com/mcmero/SVclone) - [paper](https://www.biorxiv.org/content/early/2017/08/04/172486)

### Other
* TrackSig - [code](https://github.com/YuliaRubanova/TrackSig) - [paper](https://www.biorxiv.org/content/early/2018/03/29/260471)
* MutationTime.R - [code](https://github.com/gerstung-lab/MutationTimeR) - [paper](https://www.biorxiv.org/content/early/2017/08/30/161562)
