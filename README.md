# Spatially resolved multiomics of human cardiac niches

Processed data of sc/snRNAseq and Visium data are available for browsing and download via [heartcellatlas.org](https://www.heartcellatlas.org/)

## Contents
### RNA
* QC
* Integration, Celltype annotation
* CCS

### ATAC
* QC
* Peak calling
* Peak-to-Gene (ArchR)
* Dimensional reduction, integration (PeakVI)

### Visium
* QC, Cell2location mapping
* Cell2location NMF
* Cell state spatial enrichment analysis
* Matching NMF-factors and histological annotations 
* CCS niche
* Immune niche
* Myocardial-stress niche
* SAN-RAGP

## Citation
Spatially resolved multiomics of human cardiac niches
[https://doi.org/10.1101/2023.01.30.526202](https://www.biorxiv.org/content/10.1101/2023.01.30.526202v2)

```
@article {Kanemaru2023.01.30.526202,
	author = {Kanemaru, Kazumasa and Cranley, James and Muraro, Daniele and Miranda, Antonio M.A. and Pett, Jan Patrick and Litvinukova, Monika and Kumasaka, Natsuhiko and Ho, Siew Yen and Polanski, Krzysztof and Richardson, Laura and Mach, Lukas and Dabrowska, Monika and Richoz, Nathan and Barnett, Sam N. and Perera, Shani and Wilbrey-Clark, Anna L and Talavera-L{\'o}pez, Carlos and Mulas, Ilaria and Mahbubani, Krishnaa T. and Bolt, Liam and Mamanova, Lira and Tuck, Liz and Wang, Lu and Huang, Margaret M. and Prete, Martin and Pritchard, Sophie and Dark, John and Saeb-Parsy, Kourosh and Patel, Minal and Clatworthy, Menna R. and Chowdhury, Rasheda A. and Noseda, Michela and Teichmann, Sarah A.},
	title = {Spatially resolved multiomics of human cardiac niches},
	elocation-id = {2023.01.30.526202},
	year = {2023},
	doi = {10.1101/2023.01.30.526202},
	publisher = {Cold Spring Harbor Laboratory},
	abstract = {A cell{\textquoteright}s function is defined by its intrinsic characteristics and its niche: the tissue microenvironment in which it dwells. Here, we combine single-cell and spatial transcriptomic data to discover cellular niches within eight regions of the human heart. We map cells to micro-anatomic locations and integrate knowledge-based and unsupervised structural annotations. For the first time, we profile the cells of the human cardiac conduction system, revealing their distinctive repertoire of ion channels, G-protein coupled receptors, and cell interactions using a custom CellPhoneDB.org module. We show that the sinoatrial node is compartmentalised, with a core of pacemaker cells, fibroblasts and glial cells supporting paracrine glutamatergic signalling. We introduce a druggable target prediction tool, drug2cell, which leverages single-cell profiles and drug-target interactions, providing unexpected mechanistic insights into the chronotropic effects of drugs, including GLP-1 analogues. In the epicardium, we show enrichment of both IgG+ and IgA+ plasma cells forming immune niches which may contribute to infection defence. We define a ventricular myocardial-stress niche enriched for activated fibroblasts and stressed cardiomyocytes, cell states that are expanded in cardiomyopathies. Overall, we provide new clarity to cardiac electro-anatomy and immunology, and our suite of computational approaches can be deployed to other tissues and organs.Competing Interest StatementIn the past three years, S.A.T. has consulted or been a member of scientific advisory boards at Roche, Genentech, Biogen, GlaxoSmithKline, Qiagen and ForeSite Labs and is an equity holder of Transition Bio. The remaining authors declare no competing interests.},
	URL = {https://www.biorxiv.org/content/early/2023/02/01/2023.01.30.526202},
	eprint = {https://www.biorxiv.org/content/early/2023/02/01/2023.01.30.526202.full.pdf},
	journal = {bioRxiv}
}
```

## Acknowledgments
This repository is dedicated to the memory of our dear friend and colleague Dr. Daniele Muraro who contributed to this analysis.
