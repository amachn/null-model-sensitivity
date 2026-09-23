# Project Papers

This directory contains metadata for the studies selected for reproduction and sensitivity analysis.

Each paper is assigned a stable ID (`P01`-`P12`) and stored as a YAML file. These IDs should ideally be used in reproduction directories, documentation, results, and project tracking. 

The table below is generated automatically from the YAML metadata files. Do not edit the table manually.

<!-- PAPER_TABLE_START -->

| ID | Paper | Owner | DOI | PMID | Data | Code | Reproduction |
|---|---|---|---|---|---|---|---|
| [P01](/reproductions/P01/) | Short-term diet intervention comprising of olive oil, vitamin D, and omega-3 fatty acids alters the small non-coding RNA (sncRNA) landscape of human sperm | Addison | [Link](https://doi.org/10.1038/s41598-024-83653-4) |[Link](https://pubmed.ncbi.nlm.nih.gov/40044751) | [Link](https://github.com/Candida05/Diet_Intervention_Sperm/blob/main/Genomic-Annot-Enrichment.zip) | [Link](https://github.com/Candida05/Diet_Intervention_Sperm) | QC |
| P02 | TADA-a machine learning tool for functional annotation-based prioritisation of pathogenic CNVs | Addison | [Link](https://doi.org/10.1186/s13059-022-02631-z) |[Link](https://pubmed.ncbi.nlm.nih.gov/35232478) | [Link](https://github.com/jakob-he/TADA/) | [Link](https://github.com/jakob-he/TADA/) | Not Started |
| [P03](/reproductions/P03/) | Novel insights from a multiomics dissection of the Hayflick limit | Addison | [Link](https://doi.org/10.7554/eLife.70283) |[Link](https://pubmed.ncbi.nlm.nih.gov/35119359) | [Link](https://github.com/dghendrickson/hayflick/) | [Link](https://github.com/dghendrickson/hayflick/) | QC |
| P04 | LAP2alpha facilitates myogenic gene expression by preventing nucleoplasmic lamin A/C from spreading to active chromatin regions | Anugraka | [Link](https://doi.org/10.1093) |[Link](https://pubmed.ncbi.nlm.nih.gov/39228367) | [Link](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE247774) |  | In Progress |
| P05 | Integrative analysis identifies bHLH transcription factors as contributors to Parkinson's disease risk mechanisms. | Anugraka | [Link](https://doi.org/10.1038) |[Link](https://pubmed.ncbi.nlm.nih.gov/33568722) | [Link](https://bendlj01.u.hpc.mssm.edu/multireg/resources/boca_peaks.zip) |  | In Progress |
| P06 | Integration of epigenetic and genetic profiles identifies multiple sclerosis disease-critical cell types and genes | Anugraka | [Link](https://doi.org/10.1038) |[Link](https://pubmed.ncbi.nlm.nih.gov/36997638) | [Link](https://screen.wenglab.org/) |  | Not Started |
| P07 | Archaic Introgression Shaped Human Circadian Traits | Sheridan | [Link](https://doi.org/10.1093/gbe/evad203) |[Link](https://pubmed.ncbi.nlm.nih.gov/38095367) | [Link](will be in folder for use) |  | In Progress (cleanup) |
| P08 | Genome-wide compendium and functional assessment of in vivo heart enhancers | Sheridan | [Link](https://doi.org/10.1038/ncomms12923) |[Link](https://pubmed.ncbi.nlm.nih.gov/27703156) |  |  | Not Started |
| P09 | Genome Sequencing of Autism-Affected Families Reveals Disruption of Putative Noncoding Regulatory DNA | Sheridan | [Link](https://doi.org/10.1016/j.ajhg.2015.11.23) |[Link](https://pubmed.ncbi.nlm.nih.gov/26749308) | [Link](will be in folder) |  | Not Started |
| [P10](/reproductions/P10/) | Expanded encyclopaedias of DNA elements in the human and mouse genomes | Advika | [Link](https://doi.org/10.1038/s41586-020-2493-4) |[Link](https://pubmed.ncbi.nlm.nih.gov/32728249) |  |  | QC |
| P11 | GARFIELD classifies disease-relevant genomic features through integration of functional annotations with association signals | Advika | [Link](https://doi.org/10.1038/s41588-018-0322-6) |[Link](https://pubmed.ncbi.nlm.nih.gov/30692680) |  |  | Not Started |
| P12 | Chromatin marks identify critical cell types for fine mapping complex trait variants | Advika | [Link](https://doi.org/10.1038/ng.2504) |[Link](https://pubmed.ncbi.nlm.nih.gov/23263488) |  |  | Not Started |

<!-- PAPER_TABLE_END -->

To regenerate the table:

```bash
nullsens papers refresh
```