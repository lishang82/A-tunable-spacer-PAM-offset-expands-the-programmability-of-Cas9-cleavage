# A-tunable-spacer-PAM-offset-expands-the-programmability-of-Cas9-cleavage
This project focuses on the study of CRISPR/Cas9 reprogramming with spacer and PAM separated. It examines the cutting activity and the changes in cutting sites, and uses smFRET to explore the underlying mechanism and its potential application in cells.

## 📂 Repository Structure

Below is an overview of the data processing and statistical tracking scripts compiled across the respective figure panels:

```text
├── README.md                           # This global navigation map
├── 01_upstream_pipeline/               # Linux server-side raw data processing
│   └── 01_alignment_and_analysis.sh   # Integrated pipeline (fastp, novoalign, samtools, castool)
│
├── 02_fig1_cleavage_sites_statistical/ # In vitro cleavage site alignment (Figure 1)
│   ├── findcutsiterepairCigarCdh23FBatch.R
│   ├── findcutsiterepairCigarCdh23RBatch.R
│   ├── findcutsiterepairCigarMyo7aRBatch.R
│   ├── findcutsiterepairCigarOtofRBatch.R
│   ├── findcutsiterepairCigarPcdh15RBatch.R
│   └── findcutsiterepairCigarTmc1RBatch.R
│
├── 03_fig3_random_library/             # High-throughput substrate random library profiling (Figure 3)
│   ├── 01_library.R                    # Decoding substrate random libraries to map barcode and random region sequences
│   ├── 02_NSHenrichmentfactor.R        # Target enrichment computations for non-offset gRNA library-mediated cleavage
│   ├── 03_PAMclassify.R                # Substrate PAM classifications
│   ├── 04_OSHenrichmentfactor.R        # Target enrichment computations for out1-gRNA library-mediated cleavage
│   └── 05_OSH1fixbasefrequencycut.R    # Characterization of high-efficiency cleavage sequence motifs

└── 04_cellassay_statistical/           # Downstream validation in cellular models (Figure 3)
    └── 01_InDel_pattern.R              # Multi-batch CRISPResso2 integration & 2x5 programmatic plotting

## 🛠️ System Requirements & Runtime Dependencies
The downstream computational pipelines and statistical scaling models are implemented and validated using R (version 4.2.0 or higher).
Required Packages
Run the following initialization block within your R console to establish the environment dependencies:
install.packages(c("tidyverse", "VennDiagram", "ggVennDiagram", "patchwork", "scatterplot3d"))
tidyverse (dplyr, ggplot2, tidyr, stringr) - Core data manipulation and matrix operations.
ggVennDiagram / VennDiagram - Cross-batch intersection quality control for library replicability.
scatterplot3d - Projecting duplicate trends to inspect chimera noise convergence.
patchwork - Multi-panel alignment for multi-gRNA array configurations.

## 🚀 Quick Start Examples
1. Upstream Alignment Pipeline
cd 01_upstream_pipeline/
bash 01_alignment_and_analysis.sh <forward_reads.fastq.gz> <reverse_reads.fastq.gz>

2. Random Library Decoding
cd 03_fig3_random_library/
Rscript 01_library.R <path_to_replicate_1.txt> <path_to_replicate_2.txt> <path_to_replicate_3.txt>

3. Cellular InDel pattern Aggregation
cd 04_cellassay_statistical/
Rscript 01_InDel_pattern.R

## 📄 Data Availability & Citation
Sequencing Depositions: Sequencing Depositions: The raw deep sequencing fastq files and processed tracking datasets generated in this study have been deposited in the Genome Sequence Archive (GSA) under accession number PRJCA067153. (Note: Data access will be fully released publicly upon formal manuscript acceptance).

Academic Citation: If these custom-tailored analysis workflows or figure templates support your study, please cite our core paper:
[Authors]. (2026). A tunable spacer-PAM offset expands the programmability of Cas9 cleavage. [Manuscript under submission].
