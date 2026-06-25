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
