# Transcriptomic Profiling of Pediatric Glioblastoma Multiforme: Mapping the Molecular Pathways of Tumor Aggressiveness

**Submitted by:** Nazmi Soufi Alnifari  
**Course:** Omics Lite - Week 03  

### I. Introduction
Pediatric Glioblastoma Multiforme (GBM) is not merely a localized tumor but a systemic molecular hijacking of the developing brain's architecture. As the most aggressive and lethal primary brain cancer (WHO Grade IV), GBM is characterized biologically by rapid cellular proliferation, extensive tissue invasion, and the induction of severe intracranial pressure that ultimately leads to fatal brainstem compression. 

Despite originating within the central nervous system, Glioblastoma tissue exhibits a radically different molecular phenotype compared to healthy neural tissue. This study aims to decipher the transcriptomic landscape driving this functional dichotomy. Using a differential gene expression (DEG) approach on pediatric samples, we seek to identify the precise molecular machinery and biological pathways that the tumor exploits to execute its uncontrolled expansion and suppress normal neurological functions.

### II. Methods

**Dataset & Study Design**
The transcriptomic dataset was obtained from the NCBI Gene Expression Omnibus (GEO) database under accession number GSE50161. The data was generated using the Affymetrix Human Genome U133 Plus 2.0 Array. 
To ensure strict comparative purity, the analysis utilized a Case vs. Control design comparing Glioblastoma samples against Normal Brain Tissue. Data from other low-grade brain tumors present in the original dataset were programmatically filtered out prior to analysis.

**Analysis Parameters & Reproducibility Workflow**
To ensure reproducibility, the analysis was conducted using a standardized R/Bioconductor pipeline with the following steps:
1. **Data Retrieval:** Raw matrix data and platform annotations were fetched directly via the `GEOquery` package.
2. **Pre-processing & Exploratory Data Analysis (EDA):** Expression values were evaluated for log2 transformation. Data normalization and global sample clustering were validated using Boxplot, Density Plot, and UMAP dimensionality reduction (`umap` and `ggplot2` packages).
3. **DEG Statistics:** Differential expression analysis was executed using the `limma` linear modeling framework. The statistical thresholds were strictly set at an Adjusted P-value < 0.01 and |LogFC| > 1.
4. **Functional Enrichment:** Over-representation analysis for Gene Ontology (GO - Biological Process) and KEGG pathways was performed using `clusterProfiler` to translate significant DEGs into biological insights.

### III. Results and Biological Interpretation

**A. Exploratory Data Analysis (EDA): The Global Shift**
Before gene-level analysis, the quality and distribution of the data were validated to prevent batch effects and ensure statistical reliability.

![Boxplot](https://github.com/user-attachments/assets/5476568a-ba51-410c-95b1-8daf6f29e7c5)
 **Figure 1. Boxplot of Expression Value Distribution.** The median expression values across all samples align horizontally, indicating that the data is well-normalized and comparable across arrays.

![Density Plot](https://github.com/user-attachments/assets/a44133d6-acff-462f-860b-93a31570a782)
 **Figure 2. Density Plot of Gene Expression.** The density curves for both Glioblastoma and Normal groups overlap almost perfectly. This confirms the absence of massive global technical biases, ensuring that downstream DEG analysis captures true biological variance.

![UMAP Plot](https://github.com/user-attachments/assets/ea85c20b-d712-4d0e-b941-c99108d73f39)
 **Figure 3. UMAP Dimensionality Reduction.** The UMAP plot reveals a perfect, distinct spatial separation between Glioblastoma samples (Salmon) and Normal Brain samples (Turquoise). This confirms that the cancer induces a fundamental, global shift in the brain's transcriptomic identity prior to strict DEG filtering.

**B. Transcriptomic Divergence (Volcano Plot)**
The analysis revealed a stark and massive molecular separation between the cancerous and normal tissues. 

![Volcano Plot](https://github.com/user-attachments/assets/1c810f66-2480-4191-8282-bd9fcdc28364)
 **Figure 4. Volcano Plot of Glioblastoma vs. Normal Brain Tissue.** * **Up-regulated (Red):** Represents hyperactive oncogenes driving rapid cell division and tissue invasion.
* **Down-regulated (Blue):** Represents normal neurological genes that are forcefully suppressed by the tumor, reflecting the loss of normal cognitive and sensory functions in the affected brain regions.

**C. Core Tumor Signature (Top 50 DEGs Heatmap)**

![Heatmap Glioblastoma](https://github.com/user-attachments/assets/a79b73af-7746-4b1a-851e-0e452e68c455)
 **Figure 5. Heatmap of the Top 50 DEGs.** The heatmap acts as the definitive molecular fingerprint of the cancer. A stark color demarcation exists between the Normal tissue block and the Glioblastoma block. The hierarchical clustering at the top confirms that these 50 genes alone are sufficient to perfectly distinguish a healthy pediatric brain from a terminally ill one.

**D. Functional Enrichment: The Modus Operandi**
To understand the functional consequence of these transcriptomic shifts, GO and KEGG pathway analyses were conducted on the significant DEGs.

![Gene Ontology Glioblastoma](https://github.com/user-attachments/assets/b88cf2e2-1ea7-424f-9519-1cdb7b934c86)
 **Figure 6. Gene Ontology (GO) Enrichment.** The biological processes represent a dual narrative: the hyperactivation of *Chromosome Segregation* and *Ribonucleoprotein Complex Biogenesis* reflects the tumor's unchecked cellular replication, while terms like *Regulation of Neuron Projection Development* and *Synaptic Vesicle Cycle* correspond to the massive disruption of normal neural pathways.

![KEGG Pathway Glioblastoma](https://github.com/user-attachments/assets/d5f5d109-e9c7-49f4-961f-2def8ac49ae3)
 **Figure 7. KEGG Pathway Enrichment.** The primary hijacked biochemical routes include the **Cell Cycle** and **MAPK signaling pathway**, allowing the tumor to grow relentlessly. Furthermore, the activation of the **Focal Adhesion** and **Axon Guidance** pathways highlights the cancer's invasive capability to navigate, anchor into, and degrade the surrounding healthy brain matrix.

### IV. Conclusion
This study concludes that Pediatric Glioblastoma Multiforme (Grade IV) possesses a transcriptomic profile that completely overhauls normal neurological functions. The genetic expression shifts in GBM are heavily concentrated on hijacking the *Cell Cycle*, *MAPK signaling*, and *Chromosome Segregation* pathways to facilitate lethal, unchecked proliferation, while simultaneously suppressing normal synaptic transmissions. Mapping these specific pathways provides a clearer understanding of the tumor's aggressiveness and highlights critical vulnerabilities for future targeted therapies.

### V. References
[1] Griesinger, A. M., et al. (2013). *Dataset GSE50161: Gene expression profiling of human brain tumors*. Gene Expression Omnibus (GEO).  
[2] Davis, S., & Meltzer, P. S. (2007). *GEOquery: a bridge between the Gene Expression Omnibus (GEO) and BioConductor*. Bioinformatics, 23(14), 1846-1847.  
[3] Ritchie, M. E., et al. (2015). *limma powers differential expression analyses for RNA-sequencing and microarray studies*. Nucleic Acids Research, 43(7), e47.  
[4] Wu, T., et al. (2021). *clusterProfiler 4.0: A universal enrichment tool for interpreting omics data*. The Innovation, 2(3), 100141.
