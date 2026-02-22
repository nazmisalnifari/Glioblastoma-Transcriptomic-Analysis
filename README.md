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

### III. Results and Interpretations

**A. Exploratory Data Analysis (EDA): The Global Shift**
Before gene-level analysis, the quality and distribution of the data were validated to prevent batch effects and ensure statistical reliability.

![Boxplot](https://github.com/user-attachments/assets/615763e4-6336-43e4-ae5d-5606d891d240)

 **Figure 1. Boxplot of Expression Value Distribution.** 
 
 ![Density Plot](https://github.com/user-attachments/assets/26f182d4-4ecc-4285-bf0a-db082c2b7bdc)
 
 **Figure 2. Density Plot of Gene Expression.**
 
 ![UMAP Plot](https://github.com/user-attachments/assets/a700ed9d-251f-4659-8c57-956de0243897)
 
 **Figure 3. UMAP Dimensionality Reduction.**  
 
 **Results:** The median expression values (Boxplot) and density curves (Density Plot) across all samples align perfectly, confirming successful data normalization. Crucially, the UMAP plot reveals a distinct spatial separation between Glioblastoma samples (Salmon) and Normal Brain samples (Turquoise).
* **Biological Interpretation & Impact:** This global separation proves that the cancer induces a fundamental reprogramming of the brain's transcriptomic identity. From a research perspective, this validates that the downstream DEG analysis will capture true pathological variance rather than technical noise.

**B. Transcriptomic Divergence (Volcano Plot & Heatmap)**

![Volcano Plot](https://github.com/user-attachments/assets/6aa2ac1d-a2b9-43de-98a1-363412b61ca1)

 **Figure 4. Volcano Plot of Glioblastoma vs. Normal Brain Tissue.**
 
 ![Heatmap Glioblastoma](https://github.com/user-attachments/assets/64e444e7-26e0-42ad-b1f8-1034e4d84f1b)
 
 **Figure 5. Heatmap of the Top 50 DEGs.** 
 
 **Results:** The Volcano plot displays a massive molecular divergence, with red dots indicating hyperactive oncogenes and blue dots representing suppressed normal genes. The Heatmap of the top 50 DEGs demonstrates a stark, absolute color demarcation between the Normal tissue block and the Glioblastoma block.
* **Biological Interpretation & Impact:** The severe downregulation of normal neural genes illustrates how the tumor forcefully shuts down healthy cognitive and sensory functions to prioritize its own survival. Conversely, the top upregulated genes represent the core drivers of tumor malignancy. For future clinical applications, these top 50 DEGs serve as highly promising diagnostic biomarkers. Identifying these specific gene signatures in clinical biopsies could accelerate the accurate grading of pediatric brain tumors.

**C. Functional Enrichment: The Modus Operandi (GO & KEGG)**

![Gene Ontology Glioblastoma](https://github.com/user-attachments/assets/537a9a87-959a-4fc7-97d3-f1e86f2b63d5)
 
 **Figure 6. Gene Ontology (GO) Enrichment.**
 
 ![KEGG Pathway Glioblastoma](https://github.com/user-attachments/assets/91f70e29-9afe-42b1-a324-0e41323c7ca0)
 
 **Figure 7. KEGG Pathway Enrichment.** 
 
 **Results:** GO analysis shows a hyperactivation in *Chromosome Segregation* and *Ribonucleoprotein Complex Biogenesis*, coupled with a suppression in *Synaptic Vesicle Cycle*. KEGG analysis reveals that the primary hijacked biochemical routes are the **Cell Cycle**, **MAPK signaling pathway**, and **Focal Adhesion**.
* **Biological Interpretation & Impact:** The enrichment data perfectly decodes the tumor's lethal behavior. The hijacking of the *Cell Cycle* and *MAPK signaling* pathways explains the tumor's relentless, unconstrained growth, while *Focal Adhesion* activation explains its aggressive ability to infiltrate and degrade surrounding healthy brain matter. 
* **Future Research Implications:** This finding is critical for precision medicine. Because we now know the exact pathways driving the cancer, future therapeutic developments can pivot away from generalized chemotherapy and focus on targeted therapies—such as CDK inhibitors to block the hijacked cell cycle, or specific MAPK pathway antagonists, directly disarming the tumor's primary weapons.

### IV. Conclusion
This study concludes that Pediatric Glioblastoma Multiforme (Grade IV) possesses a transcriptomic profile that completely overhauls normal neurological functions. The genetic expression shifts in GBM are heavily concentrated on hijacking the *Cell Cycle*, *MAPK signaling*, and *Chromosome Segregation* pathways to facilitate lethal, unchecked proliferation, while simultaneously suppressing normal synaptic transmissions. Mapping these specific pathways not only provides a molecular explanation for the tumor's extreme aggressiveness but also highlights critical vulnerabilities (such as the MAPK and Cell Cycle pathways) that can be exploited for future targeted therapies in pediatric neuro-oncology.

### V. References
[1] Griesinger, A. M., et al. (2013). *Dataset GSE50161: Gene expression profiling of human brain tumors*. Gene Expression Omnibus (GEO).  
[2] Davis, S., & Meltzer, P. S. (2007). *GEOquery: a bridge between the Gene Expression Omnibus (GEO) and BioConductor*. Bioinformatics, 23(14), 1846-1847.  
[3] Ritchie, M. E., et al. (2015). *limma powers differential expression analyses for RNA-sequencing and microarray studies*. Nucleic Acids Research, 43(7), e47.  
[4] Wu, T., et al. (2021). *clusterProfiler 4.0: A universal enrichment tool for interpreting omics data*. The Innovation, 2(3), 100141.
