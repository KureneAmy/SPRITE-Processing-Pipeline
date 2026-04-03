# {{ report_title }}

---

| Field | Value |
|-------|-------|
| **Project ID** | {{ project_id }} |
| **Reference Genome** | {{ assembly }} |
| **Report Date** | {{ analysis_date }} |
| **Institution** | {{ institution }} |
| **PI / Investigator** | {{ pi_name }} |
| **Pipeline** | SPRITE v1.0 |
| **Number of Samples** | {{ num_samples }} |
| **Samples** | {{ samples_list }} |

---

## Table of Contents

1. [Executive Summary](#executive-summary)
2. [Data Quality Assessment](#data-quality-assessment)
3. [Barcode Analysis](#barcode-analysis)
4. [Cluster Analysis](#cluster-analysis)
5. [3D Chromatin Interactions](#3d-chromatin-interactions)
6. [Genome-wide Contact Map](#genome-wide-contact-map)
7. [QC Metrics](#qc-metrics)
8. [Methods](#methods)

---

## Executive Summary

This report summarises the results of SPRITE (Split-Pool Recognition of
Interactions by Tag Extension) sequencing data processed through the EasyOmics
SPRITE pipeline. The analysis detects and quantifies 3D chromatin interactions
across the genome using multi-way contact information encoded in SPRITE clusters.

**Key highlights:**

- **{{ num_samples }}** sample(s) analysed
- Reference genome: **{{ assembly }}**
- Interaction type: **DNA-DNA multi-way contacts**
- Analysis pipeline: **EasyOmics SPRITE v1.0**

---

## Data Quality Assessment

### Ligation Efficiency

{{ ligation_efficiency_table }}

### Sample Summary

{{ sample_summary_table }}

> **Note:** Total clusters and contact counts depend on sequencing depth and
> library complexity. Mean and median cluster sizes reflect the number of
> genomic loci co-captured in each SPRITE complex.

---

## Barcode Analysis

Barcodes were identified using `BarcodeIdentification_v1.2.0.jar`. Reads with
incomplete barcodes were removed before alignment and cluster construction.

Key steps:

- **Full-barcode filtering:** Reads lacking complete tag sequences excluded
  via `get_full_barcodes.py`
- **DPM adapter trimming:** Applied with Cutadapt (DPM sequence `GATCGGAAGAG`)
- **Anchor DPM removal:** Left-anchored `GGTGGTCTT` sequences trimmed

---

## Cluster Analysis

SPRITE clusters represent groups of genomic loci that were co-captured in the
same split-pool reaction, indicating spatial co-localisation in the nucleus.

### Cluster Size Distribution

Cluster size distribution plots are available at:

```
workup/clusters/cluster_sizes.pdf
workup/clusters/cluster_sizes.png
```

Clusters were filtered to retain sizes between the configured minimum and
maximum thresholds before contact matrix construction.

---

## 3D Chromatin Interactions

### Contact Matrix Construction

DNA-DNA contacts were derived from SPRITE clusters using
`get_sprite_contacts.py`. Each pair of genomic loci co-occurring within a
cluster contributes to the contact frequency matrix.

- **Resolution:** 1 Mb bins (configurable)
- **Downweighting:** `two_over_n` (configurable)
- **ICE normalisation:** HiCorrector v1.2 (`{{ ice_iterations | default(100) }}` iterations)

### Available Matrix Files (per sample)

| File | Description |
|------|-------------|
| `{sample}.DNA.raw.txt` | Raw contact frequency matrix |
| `{sample}.DNA.iced.txt` | ICE-normalised contact matrix |
| `{sample}.DNA.bias.txt` | ICE bias correction values |
| `{sample}.DNA.final.txt` | Final processed contact matrix |
| `{sample}.DNA.final.png` | Contact heatmap image |

---

## Genome-wide Contact Map

Genome-wide contact heatmaps are generated for each sample at 1 Mb resolution.
High-contact regions correspond to TADs (Topologically Associating Domains),
compartments, and known chromatin hubs.

Heatmap images are located at:

```
workup/heatmap/{sample}.DNA.final.pdf
workup/heatmap/{sample}.DNA.final.png
```

---

## QC Metrics

A comprehensive MultiQC report covering all pipeline steps is available at:

```
workup/qc/multiqc_report.html
```

This report includes:

- **Trim Galore** — adapter content and read quality before/after trimming
- **FastQC** — per-base quality scores, GC content, duplication levels
- **Bowtie2** — overall alignment rates and multi-mapping statistics
- **Cutadapt** — DPM adapter trimming statistics

---

## Methods

### Pipeline Overview

The SPRITE pipeline processes paired-end sequencing reads as follows:

1. **Adapter trimming** with Trim Galore (quality ≥ 20)
2. **Barcode identification** with BarcodeIdentification v1.2.0
3. **Full-barcode filtering** with `get_full_barcodes.py`
4. **DPM trimming** with Cutadapt
5. **Alignment** to `{{ assembly }}` with Bowtie2 (MAPQ ≥ 20)
6. **Chromosome annotation** converted from Ensembl to UCSC format
7. **Repeat masking** with bedtools (blacklist + RepeatMasker elements)
8. **Cluster construction** with `get_clusters.py`
9. **Contact matrix** generation with `get_sprite_contacts.py`
10. **ICE normalisation** with HiCorrector v1.2
11. **Heatmap visualisation** with R / ggplot2

### Software Versions

| Software | Version | Purpose |
|----------|---------|---------|
| Snakemake | ≥ 5.x | Workflow management |
| Trim Galore | ≥ 0.6 | Adapter trimming |
| Cutadapt | ≥ 3.x | DPM adapter removal |
| Bowtie2 | ≥ 2.4 | Short-read alignment |
| SAMtools | ≥ 1.12 | BAM file processing |
| bedtools | ≥ 2.29 | Repeat element masking |
| HiCorrector | 1.2 | ICE normalisation |
| MultiQC | ≥ 1.11 | QC report aggregation |
| R | ≥ 4.0 | Visualisation |
| Python | ≥ 3.6 | Data processing |

### References

1. Quinodoz SA, et al. (2018). Higher-Order Inter-chromosomal Hubs Shape 3D
   Genome Organization in the Nucleus. *Cell* 174(3):744–757.e24.

2. Langmead B & Salzberg SL (2012). Fast gapped-read alignment with Bowtie 2.
   *Nature Methods* 9:357–359.

3. Servant N, et al. (2015). HiC-Pro: an optimized and flexible pipeline for
   Hi-C data processing. *Genome Biology* 16:259.

---

## Contact Information

For questions regarding this report or the EasyOmics SPRITE pipeline, please
contact:

- **Institution:** {{ institution }}
- **PI / Investigator:** {{ pi_name }}
- **Project ID:** {{ project_id }}

---

*Report generated by EasyOmics SPRITE Pipeline on {{ analysis_date }}.*
