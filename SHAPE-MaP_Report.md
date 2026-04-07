# {{ report_title }}

---

**Report Date:** {{ report_date }}  
**Analysis Date:** {{ analysis_date }}  
{% if institution %}**Institution:** {{ institution }}  
{% endif %}{% if pi_name %}**Principal Investigator:** {{ pi_name }}  
{% endif %}{% if project_id %}**Project ID:** {{ project_id }}  
{% endif %}**Output Directory:** `{{ output_dir }}`

---

## Table of Contents

1. [Project Information](#1-project-information)
2. [Analysis Overview](#2-analysis-overview)
3. [Data Quality Assessment](#3-data-quality-assessment)
4. [SHAPE Reactivity Analysis](#4-shape-reactivity-analysis)
5. [ShapeMapper Results Summary](#5-shapemapper-results-summary)
6. [RNA Structure Prediction](#6-rna-structure-prediction)
7. [QC Metrics](#7-qc-metrics)
8. [Methods](#8-methods)
9. [Contact Information](#9-contact-information)

---

## 1. Project Information

| Field              | Value                          |
|--------------------|--------------------------------|
| Report Title       | {{ report_title }}             |
| Report Date        | {{ report_date }}              |
| Analysis Date      | {{ analysis_date }}            |
| Institution        | {{ institution or "N/A" }}     |
| Principal Investigator | {{ pi_name or "N/A" }}     |
| Project ID         | {{ project_id or "N/A" }}      |
| Number of Samples  | {{ sample_count }}             |
| Output Directory   | `{{ output_dir }}`             |
| MultiQC Report     | `{{ multiqc_report }}`         |

---

## 2. Analysis Overview

This report summarises the results of a **SHAPE-MaP** (Selective 2′-Hydroxyl Acylation
analysed by Primer Extension and Mutational Profiling) experiment processed with
**ShapeMapper2** and the **EasyOmics SHAPE-MaP Processing Pipeline**.

A total of **{{ sample_count }}** sample(s) were processed. SHAPE reactivity values
were computed and used for RNA secondary-structure prediction with **RNAstructure**.

### Samples Analysed

{% for sample_name, sample_data in samples.items() %}
- **{{ sample_name }}** — {{ sample_data.targets | length }} target(s):
  {% for seq_name in sample_data.targets %}  - `{{ seq_name }}`
  {% endfor %}
{% endfor %}

---

## 3. Data Quality Assessment

{% for sample_name, sample_data in samples.items() %}
### Sample: {{ sample_name }}

{% for seq_name, tgt in sample_data.targets.items() %}
#### Target: {{ seq_name }}

**Read Processing Summary**

| Metric | Value |
|--------|-------|
{% if tgt.log_stats.paired_reads_total is defined %}| Total Input Pairs | {{ "{:,}".format(tgt.log_stats.paired_reads_total) }} |
{% endif %}{% if tgt.log_stats.paired_reads_merged is defined %}| Pairs Merged (BBmerge) | {{ "{:,}".format(tgt.log_stats.paired_reads_merged) }}{% if tgt.log_stats.merge_rate is defined %} ({{ "%.1f"|format(tgt.log_stats.merge_rate) }}%){% endif %} |
{% endif %}{% if tgt.log_stats.alignment_rate is defined %}| Overall Alignment Rate | {{ "%.2f"|format(tgt.log_stats.alignment_rate) }}% |
{% endif %}{% if tgt.log_stats.mapped_reads_concordant_1 is defined %}| Concordantly Mapped Reads (1×) | {{ "{:,}".format(tgt.log_stats.mapped_reads_concordant_1) }} |
{% endif %}{% if tgt.log_stats.mean_read_depth is defined %}| Mean Read Depth | {{ "%.1f"|format(tgt.log_stats.mean_read_depth) }} |
{% endif %}{% if tgt.log_stats.quality_filtered_reads is defined %}| Quality-Filtered Reads | {{ "{:,}".format(tgt.log_stats.quality_filtered_reads) }} |
{% endif %}{% if not tgt.log_stats %}| (no log data available) | — |
{% endif %}

**Per-Channel Quality Control**

| Channel | Total Read Depth | Total Mutations | Median Mutation Rate | Mean Read Depth |
|---------|-----------------|-----------------|---------------------|-----------------|
{% for label in ["modified", "untreated", "denatured"] %}{% set depth_key = "total_read_depth_" ~ label %}{% set mut_key = "total_mutations_" ~ label %}{% set rate_key = "mutation_rate_" ~ label %}{% set mean_key = "mean_read_depth_" ~ label %}{% if tgt.profile_stats[depth_key] is defined or tgt.log_stats[rate_key] is defined %}| {{ label | capitalize }} | {{ "{:,}".format(tgt.profile_stats[depth_key]) if tgt.profile_stats[depth_key] is defined else "N/A" }} | {{ "{:,}".format(tgt.profile_stats[mut_key]) if tgt.profile_stats[mut_key] is defined else "N/A" }} | {{ "%.4f"|format(tgt.log_stats[rate_key]) if tgt.log_stats[rate_key] is defined else ("%.4f"|format(tgt.profile_stats[rate_key]) if tgt.profile_stats[rate_key] is defined else "N/A") }} | {{ "%.1f"|format(tgt.profile_stats[mean_key]) if tgt.profile_stats[mean_key] is defined else "N/A" }} |
{% endif %}{% endfor %}{% if not tgt.profile_stats %}| (profile data not available) | — | — | — | — |
{% endif %}

{% endfor %}
{% endfor %}

---

## 4. SHAPE Reactivity Analysis

{% for sample_name, sample_data in samples.items() %}
### Sample: {{ sample_name }}

{% for seq_name, tgt in sample_data.targets.items() %}
#### Target: {{ seq_name }}

**Mutation Rates by Channel**

| Channel | Median Mutation Rate | Total Mutations | Total Read Depth |
|---------|---------------------|-----------------|-----------------|
{% for label in ["modified", "untreated", "denatured"] %}{% set rate_key = "mutation_rate_" ~ label %}{% set mut_key = "total_mutations_" ~ label %}{% set depth_key = "total_read_depth_" ~ label %}{% set rate_val = tgt.log_stats[rate_key] if tgt.log_stats[rate_key] is defined else tgt.profile_stats[rate_key] %}{% if rate_val is defined or tgt.profile_stats[depth_key] is defined %}| {{ label | capitalize }} | {{ "%.4f"|format(rate_val) if rate_val is defined else "N/A" }} | {{ "{:,}".format(tgt.profile_stats[mut_key]) if tgt.profile_stats[mut_key] is defined else "N/A" }} | {{ "{:,}".format(tgt.profile_stats[depth_key]) if tgt.profile_stats[depth_key] is defined else "N/A" }} |
{% endif %}{% endfor %}{% if tgt.log_stats.high_quality_profile_pct is defined %}| High-Quality Profile (%) | {{ "%.1f"|format(tgt.log_stats.high_quality_profile_pct) }}% | — | — |
{% elif tgt.profile_stats.high_quality_profile_pct is defined %}| High-Quality Profile (%) | {{ "%.1f"|format(tgt.profile_stats.high_quality_profile_pct) }}% | — | — |
{% endif %}{% if not (tgt.log_stats.mutation_rate_modified is defined or tgt.log_stats.mutation_rate_untreated is defined or tgt.profile_stats.mutation_rate_modified is defined or tgt.profile_stats.mutation_rate_untreated is defined) %}| (no mutation rate data available) | — | — | — |
{% endif %}

**Reactivity Summary**

{% if tgt.profile_stats %}
| Metric | Value |
|--------|-------|
| Nucleotides with Valid Reactivity | {{ "{:,}".format(tgt.profile_stats.nucleotide_count) }} |
| Mean Reactivity | {{ "%.4f"|format(tgt.profile_stats.reactivity_mean) }} |
| Median Reactivity | {{ "%.4f"|format(tgt.profile_stats.reactivity_median) }} |
| Std Dev | {{ "%.4f"|format(tgt.profile_stats.reactivity_stdev) }} |
| Maximum Reactivity | {{ "%.4f"|format(tgt.profile_stats.reactivity_max) }} |
| Reactive Fraction (> 0.4) | {{ "%.1f"|format(tgt.profile_stats.reactive_fraction * 100) }}% |
{% else %}
*SHAPE profile file not available.*
{% endif %}

{% endfor %}
{% endfor %}

---

## 5. ShapeMapper Results Summary

| Sample | Target | Alignment Rate | Mean Depth | Mut. Rate (mod.) | Mut. Rate (untr.) | Reactive Fraction |
|--------|--------|---------------|------------|------------------|-------------------|-------------------|
{% for sample_name, sample_data in samples.items() %}{% for seq_name, tgt in sample_data.targets.items() %}| {{ sample_name }} | {{ seq_name }} | {{ "%.2f"|format(tgt.log_stats.alignment_rate) ~ "%" if tgt.log_stats.alignment_rate is defined else "N/A" }} | {{ "%.1f"|format(tgt.log_stats.mean_read_depth) if tgt.log_stats.mean_read_depth is defined else "N/A" }} | {{ "%.4f"|format(tgt.log_stats.mutation_rate_modified) if tgt.log_stats.mutation_rate_modified is defined else "N/A" }} | {{ "%.4f"|format(tgt.log_stats.mutation_rate_untreated) if tgt.log_stats.mutation_rate_untreated is defined else "N/A" }} | {{ "%.1f"|format(tgt.profile_stats.reactive_fraction * 100) ~ "%" if tgt.profile_stats.reactive_fraction is defined else "N/A" }} |
{% endfor %}{% endfor %}{% if not samples %}| (no data) | — | — | — | — | — | — |
{% endif %}

---

## 6. RNA Structure Prediction

RNA secondary structures were predicted using **RNAstructure Fold** with
SHAPE-directed folding constraints.

### Method

SHAPE reactivity values from each sample were incorporated as pseudo-free-energy
constraints to guide the thermodynamic folding algorithm. The minimum free energy (MFE)
structure is reported.

### Output Files

{% for sample_name, sample_data in samples.items() %}
#### Sample: {{ sample_name }}

{% for seq_name, tgt in sample_data.targets.items() %}
**Target: {{ seq_name }}**

- CT file: `{{ output_dir }}/{{ sample_name }}/{{ seq_name }}/structure/{{ sample_name }}_{{ seq_name }}.ct`
- DBN file: `{{ output_dir }}/{{ sample_name }}/{{ seq_name }}/structure/{{ sample_name }}_{{ seq_name }}.dbn`
- SVG visualisation: `{{ output_dir }}/{{ sample_name }}/{{ seq_name }}/structure/{{ sample_name }}_{{ seq_name }}_folding.svg`

{% endfor %}
{% endfor %}

---

## 7. QC Metrics

### MultiQC Report

A comprehensive MultiQC report aggregating FastQC results for all samples is available at:

```
{{ multiqc_report }}
```

### Sample QC Summary

| Sample | Target | Total Input Pairs | Pairs Merged | Alignment Rate | Mean Read Depth | Status |
|--------|--------|-------------------|--------------|----------------|-----------------|--------|
{% for sample_name, sample_data in samples.items() %}{% for seq_name, tgt in sample_data.targets.items() %}{% set align = tgt.log_stats.alignment_rate %}| {{ sample_name }} | {{ seq_name }} | {{ "{:,}".format(tgt.log_stats.paired_reads_total) if tgt.log_stats.paired_reads_total is defined else "N/A" }} | {{ "{:,}".format(tgt.log_stats.paired_reads_merged) if tgt.log_stats.paired_reads_merged is defined else "N/A" }}{% if tgt.log_stats.merge_rate is defined %} ({{ "%.1f"|format(tgt.log_stats.merge_rate) }}%){% endif %} | {{ "%.2f"|format(align) ~ "%" if align is defined else "N/A" }} | {{ "%.1f"|format(tgt.log_stats.mean_read_depth) if tgt.log_stats.mean_read_depth is defined else "N/A" }} | {% if align is defined %}{% if align >= 70 %}PASS{% elif align >= 40 %}WARN{% else %}LOW{% endif %}{% else %}N/A{% endif %} |
{% endfor %}{% endfor %}{% if not samples %}| (no data) | — | — | — | — | — | — |
{% endif %}

---

## 8. Methods

### Analysis Pipeline

SHAPE-MaP experiments were processed using the **EasyOmics SHAPE-MaP Processing Pipeline**
based on Snakemake. The pipeline encompasses the following major steps:

1. **Quality Control** – Raw reads were assessed with FastQC and aggregated with MultiQC.
2. **Reference Splitting** – Multi-target reference FASTA files were split into individual sequences.
3. **ShapeMapper2** – Modified, untreated (and optionally denatured) read pairs were jointly
   processed to derive per-nucleotide SHAPE reactivity values.
4. **RNA Structure Prediction** – SHAPE-constrained secondary-structure prediction was performed
   with RNAstructure `Fold`.
5. **Structure Visualisation** – Arc diagrams in SVG format were generated with RNAstructure `draw`.

### Software

| Software | Version | Purpose |
|----------|---------|---------|
| Snakemake | — | Workflow management |
| FastQC | — | Read-level quality control |
| MultiQC | — | QC report aggregation |
| BBMerge | {{ software_versions.bbmerge if software_versions.bbmerge else "—" }} | Paired-end read merging |
| Bowtie2 | {{ software_versions.bowtie2 if software_versions.bowtie2 else "—" }} | Read alignment |
| ShapeMapper2 | {{ software_versions.shapemapper if software_versions.shapemapper else "—" }} | SHAPE reactivity computation |
| RNAstructure | — | RNA secondary-structure prediction |

*Version numbers are extracted automatically from ShapeMapper log files where available.*

### References

1. Busan S, Weeks KM. Accurate detection of chemical modifications in RNA by mutational profiling
   (MaP) with ShapeMapper 2. *RNA*. 2018;24(2):143–148.
2. Reuter JS, Mathews DH. RNAstructure: software for RNA secondary structure prediction and
   analysis. *BMC Bioinformatics*. 2010;11:129.
3. Ewels P, et al. MultiQC: summarize analysis results for multiple tools and samples in a single
   report. *Bioinformatics*. 2016;32(19):3047–3048.

---

## 9. Contact Information

For questions regarding this analysis, please contact:

- **Institution:** {{ institution or "Please update in config.yaml" }}
- **Principal Investigator:** {{ pi_name or "Please update in config.yaml" }}
- **Pipeline:** EasyOmics SHAPE-MaP Processing Pipeline

---

*This report was automatically generated by `scripts/compile_report.py` on {{ report_date }}.*
