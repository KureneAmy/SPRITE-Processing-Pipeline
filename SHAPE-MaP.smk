configfile : "config.yaml"

import os

def get_seq_name(fasta_file):
    all_seq_names = []
    with open(fasta_file, 'r') as f:
        for line in f:
            if line.startswith('>'):
                seq_name = line[1:].strip().split()[0]
                all_seq_names.append(seq_name)
    return all_seq_names

TARGET_FA = config["target"]
ALL_SEQ_NAMES = get_seq_name(TARGET_FA)

SPLIT_DNA_DIR = os.path.join(config["output_dir"], "reference_split/DNA")
SPLIT_RNA_DIR = os.path.join(config["output_dir"], "reference_split/RNA")
os.makedirs(SPLIT_DNA_DIR, exist_ok=True)
os.makedirs(SPLIT_RNA_DIR, exist_ok=True)


rule all:
    input:
        expand("{output_dir}/{sample}/{seq_name}/{sample}_{seq_name}.shape", sample=config["samples"], output_dir=config["output_dir"], seq_name=ALL_SEQ_NAMES),
        expand("{output_dir}/{sample}/{seq_name}/{sample}_{seq_name}_profiles.pdf", sample=config["samples"], output_dir=config["output_dir"], seq_name=ALL_SEQ_NAMES),
        expand("{output_dir}/{sample}/{seq_name}/structure/{sample}_{seq_name}.dbn", sample=config["samples"], output_dir=config["output_dir"], seq_name=ALL_SEQ_NAMES),
        expand("{output_dir}/{sample}/{seq_name}/structure/{sample}_{seq_name}.ct", sample=config["samples"], output_dir=config["output_dir"], seq_name=ALL_SEQ_NAMES),
        expand("{output_dir}/{sample}/{seq_name}/structure/{sample}_{seq_name}_folding.svg", sample=config["samples"], output_dir=config["output_dir"], seq_name=ALL_SEQ_NAMES),
        expand("{output_dir}/multiqc_report.html", output_dir=config["output_dir"]),
        os.path.join(config["output_dir"], "reports", "report_data.json"),
        os.path.join(config["output_dir"], "reports", "SHAPE-MaP_Analysis_Report.html"),
        os.path.join(config["output_dir"], "reports", "SHAPE-MaP_Analysis_Report.md"),


rule split_reference_sequences:
    input:
        target = config["target"]
    output:
        dna_dir = directory(SPLIT_DNA_DIR),
        rna_dir = directory(SPLIT_RNA_DIR),
        dna_files = expand("{dna_dir}/{seq_name}.fa", dna_dir=SPLIT_DNA_DIR, seq_name=ALL_SEQ_NAMES),
        rna_files = expand("{rna_dir}/{seq_name}.fa", rna_dir=SPLIT_RNA_DIR, seq_name=ALL_SEQ_NAMES)
    params:
        script = config["scripts"]["split"]
    container:
        config["container"]
    log:
        os.path.join(config["output_dir"], "logs", "split_reference_sequences.log")
    shell:
        "python {params.script} -i {input.target} -d {output.dna_dir} -r {output.rna_dir} --overwrite > {log} 2>&1"


rule qc_modified:
    output:
        html_R1 = os.path.join(config["output_dir"], "qc", "{sample}", "modified", "{sample}_modified_R1_fastqc.html"),
        html_R2 = os.path.join(config["output_dir"], "qc", "{sample}", "modified", "{sample}_modified_R2_fastqc.html"),
    params:
        FQ1 = lambda wildcards: config["samples"][wildcards.sample]["modified"]["R1"],
        FQ2 = lambda wildcards: config["samples"][wildcards.sample]["modified"]["R2"],
        threads = 8
    container:
        config["container"]
    log:
        os.path.join(config["output_dir"], "logs", "qc_modified_{sample}.log")
    shell:
        """
        output_dir=$(dirname {output.html_R1})
        mkdir -p "$output_dir"
        
        fastqc -t {params.threads} {params.FQ1} {params.FQ2} -o "$output_dir" > {log} 2>&1
        
        base_fq1=$(basename {params.FQ1})
        base_fq1_noext="${{base_fq1%.*}}"
        base_fq1_noext="${{base_fq1_noext%.*}}"
        
        base_fq2=$(basename {params.FQ2})
        base_fq2_noext="${{base_fq2%.*}}"
        base_fq2_noext="${{base_fq2_noext%.*}}"

        mv "$output_dir/${{base_fq1_noext}}_fastqc.html" {output.html_R1} 2>/dev/null || true
        mv "$output_dir/${{base_fq2_noext}}_fastqc.html" {output.html_R2} 2>/dev/null || true
        """

rule qc_untreated:
    output:
        html_R1 = os.path.join(config["output_dir"], "qc", "{sample}", "untreated", "{sample}_untreated_R1_fastqc.html"),
        html_R2 = os.path.join(config["output_dir"], "qc", "{sample}", "untreated", "{sample}_untreated_R2_fastqc.html"),
    params:
        FQ1 = lambda wildcards: config["samples"][wildcards.sample]["untreated"]["R1"],
        FQ2 = lambda wildcards: config["samples"][wildcards.sample]["untreated"]["R2"],
        threads = 8
    container:
        config["container"]
    log:
        os.path.join(config["output_dir"], "logs", "qc_untreated_{sample}.log")
    shell:
        """
        output_dir=$(dirname {output.html_R1})
        mkdir -p "$output_dir"
        
        fastqc -t {params.threads} {params.FQ1} {params.FQ2} -o "$output_dir" > {log} 2>&1
        
        base_fq1=$(basename {params.FQ1})
        base_fq1_noext="${{base_fq1%.*}}"
        base_fq1_noext="${{base_fq1_noext%.*}}"
        
        base_fq2=$(basename {params.FQ2})
        base_fq2_noext="${{base_fq2%.*}}"
        base_fq2_noext="${{base_fq2_noext%.*}}"

        mv "$output_dir/${{base_fq1_noext}}_fastqc.html" {output.html_R1} 2>/dev/null || true
        mv "$output_dir/${{base_fq2_noext}}_fastqc.html" {output.html_R2} 2>/dev/null || true
        """

if config.get("denatured"):
    rule qc_denatured:
        output:
            html_R1 = os.path.join(config["output_dir"], "qc", "{sample}", "denatured", "{sample}_denatured_R1_fastqc.html"),
            html_R2 = os.path.join(config["output_dir"], "qc", "{sample}", "denatured", "{sample}_denatured_R2_fastqc.html"),
        params:
            FQ1 = lambda wildcards: config["samples"][wildcards.sample]["denatured"]["R1"],
            FQ2 = lambda wildcards: config["samples"][wildcards.sample]["denatured"]["R2"],
            threads = 8
        container:
            config["container"]
        log:
            os.path.join(config["output_dir"], "logs", "qc_denatured_{sample}.log")
        shell:
            """
            output_dir=$(dirname {output.html_R1})
            mkdir -p "$output_dir"
            
            fastqc -t {params.threads} {params.FQ1} {params.FQ2} -o "$output_dir" > {log} 2>&1
            
            base_fq1=$(basename {params.FQ1})
            base_fq1_noext="${{base_fq1%.*}}"
            base_fq1_noext="${{base_fq1_noext%.*}}"
            
            base_fq2=$(basename {params.FQ2})
            base_fq2_noext="${{base_fq2%.*}}"
            base_fq2_noext="${{base_fq2_noext%.*}}"

            mv "$output_dir/${{base_fq1_noext}}_fastqc.html" {output.html_R1} 2>/dev/null || true
            mv "$output_dir/${{base_fq2_noext}}_fastqc.html" {output.html_R2} 2>/dev/null || true
            """


def get_shapemapper_params(wildcards):
    params = [
        "--target", os.path.join(SPLIT_DNA_DIR, f"{wildcards.seq_name}.fa"),
        "--name", wildcards.sample,
        "--out", os.path.join(config["output_dir"], wildcards.sample, wildcards.seq_name),
        "--log", os.path.join(config["output_dir"], wildcards.sample, wildcards.seq_name, f"{wildcards.sample}_{wildcards.seq_name}_shapemapper_log.txt"),
        "--verbose",
        "--modified", "--R1", config["samples"][wildcards.sample]["modified"]["R1"], "--R2", config["samples"][wildcards.sample]["modified"]["R2"],
        "--untreated", "--R1", config["samples"][wildcards.sample]["untreated"]["R1"], "--R2", config["samples"][wildcards.sample]["untreated"]["R2"],
        "--output-parsed-mutations",
        "--output-counted-mutations",
        "--per-read-histograms",
    ]
    
    if config["denatured"]:
        params.extend(["--denatured", "--R1", config["samples"][wildcards.sample]["denatured"]["R1"], "--R2", config["samples"][wildcards.sample]["denatured"]["R2"]])
    if config["overwrite"]:
        params.extend(["--overwrite"])
    if config.get("min_depth"):
        params.extend(["--min-depth", str(config["min_depth"])])
    if config.get("max_bg"):
        params.extend(["--max-bg", str(config["max_bg"])])
    if config.get("min_qual_to_trim"):
        params.extend(["--min-qual-to-trim", str(config["min_qual_to_trim"])])
    if config.get("window_to_trim"):
        params.extend(["--window-to-trim", str(config["window_to_trim"])])
    if config.get("min_length_to_trim"):
        params.extend(["--min-length-to-trim", str(config["min_length_to_trim"])])
    if config.get("min_qual_to_count"):
        params.extend(["--min-qual-to-count", str(config["min_qual_to_count"])])
    if config["indiv_norm"]:
        params.extend(["--indiv-norm"])
    
    if config["amplicon"]:
        params.extend(["--amplicon"])
    if config["certain_primer"]:
        params.extend(["--primers", str(config["primers_file"])])
    if config["random_primer"] and config["amplicon"] == False and config["certain_primer"] == False:
        params.extend(["--random-primer-len", str(config["random_primer_length"])])
    if config["amplicon"] or config["certain_primer"]:
        if config.get("max_primer_offset"):
            params.extend(["--max-primer-offset", str(config["max_primer_offset"])])
    return " ".join(params)

rule run_shapemapper:
    input:
        target = lambda wildcards: os.path.join(SPLIT_DNA_DIR, f"{wildcards.seq_name}.fa"),
        modified_R1 = lambda wildcards: config["samples"][wildcards.sample]["modified"]["R1"],
        modified_R2 = lambda wildcards: config["samples"][wildcards.sample]["modified"]["R2"],
        untreated_R1 = lambda wildcards: config["samples"][wildcards.sample]["untreated"]["R1"],
        untreated_R2 = lambda wildcards: config["samples"][wildcards.sample]["untreated"]["R2"],
    output:
        shape = os.path.join(config["output_dir"], "{sample}", "{seq_name}", "{sample}_{seq_name}.shape"),
        profile = os.path.join(config["output_dir"], "{sample}", "{seq_name}", "{sample}_{seq_name}_profiles.pdf"),
        log = os.path.join(config["output_dir"], "{sample}", "{seq_name}", "{sample}_{seq_name}_shapemapper_log.txt"),
    params:
        cmd_params = get_shapemapper_params,
    threads: 8
    container:
        config["container"]
    log:
        os.path.join(config["output_dir"], "logs", "shapemapper_{sample}_{seq_name}.log")
    shell:
        "shapemapper {params.cmd_params} > {log} 2>&1"

def get_rnastructure_params(wildcards): 
    params = []
    if config.get("temperature"):
        params.extend(["--temperature", str(config["temperature"])])
    if config.get("maximum"):
        params.extend(["--maximum", str(config["maximum"])])
    if config.get("loop"):
        params.extend(["--loop", str(config["loop"])])
    return " ".join(params)

rule predict_structure:
    input:
        shape = os.path.join(config["output_dir"], "{sample}", "{seq_name}", "{sample}_{seq_name}.shape"),
        sequence = lambda wildcards: os.path.join(SPLIT_RNA_DIR, f"{wildcards.seq_name}.fa"),
    output:
        ct = os.path.join(config["output_dir"], "{sample}", "{seq_name}", "structure", "{sample}_{seq_name}.ct"),
    params:
        program_params = get_rnastructure_params
    container:
        config["container"]
    log:
        os.path.join(config["output_dir"], "logs", "rnastructure_{sample}_{seq_name}_ct.log")
    shell:
        "Fold --SHAPE {input.shape} {input.sequence} {output.ct} {params.program_params} > {log} 2>&1"

rule predict_structure_dbn:
    input:
        shape = os.path.join(config["output_dir"], "{sample}", "{seq_name}", "{sample}_{seq_name}.shape"),
        sequence = lambda wildcards: os.path.join(SPLIT_RNA_DIR, f"{wildcards.seq_name}.fa"),
    output:
        dbn = os.path.join(config["output_dir"], "{sample}", "{seq_name}", "structure", "{sample}_{seq_name}.dbn"),
    params:
        program_params = get_rnastructure_params
    container:
        config["container"]
    log:
        os.path.join(config["output_dir"], "logs", "rnastructure_{sample}_{seq_name}_dbn.log")
    shell:
        "Fold --SHAPE {input.shape} -k {input.sequence} {output.dbn} {params.program_params} > {log} 2>&1"

rule visualize_structure:
    input:
        shape = os.path.join(config["output_dir"], "{sample}", "{seq_name}", "{sample}_{seq_name}.shape"),
        ct = os.path.join(config["output_dir"], "{sample}", "{seq_name}", "structure", "{sample}_{seq_name}.ct"),
    output:
        svg = os.path.join(config["output_dir"], "{sample}", "{seq_name}", "structure", "{sample}_{seq_name}_folding.svg"),
    container:
        config["container"]
    log:
        os.path.join(config["output_dir"], "logs", "visualize_structure_{sample}_{seq_name}.log")
    shell:
        "draw --SHAPE {input.shape} --svg -n 1 {input.ct} {output.svg} > {log} 2>&1"


def get_all_fastqc_outputs(wildcards):
    all_reports = []

    for sample in config["samples"]:
        # modified
        r1_path = os.path.join(config["output_dir"], "qc", sample, "modified", f"{sample}_modified_R1_fastqc.html")
        r2_path = os.path.join(config["output_dir"], "qc", sample, "modified", f"{sample}_modified_R2_fastqc.html")
        all_reports.extend([r1_path, r2_path])

        # untreated
        r1_path = os.path.join(config["output_dir"], "qc", sample, "untreated", f"{sample}_untreated_R1_fastqc.html")
        r2_path = os.path.join(config["output_dir"], "qc", sample, "untreated", f"{sample}_untreated_R2_fastqc.html")
        all_reports.extend([r1_path, r2_path])
        
        # denatured
        if config.get("denatured"):
            r1_path = os.path.join(config["output_dir"], "qc", sample, "denatured", f"{sample}_denatured_R1_fastqc.html")
            r2_path = os.path.join(config["output_dir"], "qc", sample, "denatured", f"{sample}_denatured_R2_fastqc.html")
            all_reports.extend([r1_path, r2_path])

    for sample in config["samples"]:
        for seq_name in ALL_SEQ_NAMES:
            log_path = os.path.join(
                config["output_dir"], 
                sample, 
                seq_name, 
                f"{sample}_{seq_name}_shapemapper_log.txt"
            )
            all_reports.append(log_path)
    
    return all_reports


rule multiqc:
    input:
        get_all_fastqc_outputs
    output:
        html = os.path.join(config["output_dir"], "multiqc_report.html"),
    params:
        analysis_dir = config["output_dir"]
    container:
        config["container"]
    log:
        os.path.join(config["output_dir"], "logs", "multiqc.log")
    shell:
        """
        multiqc {params.analysis_dir} -o {params.analysis_dir} --force > {log} 2>&1
        """

REPORT_CONFIG_FILE = workflow.configfiles[0] if workflow.configfiles else "config.yaml"

rule collect_report_stats:
    input:
        multiqc = os.path.join(config["output_dir"], "multiqc_report.html"),
        svg = expand(
            os.path.join(config["output_dir"], "{sample}", "{seq_name}", "structure", "{sample}_{seq_name}_folding.svg"),
            sample=config["samples"],
            seq_name=ALL_SEQ_NAMES,
        ),
    output:
        json = os.path.join(config["output_dir"], "reports", "report_data.json"),
    params:
        config_file = REPORT_CONFIG_FILE,
        script = config["scripts"]["compile_report"],
        reports_dir = os.path.join(config["output_dir"], "reports"),
    log:
        os.path.join(config["output_dir"], "logs", "compile_report_stats.log"),
    shell:
        """
        python {params.script} \
            --config {params.config_file} \
            --output-dir {params.reports_dir} \
            --stats-json {output.json} \
            --no-html \
            --no-markdown > {log} 2>&1
        """


rule generate_html_report:
    input:
        json = os.path.join(config["output_dir"], "reports", "report_data.json"),
        template = config["report_config"]["html_template"],
    output:
        html = os.path.join(config["output_dir"], "reports", "SHAPE-MaP_Analysis_Report.html"),
    params:
        config_file = REPORT_CONFIG_FILE,
        script = config["scripts"]["compile_report"],
        reports_dir = os.path.join(config["output_dir"], "reports"),
    log:
        os.path.join(config["output_dir"], "logs", "compile_report_html.log"),
    shell:
        """
        python {params.script} \
            --config {params.config_file} \
            --output-dir {params.reports_dir} \
            --html-template {input.template} \
            --no-markdown > {log} 2>&1
        """


rule generate_markdown_report:
    input:
        json = os.path.join(config["output_dir"], "reports", "report_data.json"),
        template = config["report_config"]["markdown_template"],
    output:
        md = os.path.join(config["output_dir"], "reports", "SHAPE-MaP_Analysis_Report.md"),
    params:
        config_file = REPORT_CONFIG_FILE,
        script = config["scripts"]["compile_report"],
        reports_dir = os.path.join(config["output_dir"], "reports"),
    log:
        os.path.join(config["output_dir"], "logs", "compile_report_markdown.log"),
    shell:
        """
        python {params.script} \
            --config {params.config_file} \
            --output-dir {params.reports_dir} \
            --md-template {input.template} \
            --no-html > {log} 2>&1
        """