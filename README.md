# SPRITE Processing Pipeline
The analysis pipeline for SPRITE (Split-Pool Recognition of Interactions by Tag Extension) profiling assays, processes raw FASTQ data through sequential steps such as identifying barcodes, genome alignment, discard alignments that don't meet certain criteria, group alignments into clusters, and create heatmaps from clusters, and it supports multiple samples. Also, we provide a fully containerized Singularity environment that bundles all required tools and dependencies, and with a single command, the entire workflow can be executed reproducibly on any compatible system, supporting multiple samples.

# Part I Workflow
Here stands an throughout workflow of data analysis.

<img width="1683" height="347" alt="1 workflow" src="https://github.com/user-attachments/assets/f2e44011-9d84-4efe-8faa-8c34cdc2505d" />

# Part II Requirements
1.  **Recommended System Configuration**:

      * 8-core CPU
      * 24 GB RAM

2.  **Singularity**: Must be installed on your system. Below are the detailed steps for installing on an Ubuntu 22.0.4 system. For other operating systems, please refer to the official installation guide: [https://docs.sylabs.io/guides/3.0/user-guide/installation.html](https://docs.sylabs.io/guides/3.0/user-guide/installation.html)

      * **Step 1: Install System Dependencies**

        ```bash
        # Update package lists and install dependencies
        sudo apt-get update
        sudo apt-get install -y \
            build-essential \
            libseccomp-dev \
			libfuse3-dev \
            pkg-config \
            squashfs-tools \
            cryptsetup \
            curl wget git
        ```

      * **Step 2: Install Go Language**

        ```bash
        # Download and install Go
        wget https://go.dev/dl/go1.21.3.linux-amd64.tar.gz
        sudo tar -C /usr/local -xzvf go1.21.3.linux-amd64.tar.gz
        rm go1.21.3.linux-amd64.tar.gz

        # Configure Go environment variables and apply them
        echo 'export GOPATH=${HOME}/go' >> ~/.bashrc
        echo 'export PATH=/usr/local/go/bin:${PATH}:${GOPATH}/bin' >> ~/.bashrc
        source ~/.bashrc
        ```

      * **Step 3: Download, Build, and Install Singularity**

        ```bash
        # Note: The script navigates to /mnt/share/software. 
        # You can change this to your preferred directory for source code.
        cd /mnt/share/software

        # Download the Singularity CE source code
        wget https://github.com/sylabs/singularity/releases/download/v4.0.1/singularity-ce-4.0.1.tar.gz

        # Extract the archive and clean up
        tar -xvzf singularity-ce-4.0.1.tar.gz
        rm singularity-ce-4.0.1.tar.gz
        cd singularity-ce-4.0.1

        # Configure the build
        ./mconfig

        # Build Singularity (this can be time-consuming)
        cd builddir
        make

        # Install Singularity to the system
        sudo make install
        ```

      * **Step 4: Verify the Installation**

        ```bash
        # Check the installed version
        singularity --version

        # Display help information
        singularity -h
        ```


3.  **snakemake**: Snakemake must be installed on your system and requires a Python 3 distribution.

      ```bash
      pip install snakemake
      ```

4.  **Reference Data**: A directory containing bowtie2 index (Below are the detailed steps for the human hg38 genome. For other reference genomes, please download the corresponding files and replace them as needed).

      ```bash
      mkdir index
      cd index

      # Download Genome FASTA
      wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_48/GRCh38.primary_assembly.genome.fa.gz

      # Unzip the files
      gunzip GRCh38.primary_assembly.genome.fa.gz
      
      # Build the main Genome Index
      singularity exec ../SPRITE.sif bowtie2-build --threads 16 GRCh38.primary_assembly.genome.fa GRCh38.primary_assembly.genome

      # Remove unnecessary files
      rm GRCh38.primary_assembly.genome.fa
      ```

5.  **Data Preparation**: The data run by this pipeline is from SRR7216005 in the SRA database.The specific processing method is as follows

    ```bash
    # Download the test sra data
    mkdir -p samples
    cd samples
    prefetch SRR7216005

    # Convert sra data to fastq data
    fastq-dump --split-files SRR7216005\SRR7216005.sra
    
    # Randomly sample fastq data
    singularity exec ../SPRITE.sif seqtk sample -s100 SRR7216005_1.fastq 53086 > SRR7216005_R1.fastq
    singularity exec ../SPRITE.sif seqtk sample -s100 SRR7216005_2.fastq 53086 > SRR7216005_R2.fastq
    pigz -p 8 SRR7216005_R1.fastq
    pigz -p 8 SRR7216005_R2.fastq
    ```

6.  **Sample Summary**: The location, file name and file pairs are summarised in json format by running this.

      ```bash
      singularity exec SPRITE.sif python fastq2json.py --fastq_dir samples
      ```

7.   **Required File Structure**

      ```
      root/
        ├── config.txt
        ├── config.yaml
        ├── dpm96.fasta
        ├── fastq2json.py
        ├── hg38_blacklist_rmsk.milliDivLessThan140.bed.gz
        ├── index
        │   ├── GRCh38.primary_assembly.genome.1.bt2
        │   ├── GRCh38.primary_assembly.genome.2.bt2
        │   ├── GRCh38.primary_assembly.genome.3.bt2
        │   ├── GRCh38.primary_assembly.genome.4.bt2
        │   ├── GRCh38.primary_assembly.genome.rev.1.bt2
        │   └── GRCh38.primary_assembly.genome.rev.2.bt2
        ├── samples
        │   ├── SRR7216005_R1.fastq.gz
        │   └── SRR7216005_R2.fastq.gz
        ├── samples.json
        ├── scripts
        │   ├── HiCorrector_1.2
        │   │   ├── bin
        │   │   │   ├── export_norm_data
        │   │   │   ├── ic
        │   │   │   ├── ic_mep
        │   │   │   ├── ic_mes
        │   │   │   ├── split_data
        │   │   │   └── split_data_parallel
        │   │   ├── example
        │   │   │   ├── contact.matrix
        │   │   │   ├── run_example_ic_mep.sh
        │   │   │   ├── run_example_ic_mes.sh
        │   │   │   ├── run_example_ic.sh
        │   │   │   └── run_export_norm_data.sh
        │   │   ├── Manual_HiCorrector_1.2.pdf
        │   │   └── src
        │   │       ├── export_norm_data.c
        │   │       ├── ic.c
        │   │       ├── ic_mep.c
        │   │       ├── ic_mes.c
        │   │       ├── Makefile
        │   │       ├── matutils.c
        │   │       ├── matutils.h
        │   │       ├── my_getopt-1.5
        │   │       │   ├── ChangeLog
        │   │       │   ├── getopt.3
        │   │       │   ├── getopt.h
        │   │       │   ├── getopt.txt
        │   │       │   ├── LICENSE
        │   │       │   ├── main.c
        │   │       │   ├── Makefile
        │   │       │   ├── my_getopt.c
        │   │       │   ├── my_getopt.h
        │   │       │   └── README
        │   │       ├── split_data.c
        │   │       └── split_data_parallel.c
        │   ├── images
        │   │   └── tag_layout.png
        │   ├── java
        │   │   ├── BarcodeIdentification_v1.2.0.jar
        │   │   ├── example_config.txt
        │   │   └── src
        │   │       └── edu
        │   │           └── caltech
        │   │               └── lncrna
        │   │                   └── barcode
        │   │                       └── core
        │   │                           ├── BarcodeIdentification.java
        │   │                           ├── Origin.java
        │   │                           ├── Read.java
        │   │                           ├── TagCategory.java
        │   │                           └── Tag.java
        │   ├── python
        │   │   ├── assembly.py
        │   │   ├── cluster.py
        │   │   ├── contact.py
        │   │   ├── convert_stripmask_to_bed.py
        │   │   ├── count_mouse_human.py
        │   │   ├── ensembl2ucsc.py
        │   │   ├── filter_bam_by_edit_distance.py
        │   │   ├── get_aiden_hic_contacts.py
        │   │   ├── get_clusters.py
        │   │   ├── get_full_barcodes.py
        │   │   ├── get_hub_contacts.py
        │   │   ├── get_ligation_efficiency.py
        │   │   ├── get_ren_hic_contacts.py
        │   │   ├── get_sprite_contacts_johnbot.py
        │   │   ├── get_sprite_contacts.py
        │   │   ├── merge_two_clusters.py
        │   │   └── __pycache__
        │   │       ├── assembly.cpython-39.pyc
        │   │       ├── cluster.cpython-39.pyc
        │   │       └── contact.cpython-39.pyc
        │   └── r
        │       ├── get_cluster_size_distribution.r
        │       ├── plot_chromosome_coverage.r
        │       └── plot_heatmap.R
        ├── SPRITE.sif
        └── SPRITE.smk

      ```
      
      - **SPRITE.smk** — The main Snakemake workflow script.  
      - **config.yaml** — Configuration file containing paths, parameters, and sample information.  
        ⚠️ Must be located in the same directory as `SPRITE.smk`.
      - **SPRITE.sif** — Singularity container image with all required software and dependencies pre-installed.
      - **dpm96.fasta** — FASTA file containing barcode sequences; replace with your own if needed. 
      - **config.txt** — Config file for barcode id; replace with your own if needed. 
      - **hg38_blacklist_rmsk.milliDivLessThan140.bed.gz** — Blacklist of sequences for filtering DNA contacts; replace with your own if needed.
      - **samples.json** — Json file that summerizes the information of samples; run fastq2json.py to get this. 
      - **index** — Reference genome index for Bowtie; replace with your preferred reference.

# Part III Running

   * **Example code**

      * **Step 1: Edit `config.yaml`**

        ```bash
        #Location of the container image SPRITE.sif
        container: "SPRITE.sif"
        #Location of the config file for barcodeIdentification
        bID: "./config.txt"
        #Location of the samples json file produced with fastq2json.py script
        samples: "./samples.json"
        #output directory
        output_dir: ""
        #root directory for the scripts
        scripts_dir: "./scripts"
        #Bowtie2 assembly type
        assembly: "hg38"
        #Number of barcodes used
        num_tags: "5"
        #Repeat mask used for filtering DNA contacts
        mask:
            hg38: "hg38_blacklist_rmsk.milliDivLessThan140.bed.gz"
        #Bowtie2 indexes location with prefix
        bowtie2_index:
            hg38: "./index/GRCh38.primary_assembly.genome"
        #Setting for mating heatmap matrix
        #Plot a chromosome 'chr1' or plot genome wide 'genome'
        chromosome:
            - genome
        # Options 'none', 'n_minus_one', 'two_over_n'
        downweighting:
            - two_over_n
        #Number of ICE iterations 
        ice_iterations:
            - 100
        #Contact matrix resolution in bp
        resolution:
            - 1000000
        min_cluster_size:
            - 2
        max_cluster_size:
            - 1000
        #Max value for heatmap plotting
        max_value:
            - 255
        ```

      * **Step 2: run snakemake**
      Here /mnt/zhangam/SPRITE/ represents the root directory.

        ```bash
        snakemake \
        --snakefile SPRITE.smk \
        --use-singularity \
        --cores 8 \
        --singularity-args "--bind /mnt/zhangam/SPRITE/:/root/" \
        --configfile config.yaml
        ```
      
   * **Command Parameters**

      **edit `config.yaml`**
      - `container`:(required) Path to the Singularity container image (`SPRITE.sif`) containing all required software and dependencies.
      - `bID`:        (required)Location of the config file for barcode identification.
      - `samples`:    (required) Location of the samples json file produced with fastq2json.py.
      - `output_dir`:          (optional) Path to the directory where the output will be stored. After setting your output can be found in "{output_dir}/workup/". Default to working directory as output directory.
      - `scripts_dir`:          (optional) Path to the directory where the scripts will be stored. After setting your scripts should be found in "{scripts_dir}". Default to working directory as scripts directory.
      - `assembly`:       (required) Type of assembly you use. Must correspond to the mask and bowtie2_index subsequently.
      - `num_tags`:          (optional) Number of barcodes used. Default: 5.
      - `mask`:                (required)Repeat mask used for filtering DNA contacts.
      - `bowtie2_index`:            (required) Bowtie2 indexes location with prefix.
      - `chromosome`: (optional) For intrachromosomal heatmaps, one of "chr1", "chr2", ..., "chrX". For interchromosomal heatmaps, "genome". Default: genome.
      - `downweighting`: The downweighting strategy, one of "none", "n_minus_one", and "two_over_n". Default: two_over_n
        * none: Each contact contributes a value of 1 to the contact matrix in step 2. For example, a contact from a 5-cluster will contribute 1.
        * n_minus_one: Each contact contributes a value of 1/(N - 1), where N is the size of the corresponding cluster. For example, a contact from a 5-cluster will contribute 1/4.
        * two_over_n: Each contact contributes a value of 2/N, where N is the size of the corresponding cluster. For example, a contact from a 5-cluster will contribute 2/5.
      - `ice_iterations`: The number of iterations to perform when running Hi-Corrector. Default: 100.
      - `resolution`: The binning resolution in nt. Default: 1000000 (1 Mb).
      - `min_cluster_size` and `max_cluster_size`: Ignore clusters that fall outside these parameters, i.e, clusters that are too big or too small. Default: 2 to 1000.
      - `max_value`: Max value for heatmap plotting that defines the cutoff of the heatmap. Default: 255.


      **run snakemake**
      - `--use-singularity`: Enables execution of rules within a Singularity container to ensure a fully reproducible environment.
      - `--singularity-args`: Allows passing additional arguments to the Singularity runtime (e.g., `--bind`, `--nv`, or custom options).
      - `--cores`: Specifies the maximum number of CPU cores (threads) that Snakemake can use in parallel when executing workflow rules.
      - `--bind`: Specifies the directories to be mounted within the Singularity container. Include all required paths such as raw data, scripts, container images, and references. The format is `/project_directory:/project_directory`. Multiple directories can be mounted by separating them with commas, for example: `/path1:/path1,/path2:/path2`. If you're setting `/project_directory:/root/`, note that the path in the config file must be set as a relative path to `/project_directory`.(e.g., `./data` in `config.yaml` means `/project_directory/data`.) (required)

# Part IV Output

   * **Output Structure**
      ```bash
      workup/
        ├── alignments
        │   ├── SRR7216005.DNA.bowtie2.mapq20.bam
        │   ├── SRR7216005.DNA.chr.bam
        │   └── SRR7216005.DNA.chr.masked.bam
        ├── clusters
        │   ├── cluster_sizes.pdf
        │   ├── cluster_sizes.png
        │   ├── SRR7216005.DNA.clusters
        │   └── SRR7216005.DNA.matrix.log
        ├── fastqs
        │   ├── SRR7216005_R1.barcoded.fastq.gz
        │   ├── SRR7216005_R1.barcoded_full.fastq.gz
        │   ├── SRR7216005_R1.barcoded_short.fastq.gz
        │   └── SRR7216005_R2.barcoded.fastq.gz
        ├── heatmap
        │   ├── SRR7216005.DNA.bias.txt
        │   ├── SRR7216005.DNA.final.pdf
        │   ├── SRR7216005.DNA.final.png
        │   ├── SRR7216005.DNA.final.txt
        │   ├── SRR7216005.DNA.iced.txt
        │   └── SRR7216005.DNA.raw.txt
        ├── ligation_efficiency.txt
        ├── logs
        │   ├── cluster
        │   ├── cluster_sizes.log
        │   ├── config_2025.11.24.yaml
        │   ├── multiqc.log
        │   ├── SRR7216005.bID.log
        │   ├── SRR7216005.bowtie2.log
        │   ├── SRR7216005.DNA_chr.log
        │   ├── SRR7216005_DPM.log
        │   ├── SRR7216005.heatmap.log
        │   ├── SRR7216005.make_clusters.log
        │   └── SRR7216005.trim_galore.logs
        ├── qc
        │   ├── multiqc_data
        │   │   ├── bowtie2_se_plot.txt
        │   │   ├── cutadapt_filtered_reads_plot.txt
        │   │   ├── cutadapt_trimmed_sequences_plot_3_Counts.txt
        │   │   ├── cutadapt_trimmed_sequences_plot_3_Obs_Exp.txt
        │   │   ├── cutadapt_trimmed_sequences_plot_5_Counts.txt
        │   │   ├── cutadapt_trimmed_sequences_plot_5_Obs_Exp.txt
        │   │   ├── llms-full.txt
        │   │   ├── multiqc_bowtie2.txt
        │   │   ├── multiqc_citations.txt
        │   │   ├── multiqc_cutadapt.txt
        │   │   ├── multiqc_data.json
        │   │   ├── multiqc_general_stats.txt
        │   │   ├── multiqc.log
        │   │   ├── multiqc.parquet
        │   │   ├── multiqc_software_versions.txt
        │   │   └── multiqc_sources.txt
        │   └── multiqc_report.html
        └── trimmed
            ├── SRR7216005_R1.barcoded.RDtrim.fastq.gz
            ├── SRR7216005_R1.barcoded.RDtrim.qc.txt
            ├── SRR7216005_R1.fastq.gz_trimming_report.txt
            ├── SRR7216005_R1_val_1.fq.gz
            ├── SRR7216005_R2.fastq.gz_trimming_report.txt
            └── SRR7216005_R2_val_2.fq.gz

      ```
    
   * **Output Interpretation**

      - **`multiqc_report.html`**: Open multiqc_report.html in a web browser to explore all sections interactively.

        - **General Statistics**: A combined table summarizing important metrics for each sample:

        <img width="1356" height="225" alt="2 multiqc general" src="https://github.com/user-attachments/assets/214eb743-7140-46b7-b244-f869cb624390" />

        - **Bowtie2**: The plot shows the number of reads aligning to the reference in different ways.

        <img width="1395" height="402" alt="3 multiqc bowtie2" src="https://github.com/user-attachments/assets/696106ff-d4a0-4fd2-abf4-b48d6eea13a4" />

        - **Cutadapt**: The statistics of Cutadapt cuts and removes adapters, including 'Filtered Reads' and 'Trimmed Sequence Lengths (5'/3')'.

        <img width="1338" height="441" alt="4 multiqc cutadapt" src="https://github.com/user-attachments/assets/c6be6345-2cd5-4822-98e5-41738fee63e2" />

      - **`*.DNA.chr.masked.bam`**

        - **Content**: This is the main alignment file in Binary Alignment Map (BAM) format. It contains all the sequencing reads and their mapping coordinates on the reference genome. This version has excluded reads overlapping with the annotations in the mask file and converted chromosome identifiers to UCSC format (e.g., chr1, chr2, etc.). For more information please refer to: https://genome.ucsc.edu/goldenpath/help/bam.html.
        - **Application**: It's the primary evidence for read alignment and can be used for detailed inspection in genome browsers or for downstream analyses.

      - **`*.DNA.clusters`**

        - **Content**: This is the SPRITE cluster file generated by grouping reads with the same barcode sequence into a cluster, which characterizes read clusters that are considered to have spatially adjacent features.
        - **Application**: Primarily used to construct contacts matrices.It can also be used for simple statistics and visualization (see `cluster_sizes.pdf/png`)

        <img width="480" height="480" alt="5 cluster_sizes" src="https://github.com/user-attachments/assets/97de6613-37b4-40d8-bf4a-ae3dc1781095" />

      - **`*.DNA.final.txt`**

        - **Content**: A contacts file containing a simple square matrix of values representing the contact strength or contact frequency between any two points on the genome, which is similar to an adjacency matrix.
        - **Application**: It can be applied to visualization (see `*.DNA.final.pdf/png`) as well as downstream three-dimensional genomic structure or chromatin interaction identification.
       
          <img width="2400" height="2400" alt="6 DNA final" src="https://github.com/user-attachments/assets/fcc3ca67-7a06-4ac0-86da-0c9ce937b6c8" />


# Part V Reference

Quinodoz, S.A., Bhat, P., Chovanec, P. et al. SPRITE: a genome-wide method for mapping higher-order 3D interactions in the nucleus using combinatorial split-and-pool barcoding. Nat Protoc 17, 36–75 (2022).


