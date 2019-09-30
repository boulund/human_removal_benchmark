# vim: syntax=python expandtab
# Evaluate best human removal tool
# xxx_MGISEQ_20190817_covaris_vs_enzyme
# Fredrik Boulund 2019

SAMPLES = [
    "fecal_200__V300019092_L01_1",
    "vag_200__V300019092_L01_2",
]
MAPPERS = [
    "bbmap",
    "bowtie2",
    "bmtagger",
]


rule all:
    input:
        expand("bbmap/{sample}.human.fq.gz", sample=SAMPLES),
        expand("bowtie2/{sample}.human_1.fq.gz", sample=SAMPLES),
        expand("bmtagger/{sample}_1.fastq", sample=SAMPLES),
        expand("kraken2/{mapper}.{sample}.human.kreport", mapper=MAPPERS[:-1], sample=SAMPLES),  #BMTagger doesn't give a file with human seqs
        expand("kraken2/{mapper}.{sample}.kreport", mapper=MAPPERS, sample=SAMPLES),
        #expand("summary/bowtie2_reads_{sample}.txt", sample=SAMPLES),
        #expand("summary/bowtie2_human_{sample}.txt", sample=SAMPLES),
        #expand("summary/bbmap_reads_{sample}.txt", sample=SAMPLES),
        #expand("summary/bbmap_human_{sample}.txt", sample=SAMPLES),

rule bbmap:
    input:
        read1="input/{sample}_1.fq.gz",
        read2="input/{sample}_2.fq.gz",
    output:
        read1="bbmap/{sample}_1.fq.gz",
        read2="bbmap/{sample}_2.fq.gz",
        human="bbmap/{sample}.human.fq.gz"
    log: 
        stats="logs/bbmap/{sample}.statsfile.txt",
        stdout="logs/bbmap/{sample}.log",
    benchmark: "benchmark/bbmap.{sample}.txt"
    threads: 10
    conda: "env.yaml"
    shell:
        """
        bbmap.sh \
            threads={threads} \
            in1={input.read1} \
            in2={input.read2} \
            path=/db/hg19 \
            outu1={output.read1} \
            outu2={output.read2} \
            outm={output.human} \
            statsfile={log.stats} \
            minid=0.95 \
            maxindel=3 \
            minhits=2 \
            bandwidthratio=0.16 \
            bandwidth=12 \
            qtrim=rl \
            trimq=10 \
            quickmatch \
            fast \
            untrim \
        """

rule bowtie2:
    input:
        read1="input/{sample}_1.fq.gz",
        read2="input/{sample}_2.fq.gz",
    output:
        read1="bowtie2/{sample}_1.fq.gz",
        read2="bowtie2/{sample}_2.fq.gz",
        human1="bowtie2/{sample}.human_1.fq.gz",
        human2="bowtie2/{sample}.human_2.fq.gz",
        sam="bowtie2/{sample}.sam.gz",
    log: 
        stderr="logs/bowtie2/{sample}.log",
    benchmark: "benchmark/bowtie2.{sample}.txt"
    threads: 10
    conda: "env.yaml"
    shell:
        """
        bowtie2 \
            -x /db/hg19/hg19 \
            -1 {input.read1} \
            -2 {input.read2} \
            -p {threads} \
            --very-sensitive-local \
            --un-conc-gz bowtie2/{wildcards.sample}_%.fq.gz \
            --al-conc-gz bowtie2/{wildcards.sample}.human_%.fq.gz \
            -S {output.sam} \
            2> {log.stderr}
        """

rule summarize_readcounts:
    input:
        bowtie2_read1=rules.bowtie2.output.read1,
        bowtie2_read2=rules.bowtie2.output.read2,
        bowtie2_human1=rules.bowtie2.output.human1,
        bowtie2_human2=rules.bowtie2.output.human2,
        bbmap_read1=rules.bbmap.output.read1,
        bbmap_read2=rules.bbmap.output.read2,
        bbmap_human=rules.bbmap.output.human,
    output:
        bowtie2_reads="summary/bowtie2_reads_{sample}.txt",
        bowtie2_human="summary/bowtie2_human_{sample}.txt",
        bbmap_reads="summary/bbmap_reads_{sample}.txt",
        bbmap_human="summary/bbmap_human_{sample}.txt",
    threads: 10
    conda: "env.yaml"
    shell:
        """
        reformat.sh in1={input.bowtie2_read1} in2={input.bowtie2_read2} 2> {output.bowtie2_reads}
        reformat.sh in1={input.bowtie2_human1} in2={input.bowtie2_human2} 2> {output.bowtie2_human}
        reformat.sh in1={input.bbmap_read1} in2={input.bbmap_read1} 2> {output.bbmap_reads}
        reformat.sh in1={input.bbmap_human} 2> {output.bbmap_human}
        """


rule download_KrakenTools:
    output:
        "scripts/KrakenTools/kreport2krona.py",
    shell:
        """
        cd scripts
        rm -rvf KrakenTools  # Otherwise git complains about Snakemake autocreated output dir
        git clone https://github.com/jenniferlu717/KrakenTools.git
        """


rule classify_human:
    input:
        bbmap_human="bbmap/{sample}.human.fq.gz",
        bowtie2_human1="bowtie2/{sample}.human_1.fq.gz",
        bowtie2_human2="bowtie2/{sample}.human_2.fq.gz",
    output:
        bowtie2_kreport="kraken2/bowtie2.{sample}.human.kreport",
        bowtie2_kraken="kraken2/bowtie2.{sample}.human.kraken",
        bbmap_kreport="kraken2/bbmap.{sample}.human.kreport",
        bbmap_kraken="kraken2/bbmap.{sample}.human.kraken",
    log: "logs/kraken2/{sample}.human.log"
    threads: 20
    conda: "env.yaml"
    params:
        db="/db/kraken2/ctmr_abvhf",
        confidence=0,
    shell:
        """
        kraken2 \
            --db {params.db} \
            --threads {threads} \
            --output {output.bbmap_kraken} \
            --confidence {params.confidence} \
            --report {output.bbmap_kreport} \
            --use-names \
            {input.bbmap_human} \
            2>> {log}
        kraken2 \
            --db {params.db} \
            --threads {threads} \
            --output {output.bowtie2_kraken} \
            --confidence {params.confidence} \
            --report {output.bowtie2_kreport} \
            --paired \
            --use-names \
            {input.bowtie2_human1} {input.bowtie2_human2} \
            2>> {log}
        """

rule classify_non_human:
    input:
        read1="{mapper}/{sample}_1.fq.gz",
        read2="{mapper}/{sample}_2.fq.gz",
    output:
        kreport="kraken2/{mapper}.{sample}.kreport",
        kraken="kraken2/{mapper}.{sample}.kraken",
    log: "logs/kraken2/{mapper}.{sample}.log"
    threads: 20
    conda: "env.yaml"
    params:
        db="/db/kraken2/ctmr_abvhf",
        confidence=0,
    shell:
        """
        kraken2 \
            --db {params.db} \
            --threads {threads} \
            --output {output.kraken} \
            --confidence {params.confidence} \
            --report {output.kreport} \
            --paired \
            --use-names \
            {input.read1} {input.read2}  \
            2>> {log}
        """

rule kraken2krona:
    input:
        kreport="kraken2/{sample}.human.kreport",
        script=rules.download_KrakenTools.output,
    output:
        filtered_kreport="kraken2/{sample}.human.filtered.kreport",
        krona="kraken2/krona/{sample}.human.krona",
        filtered="kraken2/krona/{sample}.human.filtered.krona",
    log: "logs/kraken2/{sample}.kraken2krona.log"
    threads: 1
    conda: "py37.yaml"
    shell:
        """
        {input.script} \
            --report-file {input.kreport} \
            --output {output.krona} \
            > {log}

        # Filtered version
        awk '{{if ($1 > 0.05) print $0}}' {input.kreport} > {output.filtered_kreport}
        {input.script} \
            --report-file {output.filtered_kreport} \
            --output {output.filtered} \
            >> {log}
        """

rule krona_kraken:
    input:
        kronas=expand("kraken2/krona/{sample}.human.krona", sample=SAMPLES),
        filtered_kronas=expand("kraken2/krona/{sample}.human.filtered.krona", sample=SAMPLES),
    output:
        html="kraken2/krona/all_samples.kraken2.human.html",
        filtered="kraken2/krona/all_samples.kraken2.human.filtered.html",
    log: "logs/kraken2/krona.log"
    conda: "env.yaml"
    shell:
        """
        ktImportText \
            -o {output.html} \
            {input.kronas} \
            > {log}
        ktImportText \
            -o {output.filtered} \
            {input.filtered_kronas} \
            >> {log}
        """


rule bmtagger:
    input:
        read1="input/{sample}_1.fq.gz",
        read2="input/{sample}_2.fq.gz",
    output:
        temp("{sample}_1.fastq"),
        temp("{sample}_2.fastq"),
        read1=temp("bmtagger/{sample}_1.fastq"),
        read2=temp("bmtagger/{sample}_2.fastq"),
    log: 
        stdout="logs/bmtagger/{sample}.stdout.log",
        stderr="logs/bmtagger/{sample}.stderr.log",
    benchmark: "benchmark/bmtagger.{sample}.txt"
    threads: 10
    conda: "env.yaml"
    shadow: "shallow"  # Otherwise fastq files mess up the base dir
    params:
        bitmask="/db/bmtagger/hg38/hg38.bitmask",
        srprism="/db/bmtagger/hg38/hg38.srprism",
        outbasename=lambda w: f"bmtagger/{w.sample}",
    shell:
        """
        reformat.sh \
            in1={input.read1} \
            in2={input.read2} \
            out1={wildcards.sample}_1.fastq \
            out2={wildcards.sample}_2.fastq \
            > {log.stdout} \
            2> {log.stderr}
            
        bmtagger.sh \
            -b {params.bitmask} \
            -x {params.srprism} \
            -q 1 \
            -1 {wildcards.sample}_1.fastq \
            -2 {wildcards.sample}_2.fastq \
            -o {params.outbasename} \
            -X \
            > {log.stdout} \
            2> {log.stderr} 
        """


rule rename_compress_bmtagger:
    input:
        read1="bmtagger/{sample}_1.fastq",
        read2="bmtagger/{sample}_2.fastq",
    output:
        read1="bmtagger/{sample}_1.fq.gz",
        read2="bmtagger/{sample}_2.fq.gz",
    log: 
        stdout="logs/bmtagger/{sample}.stdout.log",
        stderr="logs/bmtagger/{sample}.stderr.log",
    threads: 10
    conda: "env.yaml"
    shadow: "shallow"  # Otherwise fastq files mess up the base dir
    shell:
        """
        reformat.sh \
            in1={input.read1} \
            in2={input.read2} \
            out1={output.read1} \
            out2={output.read2} \
            > {log.stdout} \
            2> {log.stderr}
        """
