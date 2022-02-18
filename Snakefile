configfile: "config.yaml"


with open(config['SAMPLES']) as fp:
    SAMPLES= fp.read().splitlines()
print(SAMPLES)


rule all:
      input:
        expand("{sample}.sorted.bam", sample=SAMPLES),
        expand("{sample}_alignments.png", sample = SAMPLES),
        expand("{sample}.coverage.histogram", sample=SAMPLES),
        expand("{sample}.alignment_metrics.txt", sample=SAMPLES),
        expand("{sample}.gc_bias_metrics.txt", sample=SAMPLES),
        expand("{sample}.insert_size_metrics.txt",sample=SAMPLES), 
        expand("{sample}.wgs_metrics.txt", sample=SAMPLES),
        #expand("{sample}_screen.txt", sample=SAMPLES)
 

if config['PAIRED']:
    rule trim:
       input:
           r1 = "{sample}.r_1.fq.gz",
           r2 = "{sample}.r_2.fq.gz"
       output:
           "galore/{sample}.r_1_val_1.fq.gz",
           "galore/{sample}.r_2_val_2.fq.gz"
       conda: 'env/env-trim.yaml'
       shell:
           """
           mkdir -p galore
           mkdir -p fastqc
           trim_galore --gzip --retain_unpaired --trim1 --fastqc --fastqc_args "--outdir fastqc" -o galore --paired {input.r1} {input.r2}
           """
    rule tosam:
        input:
             r1 = "galore/{sample}.r_1_val_1.fq.gz",
             r2 = "galore/{sample}.r_2_val_2.fq.gz"
        params:
             genome=config['GENOME'],
             mem = config['MEMORY'],
        output:
             "{sample}.bam"
        conda: 'env/env-align.yaml'
        shell:
            """
            bbmap.sh {params.mem} in={input[0]} in2={input[1]} out={output} ref={params.genome}.fa
            """
else:
     rule trim:
       input:
           "{sample}.fq.gz",

       output:
           "galore/{sample}_trimmed.fq.gz",
       conda: 'env/env-trim.yaml'
       shell:
           """
           mkdir -p galore
           mkdir -p fastqc
           trim_galore --gzip --retain_unpaired --trim1 --fastqc --fastqc_args "--outdir fastqc" -o galore {input}
           """
     rule tosam:
        input:
            "galore/{sample}_trimmed.fq.gz"
        params:
            genome=config['GENOME'],
            mem = config['MEMORY'],
        output:
           "{sample}.bam"
        conda: 'env/env-align.yaml'
        shell:
          """
          bbmap.sh {params.mem} in={input} out={output} ref={params.genome}.fa 
          """

rule sort:
       input:
            "{sample}.bam"
       output:
            "{sample}.sorted.bam"
       params:
            "{sample}.tmp.sorted"
       log:
            "{sample}.sorted.log"
       conda: 'env/env-align.yaml'
       shell:
            """
                samtools sort -T {params} -o {output} {input}
                samtools index {output}
            """

rule check_quality: 
      input:
        "{sample}.sorted.bam" 
      output:
        "{sample}_alignments.txt",
        "{sample}_alignments.png"
      shell: 
          """
          samtools view -c -F 260 {input} >> {output[0]} 
          samtools view -c -F 256 {input} >> {output[0]}
          samtools view -c -f 4 {input} >> {output[0]} 
          samtools view -c {input} >> {output[0]} 
          samtools view {input} | awk '{{if($5<5) {{print $0}}}}'  | wc -l >>  {output[0]}
          samtools view {input} | awk '{{if($5>=20) {{print $0}}}}' | wc -l >> {output[0]}
          python scripts/calc_stats.py {output[0]} {output[1]} 
          """ 

rule check_contamination: 
     input: 
         "{sample}.r_1.fq.gz"
     output: 
         "{sample}_screen.txt" 
     params: 
         CONF=config['CONF'] 
     shell: 
        """
         fastq_screen --conf {params.CONF} {input} --aligner bowtie2
        """

rule GCBias:
   input:
      "{sample}.sorted.bam",
   params:
      genome= config['GENOME']
   output:
      "{sample}.gc_bias_metrics.txt",
      "{sample}.gc_bias_metrics.pdf",
      "{sample}.summary_metrics.txt"
   shell:
       """
       picard CollectGcBiasMetrics I={input} O={output[0]} CHART={output[1]} \
       S={output[2]} R={params}.fa
       """

rule CollectAlignmentSummaryMetrics: 
   input:
      "{sample}.sorted.bam",
   params:
      genome=config['GENOME']
   output:
      "{sample}.alignment_metrics.txt"
   shell:
       """
       picard CollectAlignmentSummaryMetrics I={input} O={output[0]} R={params}.fa
       """

rule CoverageHistogram: 
    input: 
      "{sample}.sorted.bam",
    output: 
       "{sample}.coverage.histogram"
    shell: 
      """
      samtools coverage {input}  --histogram  -o {output}
      """

rule InsertSize:
    input:
       "{sample}.sorted.bam"
    output:
       "{sample}.insert_size_metrics.txt",
       "{sample}.insert_size_histogram.pdf"
    shell:
       """
        picard CollectInsertSizeMetrics I={input} O={output[0]} H={output[1]} M=0.5
       """

rule WGSMetrics:
    input:
      "{sample}.sorted.bam"
    params:
       genome=config['GENOME'],
       min_base_quality = config['MIN_BASE_QUALITY']
    output:
       "{sample}.wgs_metrics.txt"
    shell:
        """
        picard CollectWgsMetrics I={input} O={output} R={params.genome}.fa Q={params.min_base_quality}
        """
