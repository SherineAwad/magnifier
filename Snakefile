configfile: "config.yaml"


with open(config['SAMPLES']) as fp:
    SAMPLES= fp.read().splitlines()
print(SAMPLES)


rule all:
      input:
        expand("{sample}.sam", sample=SAMPLES),
        expand("{sample}.coverage.histogram", sample=SAMPLES),
        expand("{sample}.alignment_metrics.txt", sample=SAMPLES),
        expand("{sample}.gc_bias_metrics.txt", sample=SAMPLES),
        expand("{sample}_quality.txt", sample=SAMPLES),
        expand("{sample}.insert_size_metrics.txt",sample=SAMPLES), 
        expand("{sample}.wgs_metrics.txt", sample=SAMPLES),
        expand("{sample}_screen.txt", sample=SAMPLES)
 


rule check_quality: 
      input:
        "{sample}.sam" 
      output:
        "{sample}_quality.txt"
      shell: 
          """ 
          samtools view {input} | awk '{{if($5<5) {{print $0}}}}'  | wc -l >  {output}
          samtools view {input} | awk '{{if($5>=5) {{print $0}}}}' | wc -l >> {output}
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
         fastq_screen --conf {params.CONF} {input}
        """

rule GCBias:
   input:
      "{sample}.sam",
   params:
      genome=config['GENOME']
   output:
      "{sample}.gc_bias_metrics.txt",
      "{sample}.gc_bias_metrics.pdf",
      "{sample}.summary_metrics.txt"
   shell:
       """
       picard CollectGcBiasMetrics I={input} O={output[0]} CHART={output[1]} \
       S={output[2]} R={params}.fasta
       """

rule CollectAlignmentSummaryMetrics: 
   input:
      "{sample}.sam",
   params:
      genome=config['GENOME']
   output:
      "{sample}.alignment_metrics.txt"
   shell:
       """
       picard CollectAlignmentSummaryMetrics I={input} O={output[0]} R={params}.fasta
       """

rule CoverageHistogram: 
    input: 
      "{sample}.sam",
    output: 
       "{sample}.coverage.histogram"
    shell: 
      """
      samtools coverage {input}  --histogram  -o {output}
      """

rule InsertSize:
    input:
       "{sample}.sam"
    output:
       "{sample}.insert_size_metrics.txt",
       "{sample}_insert_size_histogram.pdf"
    shell:
       """
        picard CollectInsertSizeMetrics I={input} O={output[0]} H={output[1]} M=0.5
      """

rule WGSMetrics:
    input:
      "{sample}.sam"
    params:
       genome=config['GENOME'],
       min_base_quality = config['MIN_BASE_QUALITY']
    output:
       "{sample}.wgs_metrics.txt"
    shell:
        """
        picard CollectWgsMetrics I={input} O={output} R={params.genome}.fasta Q={params.min_base_quality}
        """
