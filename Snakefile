configfile: "config.yaml"

rule all:
      input:
        expand("{sample}.sorted.bam", sample=config['SAMPLES']),
        expand("{sample}.coverage.histogram", sample=config['SAMPLES']),
        expand("{sample}.alignment_metrics.txt", sample=config['SAMPLES']),
        expand("{sample}.gc_bias_metrics.txt", sample=config['SAMPLES']),
        expand("{sample}_quality.txt", sample=config['SAMPLES']),
        #expand("{sample}_screen.txt", sample=config['SAMPLES'])
 

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
        "{sample}_quality.txt"
      shell: 
          """ 
          samtools view {input} | awk '{{if($5<5) {{print $0}}}}'  | wc -l >  {output}
          samtools view {input} | awk '{{if($5>=5) {{print $0}}}}' | wc -l >> {output}
          """ 

rule check_contaminate: 
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
      "{sample}.sorted.bam",
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
      "{sample}.sorted.bam",
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
      "{sample}.sorted.bam",
    output: 
       "{sample}.coverage.histogram"
    shell: 
      """
      samtools coverage {input}  --histogram  -o {output}
      """
