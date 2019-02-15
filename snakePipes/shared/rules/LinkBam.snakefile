rule link_bam:
    input:
        indir+"/{sample}"+bam_ext
    output:
        bam=mapping_prg+"/{sample}.bam",
        bai=mapping_prg+"/{sample}.bam.bai"
    conda: CONDA_SHARED_ENV
    shell:"""
        ( [ -f {output.bam} ] || ln -s -r {input} {output.bam} ); 
        if [ -f {input}.bai ]; then ln -s -r {input}.bai {output.bai}; else samtools index {output.bam}; fi
        touch -h -m {output.bai}
        """

#rule samtools_index_external:
#    input:
#        mapping_prg+"/{sample}.bam"
#    output:
#        mapping_prg+"/{sample}.bam.bai"
#    conda: CONDA_SHARED_ENV
#    shell: "samtools index {input}"


rule link_bam_bai_external:
    input:
        bam = mapping_prg+"/{sample}.bam",
        bai = mapping_prg+"/{sample}.bam.bai"
    output:
        bam_out = "filtered_bam/{sample}.filtered.bam",
        bai_out = "filtered_bam/{sample}.filtered.bam.bai",
    shell:
        "( [ -f {output.bam_out} ] || ( ln -s -r {input.bam} {output.bam_out} && ln -s -r {input.bai} {output.bai_out} ) )"

# rule samtools_index_filtered:
#     input:
#         "filtered_bam/{sample}.filtered.bam"
#     output:
#         "filtered_bam/{sample}.filtered.bam.bai"
#     conda: CONDA_SHARED_ENV
#     shell: "samtools index {input}"

rule sambamba_flagstat:
       input:
           mapping_prg+"/{sample}.bam"
       output:
           "Sambamba/{sample}.markdup.txt"
       conda: CONDA_SHARED_ENV
       shell: """
           sambamba flagstat -p {input} > {output}
           """