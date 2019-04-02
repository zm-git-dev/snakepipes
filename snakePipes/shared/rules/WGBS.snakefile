import os
import re
from operator import is_not
import tempfile
import pandas


## function to get the name of the samplesheet and extend the name of the folder for all analyses relying on sample_info
def get_outdir(folder_name):
    sample_name = re.sub('_sampleSheet.[a-z]{3}$','',os.path.basename(sampleSheet))
    return("{}_{}".format(folder_name, sample_name))

## count the number of fields in the chromosome name and generate awk string
def get_awk_cmd(fasta):
    with open(fasta) as f:
        line = f.readline()
    nF=len(line.split(' '))
    return ('\'{{print $1, ${}, ${}+1, ${}, ${}}}\''.format(nF+1,nF+1,nF+2,nF+4))

###symlink bams if this is the starting point
if fromBam:
    rule link_bam:
        input:
            indir+"/{sample}"+bam_ext
        output:
            "bams/{sample}.bam"
        shell:
            "( [ -f {output} ] || ln -s -r {input} {output} ) " #&& touch -h {output}"


###get automatic cut threshold for hard-trimming of 5' ends
if trimReads=='auto':
    if fqcin:
        rule get_cut_thd:
            input:
                R1zip = fqcin+"/{sample}"+reads[0]+"_fastqc.zip",
                R2zip = fqcin+"/{sample}"+reads[1]+"_fastqc.zip"
            output:
                R12ct= "FastQC_In/{sample}.R12.ct.txt"
            log:"FastQC_In/logs/{sample}.get_cut_thd.log"
            threads: 1
            run:
                for f,g,z,l in zip ({input.R1zip},{input.R2zip},output,log):
                    with open(z, "w") as oo:
                        cutThdRes=calc_cutThd([f,g],fqcin,l,outdir)
                        oo.write('\n'.join('%s\t%s\n' % x for x in cutThdRes))
                os.chdir(outdir)


    rule trimReads:
        input:
            R1 = "FASTQ/{sample}"+reads[0]+".fastq.gz",
            R2 = "FASTQ/{sample}"+reads[1]+".fastq.gz",
            R12ct= "FastQC_In/{sample}.R12.ct.txt"
        output:
            R1cut=temp("FASTQ_Cutadapt/{sample}"+reads[0]+".fastq.gz"),
            R2cut=temp("FASTQ_Cutadapt/{sample}"+reads[1]+".fastq.gz")
        log:
            err="FASTQ_Cutadapt/logs/{sample}.trimReads.err",
            out="FASTQ_Cutadapt/logs/{sample}.trimReads.out"
        threads: nthreads
        conda: CONDA_SHARED_ENV
        shell: "ct=($(cat {input.R12ct})); cutadapt -a AGATCGGAAGAGC -A AGATCGGAAGAGC --minimum-length 30  -n 5 -j {threads} -u ${{ct[0]}}  -U ${{ct[1]}} -o {output.R1cut} -p {output.R2cut} {input.R1} {input.R2} 1>{log.out} 2>{log.err}"


elif trimReads=='user':
    rule trimReads:
        input:
            R1 = "FASTQ/{sample}"+reads[0]+".fastq.gz",
            R2 = "FASTQ/{sample}"+reads[1]+".fastq.gz"
        output:
            R1cut="FASTQ_Cutadapt/{sample}"+reads[0]+".fastq.gz",
            R2cut="FASTQ_Cutadapt/{sample}"+reads[1]+".fastq.gz"
        log:
            err="FASTQ_Cutadapt/logs/{sample}.trimReads.err",
            out="FASTQ_Cutadapt/logs/{sample}.trimReads.out"
        params:
            adapterSeq=adapterSeq,
            trimThreshold=trimThreshold,
            trimOtherArgs=lambda wildcards: '' if trimOtherArgs is None else str(trimOtherArgs)
        threads: nthreads
        conda: CONDA_SHARED_ENV
        shell: "cutadapt -a {params.adapterSeq} -A {params.adapterSeq} -q {params.trimThreshold} -m 30 -j {threads} {params.trimOtherArgs} -o {output.R1cut} -p {output.R2cut} {input.R1} {input.R2} 1>{log.out} 2>{log.err}"



#if not trimReads is None:

rule preTrimFastQC:
    input:
        R1="FASTQ/{sample}"+reads[0]+".fastq.gz",
        R2="FASTQ/{sample}"+reads[1]+".fastq.gz"
    output:
        R1fqc="FastQC/{sample}"+reads[0]+"_fastqc.html",
        R2fqc="FastQC/{sample}"+reads[1]+"_fastqc.html"
    log:
        err="FastQC/logs/{sample}.preTrimFastQC.err",
        out="FastQC/logs/{sample}.preTrimFastQC.out"
    params:
        fqcout=os.path.join(outdir,'FastQC')
    threads: nthreads
    conda: CONDA_SHARED_ENV
    shell: "fastqc --outdir {params.fqcout} -t  {threads} {input.R1} {input.R2} 1>{log.out} 2>{log.err}"


rule postTrimFastQC:
    input:
        R1cut="FASTQ_Cutadapt/{sample}"+reads[0]+".fastq.gz",
        R2cut="FASTQ_Cutadapt/{sample}"+reads[1]+".fastq.gz"
    output:
        R1fqc="FastQC_Cutadapt/{sample}"+reads[0]+"_fastqc.html",
        R2fqc="FastQC_Cutadapt/{sample}"+reads[1]+"_fastqc.html"
    log:
        err="FastQC_Cutadapt/logs/{sample}.postTrimFastQC.err",
        out="FastQC_Cutadapt/logs/{sample}.postTrimFastQC.out"
    params:
        fqcout=os.path.join(outdir,'FastQC_Cutadapt')
    threads: nthreads
    conda: CONDA_SHARED_ENV
    shell: "fastqc --outdir {params.fqcout} -t  {threads} {input.R1cut} {input.R2cut} 1>{log.out} 2>{log.err}"

if convRef:
    rule conv_ref:
        input:
            refG=refG
        output:
            cref_sa=os.path.join("aux_files",re.sub('.fa','.fa.bwameth.c2t.sa',os.path.basename(refG))),
            cref_amb=os.path.join("aux_files",re.sub('.fa','.fa.bwameth.c2t.amb',os.path.basename(refG))),
            locrefG=os.path.join("aux_files",os.path.basename(refG))
        params:
            locdict=os.path.join("aux_files",re.sub('.fa','.dict',os.path.basename(refG)))
        log:
            err="aux_files/logs/conv_ref.err",
            out="aux_files/logs/conv_ref.out"
        threads: 1
        conda: CONDA_WGBS_ENV
        shell:"ln -s {input.refG} {output.locrefG}; bwameth.py index {output.locrefG}; samtools faidx {output.locrefG}; picard CreateSequenceDictionary R={output.locrefG} O={params.locdict}  1>{log.out} 2>{log.err}"

if not trimReads is None:
    rule map_reads:
        input:
            lambda convRef: os.path.join("aux_files",re.sub('.fa','.fa.bwameth.c2t.sa',os.path.basename(refG))) if convRef is True else [],
            lambda convRef: os.path.join("aux_files",re.sub('.fa','.fa.bwameth.c2t.amb',os.path.basename(refG))) if convRef is True  else [],
            R1cut="FASTQ_Cutadapt/{sample}"+reads[0]+".fastq.gz",
            R2cut="FASTQ_Cutadapt/{sample}"+reads[1]+".fastq.gz",
            crefG=crefG
        output:
            sbam="bams/{sample}.sorted.bam"
        log:
            err="bams/logs/{sample}.map_reads.err",
            out="bams/logs/{sample}.map_reads.out"
        params:
            #tempdir=tempfile.mkdtemp(suffix='',prefix="{sample}",dir=tempdir),
            sortThreads=min(nthreads,4),
            RG=lambda wildcards: RG_dict[wildcards.sample]
        threads: nthreads
        conda: CONDA_WGBS_ENV
        shell: "tmp_map=$(mktemp -d -p $TMPDIR -t XXXXX.{wildcards.sample});echo $tmp_map;  bwameth.py --threads  {threads}  --read-group {params.RG} --reference {input.crefG} {input.R1cut} {input.R2cut} | samtools sort -T $tmp_map -m 3G -@ {params.sortThreads} -o {output.sbam} 1>{log.out} 2>{log.err}"

if trimReads is None and not fromBam:
    rule map_reads:
        input:
            lambda convRef: os.path.join("aux_files",re.sub('.fa','.fa.bwameth.c2t.sa',os.path.basename(refG))) if convRef is True else [],
            lambda convRef: os.path.join("aux_files",re.sub('.fa','.fa.bwameth.c2t.amb',os.path.basename(refG))) if convRef is True  else [],
            R1="FASTQ/{sample}"+reads[0]+".fastq.gz",
            R2="FASTQ/{sample}"+reads[1]+".fastq.gz",
            crefG=crefG
        output:
            sbam=temp("bams/{sample}.sorted.bam")
        log:
            err="bams/logs/{sample}.map_reads.err",
            out="bams/logs/{sample}.map_reads.out"
        params:
            #tempdir=tempfile.mkdtemp(suffix='',prefix="{sample}",dir=tempdir),
            sortThreads=min(nthreads,4),
            RG=lambda wildcards: RG_dict[wildcards.sample]
        threads: nthreads
        conda: CONDA_WGBS_ENV
        shell: "tmp_map=$(mktemp -d -p $TMPDIR -t XXXXX.{wildcards.sample});echo $tmp_map; bwameth.py --threads  {threads}  --read-group {params.RG} --reference {input.crefG} {input.R1} {input.R2} | samtools sort -T $tmp_map -m 3G -@ {params.sortThreads} -o {output.sbam} 1>{log.out} 2>{log.err}"

if not fromBam:
    rule index_bam:
        input:
            sbam="bams/{sample}.sorted.bam"
        output:
            sbami="bams/{sample}.sorted.bam.bai"
        log:
            err="bams/logs/{sample}.index_bam.err",
            out="bams/logs/{sample}.index_bam.out"
        conda: CONDA_SHARED_ENV
        shell: "samtools index {input.sbam} 1>{log.out} 2>{log.err}"

    rule rm_dupes:
        input:
            sbami="bams/{sample}.sorted.bam.bai",
            sbam="bams/{sample}.sorted.bam"
        output:
            rmDupbam="bams/{sample}.PCRrm.bam"
        log:
            err="bams/logs/{sample}.rm_dupes.err",
            out="bams/logs/{sample}.rm_dupes.out"
        params:
            #tempdir=tempfile.mkdtemp(suffix='',prefix='',dir=tempdir)
        threads: nthreads
        conda: CONDA_SHARED_ENV
        shell: "tmp_dupes=$(mktemp -d -p $TMPDIR -t XXXXX.{wildcards.sample}); echo $tmp_dupes; sambamba markdup --hash-table-size=4194304 --remove-duplicates --tmpdir $tmp_dupes -t {threads} {input.sbam} {output.rmDupbam} 1>{log.out} 2>{log.err}"

rule index_PCRrm_bam:
    input:
        sbam="bams/{sample}.PCRrm.bam"
    output:
        sbami="bams/{sample}.PCRrm.bam.bai"
    params:
    log:
        err="bams/logs/{sample}.index_PCRrm_bam.err",
        out="bams/logs/{sample}.index_PCRrm_bam.out"
    threads: 1
    conda: CONDA_SHARED_ENV
    shell: "samtools index {input.sbam} 1>{log.out} 2>{log.err}"


rule get_ran_CG:
    input:
        refG=refG
    output:
        pozF="aux_files/"+re.sub('.fa*','.poz.gz',os.path.basename(refG)),
        ranCG=os.path.join("aux_files",re.sub('.fa','.poz.ran1M.sorted.bed',os.path.basename(refG)))
    params:
        awkCmd=get_awk_cmd(refG)
    log:
        err="aux_files/logs/get_ran_CG.err"
    threads: 1
    conda: CONDA_SHARED_ENV
    shell: 'set +o pipefail; ' + os.path.join(workflow_tools,'methylCtools') + " fapos {input.refG}  " + re.sub('.gz','',"{output.pozF}") + ';cat '+ re.sub('.gz','',"{output.pozF}") +' | grep "+" -' + " | shuf | head -n 1000000 | awk {params.awkCmd}" + ' - | tr " " "\\t" | sort -k 1,1 -k2,2n - > ' + "{output.ranCG} 2>{log.err}"


rule calc_Mbias:
    input:
        refG=refG,
        rmDupBam="bams/{sample}.PCRrm.bam",
        sbami="bams/{sample}.PCRrm.bam.bai"
    output:
        mbiasTXT="QC_metrics/{sample}.Mbias.txt"
    log:
        out="QC_metrics/logs/{sample}.calc_Mbias.out"
    threads: nthreads
    conda: CONDA_WGBS_ENV
    shell: "MethylDackel mbias {input.refG} {input.rmDupBam} {output.mbiasTXT} -@ {threads} 1>{log.out} 2>{output.mbiasTXT}"

if convRef:
    rule calc_genome_size:
        input:
            refG=refG
        output:
            gsize="aux_files/gsize.txt"
        log:
            err="aux_files/logs/gsize.err"
        threads: 1
        conda: CONDA_WGBS_ENV
        shell: "faCount {input.refG} | awk \'END{{print $2-$7}}\'  > {output.gsize} 2>{log.err}"

    if not skipGCbias:
        rule get_twobit_genome:
            input:
                refG=refG
            output:
                twobit="aux_files/"+ re.sub(".fa",".2bit",os.path.basename(refG))
            log:
                err="aux_files/logs/fatotwobit.err"
            threads: 1
            conda: CONDA_WGBS_ENV
            shell: "faToTwoBit {input.refG} {output.twobit} 2>{log.err}"

        rule calc_GCbias:
            input:
                refG=refG,
                rmDupBam="bams/{sample}.PCRrm.bam",
                sbami="bams/{sample}.PCRrm.bai",
                gsize="aux_files/gsize.txt",
                twobit="aux_files/"+ re.sub(".fa",".2bit",os.path.basename(refG))
            output:
                GCbiasTXT="QC_metrics/{sample}.freq.txt",
                GCbiasPNG="QC_metrics/{sample}.GCbias.png"
            log:
                out="QC_metrics/logs/{sample}.calc_GCbias.out"
            threads: nthreads
            conda: CONDA_SHARED_ENV
            shell: "genomeSize=($(cat {input.gsize}));computeGCBias -b {input.rmDupBam} --effectiveGenomeSize $genomeSize -g {input.twobit} -l 300 --GCbiasFrequenciesFile {output.GCbiasTXT} -p {threads} --biasPlot {output.GCbiasPNG} --plotFileFormat png "

else:
    if not skipGCbias:
        rule calc_GCbias:
            input:
                refG=refG,
                rmDupBam="bams/{sample}.PCRrm.bam",
                sbami="bams/{sample}.PCRrm.bam.bai"
            output:
                GCbiasTXT="QC_metrics/{sample}.freq.txt",
                GCbiasPNG="QC_metrics/{sample}.GCbias.png"
            params:
                genomeSize=genome_size,
                twobitpath=genome_2bit
            log:
                out="QC_metrics/logs/{sample}.calc_GCbias.out"
            threads: nthreads
            conda: CONDA_SHARED_ENV
            shell: "computeGCBias -b {input.rmDupBam} --effectiveGenomeSize {params.genomeSize} -g {params.twobitpath} -l 300 --GCbiasFrequenciesFile {output.GCbiasTXT} -p {threads} --biasPlot {output.GCbiasPNG} --plotFileFormat png "

if not skipDOC:
    if intList:
        rule depth_of_cov:
            input:
                irefG=crefG if convRef is True else refG,
                rmDupBam="bams/{sample}.PCRrm.bam",
                sbami="bams/{sample}.PCRrm.bam.bai",
                ranCG=os.path.join("aux_files",re.sub('.fa','.poz.ran1M.sorted.bed',os.path.basename(refG))),
                intList=intList
            output:
                outFileList=calc_doc(intList,True,skipDOC)
            params:
                tempdir=tempdir,
                auxdir=os.path.join(outdir,"aux_files"),
                OUTlist=lambda wildcards,output: [w.replace('.sample_summary', '') for w in output.outFileList],
                OUTlist0=lambda wildcards,output: [w.replace('.sample_summary', '') for w in output.outFileList][0],
                OUTlist1=lambda wildcards,output: [w.replace('.sample_summary', '') for w in output.outFileList][1],
                auxshell=lambda wildcards,input,output: ';'.join("gatk3 -Xmx30g -Djava.io.tmpdir="+ tempdir +" -T DepthOfCoverage -R "+ input.irefG +" -o "+ oi +" -I " + input.rmDupBam + " -ct 0 -ct 1 -ct 2 -ct 5 -ct 10 -ct 15 -ct 20 -ct 30 -ct 50  -omitBaseOutput -mmq 10 --partitionType sample -L " + bi for oi,bi in zip([w.replace('.sample_summary', '') for w in output.outFileList][2:],input.intList))
            log:
                err="QC_metrics/logs/{sample}.depth_of_cov.err",
                out="QC_metrics/logs/{sample}.depth_of_cov.out"
            threads: 1
            conda: CONDA_WGBS_ENV
            shell: "gatk3 -Xmx30g -Djava.io.tmpdir={params.tempdir} -T DepthOfCoverage -R {input.irefG} -o {params.OUTlist0} -I {input.rmDupBam} -ct 0 -ct 1 -ct 2 -ct 5 -ct 10 -ct 15 -ct 20 -ct 30 -ct 50  -omitBaseOutput -mmq 10 --partitionType sample ; gatk3 -Xmx30g -Djava.io.tmpdir={params.tempdir}  -T DepthOfCoverage -R {input.irefG} -o {params.OUTlist1} -I {input.rmDupBam} -ct 0 -ct 1 -ct 2 -ct 5 -ct 10 -ct 15 -ct 20 -ct 30 -ct 50  -omitBaseOutput -mmq 10 --partitionType sample -L {input.ranCG}; {params.auxshell} 1>{log.out} 2>{log.err}"


    else:
        rule depth_of_cov:
            input:
                irefG=crefG if convRef is True else refG,
                rmDupBam="bams/{sample}.PCRrm.bam",
                sbami="bams/{sample}.PCRrm.bam.bai",
                ranCG=os.path.join("aux_files",re.sub('.fa','.poz.ran1M.sorted.bed',os.path.basename(refG)))
            output:
                outFileList=calc_doc(intList,True,skipDOC)
            params:
                tempdir=tempdir,
                auxdir=os.path.join(outdir,"aux_files"),
                OUTlist0=lambda wildcards,output: output.outFileList[0].replace('.sample_summary', ''),
                OUTlist1=lambda wildcards,output: output.outFileList[1].replace('.sample_summary','')
            log:
                err="QC_metrics/logs/{sample}.depth_of_cov.err",
                out="QC_metrics/logs/{sample}.depth_of_cov.out"
            threads: 1
            conda: CONDA_WGBS_ENV
            shell: "gatk3 -Xmx30g -Djava.io.tmpdir={params.tempdir} -T DepthOfCoverage -R {input.irefG} -o {params.OUTlist0} -I {input.rmDupBam} -ct 0 -ct 1 -ct 2 -ct 5 -ct 10 -ct 15 -ct 20 -ct 30 -ct 50  -omitBaseOutput -mmq 10 --partitionType sample ; gatk3 -Xmx30g -Djava.io.tmpdir={params.tempdir}  -T DepthOfCoverage -R {input.irefG} -o {params.OUTlist1} -I {input.rmDupBam} -ct 0 -ct 1 -ct 2 -ct 5 -ct 10 -ct 15 -ct 20 -ct 30 -ct 50  -omitBaseOutput -mmq 10 --partitionType sample -L {input.ranCG} 1>{log.out} 2>{log.err}"

if not trimReads is None and not fromBam:
    rule downsample_reads:
        input:
            R1cut="FASTQ_Cutadapt/{sample}"+reads[0]+".fastq.gz",
            R2cut="FASTQ_Cutadapt/{sample}"+reads[1]+".fastq.gz"
        output:
            R1downsampled="FASTQ_downsampled/{sample}"+reads[0]+".fastq.gz",
            R2downsampled="FASTQ_downsampled/{sample}"+reads[1]+".fastq.gz"
        log:
            err="FASTQ_downsampled/logs/{sample}.downsample_reads.err",
            out="FASTQ_downsampled/logs/{sample}.downsample_reads.out"
        threads: nthreads
        conda: CONDA_SHARED_ENV
        shell: """
                seqtk sample -s 100 {input.R1cut} 5000000 | pigz -p {threads} -9 > {output.R1downsampled}
                seqtk sample -s 100 {input.R2cut} 5000000 | pigz -p {threads} -9 > {output.R2downsampled}
                1>{log.out} 2>{log.err}
               """

    rule conv_rate:
        input:
            R1downsampled="FASTQ_downsampled/{sample}"+reads[0]+".fastq.gz",
            R2downsampled="FASTQ_downsampled/{sample}"+reads[1]+".fastq.gz"
        output:
            R12cr="QC_metrics/{sample}.conv.rate.txt"
        params:
            read_root=lambda wildcards: "FASTQ_downsampled/"+wildcards.sample
        log:
            err="QC_metrics/logs/{sample}.conv_rate.err",
            out="QC_metrics/logs/{sample}.conv_rate.out"
        threads: 1
        shell: os.path.join(workflow_tools,'conversionRate_KS.sh ')+ "{params.read_root} {output.R12cr} 1>{log.out} 2>{log.err}"

else:
    if not fromBam:
        rule downsample_reads:
            input:
                R1="FASTQ/{sample}"+reads[0]+".fastq.gz",
                R2="FASTQ/{sample}"+reads[1]+".fastq.gz"
            output:
                R1downsampled="FASTQ_downsampled/{sample}"+reads[0]+".fastq.gz",
                R2downsampled="FASTQ_downsampled/{sample}"+reads[1]+".fastq.gz"
            log:
                err="FASTQ_downsampled/logs/{sample}.downsample_reads.err",
                out="FASTQ_downsampled/logs/{sample}.downsample_reads.out"
            threads: nthreads
            conda: CONDA_SHARED_ENV
            shell: """
                    seqtk sample -s 100 {input.R1} 5000000 | pigz -p {threads} -9 > {output.R1downsampled}
                    seqtk sample -s 100 {input.R2} 5000000 | pigz -p {threads} -9 > {output.R2downsampled}
                    1>{log.out} 2>{log.err}
                   """

        rule conv_rate:
            input:
                R1="FASTQ_downsampled/{sample}"+reads[0]+".fastq.gz",
                R2="FASTQ_downsampled/{sample}"+reads[1]+".fastq.gz"
            output:
                R12cr="QC_metrics/{sample}.conv.rate.txt"
            params:
                pfx="FASTQ_downsampled/{sample}"
            log:
                err="QC_metrics/logs/{sample}.conv_rate.err",
                out="QC_metrics/logs/{sample}.conv_rate.out"
            threads: 1
            shell: os.path.join(workflow_tools,'conversionRate_KS.sh ')+ "{params.pfx} {output.R12cr} 1>{log.out} 2>{log.err}"


if fromBam:
    rule get_flagstat:
        input:
            bam="bams/{sample}.bam"
        output:
            fstat="QC_metrics/{sample}.flagstat"
        log:
            err="QC_metrics/logs/{sample}.get_flagstat.err"
        threads: 1
        conda: CONDA_SHARED_ENV
        shell: "samtools flagstat {input.bam} > {output.fstat} 2>{log.err}"
else:
    rule get_flagstat:
        input:
            bam="bams/{sample}.sorted.bam"
        output:
            fstat="QC_metrics/{sample}.flagstat"
        log:
            err="QC_metrics/logs/{sample}.get_flagstat.err"
        threads: 1
        conda: CONDA_SHARED_ENV
        shell: "samtools flagstat {input.bam} > {output.fstat} 2>{log.err}"

rule produce_report:
    input:
        calc_doc(intList,False,skipDOC),
        expand("QC_metrics/{sample}.conv.rate.txt",sample=samples) if not fromBam else [],
        mbiasTXT=expand("QC_metrics/{sample}.Mbias.txt",sample=samples),
        fstat=expand("QC_metrics/{sample}.flagstat",sample=samples)
    output:
        QCrep='QC_metrics/QC_report.html'
    params:
        auxdir=os.path.join(outdir,"aux_files")
    log:
        err="QC_metrics/logs/produce_report.err",
        out="QC_metrics/logs/produce_report.out"
    conda: CONDA_RMD_ENV
    threads: 1
    shell: "cp -v " + os.path.join(workflow_rscripts,"WGBS_QC_report_template.Rmd")+ " " + os.path.join("aux_files", "WGBS_QC_report_template.Rmd") + ';Rscript -e "rmarkdown::render(\''+os.path.join(outdir,"aux_files", "WGBS_QC_report_template.Rmd")+'\', params=list(QCdir=\'"' + os.path.join(outdir,"QC_metrics") +'"\' ), output_file =\'"'+ os.path.join(outdir,"QC_metrics",'QC_report.html"\'')+')"' + " 1>{log.out} 2>{log.err}"


if mbias_ignore=="auto":
    rule methyl_extract:
        input:
            rmDupbam="bams/{sample}.PCRrm.bam",
            sbami="bams/{sample}.PCRrm.bam.bai",
            refG=refG,
            mbiasTXT="QC_metrics/{sample}.Mbias.txt"
        output:
            methTab="methXT/{sample}_CpG.bedGraph"
        params:
            OUTpfx=lambda wildcards,output: os.path.join(outdir,re.sub('_CpG.bedGraph','',output.methTab))
        log:
            err="methXT/logs/{sample}.methyl_extract.err",
            out="methXT/logs/{sample}.methyl_extract.out"
        threads: nthreads
        conda: CONDA_WGBS_ENV
        shell: "mi=$(cat {input.mbiasTXT} | sed 's/Suggested inclusion options: //' );MethylDackel extract  -o {params.OUTpfx} -q 10 -p 20 $mi --minDepth 10 --mergeContext --maxVariantFrac 0.25 --minOppositeDepth 5 -@ {threads} {input.refG} " + os.path.join(outdir,"{input.rmDupbam}") + " 1>{log.out} 2>{log.err}"


else:
    rule methyl_extract:
        input:
            rmDupbam="bams/{sample}.PCRrm.bam",
            sbami="bams/{sample}.PCRrm.bam.bai",
            refG=refG
        output:
            methTab="methXT/{sample}_CpG.bedGraph"
        params:
            mbias_ignore=mbias_ignore,
            OUTpfx=lambda wildcards,output: os.path.join(outdir,re.sub('_CpG.bedGraph','',output.methTab))
        log:
            err="methXT/logs/{sample}.methyl_extract.err",
            out="methXT/logs/{sample}.methyl_extract.out"
        threads: nthreads
        conda: CONDA_WGBS_ENV
        shell: "MethylDackel extract  -o {params.OUTpfx} -q 10 -p 20 {params.mbias_ignore} --minDepth 10 --mergeContext --maxVariantFrac 0.25 --minOppositeDepth 5 -@ {threads} {input.refG} " + os.path.join(outdir,"{input.rmDupbam}") + " 1>{log.out} 2>{log.err}"

if blackList is None:
    rule CpG_filt:
        input:
            methTab="methXT/{sample}_CpG.bedGraph"
        output:
            tabFilt="methXT/{sample}.CpG.filt2.bed"
        params:
            OUTtemp=lambda wildcards,input: os.path.join(outdir,re.sub('_CpG.bedGraph','.CpG.filt.bed',input.methTab))
        log:
            err="methXT/logs/{sample}.CpG_filt.err",
            out="methXT/logs/{sample}.CpG_filt.out"
        threads: 1
        conda: CONDA_WGBS_ENV
        shell: '''awk \'(NR>1)\' {input.methTab} | awk \'{{ print $0, $5+$6, $1\"_\"$2}}\' | tr " " "\t" | sed \'1i chr\tstart\tend\tBeta\tM\tU\tCov\tms\' > {params.OUTtemp};mv -v {params.OUTtemp} {output.tabFilt} 1>{log.out} 2>{log.err}'''

else:
    rule CpG_filt:
        input:
            methTab="methXT/{sample}_CpG.bedGraph",
            blackListF=blackList
        output:
            tabFilt="methXT/{sample}.CpG.filt2.bed"
        params:
            OUTtemp=lambda wildcards,input: os.path.join(outdir,re.sub('_CpG.bedGraph','.CpG.filt.bed',input.methTab))
        log:
            err="methXT/logs/{sample}.CpG_filt.err",
            out="methXT/logs/{sample}.CpG_filt.out"
        threads: 1
        conda: CONDA_WGBS_ENV
        shell: '''awk \'(NR>1)\' {input.methTab} | awk \'{{ print $0, $5+$6, $1\"_\"$2}}\' | tr " " "\t" | sed \'1i chr\tstart\tend\tBeta\tM\tU\tCov\tms\' > {params.OUTtemp};bedtools intersect -v -a {params.OUTtemp} -b {input.blackListF} > {output.tabFilt} 1>{log.out} 2>{log.err}'''


if sampleSheet or intList:
    rule make_CG_bed:
        input:
            pozF="aux_files/"+re.sub('.fa*','.poz.gz',os.path.basename(refG))
        output:
            imdF="aux_files/"+re.sub('.fa*','.CpG.bed',os.path.basename(refG))
        params:
            awkCmd=get_awk_cmd(refG)
        log:
            err="aux_files/logs/make_CG_bed.err"
        threads: 1
        conda: CONDA_WGBS_ENV
        shell: 'grep "+"' + " {input.pozF}  | awk {params.awkCmd}" + ' - | tr " " "\\t" | sort -k 1,1 -k2,2n - > ' + "{output.imdF}"


if sampleSheet:
    rule prep_for_stats:
        input: expand("methXT/{sample}.CpG.filt2.bed",sample=samples)
        output:
            Limdat='{}/limdat.LG.RData'.format(get_outdir("merged_methylation_data")),
            MetIN='{}/metilene.IN.txt'.format(get_outdir("merged_methylation_data")),
            Gifnfo='{}/groupInfo.txt'.format(get_outdir("merged_methylation_data"))
        params:
            statdir=os.path.join(outdir,'{}'.format(get_outdir("merged_methylation_data"))),
            sampleSheet=sampleSheet,
            importfunc = os.path.join(workflow_rscripts, "WGBSstats_functions.R")
        log:
            err='{}/logs/prep_for_stats.err'.format(get_outdir("singleCpG_stats_limma")),
            out='{}/logs/prep_for_stats.out'.format(get_outdir("singleCpG_stats_limma"))
        threads: 1
        conda: CONDA_WGBS_ENV
        shell: "Rscript --no-save --no-restore " + os.path.join(workflow_rscripts,'WGBSpipe.prep_data_for_stats.R ') + "{params.statdir} {params.sampleSheet} "  + os.path.join(outdir,"methXT") + " {params.importfunc} 1>{log.out} 2>{log.err}"


    rule CpG_stats:
        input: 
            Limdat='{}/limdat.LG.RData'.format(get_outdir("merged_methylation_data"))
        output:
            RDatAll='{}/singleCpG.RData'.format(get_outdir("singleCpG_stats_limma")),
            sinfo='{}/sessionInfo.txt'.format(get_outdir("singleCpG_stats_limma"))
        params:
            statdir=os.path.join(outdir,'{}'.format(get_outdir("singleCpG_stats_limma"))),
            datdir=os.path.join(outdir,'{}'.format(get_outdir("merged_methylation_data"))),
            sampleSheet=sampleSheet,
            diff=minAbsDiff,
            fdr=FDR,
            importfunc = os.path.join(workflow_rscripts, "WGBSstats_functions.R")
        log:
            err='{}/logs/CpG_stats.err'.format(get_outdir("singleCpG_stats_limma")),
            out='{}/logs/CpG_stats.out'.format(get_outdir("singleCpG_stats_limma"))
        threads: 1
        conda: CONDA_WGBS_ENV
        shell: "Rscript --no-save --no-restore " + os.path.join(workflow_rscripts,'WGBSpipe.singleCpGstats.limma.R ') + "{params.statdir} {params.sampleSheet} {params.datdir} {params.diff} {params.fdr} {params.importfunc} 1>{log.out} 2>{log.err}"


    rule CpG_report:
        input: 
            Limdat='{}/limdat.LG.RData'.format(get_outdir("merged_methylation_data")),
            sinfo='{}/singleCpG.RData'.format(get_outdir("singleCpG_stats_limma"))
        output:
            html='{}/Stats_report.html'.format(get_outdir("singleCpG_stats_limma"))
        params:
            statdir=os.path.join(outdir,'{}'.format(get_outdir("singleCpG_stats_limma"))),
            sampleSheet=sampleSheet,
            importfunc = os.path.join(workflow_rscripts, "WGBSstats_functions.R"),
            stat_cat="single_CpGs",
            rmd_in=os.path.join(workflow_rscripts,"WGBS_stats_report_template.Rmd"),
            rmd_out=os.path.join(outdir,"aux_files", "WGBS_stats_report_template.Rmd"),
            outFull=lambda wildcards,output: os.path.join(outdir,output.html)
        log:
            err='{}/logs/stats_report.err'.format(get_outdir("singleCpG_stats_limma")),
            out='{}/logs/stats_report.out'.format(get_outdir("singleCpG_stats_limma"))
        conda: CONDA_RMD_ENV
        threads: 1
        shell: "cp -v {params.rmd_in} {params.rmd_out} ;Rscript -e 'rmarkdown::render(\"{params.rmd_out}\", params=list(outdir=\"{params.statdir}\", input_func=\"{params.importfunc}\", stat_category=\"{params.stat_cat}\",sample_sheet=\"{params.sampleSheet}\"), output_file=\"{params.outFull}\")' 1>{log.out} 2>{log.err}"


    rule run_metilene:
        input:
            MetIN='{}/metilene.IN.txt'.format(get_outdir("merged_methylation_data")),
            Ginfo='{}/groupInfo.txt'.format(get_outdir("merged_methylation_data"))
        output:
            MetBed='{}/singleCpG.metilene.bed'.format(get_outdir("metilene_out"))
        params:
            DMRout=os.path.join(outdir,'{}'.format(get_outdir("metilene_out"))),
            maxD=maxDist,
            minCG=minCpGs,
            minMD=minMethDiff
        log:
            err="{}/logs/run_metilene.err".format(get_outdir("metilene_out"))
        threads: nthreads
        conda: CONDA_WGBS_ENV
        shell: 'Gi=($(cat {input.Ginfo}));metilene -a ' + " ${{Gi[0]}} " + " -b  ${{Gi[1]}} -M {params.maxD} -m {params.minCG} -d {params.minMD} -t {threads} {input.MetIN} | sort -k 1,1 -k2,2n > {output.MetBed} 2>{log.err}"


    rule get_CG_metilene:
        input:
            refG=refG,
            MetBed='{}/singleCpG.metilene.bed'.format(get_outdir("metilene_out")),
            imdF="aux_files/"+re.sub('.fa*','.CpG.bed',os.path.basename(refG))
        output:
            MetCG=os.path.join("aux_files",re.sub('_sampleSheet.[a-z]{3}$','.metilene.CpGlist.bed',os.path.basename(sampleSheet)))
        params:
            auxdir=os.path.join(outdir,"aux_files")
        log:
            err="aux_files/logs/get_CG_metilene.err"
        threads: 1
        conda: CONDA_WGBS_ENV
        shell: "bedtools intersect -wa -a {input.imdF} -b {input.MetBed} > {output.MetCG}  2>{log.err}"


    rule cleanup_metilene:
        input:
            Limdat='{}/limdat.LG.RData'.format(get_outdir("merged_methylation_data")),
            MetBed='{}/singleCpG.metilene.bed'.format(get_outdir("metilene_out")),
            MetCG=os.path.join("aux_files",re.sub('_sampleSheet.[a-z]{3}$','.metilene.CpGlist.bed',os.path.basename(sampleSheet))),
            sampleSheet=sampleSheet
        output:
            LimBed='{}/singleCpG.metilene.limma_unfiltered.bed'.format(get_outdir("metilene_out")),
            LimAnnot='{}/metilene.limma.annotated_unfiltered.txt'.format(get_outdir("metilene_out"))
        params:
            DMRout=os.path.join(outdir,'{}'.format(get_outdir("metilene_out"))),
            gene_mod=genes_bed,
            diff=minAbsDiff,
            fdr=FDR,
            importfunc = os.path.join(workflow_rscripts, "WGBSstats_functions.R")
        log:
            err="{}/logs/cleanup_metilene.err".format(get_outdir("metilene_out")),
            out="{}/logs/cleanup_metilene.out".format(get_outdir("metilene_out"))
        threads: 1
        conda: CONDA_WGBS_ENV
        shell: 'Rscript --no-save --no-restore ' + os.path.join(workflow_rscripts,'WGBSpipe.metilene_stats.limma.R ') + "{params.DMRout} " + os.path.join(outdir,"{input.MetBed}") +' ' + os.path.join(outdir,"{input.MetCG}") + ' ' + os.path.join(outdir,"{input.Limdat}") + " {input.sampleSheet} {params.gene_mod} {params.diff} {params.fdr} {params.importfunc} 1>{log.out} 2>{log.err}"


    rule metilene_report:
        input: 
            MetBed='{}/singleCpG.metilene.bed'.format(get_outdir("metilene_out")),
            LimBed='{}/singleCpG.metilene.limma_unfiltered.bed'.format(get_outdir("metilene_out"))
        output:
            html='{}/Stats_report.html'.format(get_outdir("metilene_out"))
        params:
            statdir=os.path.join(outdir,'{}'.format(get_outdir("metilene_out"))),
            sampleSheet=sampleSheet,
            importfunc = os.path.join(workflow_rscripts, "WGBSstats_functions.R"),
            stat_cat="metilene_DMRs",
            rmd_in=os.path.join(workflow_rscripts,"WGBS_stats_report_template.Rmd"),
            rmd_out=os.path.join(outdir,"aux_files", "WGBS_stats_report_template.Rmd"),
            outFull=lambda wildcards,output: os.path.join(outdir,output.html)
        log:
            err='{}/logs/stats_report.err'.format(get_outdir("metilene_out")),
            out='{}/logs/stats_report.out'.format(get_outdir("metilene_out"))
        conda: CONDA_RMD_ENV
        threads: 1
        shell: "cp -v {params.rmd_in} {params.rmd_out} ;Rscript -e 'rmarkdown::render(\"{params.rmd_out}\", params=list(outdir=\"{params.statdir}\", input_func=\"{params.importfunc}\", stat_category=\"{params.stat_cat}\",sample_sheet=\"{params.sampleSheet}\"), output_file=\"{params.outFull}\")' 1>{log.out} 2>{log.err}"


if intList:
    rule get_CG_per_int:
        input:
            intList=intList,
            refG=refG,
            imdF="aux_files/"+re.sub('.fa*','.CpG.bed',os.path.basename(refG))
        output:
            outList=run_int_aggStats(intList,False)
        log:
            err="aux_files/logs/get_CG_per_int.err"
        params:
            auxshell=lambda wildcards,input,output: ';'.join(["bedtools intersect -wa -a "+ input.imdF + " -b " + bli + ' > ' + oli  for bli,oli in zip(input.intList,output.outList) ])
        threads: 1
        conda: CONDA_WGBS_ENV
        shell: "{params.auxshell} 2>{log.err}"


    if sampleSheet:
        rule intAgg_stats:
            input:
                Limdat='{}/limdat.LG.RData'.format(get_outdir("merged_methylation_data")),
                intList=intList,
                refG=refG,
                sampleSheet=sampleSheet,
                aux=run_int_aggStats(intList,False)
            output:
                outFiles=run_int_aggStats(intList,sampleSheet),
                sinfo='{}/sessionInfo.txt'.format(get_outdir("aggregate_stats_limma"))
            params:
                auxshell=lambda wildcards,input:';'.join(['Rscript --no-save --no-restore ' + os.path.join(workflow_rscripts,'WGBSpipe.interval_stats.limma.R ') + os.path.join(outdir,'{}'.format(get_outdir("aggregate_stats_limma"))) + ' ' + li +' '+ aui +' ' + os.path.join(outdir,input.Limdat) + ' '  + input.sampleSheet + ' ' + str(minAbsDiff) + ' ' + str(FDR) + ' ' + os.path.join(workflow_rscripts, "WGBSstats_functions.R")  for li,aui in zip(intList,[os.path.join(outdir,"aux_files",re.sub('.fa',re.sub('.bed','.CpGlist.bed',os.path.basename(x)),os.path.basename(refG))) for x in intList])])
            log:
                err="{}/logs/intAgg_stats.err".format(get_outdir("aggregate_stats_limma")),
                out="{}/logs/intAgg_stats.out".format(get_outdir("aggregate_stats_limma"))
            threads: 1
            conda: CONDA_WGBS_ENV
            shell: "{params.auxshell} 1>{log.out} 2>{log.err}"


        rule intAgg_report:
            input: 
                outFiles='{}/sessionInfo.txt'.format(get_outdir("aggregate_stats_limma"))
            output:
                html='{}/Stats_report.html'.format(get_outdir("aggregate_stats_limma"))
            params:
                statdir=os.path.join(outdir,'{}'.format(get_outdir("aggregate_stats_limma"))),
                sampleSheet=sampleSheet,
                importfunc = os.path.join(workflow_rscripts, "WGBSstats_functions.R"),
                stat_cat="user_intervals",
                rmd_in=os.path.join(workflow_rscripts,"WGBS_stats_report_template.Rmd"),
                rmd_out=os.path.join(outdir,"aux_files", "WGBS_stats_report_template.Rmd"),
                outFull=lambda wildcards,output: os.path.join(outdir,output.html)
            log:
                err='{}/logs/stats_report.err'.format(get_outdir("aggregate_stats_limma")),
                out='{}/logs/stats_report.out'.format(get_outdir("aggregate_stats_limma"))
            conda: CONDA_RMD_ENV
            threads: 1
            shell: "cp -v {params.rmd_in} {params.rmd_out} ;Rscript -e 'rmarkdown::render(\"{params.rmd_out}\", params=list(outdir=\"{params.statdir}\", input_func=\"{params.importfunc}\", stat_category=\"{params.stat_cat}\",sample_sheet=\"{params.sampleSheet}\"), output_file=\"{params.outFull}\")' 1>{log.out} 2>{log.err}"
        
    rule on_target_rate:
        input:
            bams=expand("bams/{sample}.PCRrm.bam",sample=samples)
        output:
            tab="custom_stats/on_target_stats.all_reads.txt",
            plot="custom_stats/on_target_stats.all_reads.pdf"
        params:
            targets=intList,
            labels = " ".join(samples)
        log:
            err="custom_stats/logs/on_target_stats.all_reads.err",
            out="custom_stats/logs/on_target_stats.all_reads.out"
        threads: nthreads
        shell:"""
            plotEnrichment -p {threads} \
                   -b {input.bams} \
                   --plotFile {output.plot}\
                   --BED {params.targets} \
                   --labels {params.labels} \
                   --plotTitle 'Fraction of reads in target regions' \
                   --outRawCounts {output.tab} \
                   --variableScales > {log.out} 2> {log.err}
            """

rule on_target_rate_mapq:
        input:
            bams=expand("bams/{sample}.PCRrm.bam",sample=samples)
        output:
            tab="custom_stats/on_target_stats.mapq20.txt",
            plot="custom_stats/on_target_stats.mapq20.pdf"
        params:
            targets=intList,
            labels = " ".join(samples)
        log:
            err="custom_stats/logs/on_target_stats.mapq20.err",
            out="custom_stats/logs/on_target_stats.mapq20.out"
        threads: nthreads
        conda: CONDA_SHARED_ENV
        shell:"""
            plotEnrichment -p {threads} \
                   -b {input.bams} \
                   --plotFile {output.plot}\
                   --BED {params.targets} \
                   --labels {params.labels} \
                   --plotTitle 'Fraction of reads in target regions' \
                   --outRawCounts {output.tab} \
                   --variableScales \
                   --minMappingQuality 20 1> {log.out} 2> {log.err}
            """

rule on_target_reads_region:
        input:
            bams=expand("bams/{sample}.PCRrm.bam",sample=samples)
        output:
            tabq20="custom_stats/on_target_stats.per_region.mapq20.tsv",
            tab="custom_stats/on_target_stats.per_region.tsv"
        params:
            targets=intList,
            labels = " ".join(samples)
        log:
            err="custom_stats/logs/on_target_stats.per_region.err",
            out="custom_stats/logs/on_target_stats.per_region.out"
        threads: nthreads
        conda: CONDA_SHARED_ENV
        shell:"""
            multiBamSummary BED-file \
                -b {input.bams} \
                --BED {params.targets} \
                --outRawCounts tmp.q20.tsv \
                --minMappingQuality 20 \
                --labels {params.labels} \
                -p {threads} 1> {log.out} 2> {log.err};
            sort -k1,1 -k2,2n tmp.q20.tsv > {output.tabq20}; rm tmp.q20.tsv; 
            multiBamSummary BED-file \
                -b {input.bams} \
                --BED {params.targets} \
                --outRawCounts tmp.tsv \
                --labels {params.labels} \
                -p {threads} 1>> {log.out} 2>> {log.err};
            sort -k1,1 -k2,2n tmp.tsv > {output.tab}; rm tmp.tsv;
            """

rule methyl_extract_custom:
        input:
            rmDupbam="bams/{sample}.PCRrm.bam",
            sbami="bams/{sample}.PCRrm.bam.bai",
            refG=refG
        output:
            methTab="custom_stats/{sample}_CpG.bedGraph",
            meanTab="custom_stats/{sample}.mean_methyl_per_region.tsv"
        params:
            mbias_ignore=mbias_ignore,
            targets=intList,
            OUTpfx=lambda wildcards,output: os.path.join(outdir,re.sub('_CpG.bedGraph','',output.methTab))
        log:
            err="custom_stats/logs/{sample}.methyl_extract.err",
            out="custom_stats/logs/{sample}.methyl_extract.out"
        threads: nthreads
        conda: CONDA_WGBS_ENV
        shell: """
            MethylDackel extract -o {params.OUTpfx} -l {params.targets} \
                -q 20 -p 20 --minDepth 1 --mergeContext -@ {threads} \
                {input.refG} {input.rmDupbam} 1>{log.out} 2>{log.err};
            bedtools map -a {params.targets} -b {output.methTab} \
                -c 4 -o mean -prec 4 > {output.meanTab} 2>>{log.err}
            """

rule per_base_cov_custom:
    input:
        bams=expand("bams/{sample}.PCRrm.bam",sample=samples)
    output:
        "custom_stats/coverage_per_base.targets.bed"
    params:
        targets=intList
    conda: CONDA_SHARED_ENV
    log:
        err="custom_stats/logs/coverage_per_base.targets.err",
        out="custom_stats/logs/coverage_per_base.targets.out"
    shell:"""
        cat <(echo -e '#chr\tpos\t'$(echo '{input.bams}' | tr ' ' '\n' | sed 's/.*\///' | sed 's/.PCRrm.bam//g' | tr '\n' '\t')) <(samtools depth -a -q 20 -Q 20 {input.bams} -b <(cat {params.targets} | awk '{{OFS="\t";print $1,$2-1,$3-1}}') ) > {output} 2>{log.err}
        """

rule target_cpgs:
    input:
        "aux_files/genome.CpG.bed"
    params:
        targets=intList
    output:
        "custom_stats/targets.CpG.bed"
    log:
        err="custom_stats/logs/targets.CpG.err"
    conda: CONDA_WGBS_ENV
    shell:"""
        bedtools intersect -a <(cat {params.targets} | awk '{{OFS="\t";$2=$2-1;$3=$3+2; print $0}}') -b {input} -wo | awk '{{OFS="\t"; print $4,$5,$6,$1"_"$2,0,$7,$8}}' > {output} 2>{log.err}
    """

rule target_cpg_coverage:
    input:
        bams=expand("bams/{sample}.PCRrm.bam",sample=samples),
        cpg="custom_stats/targets.CpG.bed"
    output:
        "custom_stats/targets.CpG.coverage.txt"
    log:
        err="custom_stats/logs/targets.CpG.coverage.err"
    conda: CONDA_SHARED_ENV
    shell: """
        cat <(echo -e '#chr\tpos\t'$(echo '{input.bams}' | tr ' ' '\n' | sed 's/.*\///' | sed 's/.PCRrm.bam//g' | tr '\n' '\t')) <(samtools depth -a -q 20 -Q 20 -b {input.cpg} {input.bams}) > {output} 2>{log.err}
    """
    

rule mean_target_coverage:
    input:
        "custom_stats/coverage_per_base.targets.bed"
    output:
        "custom_stats/mean_coverage_per_base.targets.bed"
    params:
        targets=intList
    log:
        err="custom_stats/logs/mean_coverage_per_base.targets.err"
    conda: CONDA_WGBS_ENV
    shell:"""
        cat <(cat {input} | head -n1 | awk '{{OFS="\t";$2="start\tend";print $0}}') \
            <(bedtools map -a {params.targets} -b <(cat {input} | \
              awk '{{OFS="\t";$2=$2-1"\t"$2; print $0}}') \
            -c $(cat {input} | awk '{{}}END{{for (i=4;i<=NF;i++){{printf i","}}; print NF+1}}') \
            -o mean -prec 5) > {output} 2>{log.err}
        """

rule mean_methyl_per_region:
    input:
        tsv=expand("custom_stats/{sample}.mean_methyl_per_region.tsv",sample=samples),
        tab="custom_stats/on_target_stats.per_region.mapq20.tsv",
        cpg=expand("custom_stats/{sample}_CpG.bedGraph",sample=samples)
    output:
        "custom_stats/mean_methyl_per_region.tsv"
    params:
        indir="custom_stats/",
        script = os.path.join(workflow_rscripts,"merge_methyl_data.R")
    log:
        err="custom_stats/logs/mean_methyl_per_region.err",
        out="custom_stats/logs/mean_methyl_per_region.out"
    conda: CONDA_WGBS_ENV
    shell:"""
            Rscript {params.script} {params.indir} {output} 2> {log.err}
        """
