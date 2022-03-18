// Run QC process for pair-end reads
process FASTQC_READS{
    tag "$sample_id"
    publishDir "$params.results/$sample_id", mode: 'copy'

    input:
        tuple val(sample_id), path(reads)
    output:
        path 'fastqc'
    script:
    """
    fastqc -t $task.cpus -o fastqc $reads
    """
}

// create multiqc report from fastqc results of all samples
process MULTIQC_READS{
    publishDir "$params.results", mode: 'copy'
    
    input:
        path fastqc
    
    script:
    def all_fastqc = fastqc.collect{"$it"}.join(' ')
    """
    multiqc -o multiqc $all_fastqc
    """
}

// trim adapters and low quality reads
process TRIM_PE_TRIMMOMATIC{
    tag "$sample_id"
    publishDir "$params.results/$sample_id/trimmomatic", mode: 'copy'

    input:
        tuple val(sample_id), path(reads)
    
    output:
        tuple \
            val(sample_id), \
            path("${sample_id}_1.fastq.gz"), \
            path("${sample_id}_2.fastq.gz")

    script:
    """
    trimmomatic.sh PE -threads $task.cpus \
        $reads \
        ${sample_id}_1.fastq.gz ${sample_id}_S1.fastq.gz \
        ${sample_id}_2.fastq.gz ${sample_id}_S2.fastq.gz \
        ILLUMINACLIP:TruSeq3-PE-2.fa:2:30:10:2:true \
        SLIDINGWINDOW:4:15 LEADING:3 TRAILING:3 MINLEN:36
    """
}

// assemble using SPAdes
process ASSEMBLE_SPADES{
    tag "$sample_id"
    publishDir "$params.results/$sample_id", mode: 'copy'

    input:
        tuple val(sample_id), path(read1), path(read2)
    
    output:
        path "spades/${sample_id}.fasta", emit: spades_assembly
    
    script:
    """
    spades.py -m $task.memory -t $task.cpus --careful -o spades -1 $read1 -2 $reads2
    mv spades/contigs.fasta spades/${sample_id}.fasta
    """
}

// Run QC process for assembly output
process ASSEMBLY_EVAL{
    tag "$sample_id"
    publishDir "$params.results/$sample_id", mode: 'copy'

    input:
        tuple val(sample_id), path(assembly)

    script:
    """
    quast.py -t $task.cpus -o quast $assembly
    """
}

// annotate genome by prokka
process ANNOTATE_PROKKA{
    tag "$sample_id"
    publishDir "$params.results/$sample_id", mode: 'copy'
    input:
        tuple val(sample_id), path(assembly)

    output:
        tuple val(sample_id), path("prokka/${sample_id}.gff")
    script:
    """
    prokka --force --cpus $task.cpus --addgenes --mincontiglen 200 \
    --prefix $sample_id --locus $sample_id --outdir prokka $assembly
    """
}

// run pan-genome analysis by Roary
process PANGENOME_ANALYSIS_ROARY{
    publishDir "$params.results", mode: 'copy'
    
    input:
        tuple val(sample_id), path(gff)
    
    script:
    def all_gff = gff.collect{"$it"}.join(' ')
    """
    roary -v -z -p $task.cpus -f roary $all_gff
    """
}
