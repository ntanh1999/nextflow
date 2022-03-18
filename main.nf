nextflow.enable.dsl = 2

params.reads = "data/*_{1,2}.fastq.gz"
params.assembly = 'data/*.fasta'
params.result = 'results'

include {
    FASTQC_READS;
    MULTIQC_READS;
    TRIM_PE_TRIMMOMATIC;
    ASSEMBLE_SPADES;
    ASSEMBLY_EVAL;
    ANNOTATE_PROKKA;
    PANGENOME_ANALYSIS_ROARY
} from './modules.nf'

workflow {
    read_pairs_ch = Channel.fromFilePairs(params.reads)
    assembly_ch = Channel.fromFilePairs(params.assembly)
    
    // Quality control
    FASTQC_READS(read_pairs_ch)
    MULTIQC_READS(FASTQC_READS.out)
    
    // Trim and assemble
    TRIM_PE_TRIMMOMATIC(read_pairs_ch)
    ASSEMBLE_SPADES(TRIM_PE_TRIMMOMATIC.out)
    assembly_ch.join(ASSEMBLE_SPADES.out.spades_assembly)
    ASSEMBLY_EVAL(assembly_ch)

    // annotation
    ANNOTATE_PROKKA(assembly_ch)

    // pan-genome analysis
    PANGENOME_ANALYSIS_ROARY(ANNOTATE_PROKKA.out)
    
}