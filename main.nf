nextflow.enable.dsl = 2

params.reads = "data/*_{1,2}.fastq.gz"
params.assembly = 'data/*.fna'
params.result = 'results'
params.adapter = '/home/ted/nextflow/adapters/TruSeq3-PE-2.fa'

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
    TRIM_PE_TRIMMOMATIC(read_pairs_ch, params.adapter)
    ASSEMBLE_SPADES(TRIM_PE_TRIMMOMATIC.out[0])
    ASSEMBLE_SPADES.out.spades_assembly.join(assembly_ch)
    ASSEMBLY_EVAL(ASSEMBLE_SPADES.out.spades_assembly)

    // annotation
    ANNOTATE_PROKKA(ASSEMBLE_SPADES.out.spades_assembly)

    // pan-genome analysis
    PANGENOME_ANALYSIS_ROARY(ANNOTATE_PROKKA.out[0])
    
}