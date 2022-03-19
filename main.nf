// Enable DSL 2 syntax
nextflow.enable.dsl = 2

// Define the default parameters
params.reads = "data/*_{1,2}.fastq.gz"
params.assembly = "data/*.fna"
params.result = "results"
params.adapter = "/home/ted/nextflow/adapters/TruSeq3-PE-2.fa"

// Import modules
include {
    FASTQC_READS;
    MULTIQC_READS;
    TRIM_PE_TRIMMOMATIC;
    ASSEMBLE_SPADES;
    ASSEMBLY_EVAL;
    ANNOTATE_PROKKA;
    PANGENOME_ANALYSIS_ROARY
} from './modules.nf'


//pipeline
workflow {
    // create input channel
    read_pairs_ch = Channel.fromFilePairs(params.reads)
    assembly_ch = Channel
        .fromPath(params.assembly)
        .map{tuple(it.baseName, it)}

    // Quality control of reads
    FASTQC_READS(read_pairs_ch)
    MULTIQC_READS(FASTQC_READS.out.collect())
    
    // Trim reads
    TRIM_PE_TRIMMOMATIC(read_pairs_ch, params.adapter)

    // Assemble reads
    ASSEMBLE_SPADES(TRIM_PE_TRIMMOMATIC.out[0])

    // mix spades assembly and input assembly
    mix_assembly_ch = ASSEMBLE_SPADES.out.spades_assembly.mix(assembly_ch)

    // evaluate the quality of assembly
    ASSEMBLY_EVAL(mix_assembly_ch)

    // genome annotation
    ANNOTATE_PROKKA(mix_assembly_ch)

    // pan-genome analysis
    PANGENOME_ANALYSIS_ROARY(ANNOTATE_PROKKA.out.gff.collect())
    
}