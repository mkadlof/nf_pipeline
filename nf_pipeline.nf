
// Pipeline configuration
params.picardPath = "/opt/picard/picard.jar"
params.gatkPath = "/opt/gatk/gatk"
params.threads = 5

// Process - we keep them in modules in ./modules dir and include them here.

// Note:
// Nextflow requires every process to have a unique name. To reuse the same process we need to
// include it multiple times, assigning a different alias each time.
// Docs: https://www.nextflow.io/docs/latest/module.html#module-aliases
include { samtoolsIndex as samtoolsIndex_1 } from './modules/samtoolsIndex.nf'
include { samtoolsIndex as samtoolsIndex_2 } from './modules/samtoolsIndex.nf'
include { samtoolsIndex as samtoolsIndex_3 } from './modules/samtoolsIndex.nf'

// Include one-time use modules
include { samtoolsSort as samtoolsSort } from './modules/samtoolsSort.nf'
include { bwaIndex } from './modules/bwaIndex.nf'
include { bwaMapping } from './modules/bwaMapping.nf'
include { samtoolsViewFilter } from './modules/samtoolsViewFilter.nf'
include { simpleFilterAmpliconMk } from './modules/simpleFilterAmpliconMk.nf'
include { fixReadGroups } from './modules/fixReadGroups.nf'
include { gatkHaplotypeCaller } from './modules/gatkHaplotypeCaller.nf'
include { samtoolsFaidx } from './modules/samtoolsFaidx.nf'
include { gatkCreateSequenceDictionary } from './modules/gatkCreateSequenceDictionary.nf'
include { vcf2fasta } from './modules/vcf2fasta.nf'
include { tabix } from './modules/tabix.nf'

// Include modules for statistics
include { genomeStats } from './modules/genomeStats.nf'
include { bamStats as bamStats_1 } from './modules/bamStats.nf'
include { bamStats as bamStats_2 } from './modules/bamStats.nf'
include { bamStats as bamStats_3 } from './modules/bamStats.nf'
include { combineJsonFiles_x4 } from './modules/combineJsonFiles.nf'

// Workflow definition

workflow {
    // Channels
    reference_genome = Channel.value(params.reference_genome as Path)
    reads = Channel.fromFilePairs(params.reads)

    // Processes
    genomeStats(reference_genome)
    index = bwaIndex(reference_genome)
    bwaMapping(reference_genome, index, reads)
    samtoolsSort(bwaMapping.out)
    samtoolsIndex_1(samtoolsSort.out)
    samtoolsViewFilter(samtoolsSort.out)
    samtoolsIndex_2(samtoolsViewFilter.out)
    bamStats_2(samtoolsViewFilter.out, samtoolsIndex_2.out)
    simpleFilterAmpliconMk(samtoolsViewFilter.out, samtoolsIndex_2.out)
    bamStats_1(samtoolsSort.out, samtoolsIndex_1.out)
    samtoolsIndex_3(simpleFilterAmpliconMk.out)
    fixReadGroups(simpleFilterAmpliconMk.out)
    samtoolsFaidx(reference_genome)
    gatkCreateSequenceDictionary(reference_genome)
    gatkHaplotypeCaller(reference_genome, samtoolsFaidx.out, gatkCreateSequenceDictionary.out, fixReadGroups.out)
    bamStats_3(simpleFilterAmpliconMk.out, samtoolsIndex_3.out)
    combineJsonFiles_x4(genomeStats.out, bamStats_1.out, bamStats_2.out, bamStats_3.out)
    tabix(gatkHaplotypeCaller.out)
    vcf2fasta(reference_genome, gatkHaplotypeCaller.out, tabix.out)
}
