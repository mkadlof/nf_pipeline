
// Pipeline configuration
params.picardPath = "/opt/picard/picard.jar"
params.gatkPath = "/opt/gatk/gatk"

// Process - we keep them in modules in ./modules dir and include them here.

// Note:
// Nextflow requires every process to have a unique name. To reuse the same process we need to
// include it multiple times, assigning a different alias each time.
// Docs: https://www.nextflow.io/docs/latest/module.html#module-aliases
include { samtoolsIndex as samtoolsIndex_1 } from './modules/samtoolsIndex.nf'
include { samtoolsIndex as samtoolsIndex_2 } from './modules/samtoolsIndex.nf'
include { samtoolsIndex as samtoolsIndex_3 } from './modules/samtoolsIndex.nf'
include { samtoolsSort as samtoolsSort_1 } from './modules/samtoolsSort.nf'
include { samtoolsSort as samtoolsSort_2 } from './modules/samtoolsSort.nf'

include { bwaIndex } from './modules/bwaIndex.nf'
include { bwaMapping } from './modules/bwaMapping.nf'
include { samtoolsViewFilter } from './modules/samtoolsViewFilter.nf'
include { simpleFilterAmpliconMk } from './modules/simpleFilterAmpliconMk.nf'
include { fixReadGroups } from './modules/fixReadGroups.nf'
include { gatkHaplotypeCaller } from './modules/gatkHaplotypeCaller.nf'
include { samtoolsFaidx } from './modules/samtoolsFaidx.nf'
include { gatkCreateSequenceDictionary } from './modules/gatkCreateSequenceDictionary.nf'

// Include modules for statistics
include { genomeStats } from './modules/genomeStats.nf'
include { bamStats as bamStats_1 } from './modules/bamStats.nf'
include { bamStats as bamStats_2 } from './modules/bamStats.nf'
include { bamStats as bamStats_3 } from './modules/bamStats.nf'
include { combineJsonFiles_x4 } from './modules/combineJsonFiles.nf'

// Workflow definition

workflow {
    reference_genome = Channel.fromPath(params.reference_genome)
    reads = Channel.fromFilePairs(params.reads)
    genomeStats(reference_genome)
    index = bwaIndex(reference_genome)
    bwaMapping(reference_genome, index, reads)
    samtoolsSort_1(bwaMapping.out)
    samtoolsIndex_1(samtoolsSort_1.out)
    samtoolsViewFilter(bwaMapping.out)
    samtoolsSort_2(samtoolsViewFilter.out)
    samtoolsIndex_2(samtoolsSort_2.out)
    bamStats_2(samtoolsSort_2.out, samtoolsIndex_2.out)
    simpleFilterAmpliconMk(samtoolsSort_2.out, samtoolsIndex_2.out)
    bamStats_1(samtoolsSort_1.out, samtoolsIndex_1.out)
    samtoolsIndex_3(simpleFilterAmpliconMk.out)
    fixReadGroups(simpleFilterAmpliconMk.out)
    samtoolsFaidx(reference_genome)
    gatkCreateSequenceDictionary(reference_genome)
    gatkHaplotypeCaller(reference_genome, samtoolsFaidx.out, gatkCreateSequenceDictionary.out, fixReadGroups.out)
    bamStats_3(simpleFilterAmpliconMk.out, samtoolsIndex_3.out)
    combineJsonFiles_x4(genomeStats.out, bamStats_1.out, bamStats_2.out, bamStats_3.out)
}
