// Process - we keep them in modules in ./modules dir and include them here.

// Note:
// Nextflow requires every process to have a unique name. To reuse the same process we need to
// include it multiple times, assigning a different alias each time.
// Docs: https://www.nextflow.io/docs/latest/module.html#module-aliases
include { samtoolsIndex as samtoolsIndex_1 } from './modules/samtoolsIndex.nf'
include { samtoolsIndex as samtoolsIndex_2 } from './modules/samtoolsIndex.nf'
include { samtoolsSort as samtoolsSort_1 } from './modules/samtoolsSort.nf'
include { samtoolsSort as samtoolsSort_2 } from './modules/samtoolsSort.nf'

include { bwaIndex } from './modules/bwaIndex.nf'
include { bwaMapping } from './modules/bwaMapping.nf'
include { samtoolsViewFilter } from './modules/samtoolsViewFilter.nf'
include { simpleFilterAmpliconMk } from './modules/simpleFilterAmpliconMk.nf'

// Include modules for statistics
include { genomeStats } from './modules/genomeStats.nf'
include { bamStats } from './modules/bamStats.nf'

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
    simpleFilterAmpliconMk(samtoolsSort_2.out, samtoolsIndex_2.out)
    bamStats(samtoolsSort_1.out, samtoolsIndex_1.out)
}
