
// Pipeline configuration
// Path to the tools used in case it is run not in a Docker
params.picardPath = "/opt/picard/picard.jar"
params.gatkPath = "/opt/gatk/gatk"

// Number of threads (used in several different modules)
params.threads = 5

// Number of cycles (used in simpleFilterAmpliconMk)
params.cycles = 30

// Coverage threshold (used in detectLowCoverage)
params.coverage_threshold = 20

// primers
params.primers = "primers.bed"

// Process - we keep them in modules in ./modules dir and include them here.

// Note:
// Nextflow requires every process to have a unique name. To reuse the same process we need to
// include it multiple times, assigning a different alias each time.
// Docs: https://www.nextflow.io/docs/latest/module.html#module-aliases
// E.g.:
//include { samtoolsIndex as samtoolsIndex_2 } from './modules/samtoolsIndex.nf'


// Include one-time use modules
include { samtoolsSort as samtoolsSort } from './modules/samtoolsSort.nf'
include { bwaIndex } from './modules/bwaIndex.nf'
include { bwaMapping } from './modules/bwaMapping.nf'
include { samtoolsViewFilter } from './modules/samtoolsViewFilter.nf'
include { simpleFilterAmpliconMk } from './modules/simpleFilterAmpliconMk.nf'
include { iVar } from './modules/ivar.nf'
include { fixReadGroups } from './modules/fixReadGroups.nf'
include { gatkHaplotypeCaller } from './modules/gatkHaplotypeCaller.nf'
include { detectLowCov } from './modules/detectLowCov.nf'
include { samtoolsFaidx } from './modules/samtoolsFaidx.nf'
include { gatkCreateSeqDict } from './modules/gatkCreateSeqDict.nf'
include { vcf2fasta } from './modules/vcf2fasta.nf'
include { tabix } from './modules/tabix.nf'
include { applyLowCovToRef } from './modules/applyLowCovToRef.nf'
include { maskLowCov } from './modules/maskLowCov.nf'

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
    primers = Channel.value(params.primers as Path)

    // Processes
    genomeStats(reference_genome)
    bwaIndex(reference_genome)
    bwaMapping(reference_genome, bwaIndex.out, reads)
    samtoolsSort(bwaMapping.out)
    samtoolsViewFilter(samtoolsSort.out)
    bamStats_2(samtoolsViewFilter.out)
    iVar(samtoolsViewFilter.out, primers)
    simpleFilterAmpliconMk(iVar.out)
    bamStats_1(samtoolsSort.out)
    fixReadGroups(simpleFilterAmpliconMk.out)
    detectLowCov(fixReadGroups.out)
    samtoolsFaidx(reference_genome)
    gatkCreateSeqDict(reference_genome)
    gatkHaplotypeCaller(reference_genome, samtoolsFaidx.out, gatkCreateSeqDict.out, fixReadGroups.out)
    bamStats_3(simpleFilterAmpliconMk.out)
    combineJsonFiles_x4(genomeStats.out, bamStats_1.out, bamStats_2.out, bamStats_3.out)
    tabix(gatkHaplotypeCaller.out)
    left = applyLowCovToRef(reference_genome, detectLowCov.out)
    right = vcf2fasta(reference_genome, gatkHaplotypeCaller.out, tabix.out)
    combined = left.join(right)
    maskLowCov(combined)
}