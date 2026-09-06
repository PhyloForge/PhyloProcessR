---
title: 'PhyloProcessR: A Modular Toolkit for Targeted Sequence-Capture Phylogenomics'
tags:
  - R
  - phylogenomics
  - sequence capture
  - target enrichment
  - bioinformatics
authors:
  - name: Carl R. Hutter
    orcid: 0000-0001-6381-6339
    corresponding: true
    affiliation: "1, 2"
affiliations:
  - index: 1
    name: Louisiana State University, Museum of Natural History, 119 Foster Hall, Baton Rouge, LA 70803, USA
  - index: 2
    name: Virus and Prion Research Unit, National Animal Disease Center, Agricultural Research Service, United States Department of Agriculture, Ames, IA 50010, USA
date: 28 August 2026
bibliography: paper.bib
---

# Summary

Targeted sequence capture allows researchers to sequence selected genomic regions across many organisms, but the resulting reads do not arrive as analysis-ready genes. Studies differ in their starting material, marker design, treatment of heterozygosity and paralogy, and requirements for downstream matrices. Because each processing decision can alter which samples and loci are retained, no single sequence of commands can represent every defensible analysis; consistent curation and explicit documentation are both essential.

PhyloProcessR functions both as a ready-to-use target-capture workflow and as a modular toolkit for constructing reproducible, project-specific phylogenomic pipelines. Complete reference workflows demonstrate common analyses, while individual operations are exposed as callable R functions. Users can combine components, omit unnecessary stages, substitute compatible approaches, or reorder operations where data dependencies permit. The package spans raw-read processing, assembly and target recovery, variant-aware consensus generation, alignment, quality control and filtering, paralog assessment, dataset construction, legacy-data integration, and recovery of novel shared loci. Recording these choices in executable R scripts supports reproducible and extensible analyses without requiring every project to follow the same path.

# Statement of need

Target-capture studies commonly include hundreds of samples and hundreds to thousands of loci. Although mature programs exist for individual tasks, researchers must still preserve sample identities, target-to-gene relationships, file formats, parameter choices, and quality-control decisions across many programs and directories. The challenge is therefore not simply to run an assembler or aligner; it is to connect heterogeneous tools while consistently handling contamination, incomplete capture, heterozygosity, paralogy, missing data, and poorly aligned sequence.

Projects also rarely require an identical chain of operations. One study may begin with assembled contigs, another may require only alignment curation, and another may combine newly captured loci with Sanger or GenBank sequences. Reproducibility therefore depends on documenting not only parameters within a pipeline, but also which operations were selected and how they were ordered. PhyloProcessR addresses this need by exposing its processing operations as interoperable R functions. The supplied workflows are complete reference implementations, not a mandatory control path: researchers can call functions directly, omit stages, and replace compatible components while preserving the analytical sequence in a project script.

This modular interface is paired with broad phylogenomic coverage. Functions support raw-read preprocessing, assembly and target recovery, mapping and variant calling, IUPAC or phased-haplotype consensus generation, alignment and filtering, capture and depth assessment, paralog evaluation, grouping multiple targets by gene, unlinked-dataset construction, legacy-data integration, and recovery of shared genomic regions from reads not assigned to designed targets. A project-specific R script can therefore reuse only the relevant components, retain intermediate products and standard sequence formats, and be modified or extended as data and research questions change.

# State of the field

PhyloProcessR enters a mature and rapidly developing software landscape. HybPiper emphasizes recovery of coding sequences and flanking introns from target-enrichment reads and flags potential paralogs [@johnson2016hybpiper]. PHYLUCE provides assembly, locus identification, alignment, and dataset preparation for conserved and ultraconserved elements [@faircloth2016phyluce]. SECAPR processes raw sequence-capture reads into filtered alignments and includes quality-control and allele-phasing functions [@andermann2018secapr]. Captus accepts target capture, genome skimming, RNA sequencing, and whole-genome data, then assembles, extracts, filters, and aligns markers [@ortiz2026captus]; direct comparisons show that pipeline architecture can affect speed, locus recovery, and downstream datasets [@raza2023pipelines]. More recently, HybSuite has integrated HybPiper-centered processing, paralog handling, and tree inference in a read-to-tree workflow [@liu2026hybsuite]. Consequently, PhyloProcessR should not be described as the first or only end-to-end target-capture pipeline.

These programs provide configurable parameters, and several support staged, restartable, or otherwise customized execution; the relevant distinction is not between flexible and rigid software. Many organize analyses primarily around named pipeline stages or command-line workflows designed for their core recovery and assembly tasks. PhyloProcessR has a different center of gravity: it deliberately exposes a broad collection of phylogenomics-specific transformations and assessments as independently callable R functions, with complete workflows supplied as examples of how those components can be assembled. Individual functions can be used together, added to outputs from another pipeline, or incorporated into a project-specific R workflow.

This architecture allows researchers to retain preferred assemblers, mappers, aligners, or variant callers while using common functions for variant-aware consensus processing, batch alignment curation, target-to-gene organization, capture and depth assessment, paralog summaries, legacy-data integration, and shared-locus discovery. PhyloProcessR is not presented as universally faster or more accurate than HybPiper, PHYLUCE, SECAPR, Captus, or HybSuite, and this draft reports no head-to-head benchmark. Its narrower contribution is an R-centered combination of ready-to-use workflows and composable phylogenomic operations. Evaluation should therefore consider recovery, data retention, reproducibility, and the consequences of curation choices on a defined dataset rather than treating function count as evidence of superiority.

# Software design

PhyloProcessR provides two complementary interfaces: staged reference workflows and an independently callable R function set. The repository includes workflow and configuration scripts for preprocessing, assembly, variant calling, alignment, trimming and dataset construction, joint genotyping, capture assessment, legacy integration, and novel-locus discovery. These scripts demonstrate complete analyses, but the underlying functions can also be assembled directly in project-specific R scripts. Users may begin at an intermediate data product, omit unneeded stages, substitute compatible methods, and reorder operations where their input-output dependencies permit.

Version 1.0.0 exports 103 user-callable R functions. These functions connect domain-specific processing with established external programs: R code manages files, parameters, sequence objects, summaries, and transformations, whereas external programs perform read cleaning, mapping, de novo assembly, similarity searching, alignment, and variant calling. Distinctive operations include `extractGenomeTarget()`, which queries NCBI by taxon name, taxon identifier, or assembly accession and recovers target loci from matching GenBank genome assemblies. `integrateLegacy()` matches Sanger or GenBank alignments to corresponding capture loci and can merge sequences from samples represented in both datasets. `discoverSharedRegions()`, `assembleSharedRegions()`, and `collectNovelContigs()` use reads not assigned to known targets to identify, assemble, and consolidate shared candidate loci. Alignment curation includes the configurable `superTrimmer()` framework; reading-frame-aware `trimExonORF()`; native coverage, similarity, and segment filters; the native R `trimSampleHMM()` profile-HMM cleaner; and an interface to TrimAl. Additional functions expose variant-aware consensus generation and paralog assessment. Standard FASTA, PHYLIP, alignment, and intermediate files provide practical boundaries for composing these operations with alternative tools.

Explicit construction in R makes the selected operations and their order part of an executable analytical record that can be reviewed, versioned, reused, and extended. Configuration files expose consequential parameters, parallel execution is available for many batch operations, and text summaries and intermediate files allow failed samples or loci to be investigated. Each workflow stage requires only the external programs it invokes, while complete Docker, Apptainer/Singularity, and Conda environments support full analyses; the R package itself can be installed and tested without every executable [@hutter2026phyloprocessr]. Core alignment, consensus, input/output, and native-trimming behavior is covered by automated tests and an `R CMD check` workflow. A compact [seven-sample laptop example](https://github.com/PhyloForge/PhyloProcessR/tree/master/examples/laptop-example) demonstrates modular workflow construction with five single-lane and two multilane samples. Its reduced FASTQ data, 40-locus target panel, and expected biological outputs are released under CC0 1.0. From 52,138 paired reads, the workflow recovered 39 loci in at least four samples and generated 38 candidate alignments. Twenty alignments were selected by occupancy, missingness, length, and locus name; all 20 passed target-region trimming and quality control, with no taxa removed. Final alignments contained five to seven taxa. In clean-install validation, processing through target filtering took approximately 17 minutes and target-region trimming and quality control took under 30 seconds using two threads and a 4 GB memory cap. These measurements demonstrate reproducibility and laptop feasibility, not comparative performance.

# Research impact statement

PhyloProcessR has been developed publicly since 2020 and has already supported research beyond a hypothetical demonstration. The current PhyloProcessR package was used for adapter filtering, contig assembly, and alignment export in the VenomCap exon-capture study of snake venom genes, and for processing raw FrogCap data in a genomic analysis of cryptic fanged-frog species [@travers2024venomcap; @chan2026genomic]. University of Florida Research Computing also distributes PhyloProcessR as an environment module on the HiPerGator high-performance computing system [@ufit2026phyloprocessr].

The software lineage also includes earlier R workflows distributed with the FrogCap resources [@hutter2022frogcap] through the FrogCap-Sequence-Capture repository and developed under the PhyloCap name. These precursor workflows supported read processing, assembly, alignment construction, and/or variant discovery in several published frog target-capture studies [@chan2020larger; @chan2020mirage; @chan2020target; @chan2022geneflow; @chan2025deforestation; @hutterduellman2023]. FrogCap data and associated precursor processing also contributed to studies of Malagasy mantellid frogs, including the *Mantidactylus ambreensis* complex, historical giant-stream-frog material, the diversity of *Brygoomantis*, and the *Guibemantis liber* complex [@rasolonjatovo2020ambreensis; @rancilhac2020historical; @scherz2022brygoomantis; @koppetsch2023guibemantis]. These papers are distinguished from direct uses of the packaged PhyloProcessR repository because they cite the antecedent bioinformatics-pipeline and variant-pipeline scripts or use FrogCap datasets, but PhyloProcessR evolved from these resources. Notably, the 2025 application was authored by Chan, Hime, and Brown and therefore documents use of the precursor workflows in a publication without Hutter as an author.

The repository provides a GPL-3.0-or-later license, citation metadata, contribution and conduct guidance, function-level documentation, tutorials, automated tests, continuous integration, and reproducible environment definitions. These materials, together with multi-year iterative development and published use, indicate that the software is maintained research infrastructure rather than a one-off analysis script. Version 1.0.0 is archived at Zenodo [@hutter2026phyloprocessr].

The immediate research value is practical: published target-capture projects can record a versioned, inspectable processing workflow instead of citing an unpublished collection of scripts. The compact review dataset, instructions, and expected outputs are distributed with [the repository example](https://github.com/PhyloForge/PhyloProcessR/tree/master/examples/laptop-example) and exercise mixed single- and multilane inputs, assembly, target recovery, alignment, trimming, and quality control without requiring a cluster.

# AI usage disclosure

PhyloProcessR, including its scientific conception, workflow design, and source code, was authored by CRH. OpenAI ChatGPT and Codex were used to assist with debugging, document formatting, software documentation, and manuscript proofreading. These tools did not originate the package's scientific ideas or independently author its code. CRH reviewed and edited all AI-assisted outputs and retains responsibility for the scientific content, software behavior, licensing, and design decisions.

# Acknowledgements

CRH was supported by the U.S. National Science Foundation Graduate Research Fellowship Program (grant numbers 1540502, 1451148, and 0907996) and a National Science Foundation Postdoctoral Research Fellowship in Biology (grant number 2010988).

# References
