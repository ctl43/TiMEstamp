# TiMEstamp: Timing Mobile ELement Insertion from Multiple Sequence Alignments
TiMEstamp is an R package for inferring evolutionary timepoints from multiple sequence alignments.

## INSTALLATION:
1. Clone from GitHub:
git clone https://github.com/ctl43/TiMEstamp.git

2. Install Using R CMD:
R CMD INSTALL TiMEstamp

## USAGE EXAMPLE:
### Get Alignment Coverage:
Calculate the clade-level alignment coverage of all annotation in RepeatMasker.

```r
get_alignment_coverage(
tree_file = "msa_464_genomes.nh",
rmsk_file = "rmsk_hg38.fa.out",
consider_missing_loci = TRUE,
folder = "timestamp",
chrom = NULL
)
```

Output: 
```
output_folder/
├── aln_coverage/
│   └── # Contains alignment coverage results at sister-clade levels
├── large_gap_len/
│   └── # Contains lengths of upstream and downstream regions around orthologous loci annotations
├── buffered_large_gap_len/
│   └── # Contains lengths of upstream and downstream regions around orthologous loci annotations, excluding regions immediately adjacent to the annotation (with a defined buffer distance upstream and downstream)
├── sister_clades.rds
│   └── # A single R object that contains a list of sister-taxon/clades relative to the reference species
└── rmsk_range.rds
    └── # A a single R object that contains a prased repeatMasker annotation
```

### Get Timepoints:
Infer the insertion timepoint of annotation in RepeatMasker.

```r
get_timepoint(
folder = "timestamp",
consider_missing_loci = TRUE,
chrom = NULL)
```

Output:
```
output_folder/
├── timepoint/
│   └── # Contains inferred timepoint of each region by chromosome.
└── timepoint.rds
    └── # A a single R object that contains all timepoint and prased repeatMasker annotation
```


## EXPLANATION OF PARAMETERS:
`get_alignment_coverage()` Parameters:
- `tree_file` (str): Path to the phylogenetic tree file (.nh format).
- `rmsk_file` (str): Path to RepeatMasker output file (.out format).
- `consider_missing_loci` (bool): If TRUE, infers the absence of orthologous loci in other species.
- `folder` (str): Folder for output results.
- `chrom` (str or NULL): Chromosome to analyze. If NULL, analyzes all chromosomes.

`get_timepoint`() Parameters:
- `folder` (str): Folder with alignment coverage results.
- `consider_missing_loci` (bool): If TRUE, includes inferred missing loci.
- `chrom` (str or NULL): Chromosome to analyze. If NULL, analyze all chromosomes.
