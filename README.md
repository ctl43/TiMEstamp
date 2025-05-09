# TiMEstamp: Timing Mobile ELement Insertion from Multiple Sequence Alignments
TiMEstamp is an R package for inferring insertion timepoints of mobile elements from multiple sequence alignments.

## INSTALLATION:
1. Clone from GitHub:
git clone https://github.com/ctl43/TiMEstamp.git

2. Install Using R CMD:
R CMD INSTALL TiMEstamp

## PREPARATION OF MULTIPLE SEQUENCE ALIGNMENT (MSA):
1. Download available MSA data in MAF format (e.g., https://hgdownload.soe.ucsc.edu/goldenPath/hg38/multiz470way/) or generate the multiple genome alignment by yourself.
2. Split the genome alignment by species using PHAST.
```mafSpeciesSubset -keepFirst```
3. Remove gaps relative to the reference genome, which should always be the first genome, using PHAST. The result files will be in FASTA format.
```msa_view --gap-strip 1```
P.S.: This prepreation step will be incorporated into TiMEstamp without relying PHAST.

## EXTRACT GAPS RELATIVE TO THE REFERENCE GENOME:
Unaligned regions, indicated by asterisks (*) or hyphens (-) in the resulting FASTA files, were extracted and annotated in BED file format.
```r
get_missing_gap(fa, out_folder)
import_missing_ranges_by_chrom(CHROM, in_folder, out_folder) # Combine all bed files into a single RDS file for R to read.
```

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
│   └── # Contains the unaligned lengths of the upstream and downstream flanking regions surrounding the region of interest (ROI) in orthologous loci.
├── buffered_large_gap_len/
│   └── # Contains the unaligned lengths of the upstream and downstream flanking regions surrounding the ROI±length buffer in orthologous loci, 
├── sister_clades.rds
│   └── # A single R object that contains a list of sister-taxon/clades relative to the reference species
└── rmsk_range.rds
    └── # A a single R object that contains prased repeatMasker annotation
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
