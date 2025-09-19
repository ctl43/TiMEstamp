# TiMEstamp: Timing Mobile Element Insertions from Multiple Sequence Alignments

TiMEstamp is an R package for inferring insertion timepoints of mobile elements from whole-genome multiple sequence alignments (MSAs). It includes utilities to:

- convert per-chromosome MAF blocks into per-species FASTA aligned to a reference genome
- extract unaligned regions (“gaps”, represented by `-`) relative to the reference as BED and `GRanges` objects
- run downstream workflows to estimate insertion timepoints and explore chimeric insertions

---

## Table of Contents

- [Requirements](#requirements)
- [Installation](#installation)
- [Input data and folder layout](#input-data-and-folder-layout)
- [Preparation of Multiple Sequence Alignment (MSA)](#preparation-of-multiple-sequence-alignment-msa)
- [Usage example 1 — RepeatMasker timepoints](#usage-example-1--repeatmasker-timepoints)
- [Usage example 2 — Chimeric insertions on chr22](#usage-example-2--chimeric-insertions-on-chr22)
- [Notes and tips](#notes-and-tips)
- [Citing TiMEstamp](#citing-timestamp)

---

## Requirements

- **R** ≥ 4.2  
- **R packages**: `Rcpp`, `BiocParallel`, `rtracklayer`, `S4Vectors`, `ape`, `IRanges`, `GenomicRanges`, `XVector`, `data.table`, `tidytree`
  (Install from CRAN/Bioconductor as needed.)
- **C++17 toolchain** (for `std::filesystem`)  
  - Linux: GCC ≥ 8 (GCC 11+ recommended)  
  - macOS: Xcode Command Line Tools
- Adequate disk space for large MAF/FASTA/derived files

> The package compiles C++ sources with **C++17**. If you see build errors mentioning `std::filesystem`, ensure your compiler is set to C++17.

---
## Installation
Clone from GitHub:
```bash
git clone https://github.com/ctl43/TiMEstamp.git
R CMD INSTALL TiMEstamp
```
---

## Input data and folder layout

### Multiple sequence alignment (MAF)

Download per-chromosome MAFs for your reference assembly (e.g., hg38 multiz):

- https://hgdownload.soe.ucsc.edu/goldenPath/hg38/multiz470way/

**Requirement:** MAF files must be stored **one per reference chromosome**.

Example layout:

```
/path/to/maf_root/
  chr1.maf
  chr2.maf
  …
  chr22.maf
  chrX.maf
  chrY.maf
```

### Species list

Plain text file with one species ID per line (prefix before the dot in MAF `s` lines):

```
hg38
panTro6
gorGor6
ponAbe3
HLhylMol2
```

### Chromosome sizes

UCSC-style `chrom.sizes` for the **reference** assembly (e.g., hg38):
- https://hgdownload.cse.ucsc.edu/goldenpath/hg38/bigZips/hg38.chrom.sizes
```
chr1    248956422
chr2    242193529
…
chr22   50818468
chrX    156040895
chrY    57227415
```

### Phylogenetic tree
Provide a Phylogenetic tree in Newick file (`.nh` or `.nwk` accepted).
Tip labels must **exactly match** your species IDs (e.g., `hg38`, `panTro6`), be **unique**, and include the **reference species**.
Branch lengths are optional. (e.g., https://hgdownload.soe.ucsc.edu/goldenPath/hg38/multiz470way/hg38.470way.nh)

### Annotation file
For example, repeatMasker annotation (https://www.repeatmasker.org/genomes/hg38/rmsk4.0.5_rb20140131/hg38.fa.out.gz)

---

## Preparation of multiple sequence alignment (MSA):

### Convert MAF to per-species FASTA
Splits the alignment by species, removes columns where the **reference** has `-`, and writes per-species FASTA sequences under `<out_folder>/fasta/<chrom>/`. FASTA headers are set to the **reference chromosome name** (e.g., `>chr22`). Inter-block gaps are filled with `N` for the reference and `-` for non-reference species. Sequences are padded to the reference chromosome length from `chrom.sizes`.

```r
convert_maf_to_fasta(
  maf_folder        = "/path/to/maf_root",
  species_file      = "/path/to/species_list.txt",
  chrom_size_file   = "/path/to/hg38.chrom.sizes",
  out_folder        = "/path/to/TiMEstamp_example/hg38_multiz470way",
  threads           = 1L
)
```

**Resulting structure:**

```
<out_folder>/
  fasta/
    chr1/
      hg38.fa
      panTro6.fa
      gorGor6.fa
      …
    chr2/
      …
    chr22/
      …
```
---

### Extract gaps as BED and RDS

Scans each per-species FASTA for contiguous `-` runs (gaps) and writes:

- **BED** per chromosome → `<out_folder>/gap_bed/<chrom>.bed`  
  (0-based, half-open; columns: `chrom`, `start`, `end`, `name`, where `name` is the FASTA basename/species tag)
- **GRanges** per chromosome → `<out_folder>/gap/<chrom>.rds`  
  (`labels` factor column contains species tags with consistent levels)

```r
extract_gaps_from_fasta(
  folder  = "/path/to/TiMEstamp_example/hg38_multiz470way",
  threads = 1L
)
```

---

## Usage example 1 — RepeatMasker timepoints

This workflow estimates insertion timepoints for all RepeatMasker annotations using an MSA from `multiz470way`. It is optimized for speed and memory efficiency and is recommended when the number of annotations exceeds ~50,000, where per-locus inspection becomes impractical. The trade-off is that it does not evaluate locus-specific availability (unlike Workflow 2); however, for family- or class-level transposable element analyses it captures overall evolutionary trends with sufficient accuracy. See Supplementary Figure 2 in the original paper.

```r
# Define a working directory
work_dir <- "/path/to/TiMEstamp_example"

# 1) Prepare RepeatMasker annotation
prepare_rmsk(
  rmsk_file = "reference/rmsk_hg38.fa.out",
  folder    = work_dir
)

# 2) (Optional) Prune the phylogenetic tree to selected species
update_tree(
  tree    = "phylogenetic/hg38.470way.nh",
  folder  = work_dir,
  species = selected_species
)

# 3) Define sister clades from the updated tree
get_sister(
  tree_file = file.path(work_dir, "updated_tree.nh"),
  folder    = work_dir
)

# 4) Identify missing loci coverage by clade
get_missing_coverage_by_clade_fast(folder = work_dir)

# 5) Clean clade data and filter loci
clean_clade_data_fast(folder = work_dir)

# 6) Estimate insertion timepoints
get_timepoint_fast(folder = work_dir)
```

---

### Explanation of output

After running the pipeline, **TiMEstamp** will generate the following key outputs inside the working directory (`work_dir`):  

- **`reference_anno/`**  
  Processed RepeatMasker annotation files for the reference genome.  

- **`fast/`**  
  Intermediate results including reference annotation, clade coverage summaries, and pre-filtered gap data.  

- **`updated_tree.nh`**: Phylogenetic tree pruned to the selected species.  

- **`sister_clades.rds`**: List of sister clade pairs used for inferring insertion timing.  

- **`predicted_tp.rds`**  
  - Final estimates of insertion timing for mobile elements, inferred from the presence/absence patterns across clades.  

---

### Example output format

The main output (`timepoints.rds`) is a **GRanges object** with annotated mobile elements and their inferred insertion timepoints.  

```
GRanges object with 9958 ranges and 5 metadata columns:
                seqnames            ranges strand |      repeat       class      family                                           members timepoint
                   <Rle>         <IRanges>  <Rle> | <character> <character> <character>                                     <GRangesList>  <factor>
       L1PA4_71     chr1       64426-64666      * |       L1PA4        LINE          L1                                chr1:64426-64666:+         5
      L1PA4_473     chr1     400706-400825      * |       L1PA4        LINE          L1                              chr1:400706-400825:-         5
      L1PA4_532     chr1     437763-438662      * |       L1PA4        LINE          L1         chr1:437763-438484:-,chr1:438485-438662:+         4
      L1PA4_764     chr1     636459-636578      * |       L1PA4        LINE          L1                              chr1:636459-636578:-         5
      L1PA4_821     chr1     672742-673641      * |       L1PA4        LINE          L1         chr1:672742-673463:-,chr1:673464-673641:+         4
```

**Column definitions:**  
- **seqnames / ranges / strand**: Representative genomic coordinates of the repeat in the reference genome.  
- **repeat / class / family**: RepeatMasker annotation (name, class, and family).  
- **members**: A `GRangesList` of the underlying genomic fragments that compose the annotated repeat.  
  - Repeats may be split into multiple segments due to **nested insertions**, repeat specific rearrangments, interruptions, or alignment artifacts.  
  - For example, `L1PA4_532` and `L1PA4_821` are represented by two consecutive fragments joined into a single annotated repeat.  
- **timepoint**: Estimated evolutionary node when the insertion occurred (stored as a factor).  

This structure preserves both the **high-level repeat annotation** (main GRanges row) and the **fragment-level details** (`members`), that in case repeats are interrupted or fragmented. 

## Usage example 2 — Chimeric insertions on chr22

The following example illustrates how to run the *chimera prediction* workflow for LINE-1 elements on chromosome 22.  

```r
# Define a working directory
work_dir <- "/path/to/TiMEstamp_example"

# 1) Define sister clades from the updated tree
get_sister(
  tree_file = file.path(work_dir, "updated_tree.nh"),
  folder    = work_dir
)

# 2) Identify portions of LINE-1 missing from the alignment
get_missing_portion(
  chrom          = "chr22",
  reference_file = "reference/basic_filtered_l1.rds",
  selected_info  = "fivep_frag",
  folder         = work_dir
)

# 3) Inspect loci for gap structure and missing data
inspect_loci(
  folder = work_dir,
  chrom  = "chr22"
)

# 4) Extract flanking segments (e.g., upstream)
extract_flanking_segment(
  folder     = work_dir,
  chrom      = "chr22",
  which_side = "upstream"
)

# 5) Predict potential chimeric insertions
predict_chimera(
  folder = work_dir,
  chrom  = "chr22"
)
```

### Example output file
The chimera step produces tab-delimited files:

```
<work_dir>/predicted_chimera_upstream.txt
```

---

A snippet from `predicted_chimera_upstream.txt`:

```
coordinate	ID	fivep	threep	fivep_frag	threep_frag	repeat	class	family	members	flanking_gcoord	full_range	timepoint
chr22:10642781-10643486:-	2468081	chr22:10643486:-	chr22:10642781:-	chr22:10642781-10643486:-	chr22:10642781-10643486:-	L1M2	LINE	L1	chr22:10642781-10643486:-	chr22:10643487-10643626:-	chr22:10642781-10643626:-	5
chr22:11580142-11584309:-	2468839	chr22:11584309:-	chr22:11580142:-	chr22:11583099-11584309:-	chr22:11580142-11580336:-	L1MEc	LINE	L1	chr22:11580142-11580336:-,chr22:11581964-11582463:-,chr22:11582766-11582796:-,chr22:11583099-11584309:-	chr22:11584310-11584383:-	chr22:11580142-11584383:-	6
```

### Column definitions
- **coordinate**: Genomic coordinate of the TE insertion.  
- **ID**: Unique identifier assigned by RepeatMasker.  
- **fivep** / **threep**: Coordinates of the 5′ and 3′ junctions of the TE.  
- **fivep_frag** / **threep_frag**: Fragment coordinates covering the 5′ and 3′ ends of the TE.  
- **repeat**: RepeatMasker subfamily (e.g., *L1M2*, *L1MEc*).  
- **class**: Repeat class (e.g., *LINE*).  
- **family**: Repeat family (e.g., *L1*).  
- **members**: Genomic fragments that make up the annotated repeat (handles splits due to nested insertions).  
- **flanking_gcoord**: Coordinates of the extracted flanking segment.  
- **full_range**: Combined span of TE plus flanking region.  
- **timepoint**: Inferred evolutionary timepoint (factor indicating lineage of insertion).  

## Notes and tips
Because MSAs spanning evolutionarily distant species can be noisy, predicted chimeric insertions may include false positives. We strongly recommend **manual review of every candidate** to: (i) verify 5′/3′ breakpoints in the MSA; (ii) confirm orthology and local synteny across species; (iii) rule out assembly gaps, paralogy, or alignment artifacts; (iv) cross-check against independent alignments (e.g., multi-assembly MSAs such as the 447-genome set of Kuderna et al., 2023); and (v) inspect molecular signatures of the TE (e.g., target site duplications [TSDs], and poly(A) tails for LINE-1).

## Citing TiMEstamp

If you use **TiMEstamp** in your research, please cite:

Law CT, Burns KH. *Comparative Genomics Reveals LINE-1 Recombination with Diverse RNAs*. **bioRxiv** [Preprint]. 2025 Feb 3:2025.02.02.635956. doi: 10.1101/2025.02.02.635956. PMID: 39975348; PMCID: PMC11838501.

