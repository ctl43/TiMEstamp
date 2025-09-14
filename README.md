# TiMEstamp: Timing Mobile ELement Insertion from Multiple Sequence Alignments
TiMEstamp is an R package for inferring insertion timepoints of mobile elements from multiple sequence alignments.

## INSTALLATION:
1. Clone from GitHub:
git clone https://github.com/ctl43/TiMEstamp.git

2. Install Using R CMD:
R CMD INSTALL TiMEstamp

## PREPARATION OF MULTIPLE SEQUENCE ALIGNMENT (MSA):
1. Download available MSA data in MAF format (e.g., https://hgdownload.soe.ucsc.edu/goldenPath/hg38/multiz470way/) or generate the multiple genome alignment by yourself. The MAF files must be stored separately by chromosome in the reference species.
2. Split the genome alignment by species using PHAST.  
```mafSpeciesSubset -keepFirst```  
3. Remove gaps relative to the reference genome, which should always be the first genome, using PHAST. The result files will be in FASTA format.  
```msa_view --gap-strip 1```  

P.S.: These prepreation steps will be incorporated into TiMEstamp later without relying PHAST.  

## EXTRACT GAPS RELATIVE TO THE REFERENCE GENOME:
Unaligned regions, indicated by asterisks (*) or hyphens (-) in the resulting FASTA files, were extracted and annotated in BED file format.
```r
# Finding all unaligned regions (*/-) and store in RDS format, an R object.
get_unaligned_regions(fa, maf, out_folder)
combined_unaligned_by_chrom(fa, input_folder, output_folder, threads = 4)
```

## USAGE EXAMPLE:

Below is a simple example of running **TiMEstamp** with a human genome–based multiple sequence alignment (MSA).  
We assume the MSA and repeatMasker annotation files are already prepared as described above.  

```r
# Define working directory for this project
work_dir <- "TiMEstamp_example/hg38_multiz470way"

# Step 1. Prepare RepeatMasker annotation
prepare_rmsk(
  rmsk_file = "reference/rmsk_hg38.fa.out",
  folder    = work_dir
)

# Step 2. (Optional) Update phylogenetic tree to include only selected species
update_tree(
  tree    = "phylogenetic/hg38.470way.nh",
  folder  = work_dir,
  species = selected_species
)

# Step 3. Define sister clades based on the updated tree
get_sister(
  tree_file = file.path(work_dir, "updated_tree.nh"),
  folder    = work_dir
)

# Step 4. Identify missing loci coverage by clade
get_missing_coverage_by_clade_fast(folder = work_dir)

# Step 5. Clean clade data and filter loci by criteria
clean_clade_data_fast(folder = work_dir)

# Step 6. Estimate insertion timepoints of mobile elements
get_timepoint_fast(folder = work_dir)
```
---

## EXPLANATION OF OUTPUT

After running the pipeline, **TiMEstamp** will generate the following key outputs inside the working directory (`work_dir`):  

- **`reference_anno/`**  
  Processed RepeatMasker annotation files for the reference genome.  

- **`phylogenetic/`**  
  - `updated_tree.nh`: Phylogenetic tree pruned to the selected species.  
  - `sister_clades.rds`: List of sister clade pairs used for inferring insertion timing.  

- **`gap/`**  
  BED or RDS files containing unaligned regions (gaps) for each chromosome relative to the reference genome.  

- **`fast/`**  
  Intermediate results including reference annotation, clade coverage summaries, and pre-filtered gap data.  

- **`results/`**  
  - `cleaned_gap_data.rds`: Processed gap annotations after filtering by size and coverage.  
  - `timepoints.rds`: Final estimates of insertion timing for mobile elements, inferred from the presence/absence patterns across clades.  

---

## EXAMPLE OUTPUT FORMAT

The main output (`timepoints.rds`) is a **GRanges object** with annotated mobile elements and their inferred insertion timepoints.  
A simplified excerpt looks like this:  

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

This structure preserves both the **high-level repeat annotation** (main GRanges row) and the **fragment-level details** (`members`), ensuring accuracy when repeats are interrupted or fragmented. 

## Example: Detecting chimeric insertions on chromosome 22

The following example illustrates how to run the *chimera prediction* workflow for LINE-1 elements on chromosome 22.  
Input includes a filtered RepeatMasker reference file (`basic_filtered_l1.rds`) and a TiMEstamp working directory (`TiMEstamp`).  

```r
# Define working directory for this project
work_dir <- "TiMEstamp_example/hg38_multiz470way"

# Step 1. Identify portions of LINE-1 that are missing from the alignment
get_missing_portion(
  "chr22",
  reference_file = "reference/basic_filtered_l1.rds",
  selected_info  = "fivep_frag",
  folder         = work_dir
)

# Step 2. Inspect loci for gap structure and missing data
inspect_loci(
  folder = "TiMEstamp",
  chrom  = "chr22"
)

# Step 3. Extract flanking segments (here, upstream of the TE)
extract_flanking_segment(
  folder     = work_dir,
  chrom      = "chr22",
  which_side = "upstream"
)

# Step 4. Predict potential chimeric insertions
predict_chimera(
  folder = work_dir,
  chrom  = "chr22"
)
```

### Example output file
The chimera prediction step produces tab-delimited files named like:

```
TiMEstamp_example/hg38_multiz470way/predicted_chimera_upstream.txt
TiMEstamp_example/hg38_multiz470way/predicted_chimera_downstream.txt
```

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
