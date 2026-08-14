# Overview

This repository contains the scripts and data used to analyze the manuscript "Rapid Speciation Characterized by Incomplete Lineage Sorting in the Globally Distributed Bacterium Sulfitobacter".

This repository focuses on **phylogenetic reconstruction, restriction-modification system identification, quartet analyses, genomic identity scanning, neutral simulations, and genomovar analysis** for Sulfitobacter genomes.

## Installation

### Python packages
```bash
pip install pandas biopython tqdm ete3 plotly numpy scipy matplotlib
```

### Required software
- **OrthoFinder** (v2.2.1): For orthologous group inference
  - Installation: [OrthoFinder GitHub](https://github.com/davidemms/OrthoFinder)
  
- **DIAMOND**: For fast protein sequence searches
  - Installation: [DIAMOND GitHub](https://github.com/bbuchfink/diamond)
  
- **MAFFT**: For multiple sequence alignment
  - Installation: [MAFFT website](https://mafft.cbrc.jp/alignment/software/)
  
- **IQ-TREE**: For phylogenetic tree construction
  - Installation: [IQ-TREE website](http://www.iqtree.org/)
  
- **trimAl**: For alignment trimming
  - Installation: [trimAl GitHub](https://github.com/scapella/trimal)
  
- **QuartetScores**: For quartet topology analysis
  - Installation: [QuartetScores GitHub](https://github.com/lutteropp/QuartetScores)
  - **Note**: The original QuartetScores needs modification to output topology counts (see `03-Quartet-analyses/README.md`)

- **GSL (GNU Scientific Library)**: Required for neutral simulation program
  - Installation: [GSL website](https://www.gnu.org/software/gsl/)
  - On macOS: `brew install gsl`
  - On Linux: `sudo apt-get install libgsl-dev` (Debian/Ubuntu) or `sudo yum install gsl-devel` (RHEL/CentOS)

### REBASE database
The restriction-modification (R-M) system reference data should be downloaded from [REBASE](http://rebase.neb.com/rebase/rebase.ftp.html), including:
- `protein_seqs.fasta`: Protein sequences of R-M enzymes
- `enzymes.txt`: Enzyme annotation file

## Main Scripts

### **01-phylogeny-construction/** - Phylogenetic Tree Construction

This directory contains scripts for building phylogenetic trees from orthologous gene families.

- **S1.orthoFinder_blast.sh**  
  Runs **OrthoFinder** to identify orthologous groups (OGs) from protein sequences using DIAMOND for sequence similarity searches. The script processes `.pep` files (amino acid sequences) from the input directory and generates orthogroup assignments.

- **S2.extract_scp_seq.pl**  
  Extracts single-copy ortholog (SCO) sequences from OrthoFinder results. The script:
  - Reads OrthoFinder output files (`SingleCopyOrthogroups.txt` and `Orthogroups.txt`)
  - Extracts gene sequences for each SCO family in both amino acid (`.faa`) and nucleotide (`.dna`) formats
  - Organizes output into `aa_seq/` and `nuc_seq/` directories

- **S3.cat_aln.pl**  
  Concatenates single-copy ortholog alignments for concatenated phylogenetic analysis.

- **S4.sbatch_iqtree.sh**  
  Submits IQ-TREE jobs for phylogenetic tree construction from concatenated alignments.

- **nuc_seq/** subdirectory scripts:
  - `step1.mafft_aa_aln.pl`: Aligns amino acid sequences using MAFFT
  - `step2.multi2oneline.sh`: Converts multi-line FASTA to single-line format
  - `step3.impose.DNA.on.pep_alignment.pl`: Imposes nucleotide sequences onto protein alignments
  - `step4.trimAl.sh`: Trims alignments using trimAl
  - `step5.batch_iqtree.pl`: Batch runs IQ-TREE for individual gene trees
  - `step6.leaves_reroot.sh`: Reroots trees using different methods (MAD_root, MV_root)
  - `step7.label_generate.pl`: Generates annotation labels for tree visualization

**Input**: `.pep` (amino acid) and `.gene` (nucleotide) sequence files for each genome  
**Output**: Phylogenetic trees for concatenated alignments and individual ortholog families

---

### **02-RMS-identification/** - Restriction-Modification System Identification

This directory contains scripts for identifying restriction-modification (R-M) systems in genomes.

- **s1_diamond_blastp.sh**  
  Performs DIAMOND BLASTP searches of query genomes against the REBASE R-M enzyme database. The script:
  - Builds a DIAMOND database from REBASE protein sequences
  - Searches each genome's protein sequences (`.faa` files) against the R-M database
  - Uses thresholds: E-value ≤ 1e-5, identity ≥ 30%, query coverage ≥ 50%

- **s2_best_hit_RM.pl**  
  Identifies the best hit R-M enzyme for each query gene. The script:
  - Processes DIAMOND BLASTP output files
  - Selects the best hit (highest scoring) for each locus
  - Annotates hits with R-M enzyme details from REBASE `enzymes.txt`
  - Outputs TSV files (`*_RM.tsv`) with locus, enzyme name, and annotation

- **s3_seek_RM_pairs.pl**  
  Identifies R-M system pairs (restriction enzymes and their cognate modification enzymes) by checking gene neighborhood relationships.

**Note**: After obtaining best hits, R-M systems should be manually verified to ensure that restriction and modification genes are neighbors in the genome.

**Input**: Genome protein sequences (`.faa`), REBASE database  
**Output**: Best hit R-M enzymes for each genome (`*_RM.tsv` files)

---

### **03-Quartet-analyses/** - Quartet Topology Analysis

This directory contains scripts for analyzing quartet topologies to assess phylogenetic signal and gene tree discordance.

- **s1.QS_uniqgnm.sh**  
  Runs **QuartetScores** to count the number of supported topologies for each quartet of genomes. This script processes gene trees and counts how many trees support each of the three possible quartet topologies.

- **s2.process_quartets.py**  
  Processes QuartetScores output and categorizes quartets based on clade membership:
  - **Same clade**: All four taxa belong to the same clade
  - **Two clade (2-2)**: Two taxa from clade A, two from clade B
  - **Two clade (3-1)**: Three taxa from one clade, one from another
  - **Three clade**: Taxa from three different clades
  - **Two-three clade**: Mixed patterns
  - Standardizes topology representations and generates summary statistics

- **s3.vis_Quartet_prop.R**  
  Visualizes quartet topology proportions using R, generating plots to compare topology distributions across different clade combinations.

**Note**: The original QuartetScores software needs modification to output topology counts. A modified `QuartetCountConverter.hpp` is provided in this directory. Users need to replace the original file in their QuartetScores installation (lines 117-120) and recompile.

**Input**: Gene trees, species tree (`45_gnm_tree.nwk`), clade assignments (`gnm_clade_id.tsv`)  
**Output**: Quartet topology counts and categorized results (`*_results.tsv` files)

---

### **04-identity-scanning/** - Genomic Identity Analysis

This directory contains scripts for calculating nucleotide identity between genomes and clades.

- **s1.cal_identity.py**  
  Calculates sliding-window nucleotide identity between clade members. The script:
  - Loads core genome alignment and population/clade assignments
  - Uses sliding windows (default: 1000 bp windows, 1000 bp steps)
  - Calculates pairwise identity between genomes from different clades
  - Outputs identity values for each window position

- **s2.vis_identity.py**  
  Visualizes identity scanning results, generating plots showing identity patterns across the genome alignment.

- **s3.check_ANI.py**  
  Calculates Average Nucleotide Identity (ANI) between clades. ANI provides a genome-wide measure of sequence similarity and is commonly used for species delineation.

**Input**: Core genome alignment (concatenated), genome-to-clade mapping file  
**Output**: Identity values per window, ANI values between clades (`Sulfito_ANI.txt`)

---

### **05-neutral-simulation/** - Neutral Evolution Simulation

This directory contains a C++ program for simulating neutral evolution and haplotype dynamics in bacterial populations, with a focus on restriction-modification (R-M) system incompatibility loci.

- **neutral_simulation.cpp**  
  A population genetics simulation program that models the evolution of haplotypes across multiple populations. The program simulates:
  - **Population structure**: Multiple populations with configurable sizes and migration rates (stepping-stone model)
  - **Incompatibility loci**: Binary loci representing R-M system haplotypes (default: 10 loci)
  - **Neutral loci**: Additional neutral sites for comparison
  - **Life cycle stages**:
    - **Migration**: Gene flow between adjacent populations
    - **Reproduction**: Fitness-based selection and recombination
    - **Recombination**: Horizontal gene transfer with configurable proportion and recombinant size
    - **Mutation**: Point mutations at incompatibility loci with Poisson-distributed mutation rates
  - **Haplotype tracking**: Monitors haplotype frequencies, fixation events, and loss across generations
  - **Output files**:
    - `haplotype.tsv`: Haplotype frequencies per generation
    - `fixed_haplotype.tsv`: Records of haplotypes that reached fixation
    - `run_<gen>.tsv`: Final generation individual genotypes and fitness
    - `prerun_haplotype_<gen>.tsv`: Haplotype metadata for restarting simulations
    - `pop_size_<gen>.tsv` and `mut_rate_<gen>.tsv`: Population parameters

  **Key parameters** (command-line options):
  - `-N`: Population size
  - `-mut`: Mutation rate per locus
  - `-mig`: Migration rate
  - `-rec_p`: Proportion of recombinants
  - `-rec_n`: Mean number of recombined loci
  - `-s`: Effect of incompatibility between parents
  - `-q`: Effect of mutations on survival rate
  - `-nl_inc`: Number of incompatibility loci
  - `-nl_neu`: Number of neutral loci
  - `-final_gen`: Number of generations to simulate
  - `-thres`: Frequency threshold for fixation (default: 0.95)
  - `-eN`: Population size change events
  - `-eU`: Mutation rate change events
  - `-Fp`, `-Fu`, `-Fh`: Input files for population sizes, mutation rates, and haplotypes

  **Compilation requirements**:
  - Requires GSL (GNU Scientific Library) for random number generation and statistical functions
  - Compile with: `g++ -o neutral_simulation neutral_simulation.cpp -lgsl -lgslcblas`

  **Usage example**:
  ```bash
  ./neutral_simulation -N 1000 -mut 1e-6 -mig 0.0001 -rec_p 0.2 -rec_n 1.0 -s 0.1 -nl_inc 10 -final_gen 10000 -o output_folder/
  ```

---

### **06-genomovar-analysis/** - Genomovar Analysis

This directory tests whether each major *Sulfitobacter* clade subdivides into finer-scale ANI-defined units ("genomovars") and whether such units, rather than the clades themselves, are the population structure at which ancestral sequence variants are fixed. 

- **genomovar_analysis.py**  
  Main analysis: within-clade ANI distribution, genomovar delineation by average-linkage clustering across a sweep of ANI thresholds (99.2–99.9%) and linkage methods, nucleotide diversity, the genomovar partition test, a size-preserving permutation control, and clade-level sorting completeness.

- **make_figure.py**  
  Builds the summary figure (within-clade ANI distribution, genomovar counts vs. threshold, and polymorphism partitioning).

**Input**: fastANI pairwise identities (`Sulfito_ANI.txt`), clade assignments (`non_redundant_gnm.txt`, `metadata.tsv`), and a joint four-clade core-LCB alignment (produced by progressiveMauve using non-redundant genomes)  
**Output**: genomovar assignments and partition-test tables, and `Fig_genomovar.{pdf,png}`

---

## Folders

- **01-phylogeny-construction/**  
  Contains scripts and intermediate files for phylogenetic tree construction, including OrthoFinder results, single-copy ortholog sequences, alignments, and gene trees. Subdirectories `with_outgroup/` and `without_outgroup/` contain trees rooted with different methods.

- **02-RMS-identification/**  
  Contains scripts for R-M system identification, DIAMOND BLASTP results, and annotated R-M enzyme tables. The reference R-M data (`KE12_RM_new.tsv`) and genome order list (`gnm_order.txt`) are included.

- **03-Quartet-analyses/**  
  Contains scripts for quartet analysis, modified QuartetScores source code (`QuartetCountConverter.hpp`), input files (species tree, clade assignments), and results tables categorizing quartet topologies by clade membership patterns.

- **04-identity-scanning/**  
  Contains scripts for identity calculation and visualization, along with ANI calculation results (`Sulfito_ANI.txt`).

- **05-neutral-simulation/**  
  Contains a C++ simulation program (`neutral_simulation.cpp`) for modeling neutral evolution and haplotype dynamics in bacterial populations. The program simulates population genetics processes including migration, recombination, mutation, and selection, with a focus on restriction-modification system incompatibility loci. Output includes haplotype frequency trajectories, fixation records, and population parameter files.

- **06-genomovar-analysis/**  
  Contains scripts testing whether *Sulfitobacter* clades subdivide into ANI-defined genomovars and whether such genomovars, rather than the clades, are the units at which ancestral sequence variants are fixed. Includes genomovar delineation, a partition test with permutation control and the summary figure. 

---

## Notes

- Most scripts contain hardcoded paths referring to the original server environment. Users need to modify paths (e.g., `/your/path/of/OrthoFinder/`, `/your/path/of/diamond/`) to match their local installation.

- For **QuartetScores** modification: Replace lines 117-120 in `/your_installation_path/QuartetScores/src/QuartetCountConverter.hpp` with the provided modified version, then recompile QuartetScores.

- **REBASE database** files need to be downloaded separately from the REBASE website and paths updated in the scripts.

- Input file formats:
  - `.pep`: Amino acid sequences in FASTA format
  - `.gene`: Nucleotide sequences in FASTA format
  - `.faa`: Protein sequences (alternative naming)
  - `.dna`: DNA sequences

---

## Publication

Manuscript: "Rapid Speciation Characterized by Incomplete Lineage Sorting in the Globally Distributed Bacterium Sulfitobacter"


---

## Data Availability

Genome sequences and associated data are available under the following NCBI BioProject IDs: PRJNA628848.


---

## Contact Us

If you have any questions or suggestions regarding these scripts, please contact:
- yoyostudents5775@gmail.com
