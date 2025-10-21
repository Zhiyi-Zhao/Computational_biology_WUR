Computational biology course note based on Wageningen University&amp;Research, course code: SSB34306
# Computational biology

# 1A Molecular biology primer: from gene to protein
### Genome
In classical genetics, the genome of a diploid organism refers to a full set of chromosomes or genes in a gamete. A regular somatic cell contains two full sets of genomes(and a mitochondrial genome)\
In haploid organisms, including bacteria, archaea and viruses, a cell contains only a single set of the genome, usually in a single circular or contiguous linear DNA.\
Mitochondrial and Chloroplast genome are also present in single copy.\
Size of the genomes is variable and does not correlate with complexity of living beings. And small and simple does not mean harmless.

### Gene
A gene is a programmable unit that can give rise to a multitude of products: protein and RNA products through (alternative) splicing and trans-splicing

- Reading Frame(RF), Open Reading Frame(ORF) and Coding Sequence(CDS)
  - A **Reading Frame** is how a genetic sequence is read during a translation process.
  - An **Open Reading Frame** is the part of a reading frame that has the potential to code for a protein or peptide and consists of a continuous stretch of codons that do not contain a stop codon (usually UAA, UAG or UGA). Note that in this case, the ORF will be interrupted with intervening sequences or introns.
  - The **Coding Sequence**, is the portion of a gene's DNA or RNA, composed **only of exons**, that codes for protein. The region is bounded at the 5' end with a start codon and at the 3' end with a stopcodon.

# 2A Genome coverage
### Pacbio and Nanopore
**PacBio** is a "sequencing by synthesis" method: fluorescent signals are recorded as DNA is synthesized by polymerase. 

**Nanopore** is a "pore-permeation sequencing" method: bases are determined by measuring changes in current as DNA/RNA passes through a pore. Nanopore offers longer reads and larger data volumes (suitable for analyzing structural variations and complex regions).

### Coverage: 
average number of times any given base in the genome is sequenced.
<center>

$$Coverage(a) = \frac{N \times L}{G}$$

</center>

  - Number of reads -- N
  - Read Length -- L
  - Genome length -- G

  Given coverage **a**, then the probability that base is     sequenced n times or **a base not being sequenced(Gap)**:
<center>

$$P(n) = \frac{a^n e^{-a}}{n!}$$

</center>

By this fomula: P(0) = the probability of the gaps = $e^{-a}$, so for a genome of size G, the number of **nucleotides in the gaps** is: $P(0)G = e^{-a}G$

The number of **contigs**:
<center>

$$Contigs = Ne^{-a}$$

</center>
where N is the number of reads.


# 2B Genome Assembly & Annotation
### The differences between mapping and assembly
  Mapping = compare to a known genome\
  Assembly = build the genome from scratch
### The general strategy for genome assembly
  1. Find overlaps among reads
  2. Build a graph to visualize the connections
  3. Make the graph simpler
  4. Walk through the graph
### The specific approaches for genome assembly
- **Greedy** overlapping based assembly: Calculate pairwise alignment of all reads vs all reads; Find the best matching read pair; Merge pair; Repeat until you cannot continue. (Heuristic algorithm, the current optimal solution)
- **Overlap–Layout–Consensus (OLC)**: Calculate the overlap of all read pairs and construct a complete overlap graph; Analyze all possible paths through a graph search or layout phase; Ultimately, find the most plausible path (one with the highest consensus support) from a **global perspective**.
- **de Bruijn graph**: decompose all reads → count **k-mers** → build graph → find paths 
### Compare the structure of contigs and scaffolds
  Scaffold is made by the contigs placed in correct **orientation** and correct **order**, with the approximately correct **distance** sometimes including repeats.
### Define the scope of structural and functional annotations
- Structural annotation: is there a gene? Three common methods:
    1. use staristical models to search for specific sequences that indicate the presence of a gene nearby, or statistial properties of the protein-coding sequence itself.
    2. Evidence-based approach: sequence mRNA
    3. Homology-based approach: align genome with known genes or proteins to find orthologs
- Functional annotation: what's the function of the sequence?
    1. From protein to Gene:Experimental function of known proteins → mapped to corresponding genes
    2. Predict function from homology protein:If a gene has a high similarity to a known protein sequence, it is inferred that it may have the same or similar function.


# 3A Transcript analysis
### Types of RNA and main characteristics
- messenger RNA(main), code for proteins
- transfer RNA(main), central to protein synthesis as adaptors between mRNA and amino acids
- ribosomal RNA(main), form the basic structure of the ribosome and catalyze protein synthesis 
  microRNA, short interfering RNA, small nucleolar RNA, long non-coding RNA, piwi-interacting RNA etc.

### Various processing steps involved between transcription and protein expression
  1. Pre-mRNA is spliced to form mature mRNA, capped and polyadenylated
  2. In eukaryotes, one gene can result in different mRNA transcripts; in bacteria, one mRNA can translate into multiple proteins
  3. Translation

### Usage and limitation of measuring mRNA expression to understand biological process:
  1. Usage: mRNA and protein often shows same trends, so we can use it for conditions comparison and time series comparison
  2. Many reason can cause difference: different genes,isoforms,tissues,developmental stages,cell cycle,circadian rhythm, individual cells,environment, and the results of mRNA synthesis and mRNA decay

### The main technologies to measure transcript concentrations and list their specific advantages and disadvantages and main use:
  1. RT-qPCR: Low throughput
  2. Microarrays: pros: highly standardized, relatively cheap, samll data size--easy to handle; cons: gene sequence should be known, no position-specific information, can't detect new isoforms, not very quantitative: low dynamic range
  3. RNA-seq: untargeted, works for species without a sequenced genome, can identify alternatively spliced transcripts, can identify SNPAs in transcripts, high dynamic range(quantitative), no strand specificity.\
  Challenges: RNA reads will span an intron on the genome; one exon can be part of multiple isoforms\
  RNA-seq mapping: pseudo-alignment based approaches(when have ref_database): use only transcript information(refrence transcriptome); de novo assembly: when genome is not well annotated

# 3B Transcriptomics Quantification and differential expression(DE) analysis
- Normalizing for sequenceing depth:\
  **CPM**\
  **C**ount **P**er **M**ilion of mapped reads: counts scaled for total number of reads, relative measure for reads count where the total amount of reads is set to 1 milion to avoid samll numbers:
  <center>
 
  $$CPM =10^{6} \frac{ReadsPerGene}{TotalReads}$$

  </center> 

  **RPKM/FPKM**\
  **R**eads **P**er **K**ilo base of transcript per **M**illion mapped reads
    <center>

    $$RPKM =10^{6} \times \frac{ReadsPerGene}{TotalReads} \times 10^{3} \times \frac{1}{GeneLength} = 10^{9} \times\frac{ReadsPerGene}{GeneLength \times TotalReads} $$

    </center> 
  As to paired-end sequence, RPKM = 2 * FPKM

  **TPM**(Recommended, we can use this to compare both inside and between samples)\
  <center>

    $$TPM =10^{6} \times \frac{ReadsPerGene}{GeneLength} \times  \frac{1}{TotalLengthCorrectedCounts} = 10^{6} \times\frac{\frac{ReadsPerGene}{GeneLength}}{\sum\frac{ReadsPerGene}{GeneLength}} $$

    </center> 

- Differential expression analysis--which gene?
  T-test
   <center>
$$\frac{Signal}{Noise}=\frac{Difference between groups means}{Variability of froups}=\frac{\bar{x}_1-\bar{x}_2}{\sqrt{\frac{s_1^2}{n_1}+\frac{s_2^2}{n_2}}}$$
   </center> 

- Differential expression analysis--significant?
  Visualizing--volcano plot

- Differential expression analysis--meaningful?
  GO analysis

# 4A MS-based proteomics  
### Why do we study proteins?
1. It can tell us what happens now, what enzymes are currently active, which signals are being transduced.
2. Transcriptomics sometimes have poor correlation with proteomics.
3. Genomics can only tell what's the potential
   
### What is proteomics?
The study of all proteins present in a sample at a given time.

### What are the challenges in analysing proteins?
  1. Proteins are diverse in functions
  2. Proteins are diverse in localizations(nucleus, secreted, mitochondria, etc.)
  3. Can work alone but mostly in complexes
  4. Can have different turn-over(e.g. 1.5h, 1-2 days, years)
  5. Post-translational modification

### How do we do it?
  1. Mass spectrometry(MS)-based proteomics
  2. MALDI(Matrix Assisted Laser Desorption Ionization)
  3. ESI(Electron Spray Ionization)

### How can we manually give meaning to the spectra?
  MS spectra containing the m/z and intensity of the charged peptides and MS/MS spectra containing the information on the fragment of the charged peptides.

# 4B Identification and quantification of proteins in complex samples
- Understanding how mass spectrometry data are analysed to produce protein indentifications:\
  MS/MS data is aligned with MS/MS database. The database is built by the true protein fragments' theoretical data. To improve the quanlity data, we also need to use the **decoy database**\
  **FDR**： measures the error rate associated with a collection of PSMs:
$$
FDR(S_T) = \frac{N_d(S_T)}{N_t(S_T)}
$$
  
- Discover how proteins are quantified by mass spectrometry-based proteomics methods\
  1. MS1 quantification:\
  Label free quantitation\
  Metabolic labeling\
  Isotopic or isobaric tags
  2. MS2 quantification:\
  Lable free\
  isobaric tagging

- Reflect on how the quantity of the entity (proteins) is inferred by quantifying its components (peptides)

# 5A Databases
- International databases:
  DNA/RNA: GenBnak for genetic sequences, SRA for next NGS data;\
  Protein: Uniprot for sequence information, PDB for structures;
  Metabolism: KEGG, BioCyc

# 5B Public resourses for genetic data & the need for data FAIRification
Nucleotide sequence database: DDBJ, EMBL, Genbank

Uniprot: is made by two part: 1. Swiss-prot: the manually curated section; 2. TrEMBL: the automatically annotated section. this databse is focused on sequence data and functional annotaions

Key features of Swiss-Prot:
- High-Quality Annotations
- Non-Redundant Data
- Cross-References
- High Accuracy
- Regular Updates

PDB: focused on 3D structural data of proteins, nucleic acid and complexes

Metadata: data that provides information about other data.\
Data provenance: data centered metadata and captures the historical record--source origin, inputs, entities, factors, processes and locations.\
FAIR data: scientific data linked to machine readable metadata: Findable, Accessible, Interoperable, Reusable.\
Data workflows can visualize the data provenance (steps)

# 6A Substitution Patterns & Matrices
### Reasons for why sequence alignments are done
   
### Compare the Jukes-Cantor and Kimura models for nucleotide substitution
   **Jukes-Cantor**:
   $$d = \frac{-3}{4} \times ln(1-\frac{4}{3} \times D)$$
   every mutation pair have same score \
   **Kimura**: take substitutions and transversions into consideration. Transitions(A-G or T-C)are more likely than transversions(A-C, A-T, G-C, G-T)

### Calculate log-adds ratios for amino acid substitutions;
   $$log_2(\frac{Pr(x,y|R)}{Pr(x,y|U)}) = \sum log_2(\frac{q_{x,y}}{p_{x}p_{y}})$$

### Why and how amino acid substitution matrices have been constructed

| Feature                        | BLOSUM62 (Default) | Other BLOSUM Matrices | PAM Matrices |
|--------------------------------|--------------------|----------------------|--------------|
| **Basis of Calculation**       | Derived from conserved blocks of protein sequences | BLOSUMX is derived from alignments with at least X% sequence identity | Based on evolutionary models of accepted mutations |

# 6B BLAST(Basic Local Alignment Search Tool)
### List the different types of BLAST
   - blastn n vs n
   - blastp p vs p
   - blastx n(trans) vs p
   - tblastn p vs n(trans)
   - tblastx n(trans) vs n(trans)  
### Summarize the four steps of the BLAST algorithm;
1. **List**: 1. split sequence into **k-mers**; 2. calculate alignment score and keep only k-mers above certain threshold
2. **Scan and extend**: 1. scan database for exact matches to each word from list step; 2. extend each hit in both directions calculate raw score until score drops a certain amount below highest score; 3. allows gaps to be filled as long as gap penalyies do not drop total score below a certain score\
The raw score (*S*) of an alignment depends on the number of matches/mismatches, the number of gaps and the length of the gaps, and can be calculated with the following formula, where $M_{i,j}$ is the sum of the scores for each match or mismatch (from the corresponding nucleotide or amino acid scoring matrices), *c* is the number of gaps, *O* is the penalty for the existence of a gap, *d* is the total length of gaps and *G* is the per-residue penalty for extending the gap.

$$
S = \sum (M_{i,j}) - cO - dG
$$

3. **Report**: High-scoring segment pairs are reported. Raw scores (S) are converted to bit scores to allow comparisons between different searches, $S' = \frac{\lambda \times S - \ln(K)}{\ln(2)}$ , K and lambda are Karlin-Altschul parameters that depend on the scoring matrix.
  
4. Explain how changing BLAST parameters influence the algorithm and output;
   We can change the databse, targeted organism and protein matrix
5. Evaluate BLAST.
   $$
   E=K*m*n*e^{-\lambda*S}
   $$
   Where **m** is the effective length of query and **n** is the effective length of database 

# 7 Protein Domains
- Domains: are distinct structural units of proteins that fold independently and have a hydrophobic core, they're the most basic functional units of proteins that are required for their activity.\
**Motifs** are found as patterns(specific or degenerate), which can define a domain, but also can be characteristic of a subset of a domain and found in DNA

- Three general methods for domain identification:
  - Patterns: highly conserved sequence motifs that are domain-unique
  - Profiles: domain-specific position specific scoring matrices
  - Hidden Markov Models(HMMs): probability-based models

- Profiles, position-specific scoring matrices (PSSMs):\
  Table of position-specific scores and gap penalties based on a multiple sequence aligment. Specific matrix for each position in multiple sequence alignment. PSSM is a dynamic, domain-specific matrix

- PSI-BLAST(Position-Specific Iterative BLAST):\
  While standard BLAST uses a fixed substitution matrix (such as BLOSUM62), PSI-BLAST automatically creates a position-specific scoring matrix (**PSSM**) for your query sequence and continuously updates it over multiple iterations. This allows starting with closely related homologs, it gradually captures more distantly related proteins (those with lower sequence similarity but related functions).

# 8 Topological signals and sequences
- How proteins are sorted in the cell
  **Three ways to move around**
  - Gated transport: cytosol and nucleus using nuclear pore (No topological changes occurred， because the nuclear membrane is two)
  - Protein translocation: use of a protein translocator to directly transport specific proteins across membrane(topologically change happened)
  - vesicular transport(no topological change)
  **How is protein location encoded?**
  - Signal sequence
  - Signal patches
  **2 types of translocation processes**
  - Co-translational translocation : signal-recognition particle(SRP) directs the ER signal sequence to a specific receptor in the rough ER membrane
  - Post-translational translocation
    
- Illustrate the bioinformatics approaches developed to solve the problem
  - PSORT： signal sequence decision tree, follow the decision tree to a cellular compartment
  - Hidden Markov Model: an HMM models statistical regularities in a set of related sequences; tolerant to insertions and deletions. A profile is encoded in position-dependent state transitions
  - TMHMM: transmembrane prediction based on HMM models
 **Summary** For sinal peptides: signalIP/DeepSig; Transmembrane domains:TMHMM/DeepTMHMM; Signal sequences: PSORT/MULocDeep/DeepLoc

==What to study for the exam== 
- What are the three main problems in protein topology?
- Which tools solve which problem
- How does it work(use ANN or pML etc.)
- How do I interpret the output?
- How does the ML work in general? How do I assess a good ML tool?

# 9 Gene Orthology
- Understand the difference between gene orthologs and paralogs\
  - Orthologs: Sequences whose last **common ancestor** diverged at a **speciation event**, they perform the same function;
  - Paralogs: Sequences diverged at a **duplication event** in the **same species**, they have a different function. 

- Use different bioinformatic approaches to infer orthologs
 - **BBH**: Gene pairs with the highest sequence identity.
   Step 1: Forward BLAST Search：Use the "human gene of interest" as the query and perform a BLAST search against the "Mouse Genome Database." Several significant matches (black dots) are obtained, with the one with the highest score being the best hit (green dot). This is the mouse gene most similar to the human gene.\
   Step 2: Reverse BLAST: Now use the mouse gene (best hit) just found as the query and perform another BLAST search against the "Human Gene Database." If its best hit happens to be the original human gene, then the two are considered bidirectional best hits (BBHs).

   Problem: 1. orthology is nit a 1-to-1 relationship; 2. Gene evolution does not always follow species evolution. BBH works well if gene duplications are rare, and the true ortholog of each gene is not missing.

   Whole-genome duplications are common in plants and limit the utility of BBH(A complementary strategy to infer orthology: Synteny)
 - **Synteny**: shared order and orientation of a string of multiple genes and their closest counterparts in a related genome. Genes in same order tend to correspond.
 - **DELTA-BLAST**: Domain-Enhanced Lookup Time Accelerated BLAST. queries the conserved domain database. Structures persist longer than sequence identity.

# 10 Multiple Sequence Alignments & Phylogeny
- Why MSA are performed
  - Allows genome annotation(domains/orthologs)
  - Good for finding important conserved amino acid and nuvleotide positions
  - Can be used to build advanced statistical models or patterns(PSI-BLAST, PSSM, HMM)
  - Difference can be used to study evolutionary relationships

- Describe multiple criteria and strategies used to select sequences for MSA and phylogenetic trees
  - Sequences must be homologous
  - Sequences must align globally
  - Should have same domains in same order

- Explain the general procedure used to generate phylogenetic trees
  1. retrieving sequences: Swissprot/Uniprot--function; Reference protein database--species of interest;
  2.  Measure sequence similarity;
  3.  Tree building: WPGMA/UPGMA/Neighbor Joining;
  4.  Testing tree accuracy: Bootstrap

- Describe a statistical process used to determine the reliability of phylogenetic trees--Bootstrap
  1. Columns of MSAs are randomly sampled to generate different sequences of the same length and composition;
  2. For each replicate, a tree is built;
  3. Repeat steps 1-2 many times
  4. Bootstrap values = macthing replicates/total replicates

# 11 Basics of protein folding
- Bonds types:
  - Covalent bonds
  - Non-covalent interactions:
    1. Electrostatic interactions
    2. Hydrogen bonds
    3. Van der Waals forces
    4. Hydrophobic interactions

- Why does a protein fold? Covalent and non-covalent bonds reduce $\Delta G$

- Dihedral angles: the only reasonable movements are around the Ca-N($\phi$) bond and the Ca-C bond($\psi$). A Ramachandran plot shows which phi and psi angles are typically found

# 12 3D protein structures

- Primary structure: the sequence
- Secondary structure:
  1. $\alpha$-helices, Often negative amino acid near N-terminus, positive near C-terminus; 
  2. $\beta$-strands, bulky aromatic side-chains fit well in strands, and beta-branched residues are very suitable for beta-sheets;
  3. turns: connect the secondary structure elements;
  4. loops;
- Tertiary structure: Interactions between the secondary structure elements to form the structured protein subunit.
- Quaternary structure: Multiple protein subunits together form a quaternary protein structure

# 13 Structural comparison
| Method | Pro | Con |
|------|------|---------|
| X-ray crystallography | - High resolution (median 2.05Å, highest 0.48Å)| - Only static image|
| | - Easy to visualize and interpret| - Limitations in which proteins can form good crystals (membrane proteins, complexes, …)  |
| | - Can be used for large molecules | - Relatively slow |
| NMR | - Allows flexibility of proteins, shows conformation average | - Relatively expensive |
| | - Can be used to study interactions | - Requires high concentration of protein |
| | | - Cannot be used for large proteins |
| | | - Data can be difficult to interpret |
| Cryo-EM | - Small sample size required | - Typical resolution 3–4Å (but best recorded at 1.22Å) |
| | - Preservation of hydrated state | |
| | - Works on large range of sample types (incl. very large proteins) |  |

- Why search for similar structures?
  1. Find homologs with low sequence similarity
  2. Explore protein evolution: similar protein folds can support different functions
  3. Identify conserved core elements to model related proteins of unlnown structure

# 14 Protein Structure Modelling and Quality Check
- Why protein structure prediction?\
  We have an experimentally determined atomic structure for only a small fraction of the known protein sequences. Modelling is a quick and cheap way to explore protein structures and functions.
  
- Why structure comparison?
  - To understand structure-function relationship;
  - To study the evolution of many key proteins(structure is more conserved)

- AlphaFold2's Deep Learning algorithm
  1. Compare input sequence to databases; MSA and align with structural database;
  2. Neural network for refinement;
  3. Construct and evaluate models

