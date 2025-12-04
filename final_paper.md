## HMMER Workflow

The effective use of HMMER relies on understanding its workflow. Unlike simple pairwise alignment tools that function as a single step (sequence A vs. sequence B), HMMER operates as a multi-stage pipeline. This pipeline transforms raw biological data into a statistical model, which is then utilized to either discover new homologs in sequence databases or to annotate functional domains by querying profile databases.

This section lists out HMMER workflow that follows a linear progression consisting of **four distinct stages**:

### 1. Input

To get started, we need to input a high-quality **Multiple Sequence Alignment (MSA)**. HMMER cannot build a profile from a single sequence; it requires the evolutionary context provided by a family of aligned sequences to determine which residues are conserved and which are variable.
        

### 2. Profile Construction

This stage uses `hmmbuild` command to convert the biological alignment into a mathematical Profile HMM. It implements the theoretical concepts of Match, Insert, and Delete states to model the evolutionary constraints of the family.

### 3. Application: The Search Strategies

Once the profile is built, the workflow can go into two directions based on the research objective:

- **HomologySearch (`hmmsearch`)**: Used for discovery. It takes the custom profile HMM and scans a target sequence database to find new homologs.
        
- **Domain Annotation (`hmmscan`):** Used for identification. It takes a single query sequence and scans it against a library of profiles HMM to identify functional domains within the sequence.

### 4. Output: Ranked Results

Both strategies generate a list of matches. HMMER filters biological signals from random noise, presenting hits ranked by bit score and e-value.
           
## HMMER Functions

While the workflow provides the roadmap, understanding the specific mechanics of each tool is essential for effective analysis. This section details the primary functions used in HMMER analysis.

### `hmmbuild`: Profile Construction

`hmmbuild` serves as the "architect" of the pipeline, responsible for transforming the evolutionary information contained in a multiple sequence alignment into a binary-compatible statistical model.

**Input** The input for this command is a Multiple Sequence Alignment (MSA). The software is designed to accept several standard bioinformatics alignment formats, including Stockholm (`.sto`), Clustal (`.aln`), or aligned FASTA.

**Process** The tool reads the alignment column-by-column and applies the Profile HMM statistical architecture (calculating Match, Insert, and Delete state probabilities) as explained in the previous section. It converts the observed biological data into a probabilistic model without requiring manual parameterization.

**Output** The command generates a Profile HMM file (ending with `.hmm` extension). This serves as a binary-compatible file that contains all information needed for later search strategies.

![hmmbuild_output](hmmbuild_output.png)

**Usage** The command requires the user to specify an output filename followed by the input alignment. For example, if you want to create a new profile HMM named `globins4.hmm` from the source alignment `globins4.sto`, you can use the command
```
# Syntax: hmmbuild [output_hmm_file] [input_msa_file]
hmmbuild globins4.hmm globins4.sto
```

