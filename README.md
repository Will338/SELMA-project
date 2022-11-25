# Methodology and code for sequence sampling and dataset refinement for phylogenetic analysis of SELMA and TIC proteins

## Curating a custom database for protein similarity searches
A custom local database was curated to include protein sequences belonging to key taxa from various sources. All protein sequences for the taxa listed in Supplementary Information were downloaded from Uniprot. Predicted protein sequences from transcriptomes reassembled(Johnson, Alexander and Brown, 2019) from the data generated in The MMETSP(Keeling et al., 2014) for the taxa listed in Supplementary Information were downloaded from an online repository(Johnson, Alexander and Brown, 2017). Protein datasets listed in Supplementary Information were downloaded from VEuPathDB(Amos et al., 2022). Predicted protein sequences from the transcriptome dataset for the Ross Sea Dinoflagellate, referred to as “RSDallcleannoPhaeo”, was obtained from a previous study(Hehenberger, Gast and Keeling, 2019). Predicted protein sequences from transcriptomes for apicomplexa, chrompodellid and squirmid taxa were obtained from previous studies(Janouškovec et al., 2019; Mathur et al., 2019).

These protein sequences were combined as a single fasta file from which a BLAST-searchable database was created.

For example:
```bash
module load blast+/2.11.0
makeblastdb -in customdatabase.fasta -dbtype prot
```

## An empty starting directory was created for the analysis of each protein

For each protein analysed (Cdc48, Uba1, Der1, Ubc4, Ubi, Hsp70, Tic110, Tic20, Tic22) a new directory was created and was given the same name as the protein.

For example, for the protein Cdc48, a new directory with the name "Cdc48" was made:
```bash
mkdir Cdc48
cd Cdc48
```

The following commands were then used to create a sub-directory structure.
```bash
#!/bin/bash

mkdir Blast_output
mkdir Querys
mkdir Tree_building
mkdir Tree_building/Separate_and_cluster_seqs
mkdir Tree_building/Separate_and_cluster_seqs/Manual_tree_refining

echo "" > Querys/queries.faa

```
## Running BLASTp search
Protein sequences in fasta format were saved in the file that was created in the previous step at the path /Query/queries.faa
These proteins are bait proteins for which there is strong evidence that these are homologues of the protein of interest. These were used as the queries for BLASTp searches against the custom protein database.

The following section of code was ran to perform the blast searches.

$1 was provided as the global path to directory that was given the same name as the protein being analysed in the previous section (and should still be the present working directory).

$2 was provided as the global path to the custom BLAST-formatted custom database 

```bash
#!/bin/bash

module load blast+/2.11.0

protein_dir=$1
custom_database=$2

blastp -query $protein_dir/Querys/queries.faa \
    -db $custom_database \
    -out $protein_dir/Blast_output/blastp_out \
    -outfmt '6 std salltitles' \
    -num_threads 25 
```
From the BLASTp results, the fasta header of each significant hit were extracted keeping just one copy of each header.
```bash
cat "$protein_dir"/Blast_output/blastp_out | sort -g -k 11 | \
sort -k2,2 -u | sort -g -k 11 | \
cut -f 2 > "$protein_dir"/Blast_output/hit_headers.txt
```

The headers were then used as identifiers to extract the sequences with significant hits to the query sequences from the custom database, using the program seqtk.

```bash
module load seqtk

seqtk subseq $custom_database "$protein_dir"/Blast_output/hit_headers.txt > \
$protein_dir/Tree_building/combined_seqs.faa

cp $protein_dir/Tree_building/combined_seqs.faa $protein_dir/seqs_from_blast_search_for_trees.faa
```

## Refining the BLASTp-generated dataset

### Removing poorly-aligned and very short sequences

The following sections of code were used as an initial coarse filter of the sequences obtained from BLASTp searches.

### 1. First aligned and trimmed the dataset 

The following script was ran to first align the sequences and select conserved sites.
$1 is given as the parent directory with the same name as the protein being investigated.

```bash
#!/bin/bash/

module load mafft/7.475
module load trimal/1.4 

Tree_dir=$1

mafft --thread 20 $Tree_dir/combined_seqs.faa > $Tree_dir/combined_seqs.mafftaln

trimal -in $Tree_dir/combined_seqs.mafftaln -out $Tree_dir/combined_seqs.automated1 -automated1
```
### 2. Remove sequences that are shorter than the cut-off

The following section of code was used to remove sequences that had retained a number of sites below a specified percentage cut-off of the number of sites in the total dataset, after the dataset was aligned and trimmed. The percentage cut-off value was provided as $1. For all proteins value "30" for 30 percent cut-off was used. This meant that any sequence that was made-up of 70% or more of gaps (empty, non-amino acid sites) were removed. However, sequences were kept for key dinoflagellate taxa *Kryptoperidinium* (called *Glenodinium* in some datasets), *Karlodinium* and *Karenia*, regardless of their length, so these were temporarily removed prior to this filtering step and added back into the dataset at a later step.

* `$1` was provided as the cut-off percentage
* `$2` was provided as the path to directory named the same as the protein being analysed.

```bash
#!/bin/bash
percent_cut=$1
Tree_dir=$2

#Take out dinos I'm interested in since I want to keep even short sequences of these
cat $Tree_dir/combined_seqs.automated1 | \
sed -E 's/(>.*)/£\1$/' | \
tr -d '\n' | \
tr '£' '\n' | \
egrep -iv 'glenodinium|ryptoperidinium|karenia|arlodinium' | \
tr '$' '\n' | \
sed '/^$/d' > $Tree_dir/combined_seqs.trimalautomated1_key_dinos_removed

#First get the length of the trimmed alignment
aln_len=$(echo | cat $Tree_dir/combined_seqs.trimalautomated1_key_dinos_removed | sed -E 's/(>.*)/£\1$/' | \
tr -d '\n' | tr '£' '\n' |  tr '$' '\n' | sed '/^$/d' | awk {'print length'} | head -n 2 | tail -n 1)

#Calculate X percent length of alignment
#inside awk you don't have direct access to shell variables, you need to pass them as an option. Hence -v AL=
len_percent_cut=$(echo | awk -v PC="$percent_cut" -v AL="$aln_len" '{ print PC*AL/100 }')

#echo $len_percent_cut

#Get headers of all sequences with length greater than percentage cut-off
cat $Tree_dir/combined_seqs.automated1 | \
sed -E 's/(>.*)/£\1$/' | \
tr -d '\n' | \
tr '£' '\n' | \
tr '$' '\t' | \
#Removes '-' only in second field, essentially removing all non-aminoacid sites from each sequence
awk '{gsub(/-/,"",$2)} 1' | \
sed '/^$/d' | \
#Looks at length of second field and keeps those longer than cut-off
awk -v LPC="$len_percent_cut" 'length($2) >= LPC { print }' | \
tr ' ' '\n' | \
grep '>' | tr -d '>'  > $Tree_dir/Separate_and_cluster_seqs/headers_short_seqs_removed.txt

#extract seqs that are longer than percentage cut-off
module load seqtk/1.3

seqtk subseq $Tree_dir/combined_seqs.faa $Tree_dir/Separate_and_cluster_seqs/headers_short_seqs_removed.txt \
> $Tree_dir/Separate_and_cluster_seqs/short_seqs_removed.faa
```

The following code recombines all of the unfiltered sequences from *Kryptoperidinium*, *Karlodinium* and *Karenia* with the dataset that includes all other taxa, for which short seqeunces have been removed.

```bash
cat $Tree_dir/combined_seqs.faa | \
sed -E 's/(>.*)/£\1$/' | \
tr -d '\n' | \
tr '£' '\n' | \
egrep -i 'glenodinium|ryptoperidinium|karenia|arlodinium' | \
tr '$' '\n' | \
sed '/^$/d' >> $Tree_dir/Separate_and_cluster_seqs/short_seqs_removed.faa
```

## Clustering of key groups separately to remove similar sequences (reducing redundancy of the dataset)

The following section of code was used to first divide indiviudal sequences into different taxonomic groups that could then by clustered to different percentage identities using the program CD-hit.

For example, the BLASTp searches often found a greater number of homologues from Stramenopiles than other groups. This could be because Stramenopiles are represented by a greater number of taxa in the custom database that was generated, are better sampled, or are naturally more diverse. This meant that to decrease the redundancy of sequences, Stramenopiles typically had to be clustered to a lower sequence identity than other groups.

In most cases this section of code was ran iteratively, adjusting the percentage identity cut-off values and then assessing the outcome in each iteration. The aim of doing this was to obtain a dataset that had approximately similar numbers of sequences for each of the major eukaryote groups that have complex red plastids, ensuring individual particular groups were not overrepresented in the final phylogenies. Likewise, having excessively numerous similar sequences from the same eukaryote clade would not have changed major outcomes of the analyses, or the relationships between the major caldes inferred in phylogenies. However, inclusion of too many of such similar sequences would make the analyses unecessarily more computationaly demanding and more time consuming. Therefore, an additional aim of clustering was to obtain a dataset with a computationally manageable number of sequences (typically several hundred sequences).

```bash
#!/bin/bash

Tree_dir=$1

module load mafft/7.475
module load trimal/1.4 
module load cdhit/4.8.1
module load seqtk/1.3
```
The following section specified the percentage identity value that was used by CD-hit for clustering. The values here are given as an example of what was used for the Cdc48 dataset, but in practice these were adjusted during each iteration of this section of code, to refine the dataset as described above.

```bash
stramenopile_cluster_value="0.80"
kareniaceae_cluster_value="1"
mainseqs_cluster_value="0.92"
diatom_cluster_value="1"
haptophyte_cluster_value="0.99"
apicomplexa_cluster_value="1"
cryptophyte_cluster_value="1"
```

## The code in the following section divided sequences into different groups based on searching for key words/taxa names in the fasta headers of the sequences, then splits the dataset into different groups of sequences so that certain sequence groups can be clustered at different identities

* Separate the *Kareniaceae*

```bash
cat $Tree_dir/Separate_and_cluster_seqs/short_seqs_removed.faa | \
sed -E 's/(>.*)/£\1$/' | \
tr -d '\n' | \
tr '£' '\n' | \
egrep -i 'karenia' | \
tr '$' '\n' | \
sed '/^$/d' > $Tree_dir/Separate_and_cluster_seqs/kareniaceae_seqs.faa
```

* Separate the *Peridiniales*

```bash
cat $Tree_dir/Separate_and_cluster_seqs/short_seqs_removed.faa | \
sed -E 's/(>.*)/£\1$/' | \
tr -d '\n' | \
tr '£' '\n' | \
egrep -i 'glenodinium|operidinium|durinskia' | \
tr '$' '\n' | \
sed '/^$/d' > $Tree_dir/Separate_and_cluster_seqs/peridiniales_seqs.faa
```

* Separate non-diatom stramenopiles since these contain a large number of sequences relative to other groups, so will need clustering at a lower percentage identity

```bash
cat $Tree_dir/Separate_and_cluster_seqs/short_seqs_removed.faa | \
sed -E 's/(>.*)/£\1$/' | \
tr -d '\n' | \
tr '£' '\n' | \
egrep -i 'stramenopile|Aureococcus|Ochromonas|Dinobryon|Pseudopedinella|Aureococcus|Heterosigma|Chattonella|Vaucheria|Ectocarpus|Nannochloropsis|Phytophthora|Albugo|Saprolegnia|Aplanochytrium|Aurantiochytrium|Schizochytrium' | \
egrep -v -i 'Thalassiosira|Phaeodactylum|Fragilariopsis|Amphora|Skeletonema|nitzschia' | \
tr '$' '\n' | \
sed '/^$/d' > $Tree_dir/Separate_and_cluster_seqs/stramenopile_seqs.faa
```

* Separate diatoms

```bash
cat $Tree_dir/Separate_and_cluster_seqs/short_seqs_removed.faa | \
sed -E 's/(>.*)/£\1$/' | \
tr -d '\n' | \
tr '£' '\n' | \
egrep -i 'Thalassiosira|Phaeodactylum|Fragilariopsis|Amphora|Skeletonema|nitzschia' | \
tr '$' '\n' | \
sed '/^$/d' > $Tree_dir/Separate_and_cluster_seqs/diatom_seqs.faa
```

* Separate haptophytes

```bash
cat $Tree_dir/Separate_and_cluster_seqs/short_seqs_removed.faa | \
sed -E 's/(>.*)/£\1$/' | \
tr -d '\n' | \
tr '£' '\n' | \
egrep -i 'haptophyt|Calcidiscus|Chrysochromulina|Chrysoculter|Coccolithus|Emiliania|Exanthemachrysis|Gephyrocapsa|Imantonia|Isochrysis|Pavlova|Phaeocystis|Pleurochrysis|Prymnesium|Scyphosphaera|unid sp' | \
tr '$' '\n' | \
sed '/^$/d' > $Tree_dir/Separate_and_cluster_seqs/haptophyte_seqs.faa
```

* Separate cryptophytes

```bash
cat $Tree_dir/Separate_and_cluster_seqs/short_seqs_removed.faa | \
sed -E 's/(>.*)/£\1$/' | \
tr -d '\n' | \
tr '£' '\n' | \
egrep -i 'cryptophyt|Guillardia' | \
tr '$' '\n' | \
sed '/^$/d' > $Tree_dir/Separate_and_cluster_seqs/cryptophyte_seqs.faa
```

* Separate apicomplexa

```bash
cat $Tree_dir/Separate_and_cluster_seqs/short_seqs_removed.faa | \
sed -E 's/(>.*)/£\1$/' | \
tr -d '\n' | \
tr '£' '\n' | \
egrep -i 'Ascogregarina|Babesia|Besnoitia|Chromera|Cryptosporidium|Cyclospora|Cytauxzoon|Eimeria|Hepatocystis|Lankesteria|Plasmodium|Sarcocystis|Theilieria|Toxoplasma|Vitrella|apicomplexa' | \
tr '$' '\n' | \
sed '/^$/d' > $Tree_dir/Separate_and_cluster_seqs/apicomplexa_seqs.faa
```

* Separate ciliates

```bash
cat $Tree_dir/Separate_and_cluster_seqs/short_seqs_removed.faa | \
sed -E 's/(>.*)/£\1$/' | \
tr -d '\n' | \
tr '£' '\n' | \
egrep -i 'ciliates_seq|tetrahymena' | \
tr '$' '\n' | \
sed '/^$/d' > $Tree_dir/Separate_and_cluster_seqs/ciliates_seqs.faa
```

* Separate remaining sequences after removing key groups

```bash
cat $Tree_dir/Separate_and_cluster_seqs/short_seqs_removed.faa | \
sed -E 's/(>.*)/£\1$/' | \
tr -d '\n' | \
tr '£' '\n' | \
egrep -iv 'ciliates_seq|Tetrahymena|Ascogregarina|Babesia|Besnoitia|Chromera|Cryptosporidium|Cyclospora|Cytauxzoon|Eimeria|Hepatocystis|Lankesteria|Plasmodium|Sarcocystis|Theilieria|Toxoplasma|Vitrella|apicomplexa|aptophyt|Calcidiscus|Chrysochromulina|Chrysoculter|Coccolithus|Emiliania|Exanthemachrysis|Gephyrocapsa|Imantonia|Isochrysis|Pavlova|Phaeocystis|Pleurochrysis|Prymnesium|Scyphosphaera|unid sp|kareni|glenodinium|operidinium|durinskia|Thalassiosira|Phaeodactylum|Fragilariopsis|Amphora|Skeletonema|nitzschia|stramenopile|Thalassiosira|Phaeodactylum|Fragilariopsis|Aureococcus|Amphora|Skeletonema|Fragilariopsis|Ochromonas|Dinobryon|Pseudopedinella|Aureococcus|Heterosigma|Chattonella|Vaucheria|Ectocarpus|Nannochloropsis|Phytophthora|Albugo|Saprolegnia|Aplanochytrium|Aurantiochytrium|Schizochytrium' | \
tr '$' '\n' | \
sed '/^$/d' > $Tree_dir/Separate_and_cluster_seqs/short_seqs_key_taxa_removed.faa
```

## The following section clusters each subset of sequences to different levels of sequence identity using the program CD-hit.

1. Before clustering, each subset was aligned using MAFFT, and then trimmed (non-conserved sites removed) with trimAl.
2. Trimmed sequences were unaligned again by removing '-', else cd-hit did not work
* > The awk part `awk '{gsub(/-/,"",$2)} 1'` removes '-' only in the second field, essentially removing all non-aminoacid sites from each sequence
3. The clustering was then performed on this aligned-trimmed-unaligned dataset, so that clustering was only based on the identity of conserved regions of the proteins. However, after clustering, the full length sequences (not aligned-trimmed-unaligned sequences) were used for downstream analyses.

* Clustering the *Kareniaceae* sequences

```bash
cd-hit \
-i $Tree_dir/Separate_and_cluster_seqs/kareniaceae_seqs.faa \
-o $Tree_dir/Separate_and_cluster_seqs/kareniaceae_seqs.clustered.faa \
-c $kareniaceae_cluster_value \
-T 30
```

* Clustering the Stramenopile (excluding diatoms) sequences

```bash
printf "\n\nCLUSTERING STRAMENOPILES\n"

mafft --thread 30 $Tree_dir/Separate_and_cluster_seqs/stramenopile_seqs.faa > $Tree_dir/Separate_and_cluster_seqs/stramenopile_seqs.mafft

trimal -in $Tree_dir/Separate_and_cluster_seqs/stramenopile_seqs.mafft -out $Tree_dir/Separate_and_cluster_seqs/stramenopile_seqs.trimalautomated1 -automated1

cat $Tree_dir/Separate_and_cluster_seqs/stramenopile_seqs.trimalautomated1 | \
sed -E 's/(>.*)/£\1$/g' | \
tr -d '\n' | \
tr '£' '\n' | \
tr '$' '\t' | \
awk '{gsub(/-/,"",$2)} 1' | \
sed '/^$/d' | \
tr ' ' '\n' > $Tree_dir/Separate_and_cluster_seqs/stramenopile_seqs.trimalautomated1_unaligned_again

printf "\n\nCLUSTERING STRAMENOPILES\n"
cd-hit \
-i $Tree_dir/Separate_and_cluster_seqs/stramenopile_seqs.trimalautomated1_unaligned_again \
-o $Tree_dir/Separate_and_cluster_seqs/stramenopile_seqs.trimalautomated1.clustered.faa \
-c $stramenopile_cluster_value \
-T 30

cat $Tree_dir/Separate_and_cluster_seqs/stramenopile_seqs.trimalautomated1.clustered.faa | grep '>' | tr -d '>' > $Tree_dir/Separate_and_cluster_seqs/stramenopile_seqs.headers

seqtk subseq $Tree_dir/combined_seqs.faa $Tree_dir/Separate_and_cluster_seqs/stramenopile_seqs.headers > $Tree_dir/Separate_and_cluster_seqs/stramenopile_seqs.aligned_then_clustered.faa


printf "\n\nSTRAMENOPILE CLUSTERING COMPLETE\n\n\n\n\n\n"
printf "###############################################################\n"

```
This section clusters diatom sequences
```bash
printf "\n\nCLUSTERING DIATOMS"

mafft --thread 30 $Tree_dir/Separate_and_cluster_seqs/diatom_seqs.faa > $Tree_dir/Separate_and_cluster_seqs/diatom_seqs.mafft

trimal -in $Tree_dir/Separate_and_cluster_seqs/diatom_seqs.mafft -out $Tree_dir/Separate_and_cluster_seqs/diatom_seqs.trimalautomated1 -automated1

cat $Tree_dir/Separate_and_cluster_seqs/diatom_seqs.trimalautomated1 | \
sed -E 's/(>.*)/£\1$/g' | \
tr -d '\n' | \
tr '£' '\n' | \
tr '$' '\t' | \
awk '{gsub(/-/,"",$2)} 1' | \
sed '/^$/d' | \
tr ' ' '\n' > $Tree_dir/Separate_and_cluster_seqs/diatom_seqs.trimalautomated1_unaligned_again

cd-hit \
-i $Tree_dir/Separate_and_cluster_seqs/diatom_seqs.trimalautomated1_unaligned_again \
-o $Tree_dir/Separate_and_cluster_seqs/diatom_seqs.trimalautomated1.clustered.faa \
-c $diatom_cluster_value \
-T 30

cat $Tree_dir/Separate_and_cluster_seqs/diatom_seqs.trimalautomated1.clustered.faa | grep '>' | tr -d '>' > $Tree_dir/Separate_and_cluster_seqs/diatom_seqs.headers

seqtk subseq $Tree_dir/combined_seqs.faa $Tree_dir/Separate_and_cluster_seqs/diatom_seqs.headers > $Tree_dir/Separate_and_cluster_seqs/diatom_seqs.aligned_then_clustered.faa

printf "\n\nDIATOMS CLUSTERING COMPLETE\n\n\n\n\n\n"
printf "###############################################################\n"
```

* Clustering cryptophyte sequences

```bash
printf "\n\n\n\n\n\n\n\n"

printf "\n\nCLUSTERING CRYPTOPHYTES\n\n"

mafft --thread 30 $Tree_dir/Separate_and_cluster_seqs/cryptophyte_seqs.faa > $Tree_dir/Separate_and_cluster_seqs/cryptophyte_seqs.mafft

trimal -in $Tree_dir/Separate_and_cluster_seqs/cryptophyte_seqs.mafft -out $Tree_dir/Separate_and_cluster_seqs/cryptophyte_seqs.trimalautomated1 -automated1

cat $Tree_dir/Separate_and_cluster_seqs/cryptophyte_seqs.trimalautomated1 | \
sed -E 's/(>.*)/£\1$/g' | \
tr -d '\n' | \
tr '£' '\n' | \
tr '$' '\t' | \
awk '{gsub(/-/,"",$2)} 1' | \
sed '/^$/d' | \
tr ' ' '\n' > $Tree_dir/Separate_and_cluster_seqs/cryptophyte_seqs.trimalautomated1_unaligned_again

cd-hit \
-i $Tree_dir/Separate_and_cluster_seqs/cryptophyte_seqs.trimalautomated1_unaligned_again \
-o $Tree_dir/Separate_and_cluster_seqs/cryptophyte_seqs.trimalautomated1.clustered.faa \
-c $cryptophyte_cluster_value \
-T 30

cat $Tree_dir/Separate_and_cluster_seqs/cryptophyte_seqs.trimalautomated1.clustered.faa | grep '>' | tr -d '>' > $Tree_dir/Separate_and_cluster_seqs/cryptophyte_seqs.headers

seqtk subseq $Tree_dir/combined_seqs.faa $Tree_dir/Separate_and_cluster_seqs/cryptophyte_seqs.headers > $Tree_dir/Separate_and_cluster_seqs/cryptophyte_seqs.aligned_then_clustered.faa

printf "\n\nCRYPTOPHYTES CLUSTERING COMPLETE\n\n\n\n\n\n"
printf "###############################################################\n"
```

* Clustering haptophyte sequences

```bash
printf "\n\n\n\n\n\n\n\n"

printf "\n\nCLUSTERING HAPTOPHYTES\n\n"

mafft --thread 30 $Tree_dir/Separate_and_cluster_seqs/haptophyte_seqs.faa > $Tree_dir/Separate_and_cluster_seqs/haptophyte_seqs.mafft

trimal -in $Tree_dir/Separate_and_cluster_seqs/haptophyte_seqs.mafft -out $Tree_dir/Separate_and_cluster_seqs/haptophyte_seqs.trimalautomated1 -automated1

cat $Tree_dir/Separate_and_cluster_seqs/haptophyte_seqs.trimalautomated1 | \
sed -E 's/(>.*)/£\1$/g' | \
tr -d '\n' | \
tr '£' '\n' | \
tr '$' '\t' | \
awk '{gsub(/-/,"",$2)} 1' | \
sed '/^$/d' | \
tr ' ' '\n' > $Tree_dir/Separate_and_cluster_seqs/haptophyte_seqs.trimalautomated1_unaligned_again


cd-hit \
-i $Tree_dir/Separate_and_cluster_seqs/haptophyte_seqs.trimalautomated1_unaligned_again \
-o $Tree_dir/Separate_and_cluster_seqs/haptophyte_seqs.trimalautomated1.clustered.faa \
-c $haptophyte_cluster_value \
-T 30

cat $Tree_dir/Separate_and_cluster_seqs/haptophyte_seqs.trimalautomated1.clustered.faa | grep '>' | tr -d '>' > $Tree_dir/Separate_and_cluster_seqs/haptophyte_seqs.headers

seqtk subseq $Tree_dir/combined_seqs.faa $Tree_dir/Separate_and_cluster_seqs/haptophyte_seqs.headers > $Tree_dir/Separate_and_cluster_seqs/haptophyte_seqs.aligned_then_clustered.faa 

printf "\n\nHAPTOPHYTES CLUSTERING COMPLETE\n\n\n\n\n\n"
printf "###############################################################\n"
```

* Clustering the apicomplexa sequences

```bash
printf "\n\nCLUSTERING APICOMPLEXA"

mafft --thread 30 $Tree_dir/Separate_and_cluster_seqs/apicomplexa_seqs.faa > $Tree_dir/Separate_and_cluster_seqs/apicomplexa_seqs.mafft

trimal -in $Tree_dir/Separate_and_cluster_seqs/apicomplexa_seqs.mafft -out $Tree_dir/Separate_and_cluster_seqs/apicomplexa_seqs.trimalautomated1 -automated1
```
Need to unalign trimmed sequences again by removing '-' else cd-hit doesn't work
```bash
cat $Tree_dir/Separate_and_cluster_seqs/apicomplexa_seqs.trimalautomated1 | \
sed -E 's/(>.*)/£\1$/g' | \
tr -d '\n' | \
tr '£' '\n' | \
tr '$' '\t' | \
awk '{gsub(/-/,"",$2)} 1' | \
sed '/^$/d' | \
tr ' ' '\n' > $Tree_dir/Separate_and_cluster_seqs/apicomplexa_seqs.trimalautomated1_unaligned_again

cd-hit \
-i $Tree_dir/Separate_and_cluster_seqs/apicomplexa_seqs.trimalautomated1_unaligned_again \
-o $Tree_dir/Separate_and_cluster_seqs/apicomplexa_seqs.trimalautomated1.clustered.faa \
-c $apicomplexa_cluster_value \
-T 30

cat $Tree_dir/Separate_and_cluster_seqs/apicomplexa_seqs.trimalautomated1.clustered.faa | grep '>' | tr -d '>' > $Tree_dir/Separate_and_cluster_seqs/apicomplexa_seqs.headers

seqtk subseq $Tree_dir/combined_seqs.faa $Tree_dir/Separate_and_cluster_seqs/apicomplexa_seqs.headers > $Tree_dir/Separate_and_cluster_seqs/apicomplexa_seqs.aligned_then_clustered.faa

printf "\n\nAPICOMPLEXA CLUSTERING COMPLETE\n\n\n\n\n\n"
printf "###############################################################\n"
```

## The code in the following section clusters the file containing the remaining sequences after sequences from specific taxa were removed.

> Clustering might be more accurate if only aligned regions are considered, therefore the procedure of align-trim-unalign-cluster described above was applied


* Align, trim and unalign the sequences that remained after sequences from specified taxa were removed

```bash
mafft --thread 30 $Tree_dir/Separate_and_cluster_seqs/short_seqs_key_taxa_removed.faa > $Tree_dir/Separate_and_cluster_seqs/short_seqs_key_taxa_removed.mafft

trimal -in $Tree_dir/Separate_and_cluster_seqs/short_seqs_key_taxa_removed.mafft -out $Tree_dir/Separate_and_cluster_seqs/short_seqs_key_taxa_removed.trimalgappyout

cat $Tree_dir/Separate_and_cluster_seqs/short_seqs_key_taxa_removed.trimalgappyout | \
sed -E 's/(>.*)/£\1$/g' | \
tr -d '\n' | \
tr '£' '\n' | \
tr '$' '\t' | \
awk '{gsub(/-/,"",$2)} 1' | \
sed '/^$/d' | \
tr ' ' '\n' > $Tree_dir/Separate_and_cluster_seqs/short_seqs_key_taxa_removed.trimalgappyout_unaligned_again
```
* Cluster the sequences with key groups removed

```bash
cd-hit \
-i $Tree_dir/Separate_and_cluster_seqs/short_seqs_key_taxa_removed.trimalgappyout_unaligned_again \
-o $Tree_dir/Separate_and_cluster_seqs/short_seqs_key_taxa_removed.clustered \
-c $mainseqs_cluster_value \
-T 30
```
* Count the number of sequences before and after clustering to get an idea of how many sequences were removed.
```bash
preclust=$(grep '>' $Tree_dir/Separate_and_cluster_seqs/short_seqs_key_taxa_removed.trimalgappyout -c)
postclust=$(grep '>' $Tree_dir/Separate_and_cluster_seqs/short_seqs_key_taxa_removed.clustered -c)

echo "Before clustering: $preclust sequences"
echo "After clustering: $postclust sequences"
```
* Retrieve the full length sequences that were retained during the clustering using their fasta header
```bash
cat $Tree_dir/Separate_and_cluster_seqs/short_seqs_key_taxa_removed.clustered | \
grep '>' | \
tr -d '>' > $Tree_dir/Separate_and_cluster_seqs/short_seqs_key_taxa_removed_headers_after_clustering.txt 

seqtk subseq $Tree_dir/combined_seqs.faa $Tree_dir/Separate_and_cluster_seqs/short_seqs_key_taxa_removed_headers_after_clustering.txt > $Tree_dir/Separate_and_cluster_seqs/clustered_short_seqs_key_taxa_removed.faa 
```
## The code in the following section was used to concatenate the different sequence subsets (those divided at the beginning of this section) after they have each been clustered separately to different levels of sequence identity.

```bash
printf "\n\n########## Concatenating sequences ##########\n\n"

cat \
$Tree_dir/Separate_and_cluster_seqs/kareniaceae_seqs.clustered.faa \
$Tree_dir/Separate_and_cluster_seqs/peridiniales_seqs.faa \
$Tree_dir/Separate_and_cluster_seqs/stramenopile_seqs.aligned_then_clustered.faa \
$Tree_dir/Separate_and_cluster_seqs/diatom_seqs.aligned_then_clustered.faa \
$Tree_dir/Separate_and_cluster_seqs/haptophyte_seqs.aligned_then_clustered.faa \
$Tree_dir/Separate_and_cluster_seqs/apicomplexa_seqs.aligned_then_clustered.faa \
$Tree_dir/Separate_and_cluster_seqs/ciliates_seqs.faa \
$Tree_dir/Separate_and_cluster_seqs/cryptophyte_seqs.aligned_then_clustered.faa \
$Tree_dir/Separate_and_cluster_seqs/clustered_short_seqs_key_taxa_removed.faa \
> $Tree_dir/Separate_and_cluster_seqs/combined_seqs.faa

printf "\n\n########## COMPLETE ##########\n\n"
```

## Iterative dataset refinement by phylogenetic inference

The clustered dataset obtained in the previous section was then further refined using a manual approach, by aligning the sequences and inferring a phylogeny, then based on this, removing sequences that aligned poorly with the rest or formed long branches, suggesting they were either not homologues, or were distant paralogues of the protein being analysed. The redundancy of the datasets was further reduced by removal of sequences from clades that were still over-represented in the dataset by large numbers of sequences, even after the clustering steps in the above sections.

```bash
#!/bin/bash/

module load mafft/7.475
module load trimal/1.4 
module load FastTree/2.1.11

#Finds sequences with duplicate headers and keeps one
cat combined_seqs.faa | sed -E 's/>(.*)/£>\1$/' | tr -d '\n' | tr '£' '\n' | tr '$' '\t' | \
sort -k 1,1 | sort -k 1,1 -u | tr '\t' '\n' > combined_seqs.faaa

mv combined_seqs.faaa combined_seqs.faa

mafft --thread 22 combined_seqs.faa > combined_seqs.mafftaln

trimal -in combined_seqs.mafftaln -out combined_seqs.mafft.trimalautomated1 -automated1

FastTreeMP combined_seqs.mafft.trimalautomated1 > combined_seqs.fasttree

```
Once a dataset had been suitably refined, the program iqtree2 was used to infer a final maximum likelihood tree and estimate support values, as described in the methods section of the manuscript. These final datasets, as well as the log files, including the commands used to run iqtree2, can be downloaded from FigShare, [here](https://doi.org/10.6084/m9.figshare.21602697).

# References
* Keeling, P.J. et al. (2014) ‘The Marine Microbial Eukaryote Transcriptome Sequencing Project (MMETSP): illuminating the functional diversity of eukaryotic life in the oceans through transcriptome sequencing’, PLoS biology, 12(6), p. e1001889. ([link](https://journals.plos.org/plosbiology/article?id=10.1371/journal.pbio.1001889))
* Johnson, L.K., Alexander, H. and Brown, C.T. (2017) ‘(all datasets) MMETSP re-assemblies’ ([link](https://doi.org/10.5281/zenodo.3247846))
* Amos, B. et al. (2022) ‘VEuPathDB: the eukaryotic pathogen, vector and host bioinformatics resource center’, Nucleic acids research, 50(D1), pp. D898–D911. ([link](https://academic.oup.com/nar/article-lookup/doi/10.1093/nar/gkab929))
* Hehenberger, E., Gast, R.J. and Keeling, P.J. (2019) ‘A kleptoplastidic dinoflagellate and the tipping point between transient and fully integrated plastid endosymbiosis’, Proceedings of the National Academy of Sciences of the United States of America, 116(36), pp. 17934–17942. ([link](https://www.pnas.org/doi/10.1073/pnas.1910121116))
* Janouškovec, J. et al. (2019) ‘Apicomplexan-like parasites are polyphyletic and widely but selectively dependent on cryptic plastid organelles’, eLife, 8. ([link](https://elifesciences.org/articles/49662))
* Mathur, V. et al. (2019) ‘Multiple Independent Origins of Apicomplexan-Like Parasites’, Current biology: CB, 29(17), pp. 2936–2941.e5. ([link](https://www.sciencedirect.com/science/article/pii/S0960982219308644))
