## Dependencies
### Software
 - GMAP
 - MAFFT
### Python Modules
 - outlier-utils

## Installation
```bash
cd /path/to/install
git clone https://github.com/sc-zhang/PanMarker.git
chmod +x PanMarker/*.py
echo 'export PATH=/path/to/install/PanMarker:$PATH' >> ~/.bash_profile
source ~/.bash_profile
```

## Usage
1. Mapping reference cds to all genomes in pan-genome
```bash
# For each sample, run command below
gmap_build -D . -d DB genomeN.fasta
gmap -D . -d DB -f 2 -n 1 -t 20 ref.cds > genomeN.gff3
```

2. Extract cds with gff3 file generated by gmap
```bash
# For each sample, run command below
./get_promoter_and_cds.py genomeN.gff3 genomeN.fasta 2000 genomeN
```

3. Multi alignment
```bash
# For each gene and each sample, run command below
# Notice for each gene_id.fa, it should be:
# >sample_id1
# XXXXXXXXX...
# >sample_id2
# XXXXXXXXX...
echo ">sample_id" >> gene_id.fa
grep -A1 gene_id genomeN.cds | tail -n1 >> gene_id.fa
mafft --adjustdirection gene_id.fa > gene_id.aln
```

4. Extract variants
```bash
# For each gene, run command below
./extract_var.py gene_id.aln gene_id.var
```

5. Extart variant sites associate with traits
```bash
# For each gene and trait
./stat_var_vs_exp.py gene_id.var trait.txt gene_id_trait.stat
# Notice: trait.txt file should contain two columns: sample_id value
```