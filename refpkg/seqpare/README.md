# Seqpare

Seqpare is a self-consistent metric and tool for comparing sequences based on the total effective overlaps of their interval sets. With this metric, the similarity of two interval sets is quantified by a single index, which directly represents the percentage of their effective overlapping: a similarity index of zero indicates totally unrelated sequences, and an index of one means that the interval sets are exactly the same. Seqpare tool provides functions for both searching and mapping large-scale interval datasets. The seqpare code is based on [ailist][ailist].

## Installing Seqpare

If zlib is not already installed, install it:
```
sudo apt-get install libpng12-0
```
Then:

```
git clone https://github.com/deepstanding/seqpare
cd seqpare
make
sudo cp bin/seqpare /usr/local/bin
```
## Running Seqpare

### Compare two interval sets

```
seqpare <BED file 1> <BED file 2> -o filename

Output: N1 (number of intervals in BED file 1), N2, teo (total effective-overlap), tc (total counts of overlaps), sm (seqpare similarity).

option:
-o: save the result to a file
```

### Search: compare a query interval set over a collection of interval sets

```
seqpare <BED file folder> <query BED file> -m 1 -o filename

Output: a list of comparison result (N1, N2, teo, tc and sm) for each file in the folder.
See the example file AffyGnf1h_ucsc100_seqpare in the Test_results folder.
```

### Similarity mapping of two different collections of interval sets  

```
seqpare <BED file folder 1> <query BED file folder 2> -m 2 -o filename

Output format: tsv. The first line: M1 (number of files in folder 1), M2 (number of files in folder 2);
Then each line lists M2 sm values for a file in folder 1.
```

### Self similarity mapping of a collection of interval sets

```
seqpare <BED file folder> -m 3 -o filename

Output format: tsv. The first line: M (number of files in folder), M (second dimension);
Then each line lists M sm values for a file in the folder.
See the example file denaseseq.sm in the Test_results folder.
```

## Data source

http://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/

## Test results

Four files in the Test_results folder.

[ailist]: https://github.com/databio/AIList

## Verify the result

The original data has 100 intervals sets of size ~5GB. Because of the space limit of Github, only 30 small datasets (~60MB) are included in the folder UCSC_30.

To run seqpare:
```
time seqpare "ucsc_30/*" "ucsc_30/affyGnf1h.bed.gz" -m 1 -o ucsc30_seqpare
```
It takes 2.99s on a laptop. The output file ucsc30_seqpare lists the similarity between affyGnf1h interval set and all 30 interval sets.

To run GIGGLE, download giggle from https://github.com/ryanlayer/giggle, then
```
ulimit -Sn 16384
giggle index -i "ucsc_30/*.*" -o ucsc_30.giggle -s
time giggle search -i ucsc_30.giggle -q "ucsc_30/affyGnf1h.bed.gz" -s > affyGnf1h_ucsc30_giggle
```
It takes 5m12.81s. The output file affyGnf1h_uscs30_giggle contains the p-values and odds-ratios.