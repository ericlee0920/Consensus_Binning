# Consensus_Binning
## Updated: December 10, 2019.
### ABOUT Version 1.0
This approach identifies pairwise MAG similarity through both average nucleotide identity (ANI) and min-wise independent permutations locality sensitive hashing methods (MinHash). The pairwise comparisons are then deduplicated, clouded, and joined to form consensus bins. Using this consensus binning technique, we have robustly and in a heuristically efficient way binned 1994 high quality MAGs from an anaerobic digester time-series dataset into ~200 consensus bins. 
## Work Flow
Make sure FastANI and Sourmash are downloaded to your user server.
### 1. 90% Completion and 5% Contamination
Ensure that your directory has MAGs of 90% completion and 5% contamination only. This ensures that it is high quality. If in doubt check your CheckM or equivalent output.
### 2. FastANI
Obtain many to many comparisons.
```
fastANI -ql all_MAGs.txt --rl all_MAGs.txt -o fastani.out
```
### 3. Sourmash
Compute a nonscaled signature from reads.
```
sourmash compute [FILENAME] -o [FILENAME].sig
```
Compare signatures within file.
```
sourmash compare MASH_SIGs/*.sig -k 31 --csv MAG_distances.csv
```
### 4. Deduplication and Clouding
In your local directory create a Rproject and run the following R script with your FastANI output.
```
run pairwise2matrix.r
```
Also, run the following R script with your Sourmash output.
NOTE that your output may have one row of column names missing, please add that before using R.
```
run matrix2pairwise.r
```
### 5. Joining
In your local directory create a Rproject and run the following R script with your csv outputs from Step 4.
```
run consensus_join.r
```
The output for this R script will be four files.
```
consensus_bin_outer.csv: the result of outer join
consensus_bin_inner.csv: the result of inner join
consensus_bin_outer_conserve.csv: the result of outer join with lineage trim
consensus_bin_inner_conserve.csv: the result of inner join with lingege trim
```
