# NucFreq plots
Script for making nucleotide frequency plots 
![clean](imgs/image.png)

# Usage 
```
usage: NucPlot.py [-h] [-d] [--legend] [--zerostart] [-a] [-r REPEATMASKER]
                  [--regions [REGIONS [REGIONS ...]]] [--bed BED] [-y YLIM]
                  [--height HEIGHT] [-w WIDTH] [-t THREADS] [--header]
                  [--psvsites PSVSITES] [-s] [-c MINCLIP]
                  infile outfile

positional arguments:
  infile                input bam file
  outfile               output plot file

optional arguments:
  -h, --help            show this help message and exit
  -d
  --legend
  --zerostart
  -a                    output all positions (default: False)
  -r REPEATMASKER, --repeatmasker REPEATMASKER
                        rm out to add to plot (default: None)
  --regions [REGIONS [REGIONS ...]]
                        regions in this format (.*):(\d+)-(\d+) (default:
                        None)
  --bed BED             bed file with regions to plot (default: None)
  -y YLIM, --ylim YLIM  max y axis limit (default: None)
  --height HEIGHT       figure height (default: 9)
  -w WIDTH, --width WIDTH
                        figure width (default: 16)
  -t THREADS, --threads THREADS
                        [8] (default: 8)
  --header
  --psvsites PSVSITES   CC/mi.gml.sites (default: None)
  -s, --soft
  -c MINCLIP, --minclip MINCLIP
                        min number of clippsed bases in order to be displayed
                        (default: 1000)
```
      
## Detecting heterozygous sites with NucFreq
In order to detect heterozygous sites in centromeric regions, 
   we first aligned CHM13 PacBio HiFi reads to the entire CHM13 v1.0 assembly using
   pbmm2 and the following parameters: 
```
pbmm2 align --log-level DEBUG --preset SUBREAD --min-length 5000 -j 8
```
   Then, we filtered the alignments to remove secondary and partial alignment using SAMtools
   flag 2308, generating a BAM file with genome-wide alignments.
   Then, we filtered the BAM to only the cenhap regions using SAMtools.
   We used `NucPlot.py` to determine the frequency of the first and second most
   common bases in the aligned PacBio HiFi reads with the following command:
```
NucPlot.py --obed {output.bed} --bed {region.bed} --minobed 2 {input.bam} {output.png}
```
   The resulting bed file was used to identify regions where the second most common base
   was present in at least 10% of reads, and this occurred in at least five positions
   within each region using the script `HetDetection.R`.

