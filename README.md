# Genomics_Sequence_Assembly
## Description: 
Write an algorithm for assembling reads from scratch, visualize the read density of contigs,
and evaluate the sequencing method by normalizing the density plot.


## Task breakdown:
1. Create a function called “assemble reads” which takes as input a file with reads and
outputs a file with a FASTA sequence of contigs labeled numerically <br>
a. Code meets specifications – the function exists and makes correct input and
output <br>
b. No pair of reads are incorrectly merged <br>
c. Code generates accurate mapping <br>

2. Evaluate the distribution of the reads across each sequence <br>
a. Create a visualization of your choice that shows the coverage of reads across each of
your sequences <br>


## Evaluation:
1. Based on the visualization, is the sequencing method biased? Explain<br>
Ans: Yes, after normalization (the second graph), my observation is that there is no consistency
pattern, and all 22 contigs are very. This insight may indicate a systematic bias in the
sequencing method. Also, there are uneven depths across contigs (the first graph), and
many contigs have sudden drop curves in the start/end position or flat lines, which are also
suspicious and might be biased.

2. Consider the assumptions you made in your algorithm – what’s the issue?<br>
Ans: This algorithm has no tolerance, and is too sensitive.
I use k ahead and behind of sequence to search the overlaps, then merge them separately
and join back the contig in the corresponding position. However, the k parameter is very
sensitive to a shorter sequence length. That's why, in the last run, the algorithm ignored the
shorter sequences and added again, which caused a repeat in the contig tail.

3. What is the relationship between how high k can go and sequencing coverage?<br>
As the k incremented, the contig number also increased significantly because the unmapped
sequence creates contigs on its own. Thus, the higher k is, the more fragmented it becomes
in sequencing coverage with shorter contigs.


## Summary: 
We used assembly algorithms such as BWA for NGS read mapping by software like Clustal Omega
and other pipeline tools we had learned from bioinformatic courses. However, due to naive
experience, we could not evaluate the algorithm and its applied data. Through this practice,
we know there are pros and cons of parameters used in an algorithm, and we have to be aware of
potentials before/while processing the data and always check the quality using visualization
tools to be mindful of the bias in advance.  
