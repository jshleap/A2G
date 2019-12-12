# Amplicons to Global Gene (A2G<sup>2</sup>)

This program implements the progressive algorithm to align a large set
of amplicons to a reference gene consensus, or a large set of sequences
to an amplicon consensus, based on a reference consensus. This program
makes use of traditional multiple aligners such as MAFFT (default), and 
muscle, and can be extended to other aligners.

## Problem
Some taxonomic assignment software require a set of align sequences, 
both in the query as in the reference. Projects such as those using
environmental DNA (eDNA) or trying to assess wide diversity using 
metagenomics often have a hard time creating such alignments, because of
memory and computational restrictions. Another observation is that 
massive alignments often introduce more gaps in the sequences, and force
alignment of segments that should not align in that region.
Here is where A2G<sup>2</sup> will use a global to local alignment to 
avoid such issues, and retained the ungapped alignment of the amplicons.