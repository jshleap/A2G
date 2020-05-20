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

## Basic usage
A2G<sup>2</sup> will give you help by:
```bash
A2G -h
```
this should give you something like this:

```bash
A2G version: 2020.0.1
Copyright 2020 Jose Sergio Hleap

usage: A2G [-h] [--cpus CPUS] [--nowrite] [--out_prefix OUT_PREFIX]
           [--remove_duplicates]
           global_consensus local_consensus fasta

positional arguments:
  global_consensus      Sequence consensus of the global region, e.g. full COI
  local_consensus       Sequence consensus of the local region, e.g. Leray
                        fragment
  fasta                 fasta file with the focal sequences

optional arguments:
  -h, --help            show this help message and exit
  --cpus CPUS           number of cpus to use
  --nowrite             return string instead of writing
  --out_prefix OUT_PREFIX
                        Prefix of outputs
  --remove_duplicates   Keep or remove duplicated sequences
```

Then to run it, you can simply type:

```bash
A2G global_target local_target query_file --cpus 10 --out_prefix prefix --remove_duplicates
```
With this command, you will use the `global_target` as the overall region, the `local_target` as the amplicon reference
 sequence to anchor the query sequences, and `query_file` contains your query sequences. Those are the required 
 arguments. The optional arguments allow you to control the execution. `--cpus` allow you to provide the number of cpus
 to use. In the example, up to 10 cpus will be used. `--out_prefix`change the prefix of the outputs generated. Finally,
  the `--remove_duplicates` option will retain only unique sequences.

If the `no_write` option is used, A2G<sup>2</sup> will output the alignment
to standard out, and other info to standard error. If you would like to pipe
only the alignment, you can redirect the standard error to a null device:

```bash
A2G global_target local_target query_file --cpus 10 --out_prefix prefix --no_write 2> /dev/null
```
