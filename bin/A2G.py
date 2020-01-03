#!/usr/bin/env python

from A2G.align2consensus import *

parser = argparse.ArgumentParser()
parser.add_argument('global_consensus', help='Sequence consensus of the '
                                             'global region, e.g. full COI'
                    )
parser.add_argument('local_consensus',
                    help='Sequence consensus of the local region, e.g. '
                         'Leray fragment')
parser.add_argument('fasta', help='fasta file with the focal sequences')
parser.add_argument('--cpus', help='number of cpus to use', default=-1,
                    type=int)
parser.add_argument('--nowrite', help='return string instead of writing',
                    action='store_false', default=False)

args = parser.parse_args()
main(args.global_consensus, args.local_consensus, args.fasta,
            no_write=args.nowrite, cpus=args.cpu)