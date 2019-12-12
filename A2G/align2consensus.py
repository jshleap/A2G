"""
**align2consensus.py** Program to align amplicon sequences to gene consensus
**Copyright** 2019  Jose Sergio Hleap

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.

E-mail: jshleap@gmail.com
"""
import argparse
import os
from subprocess import run, PIPE, CalledProcessError
from typing import List

import dill
from sklearn.ensemble import IsolationForest
from tqdm.auto import tqdm

from A2G.SeqIO import *
from justblast.utils import FastX


# TODO: add checkpoints within the alignment


def execute(args: List[str], inpt: str):
    st = run(args, input=inpt.encode('utf-8'), env=os.environ, stdout=PIPE,
             stderr=PIPE)
    try:
        st.check_returncode()
    except CalledProcessError:
        for x in inpt.split('\n'):
            if '>' in x:
                print(x)
        raise Exception(st.stderr.decode('utf-8'), st.stdout.decode('utf-8'))
    return st


def pickle_at_exit(align_instance):
    with open('intermediate.checkpoint', 'w') as d:
        dill.dump(align_instance, d)


class Align(object):
    current = None

    def __init__(self, gene_consensus: str, amplicon_consensus: str,
                 query: str, no_write: bool = False, cpus: int = -1):
        self.out_prefix = query[: query.rfind('.')]
        self.no_write = no_write
        self.cpus = cpus
        self.gene_consensus = gene_consensus
        self.amplicon_consensus = amplicon_consensus
        self.query = query
        self.clf = IsolationForest(behaviour="new", max_samples='auto',
                                   random_state=12345,  contamination='auto',
                                   n_jobs=cpus)

    @property
    def query(self):
        return self.__query

    @query.setter
    def query(self, query: str):
        self.__query = FastX(query, cpus=self.cpus)  # ,unique=True)

    @staticmethod
    def triwise(gene_consensus: str, amplicon_consensus: str, query: str,
                entropy: bool = True, long: bool = False):
        executable = os.path.abspath(os.path.join(os.path.dirname(__file__),
                                                  os.pardir, 'bin', 'mafft'))
        mafft = [executable, '--add', '-']
        # TODO: Implement plugins
        with open(amplicon_consensus) as ac, open(gene_consensus) as gec:
            fasta = '%s%s%s' % (ac.read(), gec.read(), query)
        if len(query) > 10000 and not long:
            # Very long sequences might kill the process
            return None, None, query
        st = execute(mafft, fasta)
        alignment = st.stdout.decode('utf-8')
        with open('current.aln', 'w') as a:
            a.write(alignment)
        aln = Alignment(alignment, outliers_detection=entropy, cpus=1)
        q = aln.seq.iloc[2]
        # ungap = aln.seq.iloc[[1,2]].apply(lambda x: (x != '-').all(), axis=0)
        ungap = aln.seq.iloc[[0,1]].apply(lambda x: (x != '-').all(), axis=0)
        # idx = q.where(aln.seq.iloc[2] != '-').dropna().index.values
        aln.subset_idx = (None, ungap)
        aln.entropy = True
        median_entropy = aln.entropy.median()
        aln.subset_idx = ([q.name], None)
        fasta = str(aln)
        # intermediate_checkpoint(fasta, median_entropy, longs)
        return fasta, median_entropy, None

    def run(self):
        dill_name = '%s_checkpoint.dill' % self.out_prefix
        if os.path.isfile(dill_name):
            with open(dill_name, 'rb') as d:
                results = dill.load(d)
        else:
            results = Parallel(n_jobs=self.cpus)(
                delayed(self.triwise)(self.gene_consensus,
                                      self.amplicon_consensus, sequence,
                                      entropy=True) for sequence in
                tqdm(self.query.yield_seq(), desc='Aligning',
                     total=len(self.query)))
            with open(dill_name, 'wb') as d:
                dill.dump(results, d)
        fasta, shannon, longs = zip(*results)
        longs = list(filter(lambda x: x is not None, longs))
        fasta = list(filter(lambda x: x is not None, fasta))
        shannon = list(filter(lambda x: x is not None, shannon))
        if longs:
            long_dill = '%s_long.dill' % dill_name[:dill_name.rfind('.')]
            if os.path.isfile(long_dill):
                with open(long_dill, 'rb') as d:
                    lon = dill.load(d)
            else:
                lon = [
                    self.triwise(self.gene_consensus, self.amplicon_consensus,
                                 long, entropy=True, long=True) for long in
                    tqdm(longs, desc='Aligning long sequences',
                         total=len(longs))]
                with open(long_dill, 'wb') as d:
                    dill.dump(lon, d)
            fa, sh, _ = zip(*lon)
            fasta.extend(fa)
            shannon.extend(sh)
        fasta = [x for x in fasta if x is not None]
        shannon = np.array([x for x in shannon if x is not None]).reshape(
            -1, 1)
        outl = self.clf.fit_predict(shannon)
        l = "Outliers removed in %s_aligned.withoutoutliers:" % self.out_prefix
        print(l, sum(outl < 0))
        if self.no_write:
            return list(zip(fasta, shannon)), outl

        with open('%s_aligned.fasta' % self.out_prefix, 'w') as o, open(
                '%s_aligned.withoutliers' % self.out_prefix, 'w') as w:
            o.write('\n'.join(fasta))
            w.write('\n'.join(x for i, x in enumerate(fasta) if outl[i] != -1))


if __name__ == '__main__':
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
    aln = Align(args.global_consensus, args.local_consensus, args.fasta,
                no_write=args.nowrite, cpus=args.cpus)
    aln.run()
