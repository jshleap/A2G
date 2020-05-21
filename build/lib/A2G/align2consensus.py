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
import sys
import tempfile
from subprocess import run, PIPE, CalledProcessError
from typing import List

import dill
from psutil import Process, cpu_count
from sklearn.ensemble import IsolationForest
from tqdm.auto import tqdm

from A2G.SeqIO import *
from justblast.utils import FastX


# TODO: add checkpoints within the alignment


def execute(args: List[str], inpt: Optional[str]):
    if inpt is not None:
        inpt = inpt.encode('utf-8')
    st = run(args, input=inpt, env=os.environ, stdout=PIPE,
             stderr=PIPE)
    try:
        st.check_returncode()
    except CalledProcessError:
        if inpt is not None:
            for x in inpt.split(b'\n'):
                if b'>' in x:
                    print(x, file=sys.stderr)
        raise Exception(st.stderr.decode('utf-8'), st.stdout.decode('utf-8'))
    return st


def pickle_at_exit(align_instance):
    with open('intermediate.checkpoint', 'w') as d:
        dill.dump(align_instance, d)


def exist_n_not_empty(file_path):
    return os.path.isfile(file_path) and os.path.getsize(file_path) > 0


class Align(object):
    current = None

    def __init__(self, gene_consensus: str, amplicon_consensus: str,
                 out_prefix: str = 'A2G_aln', no_write: bool = False,
                 outliers: bool = True, remove_duplicates: bool = True,
                 cpus: int = -1, quiet=False):
        self.quiet = quiet
        self.out_prefix = out_prefix
        self.entropy = outliers
        self.remove_duplicates = remove_duplicates
        self.no_write = no_write
        self.cpus = cpus
        self.gene_consensus = gene_consensus
        self.amplicon_consensus = amplicon_consensus
        self.reference = (gene_consensus, amplicon_consensus)
        self.clf = IsolationForest(max_samples='auto', random_state=12345,
                                   contamination='auto', n_jobs=cpus)

    @property
    def reference(self):
        return self.__reference

    @reference.setter
    def reference(self, sequences: Tuple[str, str]):
        self.__reference = '{}_reference.aln'.format(self.out_prefix)
        if not exist_n_not_empty(self.__reference):
            gene_consensus, amplicon_consensus = sequences
            executable = 'mafft'
            mafft = [executable, '--auto', '-']
            with open(gene_consensus) as gc, open(amplicon_consensus) as ac:
                refs = gc.read() + ac.read()
            st = execute(mafft, refs)
            alignment = st.stdout.decode('utf-8')
            with open(self.__reference, 'w') as r:
                r.write(alignment)

    @property
    def query(self):
        return self.__query

    @query.setter
    def query(self, query: str):
        self.__query = FastX(query, cpus=self.cpus,
                             unique=self.remove_duplicates)
        # if not isinstance(self.__query.handle, StringIO):
        #     self.out_prefix = os.path.basename(query[: query.rfind('.')])

    @staticmethod
    def triwise(reference: str, query: str, entropy: bool = True,
                long: bool = False, current_file: str = 'current.aln'):
        executable = 'mafft'
        fasta, median_entropy = None, None
        if len(query) > 50000 and not long:
            # Very long sequences might kill the process
            return fasta, median_entropy, query

        with tempfile.NamedTemporaryFile() as temp:
            temp.write(query.encode('utf-8'))
            temp.seek(0)
            mafft = [executable, '--thread', '1', '--keeplength', '--add',
                     temp.name, reference]
            # TODO: Implement plugins
            st = execute(mafft, None)
            alignment = st.stdout.decode('utf-8')
            with open(current_file, 'w') as a:
                a.write(alignment)
            aln = Alignment(alignment, outliers_detection=False, cpus=1)
            q = aln.seq.iloc[2]
            ungap = aln.seq.iloc[[0,1]].apply(lambda x: (x != '-').all(),
                                              axis=0)
            aln.subset_idx = (None, ungap)
            if entropy:
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
            kwargs = dict(entropy=self.entropy,
                          current_file='%s_current.aln' % self.out_prefix)
            results = Parallel(n_jobs=self.cpus)(
                delayed(self.triwise)(self.reference, sequence, **kwargs)
                for sequence in tqdm(self.query.yield_seq(), desc='Aligning',
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
                lon = [self.triwise(self.reference, long, long=True) for long
                       in tqdm(longs, desc='Aligning long sequences',
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
        if not self.quiet:
            l = "Outliers removed in %s_aligned.withoutoutliers:" % self.out_prefix
            print(l, sum(outl < 0), file=sys.stderr)
        subset = [x for i, x in enumerate(fasta) if outl[i] != -1]
        if self.no_write:
            return '\n'.join(fasta), '\n'.join(subset)
        with open('%s_aligned.fasta' % self.out_prefix, 'w') as o, open(
                '%s_aligned.withoutliers' % self.out_prefix, 'w') as w:
            o.write('\n'.join(fasta))
            w.write('\n'.join(subset))
        return None, None


def main(global_consensus: str, local_consensus: str, fasta: str,
         no_write: bool = False, cpus: int = -1, out_prefix: str = 'A2G_aln',
         remove_duplicates: bool = True):
    p = Process()
    n_cpus = cpus if cpus > 0 else cpu_count()
    all_cpus = list(range(n_cpus))
    p.cpu_affinity(all_cpus)
    kwargs = dict(out_prefix=out_prefix, no_write=no_write, cpus=cpus,
                  remove_duplicates=remove_duplicates)
    aln = Align(global_consensus, local_consensus, **kwargs)
    aln.query = fasta
    full, subset = aln.run()
    if no_write:
        return full, subset


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
                        action='store_true', default=False)
    parser.add_argument('--remove_duplicates', action='store_false',
                        help='Keep or remove duplicated sequences',
                        default=True)
    parser.add_argument('--out_prefix', action='store', default='A2G_aln',
                        help='Prefix of outputs')
    parser.add_argument('--mpi', action='store_true', default=False,
                        help='Use MPI backend in multinode executions')

    args = parser.parse_args()

    fasta = main(args.global_consensus, args.local_consensus, args.fasta,
                 no_write=args.nowrite, cpus=args.cpus,
                 remove_duplicates=args.remove_duplicates)
    if fasta is not None:
        print(fasta[0], file=sys.stdout)
