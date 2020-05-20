"""
**justblast.py** Python wrapper for multiprocessed blast
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

E-mail: jshleap@squalus.org
"""
import sys
import argparse
from io import BytesIO

import pandas as pd

from justblast.utils import *

plt.style.use('ggplot')

cols = 'qseqid sseqid pident evalue qcovs qlen length staxid stitle'
basta_dir = os.path.expanduser("~/.basta/taxonomy")
my_env = os.environ

class Blast(object):
    """
    Class wrapper for blast
    """
    blasts = None
    taxadb = basta_dir

    def __init__(self, db: str, query: str, evalue: int = 10, p_id: float = 0,
                 max_target_seqs: int = 500, cpus: int = -1,
                 query_coverage: Optional[float] = None,
                 identify: bool = False, out: str = 'hit.hits',
                 outfmt_str: str = cols, unique: bool = False) -> None:
        self.db = db
        self.query = query
        self.evalue = evalue
        self.unique = unique
        self.p_id = p_id
        self.max_target_seqs = max_target_seqs
        self.cpus = cpus
        self.out = out
        self.outfmt_str = outfmt_str
        self.columns = outfmt_str.split()
        self.run()
        self.query_coverage = query_coverage
        self.identify = identify

    @property
    def query_coverage(self) -> Optional[float]:
        return self.__query_coverage

    @query_coverage.setter
    def query_coverage(self, query_coverage:  Optional[float] = None) -> None:
        if query_coverage is None:
            self.__query_coverage = None
        else:
            self.__query_coverage = query_coverage
            self.filter()

    @property
    def identify(self) -> bool:
        return self.__identify

    @identify.setter
    def identify(self, identify: bool) -> None:
        self.__identify = identify
        if identify:
            if 'staxid' not in self.columns:
                raise Exception("staxid not requested. Please rerun blast")
            taxafile = os.path.join(self.taxadb, 'complete_taxa.db', 'CURRENT')
            gbfile = os.path.join(self.taxadb, 'gb_mapping.db', 'CURRENT')
            if not os.path.isfile(taxafile):
                args = ['basta', 'taxonomy']
                _ = stdin_run(args, None, env=os.environ)
            if not os.path.isfile(gbfile):
                args2 = ['basta', 'download', 'gb']
                _ = stdin_run(args2, None, env=os.environ)
            self.write()
            prefix = self.out[: self.out.rfind('.')]
            if not os.path.isfile('config.txt'):
                self.basta_config()
            basta_args = ['basta', 'sequence', self.out, '%s.basta' % prefix,
                          'gb', '-e', str(self.evalue), '-n', '50', '-m', '10',
                          '-c', 'confix.txt']
            _ = stdin_run(basta_args, None, env=os.environ)
            self.read_basta()

    def basta_config(self) -> None:
        line = "query_id\t%d\nsubject_id\t%d\nalign_length\t%d\nevalue\t%d\n" \
               "pident\t%d"
        line = line % (self.columns.index('qseqid'), self.columns.index(
            'sseqid'), self.columns.index('length'), self.columns.index(
            'evalue'), self.columns.index('pident'))
        with open('confix.txt', 'w') as c:
            c.write(line)

    @staticmethod
    def plot_basta(lineages: pd.core.frame.DataFrame, out: str) -> None:
        levels = lineages.shape[1]
        rows = int(np.floor(levels/2))
        cols = levels - rows
        idx_cols = iter(lineages.columns)
        fig, axes = plt.subplots(nrows=rows, ncols=cols)
        for row in range(rows):
            for col in range(cols):
                ax = axes[row, col]
                try:
                    idx = next(idx_cols)
                    data = [y for y in lineages[idx] if y is not
                            None]
                    bins = np.arange(len(set(data))) - 0.5
                    ax.hist(data, bins=bins)
                    ax.set_xticklabels(data, size=5, rotation=90)
                    ax.set_yticklabels(ax.get_yticks(), size=5)
                    ax.set_title("Level %d" % idx, fontsize=8)
                except StopIteration:
                    fig.delaxes(ax)
        plt.tight_layout()
        fig.savefig(out)
        plt.close()

    def read_basta(self) -> None:
        fn = self.out[: self.out.rfind('.')]
        basta = pd.read_csv(fn + '.basta', sep='\t', header=None, names=[
            'qseqid', 'lineage'])
        lineages = basta.lineage.str.strip().str.strip(';')
        lineages = lineages.str.split(';', expand=True)
        self.plot_basta(lineages, fn + '_taxadist.pdf')
        self.blasts = self.blasts.merge(basta, on='qseqid', how='left')
        self.out = self.out[:self.out.rfind('.')] + '_annotated.hits'
        self.write()

    def run(self) -> None:
        """
        Run blast using multiprocessing
        """
        if not os.path.isfile(self.out):
            fasta = FastX(self.query, unique=self.unique, cpus=self.cpus)
            args = ['blastn', '-db', self.db, '-query', '-', '-evalue',
                    str(self.evalue), '-perc_identity', str(self.p_id),
                    '-outfmt', '6 %s' % self.outfmt_str, '-max_target_seqs',
                    str(self.max_target_seqs)]
            blasts = Parallel(n_jobs=self.cpus)(
                delayed(stdin_run)(args, inp, env=my_env) for inp in tqdm(
                    fasta.yield_seq()))
            print(fasta.n_duplicates, 'Duplicates in', self.query,
                  file=sys.stderr)
            blasts = [pd.read_table(BytesIO(x), header=None, names=self.columns
                                    )for x in blasts]
            self.blasts = pd.concat(blasts)
        else:
            self.blasts = pd.read_csv(self.out, sep='\t', header=None,
                                      names=self.columns)

    def filter(self) -> None:
        """
        Filter the blast table
        TODO: Add more filters, for now just query coverage
        """

        if self.query_coverage is not None:
            if 'qcovs' not in self.columns:
                raise Exception('Query coverage not selected in outfmt_str. '
                                'Consider reruning the code adding qcovs to '
                                'the blast request')
            self.blasts = self.blasts[self.blasts.qcovs >= self.query_coverage]

    def write(self):
        self.blasts.to_csv(self.out, sep='\t', index=False, header=False)


def main(db: str, query: str, evalue: int = 10, p_id: float = 0,
         max_targets: int = 500, cpus: int = -1, qcovs: Optional[float] = None,
         identify: bool = False, outfile: str = 'hit.hits', outfmt: str = cols,
         unique: bool = False) -> None:
    blast = Blast(db=db, query=query, evalue=evalue, p_id=p_id,
                  max_target_seqs=max_targets, cpus=cpus, query_coverage=qcovs,
                  identify=identify, out=outfile, outfmt_str=outfmt)
    blast.write()


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        prog='justblast', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('query', help='Fasta file with query sequences')
    parser.add_argument('db', help='path to blast database', default='nt')
    parser.add_argument('-e', '--evalue', default=10, type=float,
                        help='evalue for blast search')
    parser.add_argument('-p', '--percent_id', default=0, type=float,
                        help='Minimum percent identity on blast search')
    parser.add_argument('-m', '--max_target_seqs', default=500, type=int,
                        help='Number of aligned sequences to keep')
    parser.add_argument('-q', '--query_coverage', default=None, type=float,
                        help='Minimum query coverage to retain')
    parser.add_argument('-c', '--cpus', default=-1, type=int,
                        help='Number of cpus to use')
    parser.add_argument('-i', '--identify', default=False, action='store_true',
                        help='Whether to use basta to assign taxopnomy to the '
                             'hits based on LCA. This is a rough estimate and '
                             'should be revised carefully')
    parser.add_argument('-o', '--out_filename', default='hit.hits', type=str,
                        help='name of output (filtered) file')
    parser.add_argument('-f', '--outfmt', default=cols, type=str,
                        help='Custom format for BLAST')
    parser.add_argument('-u', '--unique', default=False, action='store_true',
                        help='Process only unique sequences')

    args = parser.parse_args()
    main(db=args.db, query=args.query, evalue=args.evalue,
         p_id=args.percent_id, max_targets=args.max_target_seqs,
         cpus=args.cpus, qcovs=args.query_coverage, identify=args.identify,
         outfile=args.out_filename, outfmt=args.outfmt, unique=args.unique)

