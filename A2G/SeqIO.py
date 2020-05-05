"""
**SeqIO.py** Module to read alignment files into pandas dataframes
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
from io import StringIO, TextIOWrapper
from textwrap import wrap
from typing import Union, TextIO, Iterator, Tuple, Optional

import numpy as np
import pandas as pd
from joblib import Parallel, delayed
from scipy.stats import entropy


class Alignment(object):
    tipo = 'string'
    reader = StringIO
    full_set = None
    length = None

    def __init__(self, alignment: str, outliers_detection: bool = False,
                 subset_idx: Tuple[Optional[Union[list, np.ndarray]], Optional[
                     Union[list, np.ndarray]]] = (None, None),
                 gap_chr: str = '-', remove_gaps_only: bool = True,
                 cpus: int = -1) -> None:
        self.cpus = cpus
        self.gap = gap_chr
        self.remove_gapped = remove_gaps_only
        self.seq = self.read(alignment)
        self.subset_idx = subset_idx
        self.entropy = outliers_detection

    def __str__(self):
        return self.to_fasta()

    @property
    def subset_idx(self):
        return self.__subset_idx

    @subset_idx.setter
    def subset_idx(self, subset_idx: Tuple[Optional[Union[list, np.ndarray]],
                                           Optional[Union[list, np.ndarray]]]):
        self.__subset_idx = subset_idx
        bools = [True if x is not None else False for x in subset_idx]
        if all(bools):
            # Both aren't None
            if self.full_set is not None:
                self.full_set = self.seq.copy(deep=True)
            self.seq = self.seq.loc[self.subset_idx]
        elif any(bools):
            # at least one is None
            if self.full_set is not None:
                self.full_set = self.seq.copy(deep=True)
            if subset_idx[0] is None:
                self.seq = self.seq.loc[:, self.subset_idx[1]]
            else:
                self.seq = self.seq.loc[self.subset_idx[0], :]

    @property
    def entropy(self):
        return self.__entropy

    @entropy.setter
    def entropy(self, outliers_detection: bool):
        """
        Compute the shannon entropy per column if requested
        :param outliers_detection: Whether or not to compute the entropy
        """
        if outliers_detection:
            self.__entropy = self.seq.apply(self.shannon_entropy, axis=0)
        else:
            self.__entropy = None

    def check_type(self, string) -> None:
        """
        Check if the type of string is a fasta string or its filename
        """
        if '>' not in string:
            self.tipo = 'file'
            self.reader = open

    def yield_sequences(self, fasta_handle: Union[
        StringIO, TextIOWrapper, TextIO]) -> Iterator[Tuple]:
        """
        Form a file or StringIO handle, yield the name and sequences of a fasta
        file. This function will raise an exception if not an alignment (i.e.
        sequences not of the same length)
        :param fasta_handle:  File or StringIO handle
        """
        name, seq = None, ''
        for line in fasta_handle:
            if line.startswith(">"):
                if name:
                    if self.length is not None:
                        try:
                            assert len(seq) == self.length
                        except:
                            print(Exception("Sequences not the same length"))
                    self.length = len(seq)
                    yield name, seq
                name = line[1:].strip()
                seq = ''
            else:
                seq += line.strip()
        if name:
            if self.length is not None:
                try:
                    assert len(seq) == self.length
                except:
                    raise Exception("Sequences not the same length")
            self.length = len(seq)
            yield name, seq

    def read(self, aln: str) -> pd.core.frame.DataFrame:
        """
        Read an alignment into a dataframe
        :param aln: Alignment string or file
        :return: pandas dataframe with alignment
        """
        self.check_type(aln)
        handle = self.reader(aln)
        fasta = Parallel(n_jobs=self.cpus)(delayed(
            pd.DataFrame)([list(seqs)], [names]) for names, seqs in
                                           self.yield_sequences(handle))
        fasta = pd.concat(fasta)
        if self.remove_gapped:
            fasta = fasta.loc[:, ~(fasta == '-').all()]
        return fasta

    @staticmethod
    def shannon_entropy(input_array: np.ndarray):
        """
        Estimate shannon entropy for one column in an alignment
        """
        uniques, n_is = np.unique(input_array, return_counts=True)
        p_is = n_is/input_array.shape[0]
        return entropy(p_is)

    def to_fasta(self):
        fasta = self.seq.apply(lambda x: '>%s\n%s' % (x.name, '\n'.join(wrap(
            ''.join(x.values), 80))), axis=1).values
        return '\n'.join(fasta)