# (c) 2020 by Authors
# This file is a part of centroFlye program.
# Released under the BSD license (see LICENSE file)

import logging

import numpy as np

from utils.bio import read_bio_seqs
from utils.os_utils import expandpath

logger = logging.getLogger("centroFlye.monomers.monomer_db")


class MonomerDB:
    """
    Database of monomers. It is usually read from a fasta file

    Attributes
    ----------
    seqs : list of str
        List of all monomers in the database
    id2index : dict
        Mapping from ids (str) into indexes (int)
    index2id : dict
        Mapping from indexificators (indexes of seqs, int) into ids (str)
        For each id index2id[id2index[id]] == id
        For each index (0 <= index < len(seqs))
            id2index[index2id[index]] == index
    """

    def __init__(self, seqs, id2index, index2id):
        """
        Constuctor for MonomerDB class

        Parameters are same as attributes of the class
        """

        self.seqs = seqs
        self.id2index = id2index
        self.index2id = index2id

    @classmethod
    def from_fasta_file(cls, fn):
        """
        Class method to construct MonomerDB from a fasta file

        Sequences in MonomerDB will be sorted by indexes
        to the order in the file

        Parameters
        ----------
        fn : str
            Name of the fasta file with monomers
        """
        fn = expandpath(fn)
        logger.info(f'Creating Monomer DataBase from {fn}')
        monomers = read_bio_seqs(fn)
        id2index = {}
        index2id = {}
        seqs = []
        for i, (monomer_id, monomer_seq) in enumerate(monomers.items()):
            id2index[monomer_id] = i
            index2id[i] = monomer_id
            seqs.append(monomer_seq)
            logger.debug(f'Monomer: index = {i}, id = {monomer_id}')
            logger.debug(f'         monomer sequence = {monomer_seq}')

        monomer_db = cls(seqs=seqs,
                         id2index=id2index,
                         index2id=index2id)
        logger.info(f'Finished Creating Monomer DataBase')
        return monomer_db

    def get_seq_by_index(self, index):
        """
        Method to get sequence by index

        Parameters
        ----------
        index : int

        Returns
        -------
        str
        """
        assert 0 <= index < len(self.seqs)
        return self.seqs[index]

    def get_seq_by_id(self, id):
        """
        Method to get sequence by id

        Parameters
        ----------
        id : str

        Returns
        -------
        str
        """
        index = self.id2index[id]
        return self.seqs[index]

    def get_ids(self):
        """
        Method to get ids of all monomers in the same order as in self.seqs

        Returns
        -------
        list
        """
        ids_generator = self.id2index.keys()
        return list(ids_generator)

    def get_total_monomers(self):
        """
        Method to the total number of monomers

        Returns
        -------
        int
        """
        return len(self.seqs)

    def get_mean_monomer_len(self):
        """
        Method to return average (mean) length of monomers

        Returns
        -------
        float
        """
        sequences_lens = [len(seq) for seq in self.seqs]
        mean_sequences_len = np.mean(sequences_lens)
        return mean_sequences_len
