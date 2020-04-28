import logging

from utils.bio import read_bio_seqs

logger = logging.getLogger("centroFlye.monomers.monomers")


class MonomerDB:
    """
    Database of monomers is usually read from fasta file

    Attributes
    ----------
    seqs : list of str
        List of all monomers in the database
    names2ident : dict
        Mapping from names (str) into identificators (indexes of seqs, int)
    ident2names : dict
        Mapping from identificators (indexes of seqs, int) into names (str)
        For each name ident2names[names2ident[name]] == name
        For each index (0 <= index < len(seqs))
            names2ident[index2names[index]] == index
    """

    def __init__(self, seqs, names2ident, ident2names):
        """
        Constuctor for MonomerDB class

        Parameters are same as attributes of the class
        """

        self.seqs = seqs
        self.names2ident = names2ident
        self.ident2names = ident2names

    @classmethod
    def from_fasta_file(cls, fn):
        """
        Class method to construct MonomerDB from a fasta file

        Sequences in MonomerDB will be sorted identically
        to the order in the file

        Parameters
        ----------
        fn : str
            Name of the fasta file with monomers
        """
        logger.info(f'Creating Monomer DataBase from {fn}')
        monomers = read_bio_seqs(fn)
        names2ident = {}
        ident2names = {}
        seqs = []
        for i, (monomer_name, monomer_seq) in enumerate(monomers.items()):
            names2ident[monomer_name] = i
            ident2names[i] = monomer_name
            seqs.append(monomer_seq)
            logger.debug(f'Monomer: ident = {i}, name = {monomer_name}')
            logger.debug(f'         monomer sequence = {monomer_seq}')

        monomer_db = cls(seqs=seqs,
                         names2ident=names2ident,
                         ident2names=ident2names)
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

    def get_seq_by_name(self, name):
        """
        Method to get sequence by name

        Parameters
        ----------
        name : str

        Returns
        -------
        str
        """
        index = self.names2ident[name]
        return self.seqs[index]

    def get_names(self):
        """
        Method to get names of all monomers in the same order as in self.seqs

        Returns
        -------
        list
        """
        names_generator = self.names2ident.keys()
        return list(names_generator)

    def get_total_monomers(self):
        """
        Method to the total number of monomers

        Returns
        -------
        int
        """
        return len(self.seqs)
