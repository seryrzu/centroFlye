import logging

from utils.bio import read_bio_seqs

logger = logging.getLogger("centroFlye.monomers.monomers")


class MonomerDB:
    def __init__(self, seqs, names2ident, ident2names):
        self.seqs = seqs
        self.names2ident = names2ident
        self.ident2names = ident2names

    @classmethod
    def from_fasta_file(cls, fn):
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
        return self.seqs[index]

    def get_seq_by_name(self, name):
        index = self.names2ident[name]
        return self.seqs[index]

    def get_names(self):
        names_generator = self.names2ident.keys()
        return list(names_generator)

    def get_total_monomers(self):
        return len(self.seqs)
