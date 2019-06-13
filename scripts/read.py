#(c) 2019 by Authors
#This file is a part of centroFlye program.
#Released under the BSD license (see LICENSE file)

class Read:
    def __init__(self, id, seq=None, simulated=False):
        self.id, self.seq = id, seq
        # Parsing SimLoRD output format
        if simulated:
            spl_id = self.id.split('_')
            self.numb = int(spl_id[1])
            self.length = int(spl_id[2].split('=')[1][:-2])
            self.start_pos = int(spl_id[3].split('=')[1])
            self.n_errors = int(spl_id[6].split('=')[1])
            self.error_rate = float(spl_id[9].split('=')[1])
            self.mult = float(spl_id[-1].split('=')[1])

    @classmethod
    def FromBiopyRead(cls, biopy_read, simulated=False):
        id, seq = biopy_read.id, str(biopy_read.seq)
        return cls(id, seq, simulated)

    def __len__(self):
        return len(self.seq)

    def __getitem__(self, key):
        return self.seq[key]
