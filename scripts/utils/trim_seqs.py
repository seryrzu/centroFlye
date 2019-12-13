import argparse

from bio import read_bio_seqs, write_bio_seqs


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input", required=True)
    parser.add_argument("-o", "--output", required=True)
    parser.add_argument("--freq", help='Freq of length to remove from each side (2*freq removed total)',
                        type=float,
                        required=True)
    params = parser.parse_args()

    if not (0 <= params.freq <= 0.5):
        print(f'Freq must be in [0;0.5]')
        return

    seqs = read_bio_seqs(params.input)
    seqs = {s_id: seq[int(params.freq*len(seq)):int((1-params.freq)*len(seq))]
            for s_id, seq in seqs.items()}
    write_bio_seqs(params.output, seqs)


if __name__ == "__main__":
    main()
