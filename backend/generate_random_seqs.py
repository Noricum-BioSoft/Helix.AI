import random
import string
import pandas as pd


def generate_random_dna_csv(csv_path: str, num_sequences):
    """
    Generates a CSV file with random DNA sequences.

    Parameters:
    - csv_path: Path to save the CSV.
    - num_sequences: Number of sequences to generate.
    """

    def random_id(length=8):
        return ''.join(random.choices(string.ascii_uppercase + string.digits, k=length))

    def random_dna_sequence(length):
        return ''.join(random.choices("ACGT", k=length))

    def random_od600():
        import random
        random_number = random.random()
        return random_number

    data = []
    for _ in range(num_sequences):
        seq_id = random_id()
        seq_len = random.randint(100, 200)
        sequence = random_dna_sequence(seq_len)
        od_value = random_od600()
        data.append({"sequence_id": seq_id, "sequence": sequence, "OD600": od_value})

    df = pd.DataFrame(data)
    df.to_csv(csv_path, index=False)
    print(f"CSV with {num_sequences} sequences saved to {csv_path}")


if __name__ ==  "__main__":
    generate_random_dna_csv("../tests/data/random_sequences.csv", num_sequences=20)