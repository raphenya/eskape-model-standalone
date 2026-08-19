"""Cleans a dataset by removing molecules which cannot be parsed by RDKit."""

import csv
from rdkit import Chem
from rdkit import RDLogger
from tap import Tap

# Disable all RDKit logging messages globally
RDLogger.DisableLog('rdApp.*')


class Args(Tap):
    data_path: str  # Data CSV to sanitize
    save_path: str  # Path to CSV where sanitized data will be saved
    error_path: str  # Path to CSV where errors will be saved


def sanitize(data_path: str, save_path: str, error_path: str):
    """

    Sanitizes a CSV file by removing molecules which cannot be parsed by RDKit i.e. invalid SMILES strings. The sanitized data is saved to a new CSV file, and any errors encountered during the sanitization process are logged to a separate CSV file.

    """
    error_messages = []

    with open(data_path) as f:
        reader = csv.reader(f)
        header = next(reader)

        lines = []
        error_lines = []

        for line in reader:
            if line[0] == '':
                continue

            mol = Chem.MolFromSmiles(line[0])
            if mol is not None:
                lines.append(line)
            else:
                error_lines.append(line)
                error_messages.append(f"Error parsing SMILES: {line[0]}")

    with open(save_path, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(header)
        for line in lines:
            writer.writerow(line)

    with open(error_path, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(header)
        for error_line in error_lines:
            writer.writerow(error_line)

    return len(error_lines)


if __name__ == '__main__':
    args = Args().parse_args()

    sanitize(args.data_path, args.save_path)
