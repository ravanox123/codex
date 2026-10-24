import argparse
import json
from pathlib import Path

from .core import analyze, read_fasta


def main():
    parser = argparse.ArgumentParser(description="Analyze DNA sequences from a FASTA file")
    parser.add_argument("fasta", type=Path, help="Input DNA FASTA file")
    parser.add_argument("--motif", help="Find exact motif matches (positions start at 1)")
    parser.add_argument("--frame", type=int, choices=(0, 1, 2), default=0)
    parser.add_argument("--output", type=Path, help="Write JSON report to this file")
    args = parser.parse_args()
    try:
        report = [analyze(record, args.motif, args.frame) for record in read_fasta(args.fasta)]
        serialized = json.dumps(report, ensure_ascii=False, indent=2) + "\n"
        if args.output:
            args.output.write_text(serialized, encoding="utf-8")
        else:
            print(serialized, end="")
    except (OSError, ValueError) as error:
        parser.exit(2, f"Error: {error}\n")


if __name__ == "__main__":
    main()
