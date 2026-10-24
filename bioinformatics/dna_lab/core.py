"""FASTA parsing and DNA analysis using the standard genetic code."""

from collections import Counter
from dataclasses import dataclass
from pathlib import Path


@dataclass(frozen=True)
class Sequence:
    name: str
    dna: str


def read_fasta(path: str | Path) -> list[Sequence]:
    records = []
    name = None
    parts = []

    def append_record():
        dna = "".join(parts).upper()
        if not dna:
            raise ValueError(f"Empty sequence: {name}")
        invalid = set(dna) - set("ACGTRYSWKMBDHVN")
        if invalid:
            raise ValueError(f"Invalid DNA symbols in {name}: {', '.join(sorted(invalid))}")
        records.append(Sequence(name, dna))

    with open(path, encoding="utf-8") as handle:
        for number, raw in enumerate(handle, 1):
            line = raw.strip()
            if not line:
                continue
            if line.startswith(">"):
                if name is not None:
                    append_record()
                name = line[1:].strip()
                if not name:
                    raise ValueError(f"Missing FASTA header at line {number}")
                parts = []
            else:
                if name is None:
                    raise ValueError(f"Sequence before FASTA header at line {number}")
                parts.append("".join(line.split()))
    if name is not None:
        append_record()
    if not records:
        raise ValueError("FASTA file contains no sequences")
    return records


def reverse_complement(dna: str) -> str:
    return dna.translate(str.maketrans("ACGTRYSWKMBDHVN", "TGCAYRSWMKVHDBN"))[::-1]


BASES = "TCAG"
AMINO_ACIDS = (
    "FFLLSSSSYY**CC*W"
    "LLLLPPPPHHQQRRRR"
    "IIIMTTTTNNKKSSRR"
    "VVVVAAAADDEEGGGG"
)
CODONS = dict(zip((a + b + c for a in BASES for b in BASES for c in BASES), AMINO_ACIDS))


def translate(dna: str, frame: int = 0) -> str:
    if frame not in (0, 1, 2):
        raise ValueError("Frame must be 0, 1 or 2")
    return "".join(CODONS.get(dna[i:i + 3], "X") for i in range(frame, len(dna) - 2, 3))


def analyze(record: Sequence, motif: str | None = None, frame: int = 0) -> dict:
    counts = Counter(record.dna)
    canonical = sum(counts[base] for base in "ACGT")
    result = {
        "name": record.name,
        "length": len(record.dna),
        "base_counts": {base: counts[base] for base in "ACGT"},
        "ambiguous_bases": len(record.dna) - canonical,
        "gc_percent": round(100 * (counts["G"] + counts["C"]) / canonical, 2) if canonical else None,
        "reverse_complement": reverse_complement(record.dna),
        "protein": translate(record.dna, frame),
        "translation_frame": frame,
    }
    if motif is not None:
        motif = motif.upper()
        if not motif or set(motif) - set("ACGTRYSWKMBDHVN"):
            raise ValueError("Motif must be a nonempty DNA string")
        result["motif"] = motif
        result["motif_positions"] = [
            i + 1 for i in range(len(record.dna) - len(motif) + 1)
            if record.dna.startswith(motif, i)
        ]
    return result
