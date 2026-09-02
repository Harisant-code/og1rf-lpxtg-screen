#!/usr/bin/env python3
"""
lpxtg_screen.py — Genome-wide screen for sortase-anchored (LPXTG-type)
cell-wall surface proteins in a bacterial proteome.

Implements the tripartite cell-wall sorting signal criteria described by
Boekhorst et al. (2005) J. Bacteriol. 187:4928, as applied to
Enterococcus faecalis OG1RF in [this study].

A protein is classified as a candidate sortase substrate if ALL of the
following hold:

  1. Sortase-recognition motif: a match to [LIVF]PXT[GSAN] or NPQTN
     located within the final 60 C-terminal residues.
  2. Hydrophobic transmembrane tail: at least 10 residues follow the
     motif, of which >40% are hydrophobic (A, I, L, M, V, F, Y, W).
  3. Positively charged cytoplasmic terminus: at least one Lys or Arg
     within the last 15 residues.

Proteins shorter than 50 aa are not evaluated (insufficient C-terminal
context).

Usage:
    python lpxtg_screen.py --fasta proteome.fasta --out screen_results.csv

Input:
    --fasta   Protein FASTA file (one record per annotated protein).
    --out     Output CSV path (default: lpxtg_screen_results.csv).

Output CSV columns:
    locus_tag          FASTA record identifier (first word of header)
    description        Full FASTA header
    length             Protein length (aa)
    motif_anywhere     Motif found anywhere in the sequence (True/False)
    motif_cterm        Motif found within the final 60 residues (True/False)
    motif              Matched motif sequence (C-terminal match only)
    motif_pos_from_end Distance of motif start from the C-terminus (aa)
    tail_length        Residues following the motif
    tail_hydrophobic   Hydrophobic fraction of the tail
    charged_terminus   >=1 K/R within the last 15 residues (True/False)
    tripartite         All three criteria fulfilled (True/False)

Requirements: Python >= 3.8, Biopython, pandas
    (install: pip install biopython pandas)

The screen is deterministic; no random seeds or external services are used.
"""

import argparse
import re

import pandas as pd
from Bio import SeqIO

# Sortase-recognition motifs: canonical LPXTG-family ([LIVF]PXT[GSAN],
# e.g. LPKTG in Ace, LPETG in EbpA) plus the NPQTN variant.
SORTASE_MOTIF = re.compile(r'([LIVF]P[A-Z]T[GSAN]|NPQTN)')

CTERM_WINDOW = 60        # motif must lie within this many residues of the C-terminus
MIN_TAIL = 10            # minimum residues following the motif
HYDROPHOBIC_AA = set('AILMVFYW')
MIN_HYDROPHOBIC_FRAC = 0.4
CHARGED_WINDOW = 15      # window at the extreme C-terminus checked for K/R
MIN_PROTEIN_LEN = 50


def screen_protein(seq):
    """Apply the tripartite sorting-signal criteria to one protein sequence."""
    seq = str(seq).upper().replace('*', '')
    result = dict(
        length=len(seq),
        motif_anywhere=bool(SORTASE_MOTIF.search(seq)),
        motif_cterm=False,
        motif=None,
        motif_pos_from_end=None,
        tail_length=None,
        tail_hydrophobic=None,
        charged_terminus=False,
        tripartite=False,
    )
    if len(seq) < MIN_PROTEIN_LEN:
        return result

    cterm = seq[-CTERM_WINDOW:]
    m = SORTASE_MOTIF.search(cterm)
    if not m:
        return result

    motif = m.group()
    after = cterm[m.end():]
    result.update(
        motif_cterm=True,
        motif=motif,
        motif_pos_from_end=len(seq) - (len(seq) - CTERM_WINDOW + m.start()),
        tail_length=len(after),
    )
    if len(after) < MIN_TAIL:
        return result

    hydrophobic_frac = sum(aa in HYDROPHOBIC_AA for aa in after) / len(after)
    charged = any(aa in 'KR' for aa in seq[-CHARGED_WINDOW:])
    result.update(
        tail_hydrophobic=round(hydrophobic_frac, 3),
        charged_terminus=charged,
        tripartite=(hydrophobic_frac > MIN_HYDROPHOBIC_FRAC and charged),
    )
    return result


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--fasta', required=True, help='Protein FASTA file')
    ap.add_argument('--out', default='lpxtg_screen_results.csv', help='Output CSV')
    args = ap.parse_args()

    records = list(SeqIO.parse(args.fasta, 'fasta'))
    rows = []
    for rec in records:
        res = screen_protein(rec.seq)
        res['locus_tag'] = rec.id
        res['description'] = rec.description
        rows.append(res)

    df = pd.DataFrame(rows)[
        ['locus_tag', 'description', 'length', 'motif_anywhere', 'motif_cterm',
         'motif', 'motif_pos_from_end', 'tail_length', 'tail_hydrophobic',
         'charged_terminus', 'tripartite']
    ]
    df.to_csv(args.out, index=False)

    n_any = int(df.motif_anywhere.sum())
    n_cterm = int(df.motif_cterm.sum())
    n_tri = int(df.tripartite.sum())
    print(f'Screened {len(df)} proteins')
    print(f'  motif anywhere in sequence:      {n_any}')
    print(f'  motif within final {CTERM_WINDOW} residues: {n_cterm}')
    print(f'  full tripartite sorting signal:  {n_tri}')
    print(f'Results written to {args.out}')


if __name__ == '__main__':
    main()
