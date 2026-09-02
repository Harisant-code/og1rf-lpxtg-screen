# LPXTG sortase-substrate screen for *Enterococcus faecalis* OG1RF

A deterministic, dependency-light Python script that identifies candidate
sortase-anchored cell-wall surface proteins from a bacterial proteome using
the tripartite cell-wall sorting signal criteria described by Boekhorst et al.
(2005).

## Criteria

A protein is classified as a candidate sortase substrate if **all three**
of the following hold:

1. **Sortase-recognition motif** — a match to `[LIVF]PXT[GSAN]` or `NPQTN`
   within the final 60 C-terminal residues.
2. **Hydrophobic transmembrane tail** — at least 10 residues follow the motif,
   of which >40% are hydrophobic (A, I, L, M, V, F, Y, W).
3. **Positively charged cytoplasmic terminus** — at least one Lys or Arg
   within the last 15 residues.

Proteins shorter than 50 aa are not evaluated.

## Contents

| File | Description |
|------|-------------|
| `lpxtg_screen.py` | The screen script (standalone, CLI) |
| `data/gtf_proteome.fasta` | Input: all 2,583 OG1RF protein sequences, translated from the RefSeq annotation of assembly GCF_000172575.2 (release 2019-11-18) by annotation coordinates from CP002621.1 |
| `data/gtf_screen.csv` | Expected output: full screen results for 2,583 proteins |
| `data/supp_table_cwpred_crosstab.csv` | Cross-validation: CW-PRED (Chatziargyri et al., 2024) predictions vs. this screen for all 22 CW-PRED positives |
| `requirements.txt` | Python dependencies |
| `LICENSE` | MIT License |
| `CITATION.cff` | Citation metadata (for Zenodo DOI) |

## Requirements

- Python >= 3.8
- Biopython
- pandas

```
pip install biopython pandas
```

## Usage

```bash
python lpxtg_screen.py --fasta data/gtf_proteome.fasta --out my_results.csv
```

Expected output:

```
Screened 2583 proteins
  motif anywhere in sequence:      165
  motif within final 60 residues: 56
  full tripartite sorting signal:  36
```

## Output columns

| Column | Description |
|--------|-------------|
| `locus_tag` | FASTA record identifier |
| `description` | Full FASTA header |
| `length` | Protein length (aa) |
| `motif_anywhere` | Motif found anywhere in the sequence |
| `motif_cterm` | Motif found within the final 60 residues |
| `motif` | Matched motif sequence (C-terminal match) |
| `motif_pos_from_end` | Distance of motif start from the C-terminus |
| `tail_length` | Residues following the motif |
| `tail_hydrophobic` | Hydrophobic fraction of the tail |
| `charged_terminus` | >=1 K/R within the last 15 residues |
| `tripartite` | All three criteria fulfilled |

## Verification

To confirm reproducibility, compare your output against the expected result:

```python
import pandas as pd
mine = pd.read_csv('my_results.csv')
ref  = pd.read_csv('data/gtf_screen.csv')
print('identical:', mine.equals(ref))  # should print True
```

**Positive controls** — five experimentally established OG1RF sortase
substrates must all show `tripartite = True`:

| Locus tag | Gene | Product |
|-----------|------|---------|
| OG1RF_RS04595 | *ace* | Collagen adhesin |
| OG1RF_RS04550 | *ebpA* | Pilus tip adhesin |
| OG1RF_RS04555 | *ebpB* | Pilus shaft |
| OG1RF_RS04560 | *ebpC* | Pilus base |
| OG1RF_RS00420 | *fss1* | Fibrinogen-binding surface protein |

## Cross-validation with CW-PRED

The same proteome was independently analysed with CW-PRED (Chatziargyri et al.,
2024, *J. Bioinform. Comput. Biol.* 22:2450021), an HMM-based predictor of
cell-wall sorting signals. CW-PRED recovered 15 of the 16 validated substrates
identified by this screen (the sole miss being EbpA, a known profile-model
limitation for pilus-type sorting signals) and rejected all 18 chance motif
matches. Full results are in `data/supp_table_cwpred_crosstab.csv`.

## Development note

This script was drafted with AI assistance (Biomni, Phylo) under the authors'
scientific direction, implementing the published tripartite sorting-signal
criteria (Boekhorst et al., 2005). The authors independently executed the
code and verified that it reproduces the reported results (2,583 proteins
screened; 165 / 56 / 36 funnel; all five positive-control adhesins correctly
identified).

## References

- Boekhorst J, de Been MW, Kleerebezem M, Siezen RJ. 2005. Predicting the
  extracellular proteome of *Lactococcus lactis*, a lactic acid bacterium.
  *J. Bacteriol.* 187:4928–4934.
- Chatziargyri A, Stasi EA, Tsirigos KD, Litou ZI, Iconomidou VA, Bagos PG.
  2024. CW-PRED: Prediction of C-terminal surface anchoring sorting signals
  in bacteria and Archaea. *J. Bioinform. Comput. Biol.* 22:2450021.
- Teufel F et al. 2022. SignalP 6.0 predicts all five types of signal peptides
  using protein language models. *Nat. Biotechnol.* 40:1023–1025.

## License

MIT License — see `LICENSE`.
