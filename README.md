[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.22885794.svg)](https://doi.org/10.5281/zenodo.22885794)

# Biotech Sequence Suite: Genomic Analysis & HBB Molecular Diagnostic Engine

**Biotech Sequence Suite** is a Python bioinformatics tool for DNA sequence processing, 6-frame translation and detection of variants at residue 6 of the human beta-globin gene (*HBB*): HbA (normal), HbS (sickle cell) and HbC. It uses [Biopython](https://biopython.org/) for translation and ships with an interactive web interface built with [Gradio](https://www.gradio.app/), available in **English and Spanish**.

> **Disclaimer:** this is an educational tool. It is **not** a clinical diagnostic and must not be used to make medical decisions.

---

## Key Features

* **DNA sanitization:** removes non-genomic noise (spaces, numbers, line breaks, other characters) with regular expressions and standardizes the input to upper-case `A`, `T`, `G`, `C`.
* **6-reading-frame translation:** translates the three forward frames (+1, +2, +3) and the three reverse-complement frames (-1, -2, -3) using Biopython. Sequences are trimmed to a multiple of 3 before translation.
* **HBB variant detection (residue 6):** searches the six frames for the start of beta-globin (motif `VHLTP`, residues 1-5) and reads residue 6. Because it scans all frames, it works when the gene is preceded by extra sequence, is not in frame +1, or is given as the reverse strand.
* **Bilingual interface (EN/ES):** the language of the report is selected from the interface.
* **Sequence metrics:** length (bp) and GC content.

---

## Diagnostic Logic (*HBB* gene)

The engine locates the motif `VHLTP` (residues 1-5 of the mature beta-globin chain) and evaluates the next residue, residue 6:

| Result | Residue 6 | Representative codons | Variant (classic / HGVS) |
| :--- | :---: | :---: | :--- |
| **HbA (wild type)** | Glutamic acid (`E`) | `GAG` / `GAA` | No variant detected at residue 6 |
| **HbS (sickle hemoglobin)** | Valine (`V`) | `GTG` / `GTT` | Glu6Val / p.Glu7Val |
| **HbC (hemoglobin C)** | Lysine (`K`) | `AAG` / `AAA` | Glu6Lys / p.Glu7Lys |
| **Atypical** | Any other | - | Flagged for review |

Numbering note: the classic nomenclature (Glu6Val) does not count the initiator methionine, while HGVS does (p.Glu7Val). Both are shown in the report.

**Limitations:** only residue 6 is evaluated, so other *HBB* variants are not covered. If the `VHLTP` motif is not found in any frame, the tool reports that residue 6 cannot be evaluated.

---

## Installation and Usage

Requires Python 3.9 or newer.

```bash
git clone https://github.com/naterag/biotech-sequence-suite.git
cd biotech-sequence-suite
pip install -r requirements.txt
python app.py
```

Open the local URL that Gradio prints in the terminal, paste a DNA sequence, choose a language and read the report.

### Google Colab

The notebook `biotech_sequence_suite.ipynb` contains the original version of the project (without Biopython), designed to run in Google Colab and share a temporary public link. The up-to-date version of the code is `app.py`.

---

## Test Cases

You can test the tool with these control sequences:

* **HbA (wild type):** `ATGGTGCACCTGACTCCTGAGGAGAAGTCTGCCGTTACT`
* **HbS (sickle cell):** `ATGGTGCACCTGACTCCTGTGGAGAAGTCTGCCGTTACT`
* **HbC (hemoglobin C):** `ATGGTGCACCTGACTCCTAAGGAGAAGTCTGCCGTTACT`

### Automated tests

The same sequences (forward, with extra leading sequence and as reverse complement) are covered by automated tests:

```bash
pip install -r requirements-dev.txt
pytest
```

---

## Code Architecture

`app.py` is organized in four blocks:

1. **Translations (i18n):** `TRANSLATIONS` dictionary with every message in English and Spanish, and the helper `t()`.
2. **Processing engine:**
   * `clean_sequence()`: filters non-canonical bases with `[^ATGC]`.
   * `translate_6_frames()`: computes the six reading frames with Biopython.
3. **Diagnostic engine:**
   * `find_hbb_residue6()`: finds the `VHLTP` motif in any frame and returns residue 6.
   * `run_hbb_diagnostic()`: builds the structured report.
4. **GUI:** `analyze_dna_gui()` connects the engine with the Gradio interface.

Tests live in `tests/test_app.py`.

---

## Citation

If you use this software, please cite it using the DOI above or the metadata in `CITATION.cff`.

## License

Released under the [MIT License](LICENSE).
