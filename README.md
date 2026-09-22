[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.22885794.svg)](https://doi.org/10.5281/zenodo.22885794)
# Biotech Sequence Suite: Genomic Analysis & HBB Molecular Diagnostic Engine

**Biotech Sequence Suite** is a bioinformatics Python tool designed for DNA sequence processing, translation, and diagnostic analysis. The system is optimized to run in cloud environments such as Google Colab and features an interactive web Graphical User Interface (GUI) powered by Gradio.

---

## 🌟 Key Features

* **Rapid DNA Sanitization:** Non-genomic noise filtering using regular expressions (`regex`) to standardize sequences to upper-case canonical bases.
* **6-Reading Frame Translation:** Automatic generation of primary amino acid structures across all three forward reading frames ($+1, +2, +3$) and three reverse-complement frames ($-1, -2, -3$).
* **$HBB$ Molecular Diagnostic Engine:** Variant detection at the critical Position 6 of the Human Beta-Globin gene, featuring dynamic index offset logic to account for the initiator Methionine ($M$).
* **Web-Based Interactive GUI:** Real-time deployment using Gradio components inside Google Colab for fast sequence input and clinical report visualization.

---

## 🧬 Diagnostic Logic ($HBB$ Gene)

The diagnostic module evaluates Position 6 of the $HBB$ gene (adjusted to Position 7 when the sequence includes the start codon `ATG` / Methionine):

| Phenotype / Condition | Amino Acid at Key Position | Representative Codon | Clinical Classification |
| :--- | :---: | :---: | :--- |
| **HbA (Wild Type)** | Glutamic Acid (`E`) | `GAG` / `GAA` | **Normal:** Fully functional Beta-Globin chain. |
| **HbS (Sickle Cell Anemia)** | Valine (`V`) | `GTG` / `GTT` | **Pathogenic (Glu6Val):** Promotes HbS polymerization and erythrocyte sickling. |
| **HbC (Hemoglobin C Disease)** | Lysine (`K`) | `AAG` / `AAA` | **Pathogenic (Glu6Lys):** Leads to hemoglobin crystallization and mild-to-moderate hemolytic anemia. |

---

## 🛠️ Code Architecture

The project is structured into modular functional blocks:

1. **Genetic Mapping:** Optimized standard genetic code dictionary mapping DNA triplets to amino acid residues.
2. **Processing Engine:**
   * `clean_sequence()`: Filters non-canonical bases using `[^ATGC]`.
   * `translate_dna()`: Translates triplets into primary peptide sequences.
   * `get_reverse_complement()`: Generates reverse-complementary antiparallel DNA strands.
   * `translate_6_frames()`: Computes all $+1, +2, +3, -1, -2, -3$ reading frames in parallel.
3. **Molecular Diagnostic Engine:** `run_hbb_diagnostic()` evaluates specific biomarkers and generates structured diagnostic reports.
4. **GUI Engine:** Interactive web layout deployed seamlessly via Gradio.

---

## 🚀 Execution Guide (Google Colab)

1. Open a new notebook in **Google Colab**.
2. Paste and run the core engine script in **Cell 1**.
3. Paste and run the Gradio interface script in **Cell 2**.
4. Click the generated public URL or interact directly with the web GUI embedded inside the notebook.

---

## 📊 Test Cases

You can test the suite using the following control sequences:

* **HbA Sample (Wild Type):** `ATGGTGCACCTGACTCCTGAGGAGAAGTCTGCCGTTACT`
* **HbS Sample (Sickle Cell):** `ATGGTGCACCTGACTCCTGTGGAGAAGTCTGCCGTTACT`
* **HbC Sample (Hemoglobin C):** `ATGGTGCACCTGACTCCTAAGGAGAAGTCTGCCGTTACT`
