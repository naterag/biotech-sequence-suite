import re
import gradio as gr
from Bio.Seq import Seq

# --- PART 1: TRANSLATIONS (i18n) ---
TRANSLATIONS = {
    "English": {
        "error_no_bases": "Error: No valid DNA bases detected (A, T, G, C).",
        "forward_strand": "--- FORWARD STRAND ---",
        "reverse_strand": "--- REVERSE COMPLEMENT STRAND ---",
        "frame": "Frame",
        "length": "Length",
        "gc_content": "GC Content",
        "bp": "bp",
        "report_header": "[MOLECULAR DIAGNOSTIC REPORT: HBB GENE]",
        "diag_no_mutation": "DIAGNOSIS: NO MUTATION (Normal Phenotype - HbA)",
        "diag_no_mutation_detail": "PROTEIN: Wild Type Beta Hemoglobin",
        "detail_glu": "MOLECULAR DETAIL: Glutamic Acid ('E') at Position {pos}.",
        "diag_sickle": "DIAGNOSIS: SICKLE CELL ANEMIA (HbS)",
        "mutation_glu_val": "MUTATION DETECTED: Glu6Val substitution (Valine)",
        "detail_val": "MOLECULAR DETAIL: Pathogenic Valine ('V') at Position {pos}.",
        "diag_hbc": "DIAGNOSIS: HEMOGLOBINOPATHY C (HbC Disease)",
        "mutation_glu_lys": "MUTATION DETECTED: Glu6Lys substitution (Lysine)",
        "detail_lys": "MOLECULAR DETAIL: Pathogenic Lysine ('K') at Position {pos}.",
        "diag_atypical": "DIAGNOSIS: ATYPICAL VARIANT ('{aa}')",
        "finding_atypical": "FINDING: Non-standard amino acid detected at Position {pos}.",
        "warning_clinical": "WARNING: Clinical review is required.",
        "error_empty": "ERROR: Empty sequence.",
        "error_short": "CRITICAL ERROR: Sequence is too short to analyze the HBB gene."
    },
    "Español": {
        "error_no_bases": "Error: No se detectaron bases válidas (A, T, G, C).",
        "forward_strand": "--- CADENA DIRECTA ---",
        "reverse_strand": "--- CADENA REVERSA COMPLEMENTARIA ---",
        "frame": "Marco",
        "length": "Longitud",
        "gc_content": "Contenido GC",
        "bp": "bp",
        "report_header": "[REPORTE DE DIAGNÓSTICO MOLECULAR: GEN HBB]",
        "diag_no_mutation": "DIAGNÓSTICO: SIN MUTACIÓN (Fenotipo Normal - HbA)",
        "diag_no_mutation_detail": "PROTEÍNA: Hemoglobina Beta Salvaje (Wild Type)",
        "detail_glu": "DETALLE MOLECULAR: Ácido Glutámico ('E') en Posición {pos}.",
        "diag_sickle": "DIAGNÓSTICO: ANEMIA FALCIFORME (HbS)",
        "mutation_glu_val": "MUTACIÓN DETECTADA: Cambio Glu6Val (Sustitución por Valina)",
        "detail_val": "DETALLE MOLECULAR: Valina ('V') patogénica en Posición {pos}.",
        "diag_hbc": "DIAGNÓSTICO: HEMOGLOBINOPATÍA C (Enfermedad por HbC)",
        "mutation_glu_lys": "MUTACIÓN DETECTADA: Cambio Glu6Lys (Sustitución por Lisina)",
        "detail_lys": "DETALLE MOLECULAR: Lisina ('K') patogénica en Posición {pos}.",
        "diag_atypical": "DIAGNÓSTICO: VARIANTE ATÍPICA ('{aa}')",
        "finding_atypical": "HALLAZGO: Aminoácido no estándar detectado en Posición {pos}.",
        "warning_clinical": "ADVERTENCIA: Se requiere revisión clínica.",
        "error_empty": "ERROR: Secuencia vacía.",
        "error_short": "ERROR CRÍTICO: La secuencia es demasiado corta para analizar el gen HBB."
    }
}

def t(key, lang, **kwargs):
    """Helper function to get a translated string."""
    text = TRANSLATIONS[lang].get(key, key)
    return text.format(**kwargs) if kwargs else text


# --- PART 2: PROCESSING ENGINE ---

def clean_sequence(raw_data: str) -> str:
    """Filters genomic noise and standardizes input to uppercase."""
    return re.sub(r'[^ATGC]', '', raw_data.upper())

def translate_6_frames(dna_seq: str) -> dict:
    """Translates DNA into all 6 reading frames using Biopython."""
    seq = Seq(dna_seq)
    frames = {}
    frames['+1'] = str(seq.translate())
    frames['+2'] = str(seq[1:].translate())
    frames['+3'] = str(seq[2:].translate())
    rev_seq = seq.reverse_complement()
    frames['-1'] = str(rev_seq.translate())
    frames['-2'] = str(rev_seq[1:].translate())
    frames['-3'] = str(rev_seq[2:].translate())
    return frames


# --- PART 3: MOLECULAR DIAGNOSTIC ENGINE (HBB) ---

def run_hbb_diagnostic(protein_seq: str, lang: str) -> str:
    """Molecular biomarker detector for the HBB gene (bilingual)."""
    if not protein_seq:
        return t("error_empty", lang)

    target_idx = 6 if protein_seq.startswith('M') else 5

    if len(protein_seq) > target_idx:
        target_aa = protein_seq[target_idx]
        pos_num = target_idx + 1

        report = ["\n" + "*" * 50, f" {t('report_header', lang)} ", "*" * 50]

        if target_aa == 'E':
            report.append(t("diag_no_mutation", lang))
            report.append(t("diag_no_mutation_detail", lang))
            report.append(t("detail_glu", lang, pos=pos_num))
        elif target_aa == 'V':
            report.append(t("diag_sickle", lang))
            report.append(t("mutation_glu_val", lang))
            report.append(t("detail_val", lang, pos=pos_num))
        elif target_aa == 'K':
            report.append(t("diag_hbc", lang))
            report.append(t("mutation_glu_lys", lang))
            report.append(t("detail_lys", lang, pos=pos_num))
        else:
            report.append(t("diag_atypical", lang, aa=target_aa))
            report.append(t("finding_atypical", lang, pos=pos_num))
            report.append(t("warning_clinical", lang))

        report.append("*" * 50)
        return "\n".join(report)
    else:
        return t("error_short", lang)


# --- PART 4: GRADIO INTERFACE ---

def analyze_dna_gui(dna_sequence, language):
    """Connects the backend logic with the Gradio interface (bilingual)."""
    lang = language
    clean_dna = clean_sequence(dna_sequence)
    if not clean_dna:
        return t("error_no_bases", lang), "", ""

    frames = translate_6_frames(clean_dna)

    frames_text = (
        f"{t('forward_strand', lang)}\n"
        f"{t('frame', lang)} +1: {frames['+1']}\n"
        f"{t('frame', lang)} +2: {frames['+2']}\n"
        f"{t('frame', lang)} +3: {frames['+3']}\n\n"
        f"{t('reverse_strand', lang)}\n"
        f"{t('frame', lang)} -1: {frames['-1']}\n"
        f"{t('frame', lang)} -2: {frames['-2']}\n"
        f"{t('frame', lang)} -3: {frames['-3']}"
    )

    gc_val = (clean_dna.count("G") + clean_dna.count("C")) / len(clean_dna) * 100
    metrics = f"{t('length', lang)}: {len(clean_dna)} {t('bp', lang)} | {t('gc_content', lang)}: {gc_val:.2f}%\n"
    diagnosis = run_hbb_diagnostic(frames['+1'], lang)

    return clean_dna, frames_text, f"{metrics}\n{diagnosis}"


demo = gr.Interface(
    fn=analyze_dna_gui,
    inputs=[
        gr.Textbox(
            lines=4,
            placeholder="Paste your DNA sequence here / Pega tu secuencia aquí...",
            label="DNA Input / Entrada ADN"
        ),
        gr.Radio(
            choices=["English", "Español"],
            value="English",
            label="Language / Idioma"
        )
    ],
    outputs=[
        gr.Textbox(label="Clean Sequence / Secuencia Limpia"),
        gr.Textbox(label="6-Frame Translation / Traducción en 6 Marcos", lines=8),
        gr.Textbox(label="Metrics & Diagnosis / Métricas y Diagnóstico", lines=8)
    ],
    title="Biotech Sequence Suite - 6 Frames & Diagnostic",
    description="Advanced genomic analyzer with bilingual support (EN/ES)."
)

if __name__ == "__main__":
    demo.launch()
