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
        "found_in_frame": "Motif VHLTP (beta-globin residues 1-5) found in frame {frame}.",
        "diag_no_variant": "RESULT: No variant detected at residue 6 (Glu6, consistent with HbA)",
        "diag_no_variant_detail": "NOTE: Only residue 6 was evaluated; other HBB variants are not covered.",
        "detail_glu": "MOLECULAR DETAIL: Glutamic Acid ('E') at residue 6 (HGVS: p.Glu7).",
        "diag_sickle": "RESULT: Variant consistent with sickle hemoglobin (HbS)",
        "mutation_glu_val": "VARIANT DETECTED: Glu6Val (HGVS: p.Glu7Val)",
        "detail_val": "MOLECULAR DETAIL: Valine ('V') at residue 6 instead of Glutamic Acid.",
        "diag_hbc": "RESULT: Variant consistent with hemoglobin C (HbC)",
        "mutation_glu_lys": "VARIANT DETECTED: Glu6Lys (HGVS: p.Glu7Lys)",
        "detail_lys": "MOLECULAR DETAIL: Lysine ('K') at residue 6 instead of Glutamic Acid.",
        "diag_atypical": "RESULT: Atypical variant at residue 6 ('{aa}')",
        "warning_clinical": "WARNING: Clinical review is required.",
        "not_found": "The VHLTP motif (beta-globin residues 1-5) was not found in any of the 6 frames; residue 6 cannot be evaluated.",
        "disclaimer": "NOTE: Educational tool. It is not a clinical diagnostic.",
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
        "found_in_frame": "Motivo VHLTP (residuos 1-5 de la beta-globina) encontrado en el marco {frame}.",
        "diag_no_variant": "RESULTADO: No se detectó variante en el residuo 6 (Glu6, consistente con HbA)",
        "diag_no_variant_detail": "NOTA: Solo se evaluó el residuo 6; otras variantes de HBB no están cubiertas.",
        "detail_glu": "DETALLE MOLECULAR: Ácido Glutámico ('E') en el residuo 6 (HGVS: p.Glu7).",
        "diag_sickle": "RESULTADO: Variante consistente con hemoglobina S (HbS)",
        "mutation_glu_val": "VARIANTE DETECTADA: Glu6Val (HGVS: p.Glu7Val)",
        "detail_val": "DETALLE MOLECULAR: Valina ('V') en el residuo 6 en lugar de Ácido Glutámico.",
        "diag_hbc": "RESULTADO: Variante consistente con hemoglobina C (HbC)",
        "mutation_glu_lys": "VARIANTE DETECTADA: Glu6Lys (HGVS: p.Glu7Lys)",
        "detail_lys": "DETALLE MOLECULAR: Lisina ('K') en el residuo 6 en lugar de Ácido Glutámico.",
        "diag_atypical": "RESULTADO: Variante atípica en el residuo 6 ('{aa}')",
        "warning_clinical": "ADVERTENCIA: Se requiere revisión clínica.",
        "not_found": "No se encontró el motivo VHLTP (residuos 1-5 de la beta-globina) en ninguno de los 6 marcos; no se puede evaluar el residuo 6.",
        "disclaimer": "NOTA: Herramienta educativa. No es un diagnóstico clínico.",
    },
}


def t(key, lang, **kwargs):
    """Helper function to get a translated string."""
    text = TRANSLATIONS[lang].get(key, key)
    return text.format(**kwargs) if kwargs else text


# --- PART 2: PROCESSING ENGINE ---

def clean_sequence(raw_data: str) -> str:
    """Filters genomic noise and standardizes input to uppercase."""
    return re.sub(r'[^ATGC]', '', raw_data.upper())


def _translate(seq: Seq) -> str:
    """Translates a sequence, trimming it to a multiple of 3 to avoid Biopython warnings."""
    usable = (len(seq) // 3) * 3
    return str(seq[:usable].translate())


def translate_6_frames(dna_seq: str) -> dict:
    """Translates DNA into all 6 reading frames using Biopython."""
    seq = Seq(dna_seq)
    rev_seq = seq.reverse_complement()
    return {
        '+1': _translate(seq),
        '+2': _translate(seq[1:]),
        '+3': _translate(seq[2:]),
        '-1': _translate(rev_seq),
        '-2': _translate(rev_seq[1:]),
        '-3': _translate(rev_seq[2:]),
    }


# --- PART 3: MOLECULAR DIAGNOSTIC ENGINE (HBB) ---

def find_hbb_residue6(frames: dict):
    """Searches the 6 frames for the start of beta-globin (VHLTP = residues 1-5)
    and returns (frame_name, residue_6), or (None, None) if not found."""
    for name, prot in frames.items():
        i = prot.find("VHLTP")
        if i != -1 and len(prot) > i + 5:
            return name, prot[i + 5]
    return None, None


def run_hbb_diagnostic(frames: dict, lang: str) -> str:
    """Molecular biomarker detector for HBB residue 6 (bilingual)."""
    frame_name, aa = find_hbb_residue6(frames)

    report = ["\n" + "*" * 50, f" {t('report_header', lang)} ", "*" * 50]

    if aa is None:
        report.append(t("not_found", lang))
    else:
        report.append(t("found_in_frame", lang, frame=frame_name))
        if aa == 'E':
            report.append(t("diag_no_variant", lang))
            report.append(t("detail_glu", lang))
            report.append(t("diag_no_variant_detail", lang))
        elif aa == 'V':
            report.append(t("diag_sickle", lang))
            report.append(t("mutation_glu_val", lang))
            report.append(t("detail_val", lang))
        elif aa == 'K':
            report.append(t("diag_hbc", lang))
            report.append(t("mutation_glu_lys", lang))
            report.append(t("detail_lys", lang))
        else:
            report.append(t("diag_atypical", lang, aa=aa))
            report.append(t("warning_clinical", lang))

    report.append(t("disclaimer", lang))
    report.append("*" * 50)
    return "\n".join(report)


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
    diagnosis = run_hbb_diagnostic(frames, lang)

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
    description="Advanced genomic analyzer with bilingual support (EN/ES). Educational tool, not a clinical diagnostic."
)

if __name__ == "__main__":
    demo.launch()
