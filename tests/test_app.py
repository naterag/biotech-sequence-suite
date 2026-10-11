import warnings

import pytest
from Bio.Seq import Seq

from app import (
    analyze_dna_gui,
    clean_sequence,
    find_hbb_residue6,
    translate_6_frames,
)

# Secuencias de control (las mismas del README)
HBA = "ATGGTGCACCTGACTCCTGAGGAGAAGTCTGCCGTTACT"
HBS = "ATGGTGCACCTGACTCCTGTGGAGAAGTCTGCCGTTACT"
HBC = "ATGGTGCACCTGACTCCTAAGGAGAAGTCTGCCGTTACT"

CASOS = [("HbA", HBA, "E"), ("HbS", HBS, "V"), ("HbC", HBC, "K")]


def residuo6(seq):
    return find_hbb_residue6(translate_6_frames(seq))


# --- Limpieza de secuencia ---

def test_clean_sequence_quita_ruido_y_pasa_a_mayusculas():
    assert clean_sequence("atg gtg\n123 cac-xyz") == "ATGGTGCAC"


def test_clean_sequence_vacia():
    assert clean_sequence("12345 !!!") == ""


# --- Traducción en 6 marcos ---

def test_hay_seis_marcos():
    frames = translate_6_frames(HBA)
    assert set(frames) == {"+1", "+2", "+3", "-1", "-2", "-3"}


def test_marco_1_empieza_con_la_beta_globina():
    assert translate_6_frames(HBA)["+1"].startswith("MVHLTPEEKSAVT")


def test_traduccion_sin_warnings_si_no_es_multiplo_de_3():
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        translate_6_frames(HBA + "AG")


# --- Detección del residuo 6 ---

@pytest.mark.parametrize("nombre,seq,esperado", CASOS)
def test_secuencia_directa(nombre, seq, esperado):
    frame, aa = residuo6(seq)
    assert frame == "+1"
    assert aa == esperado


@pytest.mark.parametrize("nombre,seq,esperado", CASOS)
def test_con_texto_antes_del_gen(nombre, seq, esperado):
    _, aa = residuo6("GGGCCC" + seq)
    assert aa == esperado


@pytest.mark.parametrize("nombre,seq,esperado", CASOS)
def test_hebra_reversa(nombre, seq, esperado):
    rev = str(Seq(seq).reverse_complement())
    frame, aa = residuo6(rev)
    assert frame.startswith("-")
    assert aa == esperado


def test_sin_motivo_vhltp():
    assert residuo6("ATGAAATTTGGGCCCAAATTT") == (None, None)


# --- Interfaz (función que usa Gradio) ---

def test_gui_hbs_en_espanol():
    limpia, marcos, reporte = analyze_dna_gui(HBS, "Español")
    assert limpia == HBS
    assert "Glu6Val" in reporte
    assert "HbS" in reporte
    assert "No es un diagnóstico clínico" in reporte


def test_gui_hba_en_ingles():
    _, _, reporte = analyze_dna_gui(HBA, "English")
    assert "No variant detected" in reporte
    assert "not a clinical diagnostic" in reporte


def test_gui_hbc():
    _, _, reporte = analyze_dna_gui(HBC, "English")
    assert "Glu6Lys" in reporte


def test_gui_entrada_invalida():
    limpia, marcos, reporte = analyze_dna_gui("12345", "English")
    assert "No valid DNA bases" in limpia
    assert marcos == "" and reporte == ""


def test_gui_sin_motivo_muestra_aviso():
    _, _, reporte = analyze_dna_gui("ATGAAATTTGGGCCCAAATTT", "English")
    assert "not found" in reporte
