# -*- coding: utf-8 -*-
"""Wächter: keine .txt-Ausgaben lesen, nichts aus dem Archiv importieren, tempsim.py von niemandem importiert."""
import waechter as W


def test_bestand():
    assert W.quellen() == []
    assert W.tempsim() == []


def test_negativ_txt_und_archiv(tmp_path):
    p = tmp_path / "probe.py"
    p.write_text("import sys\nsys.path.insert(0, '../archiv/rechnungen')\nimport ftm_kern\nfrom sparc_streuung import x\n"
                 "d = open('gewicht_rueckkehr.txt').read()\n", encoding="utf-8")
    funde = W.quellen([str(p)], archiv={"ftm_kern", "sparc_streuung"})
    assert len(funde) == 4


def test_negativ_tempsim(tmp_path):
    p = tmp_path / "probe.py"
    p.write_text("import tempsim\n", encoding="utf-8")
    assert len(W.tempsim([str(p)])) == 1


def test_nur_echte_lesezugriffe(tmp_path):
    """Eine bloße Zeichenkette mit .txt (Beschreibung, Katalogtext, Wörterbuch) ist kein Fund; jeder echte Lesezugriff
    ist einer (Negativprobe: open, Variable, Path.read_text, np.loadtxt)."""
    p = tmp_path / "text.py"
    p.write_text("KATALOG = {'a': 'Ausgabe gewicht_rueckkehr.txt des Archivs'}\n"
                 "\"\"\"Die Datei bericht.txt entsteht nicht mehr.\"\"\"\n"
                 "def f(x):\n    return 'so hiess x.txt' + x\n", encoding="utf-8")
    assert W.quellen([str(p)], archiv=set()) == []
    q = tmp_path / "lesen.py"
    q.write_text("import numpy as np\nfrom pathlib import Path\n"
                 "a = open('a.txt').read()\n"
                 "pfad = 'b.txt'\nb = open(pfad, encoding='utf-8')\n"
                 "c = Path('c.txt').read_text()\n"
                 "d = np.loadtxt('d.txt')\n", encoding="utf-8")
    assert len(W.quellen([str(q)], archiv=set())) == 4
