# -*- coding: utf-8 -*-
"""Zwischenspeicher von lauf.py: der Schlüssel folgt dem Code, den Datenhashes und dem Modelltext, nicht dem Bericht."""
import lauf


def _baum(tmp_path):
    sim = tmp_path / "sim"
    (sim / "kern").mkdir(parents=True)
    (sim / "kern" / "a.py").write_text("x = 1\n", encoding="utf-8")
    (sim / "kern" / "daten_hashes.json").write_text("{}", encoding="utf-8")
    modell = tmp_path / "modell"
    modell.mkdir()
    (modell / "01-a.md").write_text("Text 1,5\n", encoding="utf-8")
    return sim, modell


def test_schluessel_folgt_code_daten_und_modell(tmp_path):
    sim, modell = _baum(tmp_path)
    k0 = lauf.schluessel(str(sim), str(modell))
    assert lauf.schluessel(str(sim), str(modell)) == k0                       # stabil
    (sim / "bericht.md").write_text("neu", encoding="utf-8")                  # Bericht und Register zählen nicht
    (sim / "register.json").write_text("{}", encoding="utf-8")
    assert lauf.schluessel(str(sim), str(modell)) == k0
    (sim / "kern" / "a.py").write_text("x = 2\n", encoding="utf-8")           # Negativproben: jede Änderung ändert ihn
    k1 = lauf.schluessel(str(sim), str(modell))
    assert k1 != k0
    (sim / "kern" / "daten_hashes.json").write_text('{"a": 1}', encoding="utf-8")
    k2 = lauf.schluessel(str(sim), str(modell))
    assert k2 != k1
    (modell / "01-a.md").write_text("Text 1,6\n", encoding="utf-8")
    assert lauf.schluessel(str(sim), str(modell)) != k2
