# -*- coding: utf-8 -*-
"""
Erzeugt die Abschnitte von doku.md zwischen <!-- erzeugt: name --> und <!-- ende --> aus den Registern
(Symbole, Formeln, Parametersätze) und aus register.json (Prüfzahlen, Gesetze, Themen, Beschlüsse).
Der Fließtext zwischen den Markern bleibt unverändert. lauf.py ruft main() am Ende eines gültigen Laufs auf.
"""
import json
import os
import re
import sys

SIM = os.path.dirname(os.path.abspath(__file__))
if SIM not in sys.path:
    sys.path.insert(0, SIM)
DOKU = os.path.join(SIM, "doku.md")
REGISTER_JSON = os.path.join(SIM, "register.json")
MARKER = re.compile(r"(<!-- erzeugt: (\w+) -->\n)(.*?)(<!-- ende -->)", re.S)
STAND = {"fest": "fest", "seit_z10": "seit z≈10", "epoche": "Epoche", "je_system": "je System"}


def z(s):
    return str(s if s is not None else "").replace("|", "\\|").replace("\n", " ").strip()


def tafel(kopf, zeilen):
    t = ["| " + " | ".join(kopf) + " |", "|" + "---|" * len(kopf)]
    return t + ["| " + " | ".join(z(c) for c in r) + " |" for r in zeilen]


def symbole(reg, _):
    from kern.symbole import SYMBOLE
    zl = [(f"`{s.zeichen}`", s.bedeutung, s.einheit_si, s.ebene, STAND[s.zustand]) for s in SYMBOLE.values()]
    return [f"{len(zl)} Zeichen, jedes einmal.", ""] + tafel(("Zeichen", "Bedeutung", "Einheit", "Ebene", "Zustand"), zl)


def formeln(reg, _):
    from kern.formeln import FORMELN
    out = []
    for e, name in (("L", "Lokal"), ("I", "Innere Richtung"), ("Z", "Zustandsraum"), ("G", "Das Ganze")):
        fs = [f for f in FORMELN.values() if f.status == "modell" and f.ebene[0] == e]
        if fs:
            out += [f"### {name} ({e}), {len(fs)} Formeln", ""]
            out += tafel(("Kennung", "Gleichung", "Herkunft"), [(f.kennung, f.gleichung, f.herkunft) for f in fs]) + [""]
    return out[:-1]


def parametersaetze(reg, _):
    from kern.parametersatz import PARAMETERSAETZE
    return [f"- **{p.kennung}**: {p.text().split(': ', 1)[-1]}" for p in PARAMETERSAETZE.values()]


def pruefstand(reg, _):
    w = reg["waechter"]
    return [f"Wächter im letzten Lauf: {len(w)} Prüfungen, Funde je Wächter: "
            + ", ".join(f"{n} {k}" for n, k in w.items()) + ".",
            "", f"Messwerte {len(reg['messwerte'])} (davon mit [G] {sum(1 for m in reg['messwerte'].values() if m.get('art') == 'G')}), "
            f"Setzungen {len(reg['setzungen'])}, Eingaben {len(reg['eingaben'])}, "
            f"Gegenproben {len(reg['gegenproben'])}."]


def ergebnisse(reg, _):
    pz, ge = reg["pruefzahlen"], reg["gesetze"]
    out = [f"Lauf {reg['commit']}, {reg['zeit']}.", "", "**Prüfzahlen je Baustein**", ""]
    zl = []
    for b in sorted({v["baustein"] for v in pz.values()}):
        v = [x for x in pz.values() if x["baustein"] == b]
        zl.append((b, len(v)) + tuple(sum(1 for x in v if x["ergebnis"] == e) for e in ("OK", "ABWEICHUNG", "OFFEN"))
                  + (sum(1 for x in v if x["g"]), sum(1 for x in v if x["w"])))
    zl.append(("alle", len(pz)) + tuple(sum(1 for x in pz.values() if x["ergebnis"] == e) for e in ("OK", "ABWEICHUNG", "OFFEN"))
              + (sum(1 for x in pz.values() if x["g"]), sum(1 for x in pz.values() if x["w"])))
    out += tafel(("Baustein", "Zahlen", "OK", "ABWEICHUNG", "OFFEN", "mit [G]", "mit [W]"), zl)
    out += ["", "**Gesetze je Klasse**", ""]
    zl = []
    for k in sorted({v["klasse"] for v in ge.values()}):
        v = [x for x in ge.values() if x["klasse"] == k]
        zl.append((k, len(v)) + tuple(sum(1 for x in v if x["ergebnis"] == e) for e in ("OK", "SPANNUNG", "WIDERSPRUCH", "HINWEIS")))
    out += tafel(("Klasse", "Gesetze", "OK", "SPANNUNG", "WIDERSPRUCH", "HINWEIS"), zl)
    ab = [(k, v["aussage"], v["ergebnis"]) for k, v in ge.items() if v["ergebnis"] != "OK"]
    if ab:
        out += ["", "**Gesetze ohne OK**", ""] + tafel(("Gesetz", "Aussage", "Ergebnis"), ab)
    ab = [(k, v["baustein"], v["modell_text"], v["wert"], v["soll"]) for k, v in pz.items() if v["ergebnis"] == "ABWEICHUNG"]
    out += ["", f"**Abweichungen der Prüfzahlen** ({len(ab)})", ""]
    out += tafel(("Prüfzahl", "Baustein", "Modelltext", "Rechnung", "Soll"),
                 [(k, b, m, f"{w:.6g}" if isinstance(w, (int, float)) else w,
                   f"{s:.6g}" if isinstance(s, (int, float)) else s) for k, b, m, w, s in ab])
    return out


def offen(reg, _):
    from forschung.beschluesse import OFFEN, VORSCHLAEGE
    from forschung.beschluss import BESCHLUESSE
    out = ["**Offene Beschlüsse** (Entscheidung der Hauptsitzung; der Modelltext bleibt bis dahin)", ""]
    out += tafel(("Kennung", "Frage", "vorgeschlagener Satz"),
                 [(k, BESCHLUESSE[k].aussage if k in BESCHLUESSE else "", f"{art}: {grund}")
                  for k, (art, grund) in OFFEN.items()]
                 + [(k, f"{v.modell} → {v.frisch} ({v.grund})", v.satz) for k, v in VORSCHLAEGE.items()])
    out += ["", "**Themen**", ""]
    out += tafel(("Thema", "Frage", "Stand"), [(k, v["frage"], v["stand"]) for k, v in reg["themen"].items()])
    wurz = sorted({w for v in list(reg["pruefzahlen"].values()) + list(reg["gesetze"].values()) for w in v["w"]})
    out += ["", "**Beschlüsse**", "",
            f"{len(reg['beschluesse'])} Beschlüsse; ins Modell: {sum(1 for b in reg['beschluesse'].values() if b['richtung'] == 'hinein')}, "
            f"heraus: {sum(1 for b in reg['beschluesse'].values() if b['richtung'] == 'heraus')}. Zur Bestätigung durch den Autor stehen die Wahlen [W] ({len(wurz)}):", ""]
    out += [f"- {w}" for w in wurz]
    return out


ABSCHNITTE = dict(symbole=symbole, formeln=formeln, parametersaetze=parametersaetze, pruefstand=pruefstand,
                  ergebnisse=ergebnisse, offen=offen)


def main():
    import laden
    laden.lade_alles()
    with open(REGISTER_JSON, encoding="utf-8") as f:
        reg = json.load(f)
    with open(DOKU, encoding="utf-8") as f:
        text = f.read()
    def ersatz(m):
        name = m.group(2)
        if name not in ABSCHNITTE:
            raise KeyError(f"doku.md: Marker {name!r} unbekannt")
        return m.group(1) + "\n".join(ABSCHNITTE[name](reg, None)) + "\n" + m.group(4)
    neu = MARKER.sub(ersatz, text)
    with open(DOKU, "w", encoding="utf-8") as f:
        f.write(neu)
    print(f"doku_erzeugen.py: {len(MARKER.findall(text))} Abschnitte → {os.path.relpath(DOKU)}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
