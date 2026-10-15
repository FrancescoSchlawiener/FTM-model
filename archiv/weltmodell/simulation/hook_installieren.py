# -*- coding: utf-8 -*-
"""
Installiert den Git-pre-commit-Hook: Vor jedem Commit laufen nur die schnellen Wächter und die Registrierung (pytest in
tests/ ohne die als "schwer" gekennzeichneten Tests, etwa 5 s; Ziel unter 10 s). Gesetze, Gegenproben, Doppelformeln
und Prüfzahlen rechnet lauf.py, nicht der Hook. Ein strenger Fund bricht den Commit ab. Aufruf: python hook_installieren.py (einmal je Klon).
"""
import os
import stat
import subprocess
import sys

SIM = os.path.dirname(os.path.abspath(__file__))
HOOK = """#!/bin/sh
# Wächter der FTM-Simulation (installiert von weltmodell/simulation/hook_installieren.py)
log="$(git rev-parse --absolute-git-dir)/ftm_waechter.log"
cd "$(git rev-parse --show-toplevel)/weltmodell/simulation" || exit 1
{python} -m pytest -q -x -p no:cacheprovider -m "not schwer" tests > "$log" 2>&1
status=$?
if [ $status -ne 0 ]; then
    tail -n 30 "$log"
    echo "Wächter rot: Commit abgebrochen (Log $log)."
fi
exit $status
"""


def main():
    git_dir = subprocess.run(["git", "rev-parse", "--git-path", "hooks"], cwd=SIM, capture_output=True, text=True,
                             check=True).stdout.strip()
    ziel = os.path.join(SIM, git_dir, "pre-commit") if not os.path.isabs(git_dir) else os.path.join(git_dir, "pre-commit")
    with open(ziel, "w", encoding="utf-8") as f:
        f.write(HOOK.format(python=sys.executable))
    os.chmod(ziel, os.stat(ziel).st_mode | stat.S_IXUSR | stat.S_IXGRP | stat.S_IXOTH)
    print(f"pre-commit installiert: {os.path.normpath(ziel)}")


if __name__ == "__main__":
    main()
