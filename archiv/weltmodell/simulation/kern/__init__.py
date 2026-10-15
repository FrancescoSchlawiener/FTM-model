# -*- coding: utf-8 -*-
"""
FTM-Kern als Paket: alle Formeln des Modells, je Ebene ein Modul, mit einheitlichen Konstanten.
Abhängigkeiten nur in eine Richtung: konstanten → faltung (L) / rueckkehr (I) → zustandsraum (Z), kosmos (G),
teilchen_stern, galaxie, naht_haufen, licht. `from kern import *` stellt alle Namen bereit.
"""
from .konstanten import *
from .faltung import *
from .rueckkehr import *
from .zustandsraum import *
from .kosmos import *
from .teilchen_stern import *
from .galaxie import *
from .naht_haufen import *
from .naht_haufen import _UMKEHR
from .licht import *
from .daten import *
from . import formeln, symbole, messwert, messwerte, setzung, setzungen, eingabe, eingaben  # noqa: F401
from . import parametersatz, parametersaetze, herkunft_g, verfolgung, datenhashes  # noqa: F401
