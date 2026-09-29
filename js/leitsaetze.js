/* ------------------------------------------------------------------------
   leitsaetze.js — über jedem Panel steht der Satz aus der Arbeit, den dieses
   Panel sichtbar macht.

   CLAUDE.md, Abschnitt 11: "Über jedem Panel steht sichtbar, ohne Klick, eine
   Zeile: der Satz aus der Arbeit, den dieses Panel sichtbar macht. Wörtlich
   zitiert, in Anführungszeichen, mit Abschnittsangabe."

   Und die Regel, die daraus folgt: findet sich kein tragender Satz, wird
   keiner erfunden. Dann steht `open: true`, das Panel zeigt sichtbar
   "Leitsatz offen", und es wird in NOTES.md unter "Panels ohne tragenden Satz
   in der Arbeit" geführt.

   Alle Zitate hier sind wörtlich aus referenz.pdf übernommen; die Seitenangabe
   in der Klammer ist die gedruckte Seite und dient nur dem Wiederfinden beim
   Abgleich. Auslassungen sind nicht vorgenommen worden — wo ein Satz eine
   Literaturangabe in eckigen Klammern trägt, endet das Zitat davor.
   ------------------------------------------------------------------------ */

var MV = window.MV || (window.MV = {});

(function (MV) {
  'use strict';

  var LEITSAETZE = {

    /* ---------------------------------------------------------- #kitaev */

    'kitaev.A': {
      section: '3.3',
      page: 12,
      text: 'In dieser Darstellung ist ν die Frage, ob die geschlossene Bahn ' +
            'den Ursprung umschließt.'
    },

    'kitaev.B': {
      section: '3.3',
      page: 11,
      text: 'Die Lücke schließt nur, wenn beide Summanden gleichzeitig ' +
            'verschwinden.'
    },

    /* ------------------------------------------------------- #nanodraht */

    'nanodraht.A': {
      section: '4.1',
      page: 17,
      text: 'Liegt μ innerhalb dieser Lücke, so ist nur noch das untere Band ' +
            'besetzt: zwei Fermi-Punkte statt vier.'
    },

    'nanodraht.B': {
      section: '4.2',
      page: 19,
      text: 'Die effektive Paarung fällt mit wachsendem E_Z ab, weil ein ' +
            'starkes Feld die Spins zunehmend ausrichtet und die helikale ' +
            'Struktur unterdrückt.'
    },

    'nanodraht.C': {
      section: '4.3',
      page: 19,
      text: 'Diese Linie trennt die beiden Phasen, und ein Vergleich mit dem ' +
            'Grenzfall Δ → 0 legt fest, welche Seite die topologische ist.'
    },

    /* --------------------------------------------------------- #auslese */

    'auslese.A': {
      section: '5.3',
      page: 28,
      text: 'Es gibt genau drei Möglichkeiten, vier Majoranas paarweise ' +
            'zusammenzufassen – längs, quer und diagonal –, und diese drei ' +
            'Möglichkeiten sind die drei Pauli-Achsen.'
    },

    'auslese.B': {
      section: '5.3',
      page: 29,
      text: 'Ohne Referenzarm, t₀ = 0, verschwindet die z-Abhängigkeit ' +
            'vollständig – man misst dann zwar einen Strom, aber keine Parität.'
    },

    /* §5.5, letzter Unterabschnitt „Alles hängt an einer einzigen Primitive",
       wörtlich gegen die Abgabefassung geprüft (S. 33). */
    'auslese.C': {
      section: '5.5',
      page: 33,
      text: 'Der Preis dafür ist, dass die gesamte experimentelle Last auf ' +
            'einer einzigen Primitive liegt – der einzelschussfähigen, ' +
            'projektiven Messung eines Majorana-Bilinearen.'
    },

    /* ---------------------------------------------------------- #gatter */

    /* §5.4, letzter Absatz vor "Der entscheidende Punkt"; woertlich gegen die
       Abgabefassung geprueft (S. 31). Der Satz beginnt im Original mitten in
       der Zeile, das Zitat setzt an seinem Anfang an. */
    'gatter.A': {
      section: '5.4',
      page: 31,
      text: 'Eine Drehung um π/4 – das T-Gatter – ist auf diesem Weg nicht ' +
            'erreichbar, und ohne sie ist die Rechnung nicht universell und ' +
            'nach dem Gottesman-Knill-Theorem klassisch effizient simulierbar.'
    },

    /* §5.5, Unterabschnitt "Alles haengt an einer einzigen Primitive";
       woertlich gegen die Abgabefassung geprueft (S. 33). */
    'gatter.B': {
      section: '5.5',
      page: 33,
      text: 'Sämtliche geschützten Operationen der Architektur reduzieren ' +
            'sich auf projektive Paritätsmessungen plus klassische ' +
            'Nachbearbeitung.'
    },

    /* ------------------------------------------------------- #pb-tetron */

    /* §6.2, Unterabschnitt „Bauelement"; wörtlich gegen die Abgabefassung
       geprüft (S. 36). */
    'pbtetron.A': {
      section: '6.2',
      page: 36,
      text: 'Sämtliche im Folgenden besprochenen Messungen erfolgen am oberen ' +
            'Nanodraht des Tetrons BA – also an einem geerdeten Bauelement.'
    },

    /* §6.2, Unterabschnitt „Charakterisierung"; wörtlich gegen die
       Abgabefassung geprüft (S. 37). */
    'pbtetron.B': {
      section: '6.2',
      page: 37,
      text: 'Das Phasendiagramm dagegen, aus dem ∆_T ≈ 70 µeV und die ' +
            'Ausdehnung der topologischen Phase von über 1,1 mV T stammen, ' +
            'wurde mit dem Topological-Gap-Protokoll an einer separaten ' +
            '3 µm-Nanodraht-Teststruktur gewonnen.'
    },

    /* §6.1, letzter Absatz; wörtlich gegen die Abgabefassung geprüft (S. 36). */
    'pbtetron.C': {
      section: '6.1',
      page: 36,
      text: 'Der Materialwechsel wirkt also weit stärker auf den ' +
            'Vergiftungskanal als auf den Anregungskanal.'
    },

    /* ---------------------------------------------------- #energieskalen */

    'energieskalen.A': {
      section: '5.1',
      page: 24,
      text: 'Ein längerer Draht drückt die Aufspaltung exponentiell nach ' +
            'unten, schwächt aber gleichzeitig den Ladungsschutz.'
    },

    'energieskalen.B': {
      section: '5.6',
      page: 34,
      text: 'Umgekehrt kann ein Bauelement eine beeindruckende Lebensdauer ' +
            'T_P aufweisen, ohne dass die Ladungsenergie überhaupt als ' +
            'Schutzmechanismus aktiv wäre.'
    },

  };

  /**
   * HTML für die Leitzeile eines Panels. Ohne Eintrag oder mit `open: true`
   * steht dort sichtbar "Leitsatz offen" — nie ein erfundener Satz.
   */
  MV.leitsatzHtml = function (key) {
    var e = LEITSAETZE[key];
    if (!e || e.open) {
      return '<span class="panel__quote-open">Leitsatz offen</span>' +
             '<span class="panel__quote-src">in der Arbeit ist noch kein ' +
             'tragender Satz für dieses Panel formuliert</span>';
    }
    return '<span class="panel__quote-text">„' + e.text + '“</span>' +
           '<span class="panel__quote-src">' + MV.secOfWork(e.section) + '</span>';
  };

  /** Schlüssel aller Panels ohne Zitat — Grundlage der Liste in NOTES.md. */
  MV.leitsatzOffen = function () {
    return Object.keys(LEITSAETZE).filter(function (k) {
      return LEITSAETZE[k].open;
    });
  };

  MV.leitsaetze = LEITSAETZE;

}(MV));
