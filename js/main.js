/* ------------------------------------------------------------------------
   main.js — wires the tab modules to the hash router.

   Canvases cannot be measured while their section is hidden, so each tab is
   drawn when it becomes active and on every resize.
   ------------------------------------------------------------------------ */

var MV = window.MV || (window.MV = {});

(function (MV) {
  'use strict';

  /* Chapter and section numbers come from the registry in equations.js; only
     the chapter word itself is literal, because chapter numbering is far more
     stable than equation numbering. */
  var CONTEXT = {
    kitaev: {
      chapter: 'Kapitel 3 — Kitaev-Kette, ' + MV.sec('3.3'),
      convention:
        'Normierung nach Alicea (2012): Hopping- und Paarungsterm mit Faktor 1/2, ' +
        'Phasenübergang daher bei |μ| = t, nicht bei 2t. ' +
        'Phase fest φ = 0. Alle Größen dimensionslos in Einheiten von t.'
    },
    nanodraht: {
      chapter: 'Kapitel 4 — Rashba-Nanodraht, ' + MV.sec('4.1', '4.2', '4.3'),
      convention:
        'ℏ = m = 1, alle Größen dimensionslos in Einheiten von Δ. ' +
        'Spin-Bahn-Achse ∥ ŷ, Zeeman-Feld senkrecht dazu. ' +
        'Topologisch für E_Z &gt; √(Δ² + μ²).'
    },
    auslese: {
      chapter: 'Kapitel 5 — Majorana-Box-Qubit, ' + MV.sec('5.2', '5.3', '5.5'),
      convention:
        'Achsenbezeichnung nach ' + MV.eq('tetron.pauli') + ', Konvention von Plugge et al.; ' +
        'Microsoft Quantum bezeichnet als Pauli-Z die Größe, die hier x̂ heißt. ' +
        'Parität P̂ = −γ₁γ₂γ₃γ₄ = +1 nach ' + MV.eq('tetron.parity') + '. ' +
        'Willkürliche Einheiten wie bei Plugge (ℏ = 1), G in Einheiten von ' +
        '(e²/h)·ν₁ν₂. φ ist hier der dimensionslose eingeschlossene Fluss — ' +
        'kein Flussquant. Die Formeln der Verfahren 2 und 3 stammen aus ' +
        'Plugge et al., nicht aus der Arbeit.'
    },
    gatter: {
      chapter: 'Kapitel 5 — Majorana-Box-Qubit, ' + MV.sec('5.4', '5.5'),
      convention:
        'Basis |0⟩ ist ẑ = +1 nach ' + MV.eq('tetron.pauli') + '. ' +
        'Alles dimensionslos; ein Zustandsraum hat keine Einheiten. ' +
        'Es gibt hier keine Zeitentwicklung — ein Gatter wird angewandt, ' +
        'nicht gefahren; die Wiedergabe in Panel B zeigt eine Schrittfolge mit ' +
        'Anzeigetempo. θ und φ sind die Winkel des Startzustands auf der ' +
        'Bloch-Kugel, weder Phase noch Fluss. Die Zusammenstellung der neun ' +
        'Gatter und das Protokoll für das T-Gatter stehen nicht in der ' +
        'Arbeit, sondern sind Standard der Quanteninformation.'
    },
    'pb-tetron': {
      chapter: 'Kapitel 6 — Experimenteller Status, ' + MV.sec('6.1', '6.2'),
      convention:
        'Quelle aller Angaben: [9] = Microsoft Quantum, arXiv:2606.03884 — ' +
        '<b>liegt nur als Preprint vor</b>. Energien in µeV, Längen in µm und nm, ' +
        'Zeiten in ms und s, Feld in mT, Temperatur in mK. „Pauli-Z" und τ_Z ' +
        'folgen der Konvention von [9]: gemeint ist die Parität eines Drahtes, ' +
        'nach ' + MV.eq('tetron.decomposition') + ' also die hier x̂ genannte ' +
        'Achse. Alle Zeichnungen sind eigene Schemata nach der Textbeschreibung, ' +
        'keine nachgezeichneten Abbildungen aus [9]. Dieser Reiter bewertet nicht.'
    },
    energieskalen: {
      chapter: 'Kapitel 5 — Majorana-Box-Qubit, ' + MV.sec('5.1'),
      convention:
        '<b>Physikalische Einheiten</b>: Energien in µeV, Längen in µm, ' +
        'Temperatur in mK, Zeiten in µs, weil ' + MV.eqn('energy.window') +
        ' absolute Skalen vergleicht. Dimensionslos bleiben die vier Reiter zu ' +
        'den Kapiteln 3 bis 5: Kitaev-Kette, Rashba-Nanodraht, Auslese und ' +
        'Gatter; der Reiter zu Kapitel 6 trägt ebenfalls Einheiten. ' +
        'Alle Vorgabewerte hier sind illustrativ; einzige Zahl aus der Arbeit ' +
        'ist T = 20 mK. φ kommt hier nicht vor.'
    }
  };

  function start() {
    var contextEl = document.getElementById('context-line');
    var tabKitaev = document.getElementById('tab-kitaev');
    var tabNano = document.getElementById('tab-nanodraht');
    var tabAuslese = document.getElementById('tab-auslese');
    var tabGatter = document.getElementById('tab-gatter');
    var tabPb = document.getElementById('tab-pb-tetron');
    var tabEnergy = document.getElementById('tab-energieskalen');

    if (MV.kitaev) MV.kitaev.init(tabKitaev);
    if (MV.nanowire) MV.nanowire.init(tabNano);
    if (MV.auslese) MV.auslese.init(tabAuslese);
    if (MV.gatter) MV.gatter.init(tabGatter);
    if (MV.pbTetron) MV.pbTetron.init(tabPb);
    if (MV.energyscales) MV.energyscales.init(tabEnergy);

    /* Einheitlich für alle sechs Reiter: Regler in eine schmale Spalte links
       neben die Abbildungen, damit beim Ziehen die Grafik sichtbar bleibt.
       CLAUDE.md, Abschnitt 9, verlangt ein konsistentes Layout — eine
       Ausnahme für einen Reiter wäre später erklärungsbedürftig. */
    [tabKitaev, tabNano, tabEnergy, tabAuslese, tabGatter, tabPb].forEach(function (t) {
      if (t) MV.ui.twoColumn(t);
    });

    /* Erst jetzt, wenn jeder Reiter steht und seine Schnittstelle gemeldet
       hat, werden die Vorführungen an den Leitsätzen registriert. */
    if (MV.demos) MV.demos.install();

    var routes = [
      {
        hash: '#kitaev',
        tabEl: tabKitaev,
        linkEl: document.getElementById('link-kitaev'),
        chapter: CONTEXT.kitaev.chapter,
        convention: CONTEXT.kitaev.convention,
        onEnter: function () { if (MV.kitaev) MV.kitaev.draw(); },
        onLeave: function () { if (MV.kitaev && MV.kitaev.stop) MV.kitaev.stop(); }
      },
      {
        hash: '#nanodraht',
        tabEl: tabNano,
        linkEl: document.getElementById('link-nanodraht'),
        chapter: CONTEXT.nanodraht.chapter,
        convention: CONTEXT.nanodraht.convention,
        onEnter: function () { if (MV.nanowire) MV.nanowire.draw(); },
        onLeave: function () { if (MV.nanowire && MV.nanowire.stop) MV.nanowire.stop(); }
      },
      {
        hash: '#energieskalen',
        tabEl: tabEnergy,
        linkEl: document.getElementById('link-energieskalen'),
        chapter: CONTEXT.energieskalen.chapter,
        convention: CONTEXT.energieskalen.convention,
        onEnter: function () { if (MV.energyscales) MV.energyscales.draw(); },
        onLeave: function () { if (MV.energyscales && MV.energyscales.stop) MV.energyscales.stop(); }
      },
      {
        hash: '#auslese',
        tabEl: tabAuslese,
        linkEl: document.getElementById('link-auslese'),
        chapter: CONTEXT.auslese.chapter,
        convention: CONTEXT.auslese.convention,
        onEnter: function () { if (MV.auslese) MV.auslese.draw(); },
        onLeave: function () { if (MV.auslese && MV.auslese.stop) MV.auslese.stop(); }
      },
      {
        hash: '#gatter',
        tabEl: tabGatter,
        linkEl: document.getElementById('link-gatter'),
        chapter: CONTEXT.gatter.chapter,
        convention: CONTEXT.gatter.convention,
        onEnter: function () { if (MV.gatter) MV.gatter.draw(); },
        onLeave: function () { if (MV.gatter && MV.gatter.stop) MV.gatter.stop(); }
      },
      {
        hash: '#pb-tetron',
        tabEl: tabPb,
        linkEl: document.getElementById('link-pb-tetron'),
        chapter: CONTEXT['pb-tetron'].chapter,
        convention: CONTEXT['pb-tetron'].convention,
        onEnter: function () { if (MV.pbTetron) MV.pbTetron.draw(); },
        onLeave: function () { if (MV.pbTetron && MV.pbTetron.stop) MV.pbTetron.stop(); }
      }
    ];

    var router = MV.ui.router(routes, contextEl);

    /* Die Vorführungen an den Leitsätzen brauchen zweierlei vom Router: einen
       Weg, einen Reiter zu öffnen, und die Handhabung der Direktlinks
       #<reiter>/<Panelbuchstabe>. Beides hier, damit demo.js den Router nicht
       kennen muss. */
    MV.demo.setTabOpener(function (hash, silent) {
      if (!silent && window.location.hash !== hash) window.location.hash = hash;
      router.show(hash);
    });

    /* Resize redraws only the visible tab. */
    MV.onResize(function () {
      if (!tabKitaev.hidden && MV.kitaev) MV.kitaev.draw();
      if (!tabNano.hidden && MV.nanowire) MV.nanowire.draw();
      if (!tabAuslese.hidden && MV.auslese) MV.auslese.draw();
      if (!tabGatter.hidden && MV.gatter) MV.gatter.draw();
      if (!tabPb.hidden && MV.pbTetron) MV.pbTetron.draw();
      if (!tabEnergy.hidden && MV.energyscales) MV.energyscales.draw();
    });
  }

  if (document.readyState === 'loading') {
    document.addEventListener('DOMContentLoaded', start);
  } else {
    start();
  }

}(MV));
