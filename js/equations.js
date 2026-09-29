/* ------------------------------------------------------------------------
   equations.js — one place for every number that refers into the thesis.

   CLAUDE.md, section 9: "Gleichungsnummern zentral. Ein Modul hält Nummer und
   Abschnittstitel je Gleichung; Info-Panels referenzieren nur Schlüssel, keine
   Zeichenketten im Fließtext. Die Nummerierung der Arbeit ist noch im Fluss."

   Whenever the thesis is renumbered, only this file changes. The generator
   tools/gleichungen.js reads the registry below and scans the tab modules for
   the keys they use, and writes GLEICHUNGEN.txt from both — so that list can
   never drift away from the code.

   Keys are spoken names, not numbers, on purpose: a key that reads
   `kitaev.spectrum` stays correct when the equation moves.
   ------------------------------------------------------------------------ */

var MV = window.MV || (window.MV = {});

(function (MV) {
  'use strict';

  /* --------------------------------------------------------------- sections */

  var SECTIONS = {
    '3.3': 'Bulk-Spektrum, Gap-Schließung und topologische Invariante',
    '3.5': 'Exponentielle Aufspaltung bei endlicher Kettenlänge',
    '4.1': 'Das Zutatenrezept',
    '4.2': 'Effektive spinlose p-Wellen-Paarung im helikalen Regime',
    '4.3': 'Topologisches Kriterium und Phasendiagramm',
    '5.1': 'Coulomb-blockierte Box: Ladungsenergie als Schutzmechanismus',
    '5.2': 'Tetron-Geometrie',
    '5.3': 'Pauli-Operatoren aus Majorana-Bilinearen',
    '5.4': 'Braiding: physikalische Prozedur und unitärer Operator',
    '5.5': 'Messbasierte Operationen',
    '5.6': 'Was das Modell vom Experiment verlangt',
    '6.1': 'Materialwahl: InAs–Al gegenüber InAs–Pb',
    '6.2': 'Das InAs–Pb-Tetron-Experiment',
    '6.5': 'Eigene Bewertung'
  };

  /* -------------------------------------------------------------- equations */

  var EQUATIONS = {
    'kitaev.hamiltonian':   { num: '3.1',  section: '3.1', desc: 'Modell-Hamiltonian der Kitaev-Kette' },
    'kitaev.majoranaForm':  { num: '3.3',  section: '3.1', desc: 'Hamiltonian in Majorana-Form' },
    'kitaev.sweetSpot':     { num: '3.4',  section: '3.2', desc: 'topologischer Grenzfall μ = 0, t = Δ' },
    'kitaev.bulk':          { num: '3.7',  section: '3.3', desc: 'Bulk-Größen ε_k und Δ̃_k' },
    'kitaev.spectrum':      { num: '3.8',  section: '3.3', desc: 'Anregungsspektrum E(k)' },
    'kitaev.criticalLine':  { num: '3.9',  section: '3.3', desc: 'kritische Linie |μ| = t' },
    'kitaev.signs':         { num: '3.10', section: '3.3', desc: 'Vorzeichen s₀ und s_π an k = 0 und k = π' },
    'kitaev.invariant':     { num: '3.11', section: '3.3', desc: 'Z₂-Invariante ν = s₀ · s_π' },
    'kitaev.splitting':     { num: '3.14', section: '3.5', desc: 'exponentielle Aufspaltung δE ~ exp(−L/ξ)' },

    'wire.model':           { num: '4.1',  section: '4.1', desc: 'Modell des Rashba-Nanodrahts' },
    'wire.bands':           { num: '4.2',  section: '4.1', desc: 'Bänder ε±(k) des Drahtes' },
    'wire.deltaEff':        { num: '4.3',  section: '4.2', desc: 'effektive p-Wellen-Paarung Δ_eff(k)' },
    'wire.criterion':       { num: '4.4',  section: '4.3', desc: 'topologisches Kriterium E_Z > √(Δ² + μ²)' },

    'energy.charging':      { num: '5.1',  section: '5.1', desc: 'Ladungsenergie H_C = E_C(N̂ − n_g)², E_C = e²/2C' },
    'energy.tradeoff':      { num: '5.3',  section: '5.1', desc: 'Zielkonflikt E_C ~ 1/L_W gegen δE ~ Δ_T exp(−L_W/ξ)' },
    'energy.window':        { num: '5.4',  section: '5.1', desc: 'Bedingungsfenster δE ≪ k_B T ≪ E_C' },

    'tetron.parity':        { num: '5.5',  section: '5.3', desc: 'Paritätskonvention P̂ = −γ₁γ₂γ₃γ₄ = +1' },
    'tetron.pauli':         { num: '5.6',  section: '5.3', desc: 'Pauli-Operatoren aus Majorana-Bilinearen' },
    'tetron.cyclic':        { num: '5.7',  section: '5.3', desc: 'zyklische Relation x̂ŷ = iẑ' },
    'tetron.decomposition': { num: '5.8',  section: '5.3', desc: 'drei Zerlegungen, jede Achse über zwei Schleifen' },
    'ops.braid':            { num: '5.11', section: '5.4', desc: 'Zopfoperator U_ij = exp((π/4) γ_i γ_j) = (1/√2)(1 + γ_i γ_j)' },
    'ops.braidZ':           { num: '5.12', section: '5.4', desc: 'U₂₃ = exp(−i(π/4) ẑ), Drehung um π/2 um die z-Achse' },
    'ops.sz':               { num: '5.13', section: '5.5', desc: 'Ŝ_z = diag(1, i) aus zwei Messungen und einer Pauli-Korrektur' },

    'tetron.readoutH':      { num: '5.9',  section: '5.3', desc: 'Hamiltonian der interferometrischen Auslese' },
    'tetron.conductance':   { num: '5.10', section: '5.3', desc: 'Leitfähigkeit G_z = (e²/h) ν₁ν₂ |t₀ + t₁z|²' },

    'pb.tauZ':              { num: '6.1',  section: '6.2', desc: 'gemessene Paritätslebensdauer τ_Z = 22 ± 1 s' }
  };

  /* ---------------------------------------------------------------- figures */

  var FIGURES = {
    'kitaev.spectrumFig': { num: '3.2a', section: '3.3', desc: 'Anregungsspektrum ±E(k) für drei Werte von μ/t' },
    'kitaev.gapFig':      { num: '3.2b', section: '3.3', desc: 'Energielücke E_gap über μ/t' },
    'kitaev.orbitFig':    { num: '3.2c', section: '3.3', desc: 'Bahn von h(k) in der (h_y, h_z)-Ebene' },
    'wire.setupFig':      { num: '4.1',  section: '4.1', desc: 'schematischer Aufbau des Nanodrahts' },
    'wire.stagesFig':     { num: '4.2',  section: '4.1', desc: 'Aufbau der topologischen Phase, Zutat für Zutat' },
    'wire.stage3Fig':     { num: '4.2c', section: '4.1', desc: 'Spektrum ±E(k) mit Paarung' },
    'wire.phaseFig':      { num: '4.3',  section: '4.3', desc: 'Phasendiagramm in den Achsen μ/Δ und E_Z/Δ' },
    'tetron.layoutFig':   { num: '5.1',  section: '5.2', desc: 'Aufbau des Majorana-Box-Qubits' },
    'ops.braidFig':       { num: '5.4',  section: '5.4', desc: 'Prozedur, Zopfstruktur und Wirkung der Vertauschung' },
    'energy.windowFig':   { num: '5.2',  section: '5.1', desc: 'Zielkonflikt in logarithmischer Auftragung, Parameterfenster' }
  };

  /* ----------------------------------------------------------------- tables */

  var TABLES = {
    'criteria':     { num: '5.1', section: '5.6', desc: 'Bedingungen K1 bis K6 mit der jeweils messbaren Größe' },
    'materials':      { num: '6.1', section: '6.1', desc: 'Energieskalen der beiden Materialplattformen, Al gegen Pb' },
    'criteriaCheck':  { num: '6.2', section: '6.5', desc: 'Prüfung der Bedingungen K1 bis K6 am InAs–Pb-Tetron' }
  };

  function need(table, key, what) {
    var e = table[key];
    if (!e) throw new Error('Unbekannter ' + what + '-Schlüssel: ' + key);
    return e;
  }

  /* ------------------------------------------------------------------- API */

  /** "(3.8)" — bare number in brackets, for use inside a running sentence. */
  MV.eq = function (key) { return '(' + need(EQUATIONS, key, 'Gleichungs').num + ')'; };

  /** "Gl. (3.8)" */
  MV.eqn = function (key) { return 'Gl. ' + MV.eq(key); };

  /** "Gl. (3.8) der Arbeit" */
  MV.eqOfWork = function (key) { return MV.eqn(key) + ' der Arbeit'; };

  /** "Abbildung 3.2c" */
  MV.fig = function (key) { return 'Abbildung ' + need(FIGURES, key, 'Abbildungs').num; };

  /** "Abb. 3.2c" — short form, where the sentence is already crowded. */
  MV.figShort = function (key) { return 'Abb. ' + need(FIGURES, key, 'Abbildungs').num; };

  /** "Tabelle 5.1" */
  MV.tab = function (key) { return 'Tabelle ' + need(TABLES, key, 'Tabellen').num; };

  /** "Abschnitt 3.3" for one, "Abschnitte 3.2 und 3.5" for several. */
  MV.sec = function () {
    var nums = Array.prototype.slice.call(arguments);
    nums.forEach(function (n) {
      if (!SECTIONS[n]) throw new Error('Unbekannter Abschnitt: ' + n);
    });
    if (nums.length === 1) return 'Abschnitt ' + nums[0];
    return 'Abschnitte ' + nums.slice(0, -1).join(', ') + ' und ' + nums[nums.length - 1];
  };

  /** "Abschnitt 3.3 der Arbeit" */
  MV.secOfWork = function () {
    return MV.sec.apply(null, arguments) + ' der Arbeit';
  };

  /** Title of a section, for the leading quotation. */
  MV.secTitle = function (n) {
    if (!SECTIONS[n]) throw new Error('Unbekannter Abschnitt: ' + n);
    return SECTIONS[n];
  };

  /* Exposed for the generator and for the self-test. */
  MV.registry = { sections: SECTIONS, equations: EQUATIONS, figures: FIGURES, tables: TABLES };

}(MV));
