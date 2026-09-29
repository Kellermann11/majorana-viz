/* ------------------------------------------------------------------------
   energyscales.js — Reiter 5: Energieskalen, Abschnitt 5.1 der Arbeit.

   Panel A: Zielkonflikt und Bedingungsfenster        (5.3), (5.4), Abb. 5.2
   Panel B: Prüfung am Arbeitspunkt gegen Tabelle 5.1

   DIESER REITER RECHNET IN PHYSIKALISCHEN EINHEITEN (CLAUDE.md, Abschnitt 2):
   (5.4) ist ein Vergleich absoluter Skalen, deshalb Energien in µeV, Längen in
   µm, Temperatur in mK, Zeiten in µs. Dimensionslos bleiben die vier Reiter zu
   den Kapiteln 3 bis 5. Der Reiter zu Kapitel 6 trägt ebenfalls Einheiten —
   dort aber gemessene Werte, hier illustrative.

   Abbildung 5.2 der Arbeit ist laut eigener Bildunterschrift schematisch. Hier
   werden die Schnittpunkte aus den eingegebenen Werten tatsächlich gerechnet —
   das ist der Zweck dieses Reiters. Alle Vorgabewerte sind illustrativ und in
   NOTES.md mit Wert und Herkunft aufgeführt; einzige Zahl aus der Arbeit ist
   T = 20 mK.
   ------------------------------------------------------------------------ */

var MV = window.MV || (window.MV = {});

(function (MV) {
  'use strict';

  /* ====================================================================
     Naturkonstanten in den Einheiten dieses Reiters
     Beide stehen so in CLAUDE.md, Abschnitt 5.
     ==================================================================== */

  var KB_UEV_PER_K = 86.17;    /* k_B  = 86,17 µeV / K    */
  var HBAR_UEV_NS = 0.6582;    /* ℏ    = 0,6582 µeV · ns  */

  /** k_B T in µeV, aus der Temperatur in mK. */
  function kBT(T_mK) { return KB_UEV_PER_K * (T_mK / 1000); }

  /** ℏ / τ_meas in µeV, aus der Messzeit in µs. */
  function hbarOverTau(tau_us) { return HBAR_UEV_NS / (tau_us * 1000); }

  /* ====================================================================
     Die beiden Skalen aus (5.3)
     ==================================================================== */

  /** δE(L_W) = Δ_T · exp(−L_W/ξ), in µeV. */
  function splitting(L, s) { return s.deltaT * Math.exp(-L / s.xi); }

  /** E_C(L_W) = E_C^ref · (1 µm / L_W), in µeV. */
  function charging(L, s) { return s.ecRef / L; }

  /* ====================================================================
     Zustand — alle Vorgabewerte illustrativ, siehe NOTES.md
     ==================================================================== */

  var state = {
    deltaT: 100,     /* µeV  topologische Lücke Δ_T                        */
    xi: 0.1,         /* µm   Kohärenzlänge                                 */
    T: 20,           /* mK   einzige Zahl aus der Arbeit                   */
    tau: 1,          /* µs   Messzeit τ_meas                               */
    ecRef: 100,      /* µeV  Ladungsenergie bei L_W = 1 µm                 */
    factor: 10,      /* Schwelle für "≪"                                   */
    LW: 3,           /* µm   Arbeitspunkt                                  */
    TP: 1,           /* s    Paritätslebensdauer — Eingabe, nicht abgeleitet */
    floating: true   /* schwebend / geerdet                                */
  };

  /* Reglerbereiche, die im Text genannt werden — eine Quelle für beides. */
  var TAU_MIN = 0.01, TAU_MAX = 100;   /* Messzeit in µs      */
  var T_MIN = 5, T_MAX = 150;          /* Temperatur in mK    */

  /* T_P läuft über zwei Regler, die denselben Wert führen: der feine deckt
     den Bereich ab, in dem die gemessenen Lebensdauern liegen (Al 1–12 ms,
     Pb 22 s, siehe Reiter zu Kapitel 6), der logarithmische reicht acht
     Dekaden hinunter. Erst dadurch ist K4 = τ_meas/T_P überhaupt verletzbar:
     mit τ_meas ≤ 100 µs und T_P ≥ 1 ms ist das Verhältnis höchstens 0,1 und
     damit bei der Vorgabeschwelle immer erfüllt. Die Werte unterhalb einer
     Millisekunde sind illustrativ und liegen weit unter allem Gemessenen. */
  var TP_FINE_MIN = 0.001, TP_MAX = 10;   /* s, feiner Regler   */
  var TP_MIN = 1e-7;                      /* s, logarithmisch   */

  var L_MAX = 8;          /* µm, Ende der x-Achse   */
  var E_FLOOR = 1e-6;     /* µeV, untere Grenze der Log-Achse */
  var N_CURVE = 400;

  var sliders = {}, panelA = null, panelB = null, readout = null;
  var groundButton = null;

  /* Was eine Vorführung hervorhebt (demo.js). Gehört nicht zum Zustand des
     Reiters: es ändert keine Zahl und keine Bedingung, sondern nur, was
     umrandet ist — und verschwindet beim ersten Handeingriff. `holdK3` hält
     die K3-Zeile einen Schritt lang stehen, nachdem die Insel geerdet wurde,
     damit zu sehen ist, dass sie entfällt und nicht bloß verschwindet. */
  var mark = { rows: {}, holdK3: false };

  /* ====================================================================
     Fenster und Kriterien

     Die Schranken folgen aus (5.4) und den Prüfgrößen zu K2 und K3 in
     Tabelle 5.1, jeweils mit dem eingestellten Faktor f für "≪":

       K2a   δE ≤ k_B T / f          →  L_W ≥ ξ · ln( f · Δ_T / k_B T )
       K2b   δE · τ_meas / ℏ ≤ 1/f   →  L_W ≥ ξ · ln( f · Δ_T / (ℏ/τ_meas) )
       K3    k_B T ≤ E_C / f         →  L_W ≤ E_C^ref / ( f · k_B T )

     K4 (τ_meas ≪ T_P) hängt nicht von L_W ab und kann das Fenster deshalb
     nicht beschneiden; es wird getrennt ausgewiesen.
     ==================================================================== */

  function derive(s) {
    var kt = kBT(s.T);
    var ht = hbarOverTau(s.tau);
    var f = s.factor;

    /* untere Schranken; ohne Logarithmus definiert, falls Δ_T bereits unter
       der Schranke liegt — dann ist die Bedingung für jedes L_W erfüllt */
    function lowerBound(limit) {
      if (!(limit > 0)) return Infinity;
      var arg = f * s.deltaT / limit;
      return arg <= 1 ? 0 : s.xi * Math.log(arg);
    }

    var lK2a = lowerBound(kt);
    var lK2b = lowerBound(ht);
    var lK3 = s.floating ? s.ecRef / (f * kt) : Infinity;

    var lo = Math.max(lK2a, lK2b);
    var hi = lK3;
    var binding = lK2b > lK2a ? 'K2b' : 'K2a';

    /* Die Kohärenzlänge, bei der sich das Fenster gerade schließt: dort fällt
       die untere Schranke aus K2b mit der oberen aus K3 zusammen.
         ξ_zu = [E_C^ref / (f · k_BT)] / ln( f · Δ_T · τ_meas / ℏ )
       Nur sinnvoll, solange die Insel schwebt und der Logarithmus positiv ist. */
    var lnArg = ht > 0 ? f * s.deltaT / ht : 0;
    var xiClose = (s.floating && lnArg > 1)
      ? (s.ecRef / (f * kt)) / Math.log(lnArg)
      : null;

    return {
      s: s,
      kBT: kt,
      hbarOverTau: ht,
      lK2a: lK2a,
      lK2b: lK2b,
      xiClose: xiClose,
      lK3: lK3,
      binding: binding,
      lo: lo,
      hi: hi,
      /* Ein Fenster besteht nur, wenn es eine obere Schranke gibt und diese
         über der unteren liegt; im geerdeten Fall gibt es keine. */
      hasWindow: s.floating && isFinite(hi) && hi > lo,
      /* Sichtbar ist nur, was in den dargestellten Bereich fällt. */
      visibleLo: Math.max(lo, 0),
      visibleHi: Math.min(hi, L_MAX)
    };
  }

  /* Verhältnisse am Arbeitspunkt, je Kriterium aus Tabelle 5.1. */
  function ratios(d) {
    var s = d.s;
    var dE = splitting(s.LW, s);
    var eC = charging(s.LW, s);
    return {
      deltaE: dE,
      eC: eC,
      K2a: dE / d.kBT,
      K2b: dE * (s.tau * 1000) / HBAR_UEV_NS,   /* δE · τ_meas / ℏ, τ in ns */
      K3: s.floating ? d.kBT / eC : null,
      K4: (s.tau * 1e-6) / s.TP                 /* τ_meas in s, T_P in s   */
    };
  }

  /** erfüllt / grenzwertig / nicht erfüllt — nie allein über Farbe. */
  function verdict(ratio, f) {
    if (ratio === null || !isFinite(ratio)) return { text: 'entfällt', kind: 'trivial' };
    if (ratio <= 1 / f) return { text: 'erfüllt', kind: 'topological' };
    if (ratio < 1) return { text: 'grenzwertig', kind: 'critical' };
    return { text: 'nicht erfüllt', kind: 'trivial' };
  }

  /* ====================================================================
     Selbstprüfung

     Zwei Aussagen dieses Reiters sind stärker als „mit diesen Vorgabewerten":
     dass K2b über den ganzen Reglerbereich bindet, und dass sich das Fenster
     bei ξ_zu schließt. Beide werden hier nachgerechnet.
     ==================================================================== */

  function selfCheck() {
    var parts = [], i;

    /* (1) K2b bindet an jeder Ecke des Reglerraums.

       K2b ist die strengere Schranke genau dann, wenn ℏ/τ_meas < k_BT ist.
       Geprüft wird an allen Ecken, nicht nur an den beiden Extremen: die
       übrigen Regler dürfen daran nichts ändern, und dass sie es nicht tun,
       soll die Prüfung zeigen und nicht ich behaupten. */
    var corners = [], notK2b = [];
    [5, 400].forEach(function (deltaT) {
      [0.02, 0.8].forEach(function (xi) {
        [5, 400].forEach(function (ecRef) {
          [T_MIN, T_MAX].forEach(function (T) {
            [TAU_MIN, TAU_MAX].forEach(function (tau) {
              [1, 100].forEach(function (factor) {
                corners.push({ deltaT: deltaT, xi: xi, ecRef: ecRef, T: T,
                               tau: tau, factor: factor, LW: 3, TP: 1,
                               floating: true });
              });
            });
          });
        });
      });
    });
    for (i = 0; i < corners.length; i++) {
      var d = derive(corners[i]);
      if (d.binding !== 'K2b') {
        notK2b.push('τ=' + corners[i].tau + ' T=' + corners[i].T +
                    ' f=' + corners[i].factor);
      }
    }
    var htMax = hbarOverTau(TAU_MIN), ktMin = kBT(T_MIN);
    parts.push({
      name: 'K2b ist die bindende untere Schranke',
      ok: notK2b.length === 0,
      detail: notK2b.length
        ? 'nicht an: ' + notK2b.slice(0, 4).join('; ')
        : 'an allen ' + corners.length + ' Ecken des Reglerraums. Grund: ' +
          'ℏ/τ_meas ist höchstens ' + MV.fmt(htMax, 4) + ' µeV, k_BT mindestens ' +
          MV.fmt(ktMin, 4) + ' µeV — die kleinere Vergleichsgröße ergibt die ' +
          'strengere Schranke, und das ist über den ganzen Bereich ℏ/τ_meas'
    });

    /* (2) Bei ξ = ξ_zu schließt sich das Fenster gerade. */
    var ref = { deltaT: 100, xi: 0.1, ecRef: 100, T: 20, tau: 1,
                LW: 3, TP: 1, factor: 10, floating: true };
    var dRef = derive(ref);
    var atClose = derive(Object.assign({}, ref, { xi: dRef.xiClose }));
    var gap = Math.abs(atClose.lo - atClose.hi);
    parts.push({
      name: 'Bei ξ_zu fallen untere und obere Schranke zusammen',
      ok: gap < 1e-9,
      detail: 'ξ_zu = ' + MV.fmt(dRef.xiClose, 4) + ' µm; dort ist die untere ' +
        'Schranke ' + MV.fmt(atClose.lo, 4) + ' µm und die obere ' +
        MV.fmt(atClose.hi, 4) + ' µm, Abstand ' + MV.fmtExp(gap, 1) + ' µm',
      dev: gap
    });

    /* (3) Die Abnahmewerte aus Auftrag 17, Teil 3. */
    var bad = [];
    function expect(name, got, soll, tol) {
      if (!(Math.abs(got - soll) <= tol)) {
        bad.push(name + ' = ' + MV.fmt(got, 3) + ' statt ' + MV.fmt(soll, 3));
      }
    }
    expect('K2a', dRef.lK2a, 0.636, 5e-3);
    expect('K2b', dRef.lK2b, 1.423, 5e-3);
    expect('K3', dRef.hi, 5.802, 5e-3);
    var d1 = derive(Object.assign({}, ref, { factor: 1 }));
    expect('untere Grenze bei f = 1', d1.lo, 1.193, 5e-3);
    expect('obere Grenze bei f = 1', d1.hi, 58.025, 5e-2);
    parts.push({
      name: 'Abnahmewerte des Auftrags',
      ok: bad.length === 0,
      detail: bad.length ? bad.join('; ')
        : 'Vorgabewerte: K2a ab ' + MV.fmt(dRef.lK2a, 3) + ' µm, K2b ab ' +
          MV.fmt(dRef.lK2b, 3) + ' µm, K3 bis ' + MV.fmt(dRef.hi, 3) +
          ' µm; bei f = 1 ' + MV.fmtRange(d1.lo, d1.hi, 3, 'µm')
    });

    /* (4) K4 muss beide Urteile erreichen können.

       Der Grund für den zweiten, logarithmischen T_P-Regler: mit τ_meas ≤ 100 µs
       und T_P ≥ 1 ms ist τ_meas/T_P höchstens 0,1 — bei der Vorgabeschwelle
       genau die Grenze, und der Vergleich ist „≤". K4 wäre damit an „erfüllt"
       festgenagelt, und ein Kriterium, das nicht scheitern kann, prüft nichts.
       Beide Aussagen werden hier nachgerechnet, statt behauptet. */
    function k4At(tau, TP) { return (tau * 1e-6) / TP; }
    var fineWorst = k4At(TAU_MAX, TP_FINE_MIN);
    var wideWorst = k4At(TAU_MAX, TP_MIN);
    var atDefault = k4At(1, TP_MIN);
    parts.push({
      name: 'K4 ist im feinen Reglerbereich allein nicht verletzbar',
      ok: verdict(fineWorst, 10).text === 'erfüllt' &&
          verdict(wideWorst, 10).text === 'nicht erfüllt',
      detail: 'größtes K4 im feinen Bereich (τ_meas = ' + MV.fmt(TAU_MAX, 0) +
        ' µs, T_P = ' + fmtTP(TP_FINE_MIN) + '): ' + MV.fmt(fineWorst, 4) +
        ' — „' + verdict(fineWorst, 10).text + '". Mit dem logarithmischen ' +
        'Regler bis ' + fmtTP(TP_MIN) + ': ' + MV.fmtExp(wideWorst, 1) + ' — „' +
        verdict(wideWorst, 10).text + '"'
    });
    parts.push({
      name: 'K4 kippt bei den Vorgabewerten innerhalb des Reglerbereichs',
      ok: verdict(atDefault, 10).text === 'nicht erfüllt' &&
          verdict(k4At(1, 1), 10).text === 'erfüllt',
      detail: 'bei τ_meas = 1 µs: T_P = ' + fmtTP(1) + ' ergibt ' +
        MV.fmtExp(k4At(1, 1), 1) + ' („erfüllt"), T_P = ' + fmtTP(TP_MIN) +
        ' ergibt ' + MV.fmt(atDefault, 1) + ' („nicht erfüllt"); die Grenze ' +
        'liegt bei T_P = τ_meas'
    });

    return {
      parts: parts,
      passed: parts.every(function (p) { return p.ok !== false; }),
      xiClose: dRef.xiClose,
      corners: corners.length
    };
  }

  MV.energyscalesSelfCheck = selfCheck;

  var check = selfCheck();

  function checkHtml() {
    return '<ul>' + check.parts.map(function (p) {
      return '<li><b>' + (p.ok === true ? '✓' : p.ok === false ? '✗' : '!') +
        '</b> ' + p.name + ' — ' + p.detail + '</li>';
    }).join('') + '</ul>';
  }

  /* ====================================================================
     Panel A — Zielkonflikt   (Abb. 5.2)
     ==================================================================== */

  function log10(v) { return Math.log(v) / Math.LN10; }

  /**
   * Bei welchem ξ sich das Fenster schließt — als Satz, Zahl aus derive().
   * Der Regler bleibt der Weg dorthin; genannt wird nur, wo es passiert.
   */
  function xiCloseText() {
    var d = derive(state);
    if (d.xiClose === null || !(d.xiClose > 0)) {
      return 'Ob und wo sich das Fenster schließt, hängt an den übrigen Reglern.';
    }
    return 'bei etwa ' + MV.fmt(d.xiClose, 2) + ' µm schließt es sich, denn dort ' +
      'fällt die untere Schranke aus K2b mit der oberen aus K3 zusammen.';
  }

  /**
   * Warum K2b strukturell bindet: innerhalb der Reglerbereiche ist ℏ/τ_meas
   * stets kleiner als k_BT, und die kleinere Vergleichsgröße ergibt die
   * strengere Schranke. Beide Grenzwerte aus den Reglerecken gerechnet.
   */
  function bindingText() {
    var htMax = hbarOverTau(TAU_MIN);
    var ktMin = kBT(T_MIN);
    return 'Innerhalb der Reglerbereiche ist ℏ/τ_meas höchstens ' +
      MV.fmt(htMax, 3) + ' µeV (bei τ_meas = ' + MV.fmt(TAU_MIN, 2) + ' µs) und ' +
      'k_B T mindestens ' + MV.fmt(ktMin, 2) + ' µeV (bei T = ' + MV.fmt(T_MIN, 0) +
      ' mK). ℏ/τ_meas ist damit <b>immer</b> die kleinere der beiden Größen, und ' +
      'K2b ist immer die strengere Schranke — das ist keine Eigenheit der ' +
      'Vorgabewerte, sondern gilt über den ganzen einstellbaren Bereich.';
  }

  function drawPanelA() {
    var d = derive(state), s = state, i, L;

    var top = Math.ceil(log10(Math.max(s.deltaT, 2 * s.ecRef, d.kBT, d.hbarOverTau)));
    var bottom = log10(E_FLOOR);

    /* x-Marken in µm, darunter dieselben Stellen in L_W/ξ */
    var xTicks = [], xTicks2 = [];
    for (i = 0; i <= 8; i++) {
      xTicks.push({ v: i, label: String(i) });
      xTicks2.push({ v: i, label: MV.fmt(i / s.xi, 0) });
    }

    var p = panelA.plot;
    p.begin({
      x: [0, L_MAX], y: [bottom, top],
      xTicks: xTicks, xTicks2: xTicks2,
      yTicks: MV.decadeTicks(bottom, top, 10),
      xLabel: 'Drahtlänge L_W — obere Marken in µm, untere in L_W/ξ',
      yLabel: 'Energie in µeV'
    });

    /* ---- Bedingungsfenster ---- */
    if (d.hasWindow && d.visibleHi > d.visibleLo) {
      p.vspan(d.visibleLo, d.visibleHi, {
        fill: MV.palette.topologicalSoft,
        edge: MV.palette.topological,
        dash: [4, 3],
        legend: 'Fenster: alle Bedingungen mit Faktor ' + MV.fmt(s.factor, 0) + ' erfüllt'
      });
    }

    /* ---- die vier Skalen ---- */
    var pts, val, clipped = false;

    pts = [];
    for (i = 0; i <= N_CURVE; i++) {
      L = L_MAX * i / N_CURVE;
      val = splitting(L, s);
      if (val < E_FLOOR) { clipped = true; continue; }
      pts.push([L, log10(val)]);
    }
    p.curve(pts, {
      color: MV.palette.topological, width: 2.2,
      legend: 'δE = Δ_T · exp(−L_W/ξ)'
    });

    pts = [];
    for (i = 1; i <= N_CURVE; i++) {
      L = L_MAX * i / N_CURVE;
      pts.push([L, log10(charging(L, s))]);
    }
    p.curve(pts, {
      color: s.floating ? MV.palette.curve1 : MV.palette.line,
      width: s.floating ? 2.2 : 1.4,
      dash: s.floating ? [] : [5, 4],
      alpha: s.floating ? 1 : 0.7,
      legend: s.floating ? 'E_C = E_C^ref · (1 µm / L_W)'
                         : 'E_C — geerdet, ohne Wirkung'
    });

    p.hline(log10(d.kBT), {
      color: MV.palette.critical, dash: [7, 4], width: 1.6,
      legend: 'k_B T'
    });
    p.hline(log10(d.hbarOverTau), {
      color: MV.palette.curve3, dash: [2, 3], width: 1.6,
      legend: 'ℏ / τ_meas'
    });

    /* ---- Schwellenlinien ----

       „≪" ist im ganzen Reiter als „um den Faktor f kleiner" gerechnet. Die
       Fenstergrenzen liegen deshalb nicht auf den Skalen selbst, sondern auf
       den um f verschobenen Linien: δE muss unter k_BT/f und unter
       (ℏ/τ_meas)/f bleiben, E_C über f·k_BT. Bei f = 1 fallen beide zusammen.
       Es sind dieselben Ungleichungen wie in derive() — keine neue Physik. */
    if (s.factor !== 1) {
      p.hline(log10(d.kBT / s.factor), {
        color: MV.palette.critical, dash: [1, 3], width: 1.2, alpha: 0.85,
        legend: 'Schwelle k_B T / ' + MV.fmt(s.factor, 0) + ' (K2a)'
      });
      p.hline(log10(d.hbarOverTau / s.factor), {
        color: MV.palette.curve3, dash: [1, 3], width: 1.2, alpha: 0.85,
        legend: 'Schwelle (ℏ/τ_meas) / ' + MV.fmt(s.factor, 0) + ' (K2b)'
      });
      if (s.floating) {
        p.hline(log10(d.kBT * s.factor), {
          color: MV.palette.curve1, dash: [1, 3], width: 1.2, alpha: 0.85,
          legend: 'Schwelle ' + MV.fmt(s.factor, 0) + ' · k_B T (K3)'
        });
      }
    }

    /* ---- Arbeitspunkt ---- */
    p.vline(s.LW, { color: MV.palette.ink, width: 1.5, legend: 'Arbeitspunkt L_W' });

    /* ---- Randbeschriftungen wie in Abb. 5.2 ---- */
    if (d.hasWindow) {
      if (d.visibleLo > 0.35) {
        p.label('zu kurz:', p.X(d.visibleLo) - 8, p.py0 + 16,
          { color: MV.palette.muted, align: 'right', box: true });
        p.label('Randmoden überlappen', p.X(d.visibleLo) - 8, p.py0 + 30,
          { color: MV.palette.muted, align: 'right', box: true });
      }
      if (d.visibleHi < L_MAX - 0.35) {
        p.label('zu lang:', p.X(d.visibleHi) + 8, p.py0 + 16,
          { color: MV.palette.muted, align: 'left', box: true });
        p.label('E_C unter der Störskala', p.X(d.visibleHi) + 8, p.py0 + 30,
          { color: MV.palette.muted, align: 'left', box: true });
      }
    }

    p.finish();

    /* ---- Zeile unter der Abbildung ---- */
    var msg, kind;
    if (!s.floating) {
      msg = 'Geerdet: ohne Ladungsenergie fehlt die obere Schranke aus ' +
            MV.eq('energy.window') + ' — K3 aus ' + MV.tab('criteria') +
            ' ist nicht erfüllbar, und es gibt kein Fenster.';
      kind = MV.palette.critical;
    } else if (!d.hasWindow) {
      msg = 'Kein Fenster: die untere Schranke (' + MV.fmt(d.lo, 3) +
            ' µm, bindend ' + d.binding + ') liegt über der oberen (' +
            MV.fmt(d.hi, 3) + ' µm). Mit diesen Werten ist ' +
            MV.eq('energy.window') + ' bei keiner Drahtlänge erfüllbar.';
      kind = MV.palette.critical;
    } else {
      msg = 'Fenster ' + MV.fmtRange(d.lo, d.hi, 3, 'µm') + ', ' +
            'entspricht L_W/ξ ' + MV.fmtRange(d.lo / s.xi, d.hi / s.xi, 1) +
            '. Untere Schranke bindend durch ' +
            d.binding + '.';
      kind = MV.palette.topological;
      if (clipped) {
        msg += ' δE fällt ab L_W = ' + MV.fmt(s.xi * Math.log(s.deltaT / E_FLOOR), 2) +
               ' µm unter die Achsengrenze von 10⁻⁶ µeV und ist dort abgeschnitten.';
      }
    }
    panelA.setStatus(msg, kind);
  }

  /* ====================================================================
     Panel B — Prüfung am Arbeitspunkt gegen Tabelle 5.1
     ==================================================================== */

  /**
   * T_P mit der Einheit, die zum Wert passt. Die Zahl läuft durch MV.fmt,
   * also mit Dezimalkomma; gewählt wird nur die Einheit.
   */
  function fmtTP(v) {
    if (!isFinite(v) || v <= 0) return '—';
    if (v < 1e-3) return MV.fmt(v * 1e6, v * 1e6 < 10 ? 2 : 1) + ' µs';
    if (v < 1) return MV.fmt(v * 1e3, v * 1e3 < 10 ? 2 : 1) + ' ms';
    return MV.fmt(v, 2) + ' s';
  }

  function fmtE(v) {
    if (v === null || !isFinite(v)) return '—';
    if (v === 0) return '0';
    if (Math.abs(v) >= 1e-3 && Math.abs(v) < 1e5) return MV.fmt(v, 4);
    return MV.fmtExp(v, 2);
  }

  function renderPanelB() {
    var d = derive(state), r = ratios(d), s = state;
    var f = s.factor;
    var html = '';

    html += '<table class="numtable"><caption>Die vier Skalen am Arbeitspunkt ' +
      'L_W = ' + MV.fmt(s.LW, 2) + ' µm</caption>' +
      '<thead><tr><th>Größe</th><th class="num">Wert in µeV</th><th>Herkunft</th>' +
      '</tr></thead><tbody>' +
      '<tr><td>δE = Δ_T · exp(−L_W/ξ)</td><td class="num">' + fmtE(r.deltaE) +
      '</td><td>' + MV.eq('energy.tradeoff') + ', aus ' + MV.eq('kitaev.splitting') + '</td></tr>' +
      '<tr><td>E_C = E_C^ref · (1 µm / L_W)</td><td class="num">' +
      (s.floating ? fmtE(r.eC) : '—') + '</td><td>' + MV.eq('energy.charging') + ', ' +
      MV.eq('energy.tradeoff') + (s.floating ? '' : ' — geerdet') + '</td></tr>' +
      '<tr><td>k_B T</td><td class="num">' + fmtE(d.kBT) +
      '</td><td>T = ' + MV.fmt(s.T, 0) + ' mK, k_B = 86,17 µeV/K</td></tr>' +
      '<tr><td>ℏ / τ_meas</td><td class="num">' + fmtE(d.hbarOverTau) +
      '</td><td>τ_meas = ' + MV.fmt(s.tau, 2) + ' µs, ℏ = 0,6582 µeV·ns</td></tr>' +
      '</tbody></table>';

    var rows = [
      {
        k: 'K2a', expr: 'δE / k_B T', val: r.K2a,
        cond: 'Die Restaufspaltung ist vernachlässigbar — thermisch.',
        bound: d.lK2a
      },
      {
        k: 'K2b', expr: 'δE · τ_meas / ℏ', val: r.K2b,
        cond: 'Die Restaufspaltung ist vernachlässigbar — dynamische Phase.',
        bound: d.lK2b
      },
      {
        k: 'K4', expr: 'τ_meas / T_P', val: r.K4,
        cond: 'Die Parität übersteht viele Operationen.',
        bound: null
      }
    ];

    /* Bei geerdeter Insel entfällt die K3-Zeile ganz, statt sie mit einem
       Ersatzwert zu füllen: ohne Ladungsenergie gibt es die Prüfgröße nicht. */
    if (s.floating) {
      rows.splice(2, 0, {
        k: 'K3', expr: 'k_B T / E_C', val: r.K3,
        cond: 'Die Ladungsenergie wirkt als Paritätsblockade.',
        bound: null
      });
    } else if (mark.holdK3) {
      /* Der Halteschritt einer Vorführung: dieselbe Zeile, aber ohne Wert.
         Wert und Status kommen aus denselben Funktionen wie sonst — val = null
         ergibt „—" und „entfällt". */
      rows.splice(2, 0, {
        k: 'K3', expr: 'k_B T / E_C', val: null,
        cond: 'Entfällt — ohne Ladungsenergie gibt es die Prüfgröße nicht.',
        bound: null
      });
    }

    html += '<table class="numtable"><caption>Die vier Verhältnisse, ' +
      'Kriterien nach ' + MV.tab('criteria') + ' — als erfüllt gilt ein ' +
      'Verhältnis unter 1/' + MV.fmt(f, 0) + ' = ' + MV.fmt(1 / f, 4) +
      (s.floating ? '' : ' — bei geerdeter Insel fehlt K3 ganz, weil es ohne ' +
        'Ladungsenergie keine Prüfgröße gibt') +
      '</caption>' +
      '<thead><tr><th>Kriterium</th><th>Prüfgröße</th><th class="num">Wert</th>' +
      '<th>Status</th></tr></thead><tbody>';

    rows.forEach(function (row) {
      var v = verdict(row.val, f);
      var isBinding = (row.k === d.binding);
      var cls = [];
      if (isBinding) cls.push('numtable__binding');
      if (mark.rows[row.k]) cls.push('numtable__marked');
      html += '<tr' + (cls.length ? ' class="' + cls.join(' ') + '"' : '') + '>' +
        '<td><b>' + row.k + '</b>' +
        (isBinding ? ' <span class="tag">bindend</span>' : '') +
        '<br><span class="numtable__sub">' + row.cond + '</span></td>' +
        '<td>' + row.expr + ' ≪ 1</td>' +
        '<td class="num">' + fmtE(row.val) + '</td>' +
        '<td><span class="verdict verdict--' + v.kind + '">' + v.text + '</span>' +
        (row.bound !== null && isFinite(row.bound)
          ? '<br><span class="numtable__sub">verlangt L_W ≥ ' + MV.fmt(row.bound, 3) + ' µm</span>'
          : '') +
        '</td></tr>';
    });
    html += '</tbody></table>';

    html += '<p class="note"><b>T_P ist eine Eingabegröße</b>, kein aus ' +
      MV.eq('energy.charging') + ' bis ' + MV.eq('energy.window') + ' abgeleiteter ' +
      'Wert. Die Paritätslebensdauer folgt nicht aus der Ladungsenergie; ' +
      MV.secOfWork('5.6') + ' hält ausdrücklich fest, dass beides ' +
      'auseinanderfallen kann. Sie wird hier gesetzt, nicht berechnet.</p>';

    panelB.body.innerHTML = html;

    var bindRow = d.binding === 'K2b' ? 'K2b, die dynamische Phase' : 'K2a, die thermische Schranke';
    panelB.setStatus(
      state.floating
        ? 'Bindend für die untere Schranke ist ' + bindRow + ': ' +
          'K2a verlangt L_W ≥ ' + MV.fmt(d.lK2a, 3) + ' µm, K2b verlangt L_W ≥ ' +
          MV.fmt(d.lK2b, 3) + ' µm.'
        : 'Geerdet: K3 entfällt, damit gibt es keine obere Schranke an L_W. ' +
          'Bindend unten bleibt ' + bindRow + '.',
      state.floating ? MV.palette.ink : MV.palette.critical
    );
  }

  /* ====================================================================
     Anzeigefeld und Aktualisierung
     ==================================================================== */

  function refresh() {
    var d = derive(state), r = ratios(d);

    readout.set('kbt', MV.fmt(d.kBT, 3) + ' µeV');
    readout.set('hbartau', MV.fmtExp(d.hbarOverTau, 2) + ' µeV');
    readout.set('window', d.hasWindow
      ? MV.fmtRange(d.lo, d.hi, 2, 'µm')
      : 'keines');
    readout.set('windowXi', d.hasWindow
      ? MV.fmtRange(d.lo / state.xi, d.hi / state.xi, 1)
      : '—');
    readout.set('binding', d.binding);
    readout.set('deltaE', fmtE(r.deltaE) + ' µeV');

    var inside = d.hasWindow && state.LW >= d.lo && state.LW <= d.hi;
    readout.status(
      !state.floating ? 'critical' : (inside ? 'topological' : 'trivial'),
      !state.floating
        ? 'Geerdet — K3 nicht erfüllbar, kein Fenster'
        : (d.hasWindow
            ? (inside ? 'Arbeitspunkt liegt im Fenster' : 'Arbeitspunkt liegt außerhalb des Fensters')
            : 'Kein Fenster bei diesen Werten')
    );

    drawPanelA();
    renderPanelB();
    panelA.info.update();
    panelB.info.update();
  }

  /* ====================================================================
     Info-Panels
     ==================================================================== */

  var UNITS_NOTE =
    'Dieser Reiter <b>rechnet in physikalischen Einheiten</b>: ' +
    MV.eq('energy.window') + ' vergleicht absolute Skalen miteinander, deshalb ' +
    'Energien in µeV, Längen in µm, Temperatur in mK und Zeiten in µs. ' +
    'Dimensionslos bleiben die vier Reiter zu den Kapiteln 3 bis 5 — ' +
    'Kitaev-Kette, Rashba-Nanodraht, Auslese und Gatter. Der Reiter zu ' +
    'Kapitel 6 trägt ebenfalls Einheiten, dort aber gemessene Werte statt ' +
    'illustrativer.';

  var VALUES_NOTE =
    '<b>Alle Vorgabewerte sind illustrativ</b>, nicht aus einer Veröffentlichung ' +
    'übernommen. ' + MV.fig('energy.windowFig') + ' der Arbeit ist laut eigener ' +
    'Bildunterschrift schematisch; hier werden die Schnittpunkte aus den ' +
    'eingegebenen Werten tatsächlich gerechnet, aber die Eingaben sind runde ' +
    'Zahlen zur Veranschaulichung. Einzige Ausnahme ist T = 20 mK, das so in der ' +
    'Arbeit steht. Jeder andere Vorgabewert ist frei gewählt und trägt keine ' +
    'Aussage über ein bestimmtes Bauelement.';

  var PHI_NOTE =
    'φ kommt in diesem Reiter nicht vor. In den Reitern zu Kapitel 3 und 4 ' +
    'bezeichnet φ die supraleitende Phase, im Tetron-Reiter den eingeschlossenen ' +
    'Fluss.';

  function paramLine(tail) {
    var d = derive(state), s = state;
    return '<p>Δ_T = <span class="num">' + MV.fmt(s.deltaT, 1) + '</span> µeV, ' +
      'ξ = <span class="num">' + MV.fmt(s.xi, 3) + '</span> µm, ' +
      'T = <span class="num">' + MV.fmt(s.T, 1) + '</span> mK, ' +
      'τ_meas = <span class="num">' + MV.fmt(s.tau, 2) + '</span> µs, ' +
      'E_C^ref = <span class="num">' + MV.fmt(s.ecRef, 1) + '</span> µeV bei L_W = 1 µm, ' +
      'T_P = <span class="num">' + fmtTP(s.TP) + '</span>, ' +
      'Schwelle für ≪ = <span class="num">' + MV.fmt(s.factor, 0) + '</span>, ' +
      'Insel <b>' + (s.floating ? 'schwebend' : 'geerdet') + '</b>.</p>' +
      '<p>Daraus: k_B T = <span class="num">' + MV.fmt(d.kBT, 4) + '</span> µeV, ' +
      'ℏ/τ_meas = <span class="num">' + MV.fmtExp(d.hbarOverTau, 3) +
      '</span> µeV, Arbeitspunkt L_W = <span class="num">' + MV.fmt(s.LW, 2) +
      '</span> µm.</p>' +
      (tail ? tail(d) : '');
  }

  var INFO_A = {
    what:
      '<p>Die vier Energieskalen des Qubits über der Drahtlänge L_W, ' +
      'logarithmisch aufgetragen. Die untere Achsenbeschriftung nennt dieselbe ' +
      'Länge in Einheiten der Kohärenzlänge, weil ' + MV.fig('energy.windowFig') +
      ' der Arbeit so skaliert ist. Der hinterlegte Streifen ist der Bereich, in ' +
      'dem alle Bedingungen gleichzeitig erfüllt sind; die senkrechte Linie ist ' +
      'der Arbeitspunkt aus Panel B.</p>',
    model:
      '<p>Coulomb-blockierte Insel, ' + MV.secOfWork('5.1') + '. Ladungsanteil des ' +
      'Hamiltonians, ' + MV.eqn('energy.charging') + ':</p>' +
      '<span class="eq">H_C = E_C ( N̂ − n_g )²,     E_C = e² / 2C</span>' +
      '<p>Die beiden maßgeblichen Skalen hängen gegenläufig von der Drahtlänge ab, ' +
      MV.eqn('energy.tradeoff') + ':</p>' +
      '<span class="eq">E_C ~ 1 / L_W          Kapazität wächst mit der Ausdehnung\n' +
      'δE  ~ Δ_T · exp( −L_W / ξ )   Aufspaltung nach ' + MV.eq('kitaev.splitting') + '</span>' +
      '<p>Beide fallen mit L_W, aber unterschiedlich schnell — daraus entsteht ein ' +
      'Fenster, kein einseitiges Optimum.</p>',
    computed:
      '<p>Aufgetragen sind die vier Größen</p>' +
      '<span class="eq">δE(L_W)  = Δ_T · exp( −L_W / ξ )\n' +
      'E_C(L_W) = E_C^ref · ( 1 µm / L_W )\n' +
      'k_B T\n' +
      'ℏ / τ_meas</span>' +
      '<p>Das Bedingungsfenster folgt aus ' + MV.eqn('energy.window') +
      ' zusammen mit den Prüfgrößen zu K2 und K3 aus ' + MV.tab('criteria') + '. ' +
      'Mit dem Faktor f für „≪" lauten die Schranken geschlossen:</p>' +
      '<span class="eq">K2a   δE ≤ k_B T / f         →  L_W ≥ ξ · ln( f · Δ_T / k_B T )\n' +
      'K2b   δE·τ_meas/ℏ ≤ 1/f     →  L_W ≥ ξ · ln( f · Δ_T / (ℏ/τ_meas) )\n' +
      'K3    k_B T ≤ E_C / f       →  L_W ≤ E_C^ref / ( f · k_B T )</span>' +
      '<p>Die Fenstergrenzen sind also nicht abgelesen, sondern aus den ' +
      'eingegebenen Werten berechnet. Sie liegen auf den Schnittpunkten mit den ' +
      '<b>gepunkteten Schwellenlinien</b> — δE mit k_B T / f und mit ' +
      '(ℏ/τ_meas) / f, E_C mit f · k_B T. Bei f = 1 fallen diese Linien mit den ' +
      'Skalen selbst zusammen und werden dann nicht eigens gezeichnet; bei ' +
      'größerem f liegen die Grenzen entsprechend weiter innen. K4 hängt nicht ' +
      'von L_W ab und kann das Fenster nicht beschneiden.</p>',
    params: function () {
      return paramLine(function (d) {
        if (!state.floating) {
          return '<p>Geerdet: ohne Ladungsenergie existiert die obere Schranke aus ' +
            MV.eq('energy.window') + ' nicht.</p>';
        }
        return '<p>Schranken: K2a ab <span class="num">' + MV.fmt(d.lK2a, 4) +
          '</span> µm, K2b ab <span class="num">' + MV.fmt(d.lK2b, 4) +
          '</span> µm, K3 bis <span class="num">' + MV.fmt(d.lK3, 4) + '</span> µm. ' +
          (d.hasWindow
            ? 'Fenster <span class="num">' + MV.fmtRange(d.lo, d.hi, 4) +
              '</span> µm, in Einheiten von ξ <span class="num">' +
              MV.fmtRange(d.lo / state.xi, d.hi / state.xi, 2) + '</span>.'
            : 'Es gibt kein Fenster.') +
          '</p>';
      });
    },
    numerics: checkHtml() +
      '<ul>' +
      '<li>' + N_CURVE + ' Stützstellen je Kurve über L_W ∈ [0, ' + L_MAX + '] µm. ' +
      'Die Fenstergrenzen werden nicht am Gitter gesucht, sondern aus den oben ' +
      'angegebenen geschlossenen Ausdrücken berechnet.</li>' +
      '<li><b>Die Schwelle für „≪" ist einstellbar und wird nicht still ' +
      'festgelegt.</b> ' + MV.eqn('energy.window') + ' schreibt „≪", nicht „&lt;" — ' +
      'ein Ungleichheitszeichen ohne Zahl. Wie viel kleiner „viel kleiner" heißt, ' +
      'ist eine Festlegung des Lesers und verschiebt die Fenstergrenzen erheblich. ' +
      'Der Regler macht diese Festlegung sichtbar; die Vorgabe 10 ist eine runde ' +
      'Wahl, keine aus der Arbeit übernommene. Auf 1 gestellt wird aus „≪" ein ' +
      'bloßes „&lt;", und das Fenster wird maximal weit.</li>' +
      '<li>Naturkonstanten in den Einheiten dieses Reiters: k_B = 86,17 µeV/K, ' +
      'also k_B T [µeV] = 0,08617 · T [mK] / 1000 · 1000; und ℏ = 0,6582 µeV·ns, ' +
      'also ℏ/τ_meas [µeV] = 6,582·10⁻⁴ / τ_meas [µs]. Bei T = 20 mK folgt ' +
      'k_B T = 1,723 µeV, was dem in der Arbeit genannten Wert von rund 1,7 µeV ' +
      'entspricht.</li>' +
      '<li>Die Höhenachse endet unten bei 10⁻⁶ µeV. Fällt δE darunter, wird die ' +
      'Kurve sichtbar abgeschnitten und der Umstand unter der Abbildung vermerkt — ' +
      'sie wird nicht auf der Achsengrenze weitergeführt.</li>' +
      '<li>Existiert kein Fenster, bleibt der Streifen weg und der Grund steht als ' +
      'Text unter der Abbildung.</li>' +
      '</ul>',
    convention:
      '<p>' + UNITS_NOTE + '</p>' +
      '<p>' + VALUES_NOTE + '</p>' +
      '<p>' + PHI_NOTE + '</p>' +
      '<p>Der Schalter schwebend/geerdet illustriert das Kriterium K3 aus ' +
      MV.tab('criteria') + ' — die Ladungsenergie wirkt nur als Paritätsblockade, ' +
      'wenn die Insel schwebend betrieben wird. Er bezeichnet kein bestimmtes ' +
      'Bauelement.</p>',
    reference:
      '<p>' + MV.secOfWork('5.1') + ', ' + MV.fig('energy.windowFig') + '. ' +
      'Verwendete Gleichungen: ' + MV.eq('energy.charging') + ', ' +
      MV.eq('energy.tradeoff') + ', ' + MV.eq('energy.window') + ', sowie ' +
      MV.eq('kitaev.splitting') + ' aus ' + MV.sec('3.5') + '. Kriterien nach ' +
      MV.tab('criteria') + '.</p>',
    reading:
      '<p><b>Ziehe ξ von ' + MV.fmt(0.1, 1) + ' auf ' + MV.fmt(0.45, 2) + ' µm</b> und ' +
      'beobachte, wie die δE-Gerade flacher wird und die linke Fenstergrenze nach ' +
      'rechts wandert: ' + xiCloseText() + ' <b>Ziehe danach E_C^ref von 100 auf ' +
      '30 µeV</b>: die E_C-Kurve sinkt, ' +
      'die rechte Grenze wandert nach links, und das Fenster wird von der anderen ' +
      'Seite zusammengedrückt. Das ist der Zielkonflikt aus ' +
      MV.eq('energy.tradeoff') + ' — ein längerer Draht hilft der einen Bedingung ' +
      'und schadet der anderen.</p>' +
      '<p><b>Stelle die Schwelle von 10 auf 1</b>: das Fenster wird sichtbar breiter, ' +
      'weil aus „viel kleiner" ein bloßes „kleiner" wird. Die Grenzen sind keine ' +
      'Naturkonstanten, sondern hängen daran, wie streng man „≪" liest.</p>' +
      '<p><b>Schalte auf „geerdet"</b>: die E_C-Kurve wird ausgegraut und der ' +
      'Streifen verschwindet. Ohne Ladungsenergie gibt es keine obere Schranke, ' +
      'und K3 ist nicht erfüllbar.</p>' +
      '<p>Beachte zuletzt, wo das Fenster liegt: mit den Vorgabewerten beginnt es ' +
      'bei L_W/ξ ≈ 14. Für L_W/ξ jenseits von etwa 15 sagt ' +
      MV.eq('energy.tradeoff') + ' eine Aufspaltung voraus, die weit unterhalb ' +
      'jeder erreichbaren Messauflösung liegt — δE ist an der linken Fenstergrenze ' +
      'bereits kleiner als 10⁻⁴ µeV. Das ist eine Aussage über die Formel, kein ' +
      'Darstellungsfehler: das Fenster ist nicht dadurch definiert, dass δE messbar ' +
      'wäre, sondern dadurch, dass es vernachlässigbar ist.</p>'
  };

  var INFO_B = {
    what:
      '<p>Die vier Energieskalen am eingestellten Arbeitspunkt L_W als Zahlen, ' +
      'darunter die vier Verhältnisse, die nach ' + MV.tab('criteria') + ' klein ' +
      'gegen eins sein müssen. Je Zeile steht das Kriterium, die Prüfgröße, der ' +
      'Wert und ein Status in Worten. Die Zeile, die die strengere untere Schranke ' +
      'an L_W erzeugt, ist als bindend gekennzeichnet.</p>',
    model:
      '<p>Dieselben Skalen wie in Panel A, ' + MV.secOfWork('5.1') + '. Die ' +
      'Kriterien und ihre Prüfgrößen stehen in ' + MV.tab('criteria') + ':</p>' +
      '<span class="eq">K2   δE ≪ k_B T  und  δE τ_meas / ℏ ≪ 1\n' +
      'K3   E_C ≫ k_B T\n' +
      'K4   T_P ≫ τ_meas</span>' +
      '<p>K2 und K3 zusammen sind ' + MV.eq('energy.window') + '; K4 stammt aus ' +
      MV.sec('5.1') + ' und ' + MV.sec('5.5') + '.</p>',
    computed:
      '<span class="eq">K2a   δE / k_B T\n' +
      'K2b   δE · τ_meas / ℏ\n' +
      'K3    k_B T / E_C\n' +
      'K4    τ_meas / T_P</span>' +
      '<p>Alle vier sind dimensionslose Verhältnisse und sollen klein gegen eins ' +
      'sein. Als erfüllt gilt ein Wert unterhalb von 1/f mit der eingestellten ' +
      'Schwelle f, als grenzwertig ein Wert zwischen 1/f und 1, als nicht erfüllt ' +
      'ein Wert ab 1.</p>' +
      '<p>K2a und K2b erzeugen beide eine untere Schranke an L_W. Welche der ' +
      'beiden bindet, hängt davon ab, ob k_B T oder ℏ/τ_meas die kleinere Größe ' +
      'ist — die kleinere Vergleichsgröße ergibt die strengere Schranke. ' +
      bindingText() + '</p>',
    params: function () {
      return paramLine(function (d) {
        var r = ratios(d);
        return '<p>Am Arbeitspunkt: δE = <span class="num">' + fmtE(r.deltaE) +
          '</span> µeV, E_C = <span class="num">' +
          (state.floating ? fmtE(r.eC) : '—') + '</span> µeV.</p>' +
          '<p>Verhältnisse: K2a = <span class="num">' + fmtE(r.K2a) +
          '</span>, K2b = <span class="num">' + fmtE(r.K2b) +
          '</span>, K3 = <span class="num">' + fmtE(r.K3) +
          '</span>, K4 = <span class="num">' + fmtE(r.K4) + '</span>.</p>';
      });
    },
    numerics: checkHtml() +
      '<ul>' +
      '<li>Alle Werte werden direkt aus den Reglerwerten gerechnet, ohne Näherung. ' +
      'τ_meas wird für K2b in Nanosekunden umgerechnet, damit es zu ℏ in µeV·ns ' +
      'passt, und für K4 in Sekunden, damit es zu T_P passt.</li>' +
      '<li>Dieselbe Schwelle f wie in Panel A trennt „erfüllt" von ' +
      '„grenzwertig"; der Übergang zu „nicht erfüllt" liegt fest bei 1, denn ab ' +
      'dort ist die Ungleichung selbst verletzt, nicht nur ihr „≪".</li>' +
      '<li><b>T_P ist eine Eingabegröße.</b> Sie folgt nicht aus ' +
      MV.eq('energy.charging') + ' bis ' + MV.eq('energy.window') + ' und wird ' +
      'hier gesetzt, nicht berechnet. ' + MV.secOfWork('5.6') + ' hält ' +
      'ausdrücklich fest, dass eine große Ladungsenergie keine lange ' +
      'Paritätslebensdauer garantiert.</li>' +
      '<li>Im geerdeten Fall entfällt die K3-Zeile ganz, statt sie mit einem ' +
      'Ersatzwert zu füllen.</li>' +
      '</ul>',
    convention:
      '<p>' + UNITS_NOTE + '</p>' +
      '<p>' + VALUES_NOTE + '</p>' +
      '<p>' + PHI_NOTE + '</p>' +
      '<p>Die Kriteriennummern K2 bis K4 sind die aus ' + MV.tab('criteria') + '; ' +
      'die Aufspaltung von K2 in K2a und K2b ist eine Benennung dieser ' +
      'Darstellung, weil die Tabelle dort zwei Prüfgrößen in einer Zeile führt.</p>',
    reference:
      '<p>' + MV.secOfWork('5.1') + ' und ' + MV.sec('5.6') + ', ' +
      MV.tab('criteria') + '. Verwendete Gleichungen: ' +
      MV.eq('energy.charging') + ', ' + MV.eq('energy.tradeoff') + ', ' +
      MV.eq('energy.window') + '.</p>',
    reading:
      '<p><b>Ziehe L_W von 1 auf 6 µm</b> und beobachte, wie K2a und K2b von ' +
      '„nicht erfüllt" über „grenzwertig" nach „erfüllt" wandern, während K3 den ' +
      'umgekehrten Weg nimmt. Der Arbeitspunkt läuft dabei durch das Fenster in ' +
      'Panel A hindurch und wieder heraus.</p>' +
      '<p><b>Ziehe τ_meas von 1 auf 100 µs</b>: K2b verschlechtert sich um zwei ' +
      'Größenordnungen, K2a bleibt unverändert. Daran sieht man, dass die beiden ' +
      'Hälften von K2 verschiedene Dinge verlangen. <b>Die Anzeige „bindend unten" ' +
      'bleibt dabei auf K2b — und zwar in jeder Reglerstellung.</b> ' +
      bindingText() + ' Das ist der eigentliche Ertrag dieses Panels: die ' +
      'naheliegende Lesart von ' + MV.eq('energy.window') + ' als thermischer ' +
      'Bedingung greift zu kurz.</p>' +
      '<p><b>Zieh den logarithmischen T_P-Regler von 1 s nach links</b>, bis K4 ' +
      'von „erfüllt" über „grenzwertig" auf „nicht erfüllt" springt — das ' +
      'geschieht, sobald T_P unter τ_meas fällt. Am Fenster in Panel A ändert ' +
      'sich dabei nichts: K4 hängt nicht von der Drahtlänge ab und lässt sich ' +
      'durch Geometrie nicht heilen. Die Insel bleibt dabei schwebend, K3 also ' +
      'erfüllt — <b>ein wirksamer Ladungsschutz und eine kurze Lebensdauer ' +
      'schließen einander nicht aus.</b></p>' +
      '<p><b>Schalte danach auf „geerdet" und zieh T_P auf 10 s</b>: jetzt ist ' +
      'es umgekehrt, K3 entfällt und K4 ist erfüllt. Beide Wege stehen in ' +
      MV.sec('5.6') + ' nebeneinander. Welcher Mechanismus eine kurze ' +
      'Lebensdauer verursacht, sagt diese Darstellung nicht — T_P wird hier ' +
      'gesetzt, nicht gerechnet.</p>'
  };

  /* ====================================================================
     Aufbau
     ==================================================================== */

  function buildControls(mount) {
    var box = MV.ui.el('div', { class: 'controls' });

    function on(key) {
      return function (v) { state[key] = v; refresh(); };
    }

    sliders.deltaT = MV.ui.slider(box, {
      label: 'Δ_T — topologische Lücke (µeV)',
      min: 5, max: 400, step: 1, value: state.deltaT, digits: 0, onInput: on('deltaT')
    });
    sliders.xi = MV.ui.slider(box, {
      label: 'ξ — Kohärenzlänge (µm)',
      min: 0.02, max: 0.8, step: 0.01, value: state.xi, digits: 2, onInput: on('xi')
    });
    sliders.ecRef = MV.ui.slider(box, {
      label: 'E_C^ref — Ladungsenergie bei L_W = 1 µm (µeV)',
      min: 5, max: 400, step: 1, value: state.ecRef, digits: 0, onInput: on('ecRef')
    });
    sliders.T = MV.ui.slider(box, {
      label: 'T — Temperatur (mK)',
      min: T_MIN, max: T_MAX, step: 1, value: state.T, digits: 0, onInput: on('T')
    });
    sliders.tau = MV.ui.slider(box, {
      label: 'τ_meas — Messzeit (µs)',
      min: TAU_MIN, max: TAU_MAX, step: 0.01, value: state.tau, digits: 2, onInput: on('tau')
    });
    sliders.LW = MV.ui.slider(box, {
      label: 'L_W — Arbeitspunkt (µm)',
      min: 0.1, max: L_MAX, step: 0.01, value: state.LW, digits: 2, onInput: on('LW')
    });
    /* Zwei Regler, ein Wert: der obere fein über den gemessenen Bereich, der
       untere logarithmisch über alle Dekaden. Wer an einem zieht, bewegt den
       anderen mit. */
    sliders.TP = MV.ui.slider(box, {
      label: 'T_P — Paritätslebensdauer, fein',
      min: TP_FINE_MIN, max: TP_MAX, step: 0.001, value: state.TP,
      format: fmtTP,
      onInput: function (v) { setTP(v); }
    });
    sliders.TPlog = MV.ui.slider(box, {
      label: 'T_P — derselbe Wert, logarithmisch',
      hint: 'Eingabegröße, nicht aus ' + MV.eq('energy.charging') + ' bis ' +
            MV.eq('energy.window') + ' abgeleitet. Acht Dekaden von ' +
            fmtTP(TP_MIN) + ' bis ' + fmtTP(TP_MAX) + '; unter einer ' +
            'Millisekunde ist der Bereich illustrativ und liegt weit unter ' +
            'allem, was gemessen wurde.',
      min: TP_MIN, max: TP_MAX, step: 0.01, value: state.TP,
      scale: 'log', format: fmtTP,
      onInput: function (v) { setTP(v); }
    });
    sliders.factor = MV.ui.slider(box, {
      label: 'Schwelle für „≪" — Faktor',
      hint: 'um diesen Faktor müssen die Ungleichungen erfüllt sein',
      min: 1, max: 100, step: 1, value: state.factor, digits: 0, onInput: on('factor')
    });

    var row = MV.ui.el('div', { class: 'controls__wide' });
    var g = MV.ui.buttonRow(row, {
      label: 'Insel:',
      buttons: [{
        text: 'schwebend',
        pressed: true,
        onClick: function () { setFloating(!state.floating); }
      }]
    });
    groundButton = g.buttons[0];

    MV.ui.buttonRow(row, {
      label: 'Voreinstellungen:',
      buttons: [
        { text: 'Vorgabewerte', onClick: function () { preset(100, 0.1, 100, 20, 1, 3, 1, 10); } },
        { text: 'Schwelle 1 (nur „<")', onClick: function () { setFactor(1); } },
        { text: 'lange Kohärenzlänge (ξ = 0,3 µm)', onClick: function () { setXi(0.3); } }
      ]
    });

    row.appendChild(MV.ui.el('p', {
      class: 'note',
      html: 'Physikalische Einheiten statt dimensionsloser Größen. <b>Alle ' +
            'Vorgabewerte hier sind illustrativ</b>; einzige Zahl aus der Arbeit ' +
            'ist T = 20 mK. Gemessene Werte stehen im Reiter zu Kapitel 6.'
    }));
    box.appendChild(row);
    mount.appendChild(box);
  }

  /* Ein Schalter, zwei Aufrufer: der Knopf und die Vorführung in demo.js.
     Beide gehen denselben Weg, damit Beschriftung und aria-pressed nie von
     der Schalterstellung abweichen. */
  function setFloating(on) {
    state.floating = !!on;
    groundButton.setAttribute('aria-pressed', state.floating ? 'true' : 'false');
    groundButton.textContent = state.floating ? 'schwebend' : 'geerdet';
    refresh();
  }

  function preset(deltaT, xi, ecRef, T, tau, LW, TP, factor) {
    state.deltaT = deltaT; state.xi = xi; state.ecRef = ecRef; state.T = T;
    state.tau = tau; state.LW = LW; state.factor = factor;
    Object.keys(sliders).forEach(function (k) {
      if (sliders[k] && state[k] !== undefined) sliders[k].set(state[k], true);
    });
    setTP(TP, true);          /* zieht beide T_P-Regler nach */
    setFloating(true);
  }

  /**
   * T_P setzen — ein Wert, zwei Regler. Der feine kann nur seinen Ausschnitt
   * zeigen; liegt der Wert darunter, steht sein Griff am Anschlag, die Zahl
   * daneben bleibt aber die richtige, und der Hinweis darunter sagt es.
   */
  function setTP(v, silent) {
    /* Nicht auf ein Raster gerundet: jeder der beiden Regler rastet beim
       Ziehen schon auf sein eigenes, und beide zeigen den gespeicherten Wert
       statt ihrer Griffstellung. Würde hier zusätzlich gerundet, käme aus
       den 22 ms des feinen Reglers eine krumme Zahl. */
    var w = Math.min(TP_MAX, Math.max(TP_MIN, v));
    state.TP = w;
    if (sliders.TP) {
      sliders.TP.set(w, true);
      sliders.TP.setHint(sliders.TP.inRange(w)
        ? ''
        : 'Wert unterhalb dieses Reglers — der Griff steht am Anschlag. ' +
          'Einstellen lässt er sich mit dem logarithmischen Regler darunter.');
    }
    if (sliders.TPlog) sliders.TPlog.set(w, true);
    if (!silent) refresh();
  }

  function setFactor(f) { state.factor = f; sliders.factor.set(f, true); refresh(); }
  function setXi(x) { state.xi = x; sliders.xi.set(x, true); refresh(); }

  function buildReadout(mount) {
    readout = MV.ui.readout(mount, {
      items: [
        { key: 'kbt',      label: 'k_B T' },
        { key: 'hbartau',  label: 'ℏ / τ_meas' },
        { key: 'deltaE',   label: 'δE am Arbeitspunkt' },
        { key: 'window',   label: 'erlaubte Drahtlänge' },
        { key: 'windowXi', label: 'Fenster in L_W/ξ' },
        { key: 'binding',  label: 'bindend unten' }
      ],
      note:
        'Die Fenstergrenzen sind aus den eingegebenen Werten berechnet, nicht ' +
        'abgelesen. ' + MV.fig('energy.windowFig') + ' der Arbeit ist schematisch; hier liegen die ' +
        'Grenzen auf den Schnittpunkten mit den gepunkteten Schwellenlinien. Bei ' +
        'Schwelle 1 fallen diese mit den Skalen selbst zusammen.'
    });
  }

  /* Die Hervorhebungen setzen und wieder wegnehmen. Sie sind reine Anzeige;
     `refresh` zeichnet Panel B ohnehin neu. */
  function markCriteria(list) {
    mark.rows = {};
    (list || []).forEach(function (k) { mark.rows[k] = true; });
    refresh();
  }

  function clearMarks() {
    if (!mark.holdK3 && !Object.keys(mark.rows).length) return;
    mark.rows = {};
    mark.holdK3 = false;
    refresh();
  }

  /* Schnittstelle für die Vorführungen (demo.js) — nur vorhandene Regler. */
  function buildDemoHandle() {
    MV.demo.provide('energieskalen', MV.demo.handle({
      state: state,
      sliders: sliders,
      refresh: refresh,
      stop: clearMarks,
      setters: { floating: setFloating, TP: function (v) { setTP(v); } },
      extra: {
        read: function () { return derive(state); },
        ratiosNow: function () { return ratios(derive(state)); },
        /* Die Kriterien, wie sie in Panel B wirklich stehen: aus der
           gezeichneten Tabelle gelesen, nicht noch einmal gerechnet. */
        criteriaShown: function () {
          if (!panelB || !panelB.body) return [];
          return Array.prototype.map.call(
            panelB.body.querySelectorAll('tbody tr td b'),
            function (e) { return e.textContent; });
        },
        /* Der Status einer Kriterienzeile, so wie er dort steht. */
        criteriaVerdict: function (key) {
          if (!panelB || !panelB.body) return null;
          var out = null;
          Array.prototype.forEach.call(
            panelB.body.querySelectorAll('tbody tr'), function (tr) {
              var b = tr.querySelector('td b'), v = tr.querySelector('.verdict');
              if (b && v && b.textContent === key) out = v.textContent;
            });
          return out;
        },
        markedCriteria: function () {
          if (!panelB || !panelB.body) return [];
          return Array.prototype.map.call(
            panelB.body.querySelectorAll('tr.numtable__marked td b'),
            function (e) { return e.textContent; });
        },
        fmtTP: fmtTP,
        tpRange: function () {
          return { min: TP_MIN, max: TP_MAX, fineMin: TP_FINE_MIN };
        },
        markCriteria: markCriteria,
        holdK3: function (on) { mark.holdK3 = !!on; refresh(); },
        clearMarks: clearMarks,
        showPanel: function (which) {
          var p = which === 'A' ? panelA : panelB;
          if (p && p.root.scrollIntoView) p.root.scrollIntoView({ block: 'center' });
        },
        defaults: function () { preset(100, 0.1, 100, 20, 1, 3, 1, 10); }
      }
    }));
  }

  MV.energyscales = {
    check: check,

    init: function (mount) {
      if (!check.passed) {
        var bad = check.parts.filter(function (p) { return p.ok === false; });
        var warn = MV.ui.el('div', {
          class: 'readout',
          html: '<b>Warnung:</b> ' + bad.length + ' Prüfung(en) dieses Reiters sind ' +
                'fehlgeschlagen: ' + bad.map(function (p) { return p.name; }).join('; ') +
                '. Die Aussagen über die bindende Schranke sind nicht belastbar.'
        });
        warn.style.borderColor = MV.palette.critical;
        warn.style.color = MV.palette.critical;
        mount.appendChild(warn);
      }

      buildControls(mount);
      buildReadout(mount);

      var panels = MV.ui.el('div', { class: 'panels' });
      mount.appendChild(panels);

      panelA = MV.ui.panel(panels, {
        title: 'Panel A — Zielkonflikt und Bedingungsfenster',
        quote: 'energieskalen.A',
        lead: 'Die vier Energieskalen über der Drahtlänge, logarithmisch. ' +
              'Interaktive Fassung von ' + MV.fig('energy.windowFig') + '.',
        wide: true,
        aspect: 0.46,
        margin: { l: 66, r: 18, t: 14, b: 62 },
        info: INFO_A
      });

      buildDemoHandle();
      MV.demo.onHandInput(clearMarks);

      panelB = MV.ui.panel(panels, {
        title: 'Panel B — Prüfung am Arbeitspunkt',
        quote: 'energieskalen.B',
        lead: 'Die vier Skalen und die vier Verhältnisse nach ' +
              MV.tab('criteria') + ', mit der bindenden Schranke.',
        wide: true,
        canvas: false,
        info: INFO_B
      });
    },

    draw: function () { refresh(); },

    /* Es läuft nichts im Hintergrund; aufzuräumen ist nur, was eine
       Vorführung hervorgehoben hat. */
    stop: clearMarks
  };

}(MV));
