/* ------------------------------------------------------------------------
   nanowire.js — tab 2: Rashba nanowire, sections 4.1 – 4.3 of the thesis.

   Panel A: bands in three stages                       (figure 4.2 a/b/c)
   Panel B: effective p-wave pairing Δ_eff(k)           (equation 4.3)
   Panel C: phase diagram                               (figure 4.3)

   Model, equation (4.1):
     H = ∫dx ψ† ( −∂ₓ²/2m − μ − i α σ_y ∂ₓ + E_Z σ_z ) ψ + ∫dx ( Δ ψ↑ψ↓ + h.c. )
   Spin-orbit axis ∥ ŷ, Zeeman field ⊥ to it. ℏ = m = 1, everything
   dimensionless.

   Axis convention, following figure 4.2 of the thesis: stages 1 and 2 plot the
   bare bands k²/2 ± sqrt(α²k² + E_Z²) and draw μ as a horizontal line, stage 3
   switches to ±E(k) where μ sits inside ξ_k and the zero line is the Fermi
   level. The switch is stated in the info panel of panel A.
   ------------------------------------------------------------------------ */

var MV = window.MV || (window.MV = {});

(function (MV) {
  'use strict';

  /* Tolerance for "gap closed" / "on the phase boundary", same value as in
     tab 1 so both tabs report a transition at the same distance from it. */
  var TOL = 0.02;

  var N_CURVE = 600;         /* samples per plotted curve (601 points)      */
  var K_PANEL_B = 3;         /* fixed momentum window of panel B            */
  var EZ_MIN_APPROX = 0.05;  /* below this the linear approximation has no
                                meaning: it diverges as Δα/E_Z             */

  /* Starting point chosen to reproduce figure 4.2b: α² > E_Z, so the lower
     band has the double well of the figure, μ lies inside the Zeeman window,
     and the criterion (4.4) is met. */
  var state = { mu: -0.3, ez: 1.0, alpha: 1.4, delta: 0.6, k: 0.8, stage: 1 };

  var sliders = {}, panelA = null, panelB = null, panelC = null;
  var readout = null, segStage = null;

  /* ====================================================================
     Physics
     ==================================================================== */

  /**
   * Bare bands of (4.2) without the −μ term, as drawn in figure 4.2 a/b:
   *   k²/2 ± sqrt( α²k² + E_Z² )
   * `ez` is passed explicitly because stage 1 shows the same expression at
   * E_Z = 0 regardless of the slider.
   */
  function bandBare(k, s, sign, ez) {
    if (ez === undefined) ez = s.ez;
    return 0.5 * k * k + sign * Math.sqrt(s.alpha * s.alpha * k * k + ez * ez);
  }

  /** Equation (4.2) itself: ε±(k) = k²/2 − μ ± sqrt(α²k² + E_Z²). */
  function bandEps(k, s, sign) {
    return bandBare(k, s, sign) - s.mu;
  }

  /**
   * Quasiparticle spectrum of (4.1) in closed form — Alicea (2012), Gl. (60).
   * Section 4.3 of the thesis quotes only its value at k = 0, Gl. (4.4);
   * that value is reproduced exactly by this expression (see selfCheck).
   *
   *   ξ_k = k²/2 − μ
   *   E±(k)² = ξ_k² + α²k² + E_Z² + Δ² ± 2 sqrt( ξ_k²(α²k² + E_Z²) + E_Z²Δ² )
   */
  function bdg(k, s) {
    var xi = 0.5 * k * k - s.mu;
    var R = s.alpha * s.alpha * k * k + s.ez * s.ez;
    var base = xi * xi + R + s.delta * s.delta;
    var root = 2 * Math.sqrt(Math.max(0, xi * xi * R + s.ez * s.ez * s.delta * s.delta));
    return {
      plus: Math.sqrt(Math.max(0, base + root)),
      minus: Math.sqrt(Math.max(0, base - root))
    };
  }

  /** Excitation energy at k = 0, Gl. (4.4). */
  function energyAtZero(s) {
    return Math.abs(s.ez - Math.sqrt(s.delta * s.delta + s.mu * s.mu));
  }

  /**
   * Effective p-wave pairing, Gl. (4.3):
   *   Δ_eff(k) = Δ · αk / sqrt( E_Z² + α²k² ),   approximation  (Δα/E_Z) · k
   * Odd in k and zero at k = 0. At E_Z = 0 the expression degenerates to
   * Δ·sgn(k), which is discontinuous at k = 0; Δ_eff(0) = 0 is kept because
   * the function has to be odd.
   */
  function deltaEff(k, s) {
    if (k === 0) return 0;
    return s.delta * s.alpha * k / Math.sqrt(s.ez * s.ez + s.alpha * s.alpha * k * k);
  }

  function deltaEffApprox(k, s) {
    return s.delta * s.alpha * k / s.ez;
  }

  /**
   * Relative deviation of the approximation from the exact expression, in
   * closed form. The quotient is
   *   Δ_approx / Δ_exact = sqrt( E_Z² + α²k² ) / E_Z,
   * so the deviation is sqrt(1 + α²k²/E_Z²) − 1, which is finite at k = 0
   * and needs no 0/0 handling.
   */
  function deltaEffDeviation(k, s) {
    if (s.ez < EZ_MIN_APPROX) return null;
    var r = s.alpha * k / s.ez;
    return Math.sqrt(1 + r * r) - 1;
  }

  /**
   * Fermi points: solutions of  k²/2 ± sqrt(α²k² + E_Z²) = μ, solved exactly.
   *
   * With x = k²:  x² − 4(μ + α²)x + 4(μ² − E_Z²) = 0,
   *   x = 2(μ + α²) ± 2 sqrt( α⁴ + 2μα² + E_Z² ).
   * Squaring can introduce spurious roots, so every candidate is put back
   * into the original equation and only kept if it satisfies it.
   *
   * Returns a sorted list of { k, sign }.
   */
  function fermiPoints(s) {
    var a2 = s.alpha * s.alpha;
    var disc = a2 * a2 + 2 * s.mu * a2 + s.ez * s.ez;
    var out = [];
    if (!(disc >= 0)) return out;

    var sq = Math.sqrt(disc);
    var tol = 1e-7 * (1 + Math.abs(s.mu) + s.ez);

    [2 * (s.mu + a2) + 2 * sq, 2 * (s.mu + a2) - 2 * sq].forEach(function (x) {
      if (x < -tol) return;
      x = Math.max(0, x);
      var k = Math.sqrt(x);
      [1, -1].forEach(function (sign) {
        var residual = 0.5 * x + sign * Math.sqrt(a2 * x + s.ez * s.ez) - s.mu;
        if (Math.abs(residual) > tol) return;      /* spurious root of squaring */
        if (k < 1e-9) out.push({ k: 0, sign: sign });
        else { out.push({ k: -k, sign: sign }); out.push({ k: k, sign: sign }); }
      });
    });

    out.sort(function (a, b) { return a.k - b.k; });
    return out.filter(function (p, i) {
      return i === 0 || Math.abs(p.k - out[i - 1].k) > 1e-7;
    });
  }

  /**
   * Was an den Fermi-Punkten des ungepaarten Systems gilt. Physikalisch zählt
   * Δ_eff genau dort — der k-Regler von Panel B misst dagegen nur, wie gut die
   * lineare Näherung an einer beliebigen Stelle ist.
   *
   * Es wird nur ausgewertet, nicht gedeutet: dass die Lücke gleich Δ_eff(k_F)
   * wäre, steht weder in der Arbeit noch in CLAUDE.md und wird hier auch nicht
   * behauptet.
   */
  function fermiInfo(s) {
    var pts = fermiPoints(s).filter(function (p) { return p.k > 0; });
    var noApprox = s.ez < EZ_MIN_APPROX;
    return {
      all: fermiPoints(s),
      positive: pts,
      count: fermiPoints(s).length,
      /* Vier Fermi-Punkte: μ liegt oberhalb des Zeeman-Fensters, beide Bänder
         sind besetzt — dann ist das System nicht effektiv spinlos. */
      spinless: fermiPoints(s).length === 2,
      rows: pts.map(function (p) {
        return {
          k: p.k,
          exact: deltaEff(p.k, s),
          approx: noApprox ? null : deltaEffApprox(p.k, s),
          dev: deltaEffDeviation(p.k, s),
          ratio: s.ez > 0 ? s.alpha * p.k / s.ez : null,
          inWindow: Math.abs(p.k) <= K_PANEL_B
        };
      })
    };
  }

  /** Everything the three panels and the read-out need. */
  function derive(s) {
    var boundary = Math.sqrt(s.delta * s.delta + s.mu * s.mu);
    var e0 = energyAtZero(s);
    var critical = e0 < TOL;
    return {
      s: s,
      e0: e0,
      boundary: boundary,
      critical: critical,
      topological: !critical && s.ez > boundary,
      fermi: fermiPoints(s)
    };
  }

  /** Der Zustand des Reiters mit einzelnen Werten überschrieben. */
  function withState(over) {
    var s = { mu: state.mu, ez: state.ez, alpha: state.alpha,
              delta: state.delta, k: state.k, stage: state.stage }, key;
    if (over) {
      for (key in over) {
        if (Object.prototype.hasOwnProperty.call(over, key)) s[key] = over[key];
      }
    }
    return s;
  }

  function phaseKind(d) {
    return d.critical ? 'critical' : (d.topological ? 'topological' : 'trivial');
  }

  function phaseName(d) {
    if (d.critical) return 'Übergang, Lücke bei k = 0 geschlossen';
    return d.topological ? 'topologisch (ν = −1)' : 'trivial (ν = +1)';
  }

  function phaseColor(d) {
    return MV.palette[phaseKind(d)];
  }

  /* ====================================================================
     Automatic consistency checks, run once when the tab is built.

       1. at k = 0 the closed form must give |E_Z − sqrt(Δ² + μ²)|, i.e. (4.4)
       2. for Δ → 0 the pair {E₊, E₋} must reduce to {|ε₊|, |ε₋|} from (4.2)
       3. Δ_eff must be odd in k and vanish at k = 0
       4. with μ inside the Zeeman window there must be exactly two Fermi
          points, above the window exactly four
     ==================================================================== */

  function selfCheck() {
    var samples = [], worst0 = 0, worstD = 0, worstOdd = 0, i, j;
    var countFail = [];
    var DTINY = 1e-6;    /* stand-in for Δ → 0 */
    var mus = [-2.3, -0.7, 0, 0.9, 2.5];
    var ezs = [0, 0.4, 1.1, 2.2, 3];
    var als = [0.15, 0.8, 1.7];
    var dls = [0.2, 1.0, 1.9];

    mus.forEach(function (mu) {
      ezs.forEach(function (ez) {
        als.forEach(function (al) {
          dls.forEach(function (dl) {
            samples.push({ mu: mu, ez: ez, alpha: al, delta: dl });
          });
        });
      });
    });

    for (i = 0; i < samples.length; i++) {
      var s = samples[i];

      /* 1 — value at k = 0 against Gl. (4.4) */
      var d0 = Math.abs(bdg(0, s).minus - energyAtZero(s));
      if (d0 > worst0) worst0 = d0;

      /* 2 — limit Δ → 0 against (4.2), tested at a small but finite Δ.
         The deviation vanishes linearly in Δ, so the bound is stated
         relative to DTINY rather than as an absolute number. */
      var sd = { mu: s.mu, ez: s.ez, alpha: s.alpha, delta: DTINY };
      for (j = -30; j <= 30; j++) {
        var k = j / 10;
        var e = bdg(k, sd);
        var want = [Math.abs(bandEps(k, sd, 1)), Math.abs(bandEps(k, sd, -1))].sort(function (a, b) { return a - b; });
        var got = [e.minus, e.plus].sort(function (a, b) { return a - b; });
        var dd = Math.max(Math.abs(got[0] - want[0]), Math.abs(got[1] - want[1]));
        if (dd > worstD) worstD = dd;
      }

      /* 3 — Δ_eff odd in k, zero at k = 0 */
      if (Math.abs(deltaEff(0, s)) > worstOdd) worstOdd = Math.abs(deltaEff(0, s));
      for (j = 1; j <= 30; j++) {
        var kk = j / 10;
        var odd = Math.abs(deltaEff(-kk, s) + deltaEff(kk, s));
        if (odd > worstOdd) worstOdd = odd;
      }

      /* 4 — number of Fermi points against the Zeeman window. Values right on
         the window edge are skipped: k = 0 is then a Fermi point itself. */
      if (s.ez > 0.05) {
        var nF = fermiPoints(s).length;
        if (Math.abs(s.mu) < s.ez - 0.05 && nF !== 2) {
          countFail.push('μ=' + s.mu + ' E_Z=' + s.ez + ' α=' + s.alpha + ': ' + nF + ' statt 2');
        }
        if (s.mu > s.ez + 0.05 && nF !== 4) {
          countFail.push('μ=' + s.mu + ' E_Z=' + s.ez + ' α=' + s.alpha + ': ' + nF + ' statt 4');
        }
      }
    }

    /* 5 — die Abnahmewerte aus Auftrag 17, Teil 2, an den Vorgabewerten.
       Sie stehen im Info-Panel als Zahlen; hier werden sie nachgerechnet. */
    var ref = { mu: -0.3, ez: 1.0, alpha: 1.4, delta: 0.6 };
    var refFail = [];
    var rfi = fermiInfo(ref);
    function expect(name, got, soll, tol) {
      if (!(Math.abs(got - soll) <= tol)) {
        refFail.push(name + ' = ' + MV.fmt(got, 4) + ' statt ' + MV.fmt(soll, 4));
      }
    }
    if (rfi.count !== 2) refFail.push('Fermi-Punkte: ' + rfi.count + ' statt 2');
    else {
      var r0 = rfi.rows[0];
      expect('k_F', r0.k, 2.674, 5e-3);
      expect('Δ_eff(k_F)', r0.exact, 0.580, 5e-3);
      expect('Näherung', r0.approx, 2.246, 5e-3);
      expect('Abweichung in %', 100 * r0.dev, 287, 1);
    }

    return {
      cases: samples.length,
      referenceFailures: refFail,
      reference: rfi.count === 2 ? rfi.rows[0] : null,
      maxDeviationAtZero: worst0,
      maxDeviationDeltaToZero: worstD,
      maxDeviationOddness: worstOdd,
      fermiCountFailures: countFail,
      deltaUsed: DTINY,
      /* (4.4) must be reproduced to machine precision; the Δ → 0 limit is
         approached linearly, so a few times Δ is the right bound; oddness is
         exact up to rounding */
      passed: worst0 < 1e-9 && worstD < 10 * DTINY &&
              worstOdd < 1e-15 && countFail.length === 0 && refFail.length === 0
    };
  }

  MV.nanowireSelfCheck = selfCheck;

  /* ====================================================================
     Panel A — bands, three stages   (figure 4.2 a/b/c)
     ==================================================================== */

  /** Momentum window: wide enough to always contain the Fermi points. */
  function windowK(s, fermi) {
    var maxKf = 0;
    fermi.forEach(function (p) { maxKf = Math.max(maxKf, Math.abs(p.k)); });
    return MV.ceilTo(Math.max(1.0, 2.2 * s.alpha, 1.25 * maxKf), 0.5);
  }

  /**
   * Parameters as panel A draws them. Stage 1 shows (4.2) at E_Z = 0
   * regardless of the slider, so the Fermi points marked there have to be
   * those of the field-free system — otherwise markers and curves would
   * disagree.
   */
  function stageState(s) {
    return {
      mu: s.mu, alpha: s.alpha, delta: s.delta,
      ez: s.stage === 0 ? 0 : s.ez
    };
  }

  /** One-line reading aid under the bands; derived, never hard coded. */
  function bannerText(stage, nFermi, s) {
    if (stage === 0) {
      return 'Stufe 1: E_Z = 0 gesetzt — ' + nFermi + ' Fermi-Punkte, spinbehaftet';
    }
    if (Math.abs(s.mu) < s.ez) {
      return 'μ im Zeeman-Fenster — ' + nFermi + ' Fermi-Punkte, effektiv spinlos';
    }
    return 'μ außerhalb des Zeeman-Fensters — ' + nFermi + ' Fermi-Punkte';
  }

  function drawPanelA(d) {
    var p = panelA.plot, s = d.s, i, k;
    var stage = s.stage;
    var sEff = stageState(s);
    var ezUsed = sEff.ez;
    /* Fermi points of the configuration actually drawn (see stageState). */
    var fermi = stage === 0 ? fermiPoints(sEff) : d.fermi;
    var K = windowK(s, fermi);
    var n = N_CURVE;

    if (stage < 2) {
      /* ---------------- stages 1 and 2: bare bands, μ as a line ---------- */
      var up = [], lo = [], yMin = Infinity;
      for (i = 0; i <= n; i++) {
        k = -K + 2 * K * i / n;
        var a = 0.5 * k * k + Math.sqrt(s.alpha * s.alpha * k * k + ezUsed * ezUsed);
        var b = 0.5 * k * k - Math.sqrt(s.alpha * s.alpha * k * k + ezUsed * ezUsed);
        up.push([k, a]); lo.push([k, b]);
        yMin = Math.min(yMin, b);
      }
      /* Scale to the structure that matters — band bottom and Zeeman window —
         and let the rising parabolas leave the frame, as in figure 4.2.
         The energy axis is deliberately independent of μ: the Zeeman window is
         2·E_Z high whatever μ does, and an axis that followed μ would change
         its apparent height while μ is dragged through it. The frame only
         widens when μ itself would otherwise leave it. */
      var Y0 = Math.min(
        -MV.ceilTo(Math.max(-yMin, ezUsed) * 1.25 + 0.2, 0.5),
        -MV.ceilTo(-s.mu + 0.4, 0.5)
      );
      var Y1 = Math.max(
        MV.ceilTo(Math.max(ezUsed, -yMin) * 2.0 + 0.3, 0.5),
        MV.ceilTo(s.mu + 0.4, 0.5)
      );

      p.begin({
        x: [-K, K], y: [Y0, Y1],
        xLabel: 'k', yLabel: 'ε(k)',
        xTickCount: 5, yTickCount: 5
      });

      /* Zeeman window: the gap of size 2·E_Z that the field opens at k = 0 */
      if (stage === 1 && s.ez > 1e-9) {
        p.hspan(-s.ez, s.ez, {
          fill: MV.palette.neutralSoft, edge: MV.palette.muted, dash: [3, 3],
          legend: 'Zeeman-Fenster'
        });
        /* bands without the field, for comparison (figure 4.2b) */
        p.fn(function (kk) { return bandBare(kk, s, 1, 0); }, -K, K, n,
          { color: MV.palette.muted, width: 1, dash: [4, 3], alpha: 0.75,
            legend: 'ohne Feld (E_Z = 0)' });
        p.fn(function (kk) { return bandBare(kk, s, -1, 0); }, -K, K, n,
          { color: MV.palette.muted, width: 1, dash: [4, 3], alpha: 0.75 });
      }

      p.curve(up, { color: MV.palette.curve1, width: 2, legend: 'ε₊ (oberes Band)' });
      p.curve(lo, { color: MV.palette.curve2, width: 2, legend: 'ε₋ (unteres Band)' });

      p.hline(s.mu, {
        color: MV.palette.ink, dash: [6, 4], width: 1.4,
        label: 'μ', labelAt: 'left'
      });

      fermi.forEach(function (pt) {
        p.marker(pt.k, s.mu, { color: MV.palette.topological, r: 4.5, shape: 'ring', width: 2 });
      });
      p.legendItems.push({
        label: 'Fermi-Punkte: ' + fermi.length,
        color: MV.palette.topological, swatch: 'dot'
      });

      if (stage === 1 && s.ez > 1e-9) {
        p.dimension(0, -s.ez, s.ez, '2·E_Z = ' + MV.fmt(2 * s.ez, 2),
          { color: MV.palette.ink, side: 'right', offset: -14 });
      }

      panelA.setStatus(bannerText(stage, fermi.length, s), MV.palette.muted);

      p.finish();
      return;
    }

    /* ------------------- stage 3: quasiparticle spectrum ---------------- */
    var eMinus = [], ePlus = [], mMax = 0;
    for (i = 0; i <= n; i++) {
      k = -K + 2 * K * i / n;
      var e = bdg(k, s);
      eMinus.push([k, e.minus]);
      ePlus.push([k, e.plus]);
      if (e.minus > mMax) mMax = e.minus;
    }
    var Y = MV.ceilTo(Math.max(mMax * 1.1, 0.5), 0.5);

    p.begin({
      x: [-K, K], y: [-Y, Y],
      xLabel: 'k', yLabel: '±E(k)',
      xTickCount: 5, yTickCount: 5
    });

    /* Fermi momenta of the unpaired system, for the link to stage 2 */
    d.fermi.forEach(function (pt) {
      p.vline(pt.k, { color: MV.palette.line, dash: [3, 3] });
    });
    if (d.fermi.length) {
      p.legendItems.push({
        label: 'Fermi-Punkte ohne Δ (' + d.fermi.length + ')',
        color: MV.palette.line, dash: [3, 3]
      });
    }

    p.curve(ePlus, { color: MV.palette.curve1, width: 1.8, legend: '±E₊(k)' });
    p.curve(ePlus.map(function (q) { return [q[0], -q[1]]; }), { color: MV.palette.curve1, width: 1.8 });
    p.curve(eMinus, { color: MV.palette.curve2, width: 2, legend: '±E₋(k)' });
    p.curve(eMinus.map(function (q) { return [q[0], -q[1]]; }), { color: MV.palette.curve2, width: 2 });

    p.hline(0, { color: MV.palette.axis, dash: [4, 4], label: 'E = 0 (Fermi-Niveau)', labelAt: 'left' });

    if (d.e0 > TOL) {
      p.dimension(0, 0, d.e0, 'E(0) = ' + MV.fmt(d.e0, 2),
        { color: MV.palette.ink, side: 'right' });
      panelA.setStatus(
        'Lücke bei k = 0: E(0) = ' + MV.fmt(d.e0, 2) + ' — ' +
        (d.topological ? 'topologische Phase' : 'triviale Phase'),
        phaseColor(d));
    } else {
      p.marker(0, 0, { color: MV.palette.critical, r: 5, shape: 'ring', width: 2 });
      panelA.setStatus('Lücke bei k = 0 geschlossen — Phasengrenze', MV.palette.critical);
    }

    p.finish();
  }

  /**
   * Der Satz für Punkt 8: wie weit die Vorgabewerte am Fermi-Punkt von der
   * Bedingung αk ≪ E_Z entfernt sind. Alle Zahlen live aus dem Zustand.
   */
  function fermiReadingText() {
    var fi = fermiInfo(state);
    if (!fi.rows.length) {
      return 'Bei der eingestellten Kombination gibt es keine Fermi-Punkte — ' +
             'μ liegt unterhalb beider Bänder.';
    }
    var r = fi.rows[fi.rows.length - 1];
    if (r.approx === null || r.ratio === null) {
      return 'Bei E_Z → 0 ist die Näherung nicht definiert; der Vergleich ' +
             'am Fermi-Punkt entfällt.';
    }
    return 'Der Regler springt auf k_F = ' + MV.fmt(r.k, 3) + ', und die ' +
      'Statuszeile zeigt dort Δ_eff exakt ' + MV.fmt(r.exact, 3) + ' gegen ' +
      MV.fmt(r.approx, 3) + ' genähert — ' + MV.fmt(100 * r.dev, 0) + ' % daneben. ' +
      'Die Näherung setzt αk ≪ E_Z voraus; am Fermi-Punkt ist αk_F/E_Z = ' +
      MV.fmt(r.ratio, 2) + ', also <b>größer</b> als eins. Genau dort, wo Δ_eff ' +
      'physikalisch zählt, taugt die lineare Näherung nicht.';
  }

  /**
   * Was an den Fermi-Punkten steht — dieselbe Zeichenkette für Statuszeile und
   * Info-Panel, damit beide nicht auseinanderlaufen können.
   */
  function fermiText(fi) {
    if (!fi.count) return 'Keine Fermi-Punkte: μ liegt unterhalb beider Bänder.';
    var parts = fi.rows.map(function (r) {
      return 'bei k_F = ±' + MV.fmt(r.k, 3) + ': Δ_eff exakt ' + MV.fmt(r.exact, 3) +
        (r.approx === null ? '' :
          ', genähert ' + MV.fmt(r.approx, 3) +
          ', Abweichung ' + MV.fmt(100 * r.dev, 1) + ' %');
    });
    return parts.join('; ') +
      (fi.spinless ? '.'
        : '. Vier Fermi-Punkte — μ liegt oberhalb des Zeeman-Fensters, beide ' +
          'Bänder sind besetzt, und das System ist dann nicht effektiv spinlos.');
  }

  /* ====================================================================
     Panel B — effective p-wave pairing   (equation 4.3)
     ==================================================================== */

  function drawPanelB(d) {
    var p = panelB.plot, s = d.s, i, k;
    var K = K_PANEL_B, n = N_CURVE;
    var noApprox = s.ez < EZ_MIN_APPROX;

    var Y = MV.ceilTo(s.delta * 1.25, 0.5);

    p.begin({
      x: [-K, K], y: [-Y, Y],
      xLabel: 'k', yLabel: 'Δ_eff(k)',
      xTickCount: 5, yTickCount: 5
    });

    p.hline(0, { color: MV.palette.line });
    p.vline(0, { color: MV.palette.line });

    /* approximation first, so the exact curve stays on top */
    if (!noApprox) {
      p.fn(function (kk) { return deltaEffApprox(kk, s); }, -K, K, 2,
        { color: MV.palette.curve3, width: 1.6, dash: [6, 4],
          legend: 'Näherung (Δα/E_Z)·k' });
    }

    /* At E_Z = 0 the exact expression jumps at k = 0; the two branches are
       drawn separately rather than joined by a vertical segment. */
    var exact = [];
    for (i = 0; i <= n; i++) {
      k = -K + 2 * K * i / n;
      if (noApprox && Math.abs(k) < 1e-9) { exact.push(null); continue; }
      exact.push([k, deltaEff(k, s)]);
    }
    if (noApprox) {
      exact = [];
      for (i = 0; i <= n; i++) {
        k = -K + 2 * K * i / n;
        if (Math.abs(k) < K / n / 2) { exact.push(null); continue; }
        exact.push([k, deltaEff(k, s)]);
      }
    }
    p.curve(exact, { color: MV.palette.curve2, width: 2, legend: 'exakt, Gl. ' + MV.eq('wire.deltaEff') + '' });

    p.marker(0, 0, { color: MV.palette.ink, r: 3.5 });

    var kk = s.k;
    var ex = deltaEff(kk, s);
    p.vline(kk, { color: MV.palette.critical, dash: [3, 3] });
    p.marker(kk, ex, { color: MV.palette.critical, r: 4.5, shape: 'ring', width: 2 });
    if (!noApprox) {
      p.marker(kk, deltaEffApprox(kk, s), { color: MV.palette.curve3, r: 3.5 });
    }

    /* Die Fermi-Punkte des ungepaarten Systems. Sie kommen aus Panel A und
       sind die Stellen, an denen Δ_eff physikalisch zählt. */
    var fi = fermiInfo(s);
    var outside = [];
    var legendDone = false;
    fi.rows.forEach(function (r) {
      if (!r.inWindow) { outside.push(r.k); return; }
      [r.k, -r.k].forEach(function (kf) {
        p.vline(kf, {
          color: MV.palette.curve1, dash: [6, 3], width: 1.4,
          legend: legendDone ? null : 'Fermi-Punkte ±k_F'
        });
        legendDone = true;
      });
    });

    var dev = deltaEffDeviation(kk, s);
    var msg = noApprox
      ? 'Für E_Z → 0 ist die Näherung nicht definiert; exakt wird Δ_eff zu Δ·sgn(k)'
      : 'bei k = ' + MV.fmt(kk, 2) + ': exakt ' + MV.fmt(ex, 3) +
        ', genähert ' + MV.fmt(deltaEffApprox(kk, s), 3) +
        ', relative Abweichung ' + MV.fmt(100 * dev, 1) + ' %';

    msg += '  ·  ' + fermiText(fi);
    if (outside.length) {
      msg += ' Außerhalb des gezeigten Bereichs |k| ≤ ' + MV.fmt(K_PANEL_B, 0) +
             ': ' + outside.map(function (k) { return '±' + MV.fmt(k, 3); }).join(', ') + '.';
    }

    panelB.setStatus(msg, noApprox ? MV.palette.critical : MV.palette.ink);

    p.finish();
  }

  /* ====================================================================
     Panel C — phase diagram   (figure 4.3)
     ==================================================================== */

  function drawPanelC(d) {
    var p = panelC.plot, s = d.s, i;
    var X = 3, Y = 3, n = 400;

    p.begin({
      x: [-X, X], y: [0, Y],
      xLabel: 'μ/Δ', yLabel: 'E_Z/Δ',
      xTickCount: 6, yTickCount: 4
    });

    /* boundary E_Z = sqrt(Δ² + μ²), in units of Δ: E_Z/Δ = sqrt(1 + (μ/Δ)²) */
    var pts = [];
    for (i = 0; i <= n; i++) {
      var x = -X + 2 * X * i / n;
      pts.push([x, Math.sqrt(1 + x * x)]);
    }
    p.areaAbove(pts, { fill: MV.palette.topologicalSoft });
    p.areaBelow(pts, { fill: MV.palette.trivialSoft });
    p.curve(pts, { color: MV.palette.topological, width: 2 });

    p.text('topologisch', 0, 2.55, { color: MV.palette.topological, align: 'center', box: true });
    p.text('ν = −1', 0, 2.25, { color: MV.palette.topological, align: 'center', box: true });
    p.text('trivial', 0, 0.62, { color: MV.palette.trivial, align: 'center', box: true });
    p.text('ν = +1', 0, 0.32, { color: MV.palette.trivial, align: 'center', box: true });

    p.marker(0, 1, { color: MV.palette.topological, r: 3.5 });
    p.label('E_Z = Δ', p.X(0) + 8, p.Y(1) + 2,
      { color: MV.palette.muted, baseline: 'top', box: true });

    /* working point */
    var mx = s.mu / s.delta, my = s.ez / s.delta;
    var outside = mx < -X || mx > X || my > Y;
    var cx = Math.min(X, Math.max(-X, mx));
    var cy = Math.min(Y, Math.max(0, my));

    p.marker(cx, cy, {
      color: phaseColor(d),
      r: outside ? 5 : 6,
      shape: outside ? 'ring' : 'dot',
      width: 2,
      label: 'Arbeitspunkt (' + MV.fmt(mx, 2) + ', ' + MV.fmt(my, 2) + ')' +
             (outside ? ' — außerhalb' : ''),
      anchor: cy > Y * 0.75 ? 'below' : 'above'
    });

    p.finish();
  }

  /* ====================================================================
     Read-out and update cycle
     ==================================================================== */

  function refresh() {
    var d = derive(state);

    readout.set('e0', MV.fmt(d.e0, 2));
    readout.set('boundary', MV.fmt(d.boundary, 2));
    readout.set('criterion', d.critical ? 'Grenzfall' : (d.topological ? 'erfüllt' : 'nicht erfüllt'));
    readout.set('fermi', String(d.fermi.length));
    readout.status(phaseKind(d), phaseName(d));

    drawPanelA(d);
    drawPanelB(d);
    drawPanelC(d);

    panelA.info.update();
    panelB.info.update();
    panelC.info.update();
    return d;
  }

  /* ====================================================================
     Info panels
     ==================================================================== */

  var ALICEA_60 =
    'Geschlossene Form des Quasiteilchenspektrums nach Alicea (2012), Gl. (60). ' +
    'Die Arbeit gibt in Abschn. 4.3 nur den Wert bei k = 0 an, Gl. ' + MV.eq('wire.criterion') + '; dieser ' +
    'wird von der hier verwendeten Form exakt reproduziert.';

  var AXIS_SWITCH =
    'Die Achsenbeschriftung wechselt zwischen den Stufen: Stufe 1 und 2 zeigen die ' +
    'nackten Bänder ε(k) nach ' + MV.eq('wire.bands') + ' ohne den Term −μ, μ erscheint als waagerechte ' +
    'Linie. Stufe 3 zeigt ±E(k); dort steckt μ in ξ_k = k²/2 − μ, und die Nulllinie ' +
    'ist das Fermi-Niveau. Ebenso ist es in ' + MV.fig('wire.stagesFig') + ' der Arbeit zwischen den ' +
    'Teilbildern (b) und (c) gezeichnet.';

  /* CLAUDE.md §2: Punkt 6 nennt die dort gültige Bedeutung von φ. In Kapitel 4
     ist φ die supraleitende Phase, fest null; im Modell (4.1) taucht sie nicht
     auf, weil Δ dort reell angesetzt ist. */
  var PHI_NOTE =
    'φ ist auch in diesem Reiter die supraleitende Phase und fest auf null ' +
    'gesetzt; im Modell wird Δ reell angesetzt, deshalb erscheint φ in keiner ' +
    'geplotteten Größe und an keinem Regler.';

  var UNITS =
    'ℏ = m = 1; alle Größen sind dimensionslos und auf Δ bezogen.';

  function baseParams() {
    var d = derive(state);
    return '<p>μ = <span class="num">' + MV.fmt(state.mu, 2) + '</span>, ' +
      'E_Z = <span class="num">' + MV.fmt(state.ez, 2) + '</span>, ' +
      'α = <span class="num">' + MV.fmt(state.alpha, 2) + '</span>, ' +
      'Δ = <span class="num">' + MV.fmt(state.delta, 2) + '</span>. ' + UNITS + '</p>' +
      '<p>Daraus: √(Δ² + μ²) = <span class="num">' + MV.fmt(d.boundary, 2) + '</span>, ' +
      'E(0) = |E_Z − √(Δ² + μ²)| = <span class="num">' + MV.fmt(d.e0, 2) + '</span>, ' +
      'Kriterium E_Z &gt; √(Δ² + μ²) ' +
      (d.critical ? '<b>im Grenzfall</b>' : (d.topological ? '<b>erfüllt</b>' : '<b>nicht erfüllt</b>')) +
      ', Fermi-Punkte: <span class="num">' + d.fermi.length + '</span>.</p>';
  }

  var INFO_A = {
    what:
      '<p>Bandstruktur des Drahtes über dem Impuls k, in drei zuschaltbaren Stufen ' +
      'entsprechend ' + MV.fig('wire.stagesFig') + ' der Arbeit: nur Rashba-Kopplung, zusätzlich ' +
      'Zeeman-Feld, zusätzlich Paarung.</p>' +
      '<p>In Stufe 1 und 2 ist senkrecht ε(k) aufgetragen und μ als waagerechte ' +
      'gestrichelte Linie eingezeichnet; die Schnittpunkte mit den Bändern sind die ' +
      'Fermi-Punkte. In Stufe 2 markiert das hinterlegte Fenster die Lücke der Größe ' +
      '2·E_Z, die das Feld bei k = 0 öffnet, und die grau gestrichelten Kurven zeigen ' +
      'die Bänder ohne Feld. In Stufe 3 ist ±E(k) aufgetragen.</p>',
    model:
      '<p>Rashba-Nanodraht, Gl. ' + MV.eq('wire.model') + ' der Arbeit:</p>' +
      '<span class="eq">H = ∫dx ψ† ( −∂ₓ²/2m − μ − i α σ_y ∂ₓ + E_Z σ_z ) ψ\n' +
      '  + ∫dx ( Δ ψ↑ψ↓ + h.c. )</span>' +
      '<p>Spin-Bahn-Achse ∥ ŷ, Zeeman-Feld senkrecht dazu; nur die Orthogonalität ' +
      'beider Achsen zählt.</p>' +
      '<p>Stufe 3 verwendet zusätzlich: ' + ALICEA_60 + '</p>',
    computed:
      '<p>Stufe 1 und 2 — nackte Bänder nach Gl. ' + MV.eq('wire.bands') + ', ohne den Term −μ:</p>' +
      '<span class="eq">ε±(k) = k²/2 ± sqrt( α²k² + E_Z² )</span>' +
      '<p>Stufe 1 setzt darin E_Z = 0, unabhängig vom Regler. Bei k = 0 öffnet das ' +
      'Feld eine Lücke der Größe 2·E_Z, nicht E_Z.</p>' +
      '<p>Fermi-Punkte sind die Lösungen von ε±(k) = μ.</p>' +
      '<p>Stufe 3 — Quasiteilchenspektrum von ' + MV.eq('wire.model') + ':</p>' +
      '<span class="eq">ξ_k = k²/2 − μ\n' +
      'E±(k)² = ξ_k² + α²k² + E_Z² + Δ²\n' +
      '         ± 2 sqrt( ξ_k² (α²k² + E_Z²) + E_Z² Δ² )</span>' +
      '<p>Bei k = 0 liefert das E₋(0) = |E_Z − √(Δ² + μ²)|, also genau Gl. ' + MV.eq('wire.criterion') + '.</p>' +
      '<p>' + AXIS_SWITCH + '</p>',
    params: function () {
      var st = ['1 — nur Rashba (E_Z = 0 gesetzt)', '2 — plus Zeeman', '3 — plus Paarung'][state.stage];
      var sEff = stageState(state);
      var fermi = state.stage === 0 ? fermiPoints(sEff) : derive(state).fermi;
      var K = MV.fmt(windowK(state, fermi), 1);
      return baseParams() +
        '<p>Angezeigte Stufe: <b>' + st + '</b>. Dargestellter Impulsbereich: ' +
        'k ∈ [−<span class="num">' + K + '</span>, <span class="num">' + K + '</span>].' +
        (fermi.length
          ? ' In dieser Stufe markierte Fermi-Punkte (bei E_Z = <span class="num">' +
            MV.fmt(sEff.ez, 2) + '</span>): k = <span class="num">' +
            fermi.map(function (q) { return MV.fmtSigned(q.k, 2); }).join(', ') + '</span>.'
          : ' Keine Fermi-Punkte: μ liegt unterhalb des Bandbodens.') +
        (state.stage === 0
          ? ' Das Anzeigefeld über den Panels nennt dagegen die Anzahl für den ' +
            'eingestellten Wert E_Z = <span class="num">' + MV.fmt(state.ez, 2) + '</span>.'
          : '') +
        '</p>';
    },
    numerics:
      '<ul>' +
      '<li>601 Stützstellen je Kurve.</li>' +
      '<li>Die Fermi-Punkte werden nicht am Gitter gesucht. Mit x = k² wird ' +
      'x² − 4(μ + α²)x + 4(μ² − E_Z²) = 0 exakt gelöst; das Quadrieren kann ' +
      'Scheinlösungen erzeugen, deshalb wird jeder Kandidat in die ' +
      'Ausgangsgleichung eingesetzt und nur bei einem Residuum unter ' +
      '10⁻⁷·(1 + |μ| + E_Z) übernommen.</li>' +
      '<li>Die Energieachse hängt nicht von μ ab, sondern nur von E_Z und vom ' +
      'Bandboden; sie weitet sich lediglich, wenn μ sonst aus dem Bild liefe. ' +
      'Damit behält das Zeeman-Fenster beim Verschieben von μ seine Höhe auf dem ' +
      'Bildschirm — es ist stets 2·E_Z hoch.</li>' +
      '<li>Der dargestellte Impulsbereich wird dagegen nachgeführt, in Stufen von ' +
      MV.fmt(0.5, 1) + ', damit alle Fermi-Punkte enthalten bleiben. Deren Lage hängt von μ ab, ' +
      'anders als die Höhe des Zeeman-Fensters; die waagerechte Skala ändert sich ' +
      'beim Verschieben von μ deshalb mit.</li>' +
      '<li>In Stufe 3 richtet sich die Energieachse nach dem unteren Zweig E₋; der ' +
      'obere Zweig E₊ verlässt den Ausschnitt, wie in ' + MV.fig('wire.stage3Fig') + '.</li>' +
      '<li>Beim Laden der Seite laufen über ein Raster aus 225 Parametersätzen vier ' +
      'automatische Prüfungen: (a) die geschlossene Form muss bei k = 0 den Wert ' +
      'aus Gl. ' + MV.eq('wire.criterion') + ' liefern; (b) im Grenzfall Δ → 0 muss |E±(k)| in |ε±(k)| aus ' +
      '' + MV.eq('wire.bands') + ' übergehen — die Abweichung fällt dabei linear mit Δ; (c) Δ_eff muss ' +
      'ungerade in k sein und bei k = 0 verschwinden; (d) mit μ innerhalb des ' +
      'Zeeman-Fensters müssen genau zwei, oberhalb genau vier Fermi-Punkte ' +
      'herauskommen. Schlägt eine Prüfung fehl, erscheint über den Reglern eine ' +
      'Warnung und die betroffene Darstellung ist nicht belastbar.</li>' +
      '</ul>',
    convention:
      '<p>' + UNITS + ' Es werden keine physikalischen Einheiten angezeigt.</p>' +
      '<p>' + PHI_NOTE + '</p>' +
      '<p>' + ALICEA_60 + '</p>' +
      '<p>' + AXIS_SWITCH + '</p>' +
      '<p>Die Bänder in Stufe 1 und 2 sind nach ihrer Energie beschriftet: ε₊ ist ' +
      'stets das obere, ε₋ stets das untere Band. Bei E_Z = 0 berühren sich beide bei ' +
      'k = 0, wo die Spinentartung nicht aufgehoben ist.</p>',
    reference:
      '<p>' + MV.sec('4.1', '4.2') + ' der Arbeit, ' + MV.fig('wire.stagesFig') + ' (a), (b), (c). ' +
      'Verwendete Gleichungen: ' + MV.eq('wire.model') + ', ' + MV.eq('wire.bands') + ', ' + MV.eq('wire.criterion') + '; für Stufe 3 zusätzlich ' +
      'Alicea (2012), Gl. (60).</p>',
    reading:
      '<p><b>Schalte die Stufe von „nur Rashba" auf „plus Zeeman" und zieh dann μ ' +
      'in das Fenster zwischen −E_Z und +E_Z.</b> Die Anzeige „Fermi-Punkte" springt ' +
      'dabei von vier auf zwei — das ist der Übergang, um den es geht.</p>' +
      '<p>Der Kern der Konstruktion ist der Übergang von vier auf zwei Fermi-Punkte. ' +
      'Ohne Feld schneidet μ oberhalb des Bandbodens stets vier Punkte, das System ' +
      'bleibt spinbehaftet. Das Zeeman-Feld öffnet bei k = 0 eine Lücke der Größe ' +
      '2·E_Z; wandert μ in dieses Fenster hinein, ist nur noch das untere Band ' +
      'besetzt und es bleiben zwei Fermi-Punkte — effektiv spinlos.</p>' +
      '<p>In Stufe 3 öffnet die Paarung an diesen Punkten eine Lücke. Die ' +
      'Darstellungen sind stationäre Gleichgewichtszustände; das Verschieben eines ' +
      'Reglers ist als quasistatische Parametervariation zu lesen.</p>'
  };

  var INFO_B = {
    what:
      '<p>Effektive Paarungsamplitude Δ_eff im unteren Band über dem Impuls k. ' +
      'Aufgetragen sind der exakte Ausdruck und seine lineare Näherung für kleine ' +
      'αk. Die senkrechte Linie markiert den mit dem Regler gewählten ' +
      'Auswertungspunkt; unten rechts steht die relative Abweichung dort.</p>',
    model:
      '<p>Projektion des Singulett-Paarungsterms aus Gl. ' + MV.eq('wire.model') + ' auf das untere Band. ' +
      'Der Spinor dieses Bandes zeigt entlang −d̂(k) mit</p>' +
      '<span class="eq">d(k) = ( 0, −α k, E_Z )</span>' +
      '<p>und dreht sich damit stetig mit k: bei k = 0 steht er entlang −ẑ, für große ' +
      '|k| kippt er in ∓ŷ-Richtung. Genau diese Drehung macht die Singulett-Paarung ' +
      'im helikalen Regime brauchbar.</p>',
    computed:
      '<p>Gl. ' + MV.eq('wire.deltaEff') + ' der Arbeit:</p>' +
      '<span class="eq">Δ_eff(k) = Δ · α k / sqrt( E_Z² + α²k² )\n\n' +
      'Näherung für α k ≪ E_Z:   Δ_eff(k) ≈ (Δ α / E_Z) · k</span>' +
      '<p>Beide sind ungerade in k und verschwinden bei k = 0. Die relative ' +
      'Abweichung der Näherung lässt sich geschlossen angeben, weil sich Δ, α und k ' +
      'im Quotienten herauskürzen:</p>' +
      '<span class="eq">Δ_approx / Δ_exakt = sqrt( E_Z² + α²k² ) / E_Z\n\n' +
      'relative Abweichung = sqrt( 1 + α²k²/E_Z² ) − 1</span>',
    params: function () {
      var dev = deltaEffDeviation(state.k, state);
      return baseParams() +
        '<p>Auswertungspunkt k = <span class="num">' + MV.fmt(state.k, 2) + '</span>. ' +
        'Dort: Δ_eff exakt = <span class="num">' + MV.fmt(deltaEff(state.k, state), 3) + '</span>, ' +
        (state.ez < EZ_MIN_APPROX
          ? 'Näherung nicht definiert (E_Z &lt; ' + MV.fmt(EZ_MIN_APPROX, 2) + ').'
          : 'genähert = <span class="num">' + MV.fmt(deltaEffApprox(state.k, state), 3) + '</span>, ' +
            'relative Abweichung <span class="num">' + MV.fmt(100 * dev, 1) + ' %</span>.') +
        '</p>' +
        '<p>Maximalwert der exakten Kurve im dargestellten Bereich: ' +
        '<span class="num">' + MV.fmt(deltaEff(K_PANEL_B, state), 3) + '</span> bei k = ' +
        K_PANEL_B + '; für |k| → ∞ strebt Δ_eff gegen Δ = ' +
        '<span class="num">' + MV.fmt(state.delta, 2) + '</span>.</p>';
    },
    numerics:
      '<ul>' +
      '<li>Fester Impulsbereich k ∈ [−3, 3] mit 601 Stützstellen; der Regler für den ' +
      'Auswertungspunkt läuft über denselben Bereich. Anders als Panel A wird hier ' +
      'nicht nachgeführt, weil die Aussage eine Aussage über kleine k ist.</li>' +
      '<li>Die exakte Kurve wird direkt aus ' + MV.eq('wire.deltaEff') + ' ausgewertet, die Näherung als ' +
      'Gerade durch zwei Punkte gezeichnet — sie ist exakt linear.</li>' +
      '<li>Die relative Abweichung wird nicht aus den beiden Zahlenwerten gebildet, ' +
      'sondern geschlossen als sqrt(1 + α²k²/E_Z²) − 1. Dadurch ist sie auch bei ' +
      'k = 0 definiert, wo beide Kurven verschwinden.</li>' +
      '<li>Die Werte an den Fermi-Punkten werden gegen die Abnahmewerte des ' +
      'Auftrags gehalten: bei μ = ' + MV.fmt(-0.3, 1) + ', E_Z = ' + MV.fmt(1, 1) +
      ', α = ' + MV.fmt(1.4, 1) + ', Δ = ' + MV.fmt(0.6, 1) + ' müssen es zwei ' +
      'Fermi-Punkte bei k_F = ±' + MV.fmt(2.674, 3) + ' sein, mit Δ_eff = ' +
      MV.fmt(0.58, 3) + ' exakt gegen ' + MV.fmt(2.246, 3) + ' genähert, also ' +
      MV.fmt(287, 0) + ' % Abweichung. ' +
      (MV.nanowire && MV.nanowire.check && MV.nanowire.check.referenceFailures &&
       MV.nanowire.check.referenceFailures.length
        ? '<b>Diese Prüfung ist fehlgeschlagen: ' +
          MV.nanowire.check.referenceFailures.join('; ') + '.</b>'
        : 'Die Prüfung läuft beim Laden und ist bestanden.') + '</li>' +
      '<li>Für E_Z &lt; ' + MV.fmt(EZ_MIN_APPROX, 2) + ' wird die Näherung nicht ' +
      'gezeichnet: sie divergiert wie ' +
      'Δα/E_Z. Der exakte Ausdruck geht dort in Δ·sgn(k) über und springt bei k = 0; ' +
      'die Kurve wird deshalb an dieser Stelle unterbrochen, während Δ_eff(0) = 0 ' +
      'als Punkt erhalten bleibt, weil die Funktion ungerade sein muss.</li>' +
      '</ul>',
    convention:
      '<p>' + UNITS + '</p>' +
      '<p>Δ ist in Kapitel 4 proximity-induziert und extern vorgegeben, nicht ' +
      'selbstkonsistent bestimmt.</p>' +
      '<p>Beachte das Minuszeichen in d(k) = (0, −αk, E_Z): der Spinor des unteren ' +
      'Bandes zeigt entlang −d̂(k).</p>' +
      '<p>' + PHI_NOTE + '</p>',
    reference:
      '<p>' + MV.sec('4.2') + ' der Arbeit, Gl. ' + MV.eq('wire.deltaEff') + '.</p>',
    reading:
      '<p><b>Zieh E_Z von 0,5 auf 3 und beobachte, wie die Kurve flacher wird.</b> ' +
      'Der angezeigte Wert von Δ_eff fällt dabei monoton — ein starkes Feld richtet ' +
      'die Spins aus und unterdrückt die helikale Struktur.</p>' +
      '<p>Δ_eff ist ungerade in k und verschwindet bei k = 0 — das ist die ' +
      'Kontinuumsform der p-Wellen-Paarung aus Kapitel 3 und keine Zufälligkeit der ' +
      'Rechnung, sondern von der Antisymmetrie der Paarwellenfunktion erzwungen. Für ' +
      'α → 0 verschwindet ' + MV.eq('wire.deltaEff') + ' vollständig: die Spin-Bahn-Kopplung ist nicht bloß ' +
      'Hilfsmittel, sondern Ursache der Paarung.</p>' +
      '<p><b>Drück „k auf k_F setzen".</b> ' + fermiReadingText() + '</p>' +
      '<p>Entscheidend ist der Zielkonflikt: Δ_eff fällt mit wachsendem E_Z, weil ein ' +
      'starkes Feld die Spins ausrichtet und die helikale Struktur unterdrückt. Ein ' +
      'großes E_Z ist zugleich Voraussetzung für die topologische Phase (Panel C). ' +
      'Die Lücke des topologischen Supraleiters ist deshalb nicht durch Δ allein ' +
      'gegeben, sondern durch ein Optimum zwischen beiden Effekten.</p>'
  };

  var INFO_C = {
    what:
      '<p>Phasendiagramm in den Achsen μ/Δ waagerecht und E_Z/Δ senkrecht. Die Kurve ' +
      'ist die Phasengrenze; oberhalb liegt die topologische Phase ν = −1, unterhalb ' +
      'die triviale ν = +1. Der Marker zeigt den aktuellen Arbeitspunkt.</p>',
    model:
      '<p>Rashba-Nanodraht, Gl. ' + MV.eq('wire.model') + '. Die Phasengrenze folgt aus dem Schließen der ' +
      'Anregungslücke bei k = 0, Gl. ' + MV.eq('wire.criterion') + '.</p>',
    computed:
      '<span class="eq">E(0) = | E_Z − sqrt( Δ² + μ² ) |\n\n' +
      'topologisch  ⇔  E_Z &gt; sqrt( Δ² + μ² )        ' + MV.eq('wire.criterion') + '</span>' +
      '<p>In den Achsen des Diagramms, also nach Division durch Δ:</p>' +
      '<span class="eq">E_Z/Δ = sqrt( 1 + (μ/Δ)² )</span>' +
      '<p>Der tiefste Punkt der Parabel liegt bei μ = 0 und verlangt E_Z &gt; Δ.</p>',
    params: function () {
      var mx = state.mu / state.delta, my = state.ez / state.delta;
      var outside = mx < -3 || mx > 3 || my > 3;
      return baseParams() +
        '<p>Arbeitspunkt: μ/Δ = <span class="num">' + MV.fmt(mx, 2) + '</span>, ' +
        'E_Z/Δ = <span class="num">' + MV.fmt(my, 2) + '</span>. Phasengrenze an ' +
        'dieser Stelle: E_Z/Δ = <span class="num">' + MV.fmt(Math.sqrt(1 + mx * mx), 2) +
        '</span>.' +
        (outside ? ' Der Punkt liegt außerhalb des dargestellten Bereichs und wird am ' +
                   'Rand als offener Kreis gezeigt.' : '') +
        '</p>';
    },
    numerics:
      '<ul>' +
      '<li>Phasengrenze aus 401 Stützstellen.</li>' +
      '<li>Achsenbereiche fest: μ/Δ von −3 bis 3, E_Z/Δ von 0 bis 3, wie in ' +
      '' + MV.fig('wire.phaseFig') + '. Da μ und Δ unabhängig einstellbar sind, kann der Quotient ' +
      'μ/Δ diesen Bereich verlassen; der Marker wird dann an den Rand geklemmt und ' +
      'als offener Kreis gezeichnet, der Zahlenwert steht daneben.</li>' +
      '<li>Als Grenzfall gilt E(0) &lt; ' + MV.fmt(TOL, 2) + ', dieselbe Toleranz wie ' +
      'in Reiter 1.</li>' +
      '</ul>',
    convention:
      '<p>' + UNITS + ' Beide Achsen sind auf Δ bezogen, deshalb verschiebt eine ' +
      'Änderung von Δ den Arbeitspunkt, nicht die Kurve.</p>' +
      '<p>ν = −1 bezeichnet die topologische, ν = +1 die triviale Phase — dieselbe ' +
      'Konvention wie in Reiter 1.</p>' +
      '<p>' + PHI_NOTE + '</p>',
    reference:
      '<p>' + MV.sec('4.3') + ' der Arbeit, ' + MV.fig('wire.phaseFig') + ', Gl. ' + MV.eq('wire.criterion') + '.</p>',
    reading:
      '<p><b>Zieh μ bei festem E_Z nach außen und beobachte, wann der Arbeitspunkt ' +
      'die Parabel überquert.</b> Die Statuszeile wechselt dort von „erfüllt" auf ' +
      '„nicht erfüllt" — und zwar genau bei E_Z = √(Δ² + μ²).</p>' +
      '<p>Aus dem Diagramm folgen zwei Bedingungen von unterschiedlichem Charakter. ' +
      'Am tiefsten Punkt der Parabel verlangt das Kriterium E_Z &gt; Δ — eine ' +
      'Materialbedingung. Und μ muss über die gesamte Drahtlänge innerhalb des ' +
      'topologischen Bereichs gehalten werden — eine Kontrollbedingung. Räumliche ' +
      'Schwankungen von μ durch Unordnung können den Draht teilweise über die ' +
      'Phasengrenze schieben.</p>'
  };

  /* ====================================================================
     Construction
     ==================================================================== */

  function buildControls(mount) {
    var box = MV.ui.el('div', { class: 'controls' });

    function onSlider(key) {
      return function (v) { state[key] = v; refresh(); };
    }

    sliders.mu = MV.ui.slider(box, {
      label: 'μ — chemisches Potential',
      min: -3, max: 3, step: 0.01, value: state.mu, digits: 2,
      onInput: onSlider('mu')
    });
    sliders.ez = MV.ui.slider(box, {
      label: 'E_Z — Zeeman-Energie',
      min: 0, max: 3, step: 0.01, value: state.ez, digits: 2,
      onInput: onSlider('ez')
    });
    sliders.alpha = MV.ui.slider(box, {
      label: 'α — Rashba-Kopplung',
      min: 0.1, max: 2.0, step: 0.01, value: state.alpha, digits: 2,
      onInput: onSlider('alpha')
    });
    sliders.delta = MV.ui.slider(box, {
      label: 'Δ — induzierte Paarung',
      min: 0.1, max: 2.0, step: 0.01, value: state.delta, digits: 2,
      onInput: onSlider('delta')
    });

    var row = MV.ui.el('div', { class: 'controls__wide' });
    segStage = MV.ui.segmented(row, {
      label: 'Panel A, Aufbau in drei Stufen:',
      options: ['1 — nur Rashba', '2 — plus Zeeman', '3 — plus Paarung'],
      value: state.stage,
      onChange: function (i) { state.stage = i; refresh(); }
    });
    row.appendChild(MV.ui.el('p', {
      class: 'note',
      html: 'Stufe 1 setzt in Panel A E_Z = 0, unabhängig vom Regler. Panel B und C ' +
            'verwenden immer die eingestellten Werte. Alle Darstellungen sind ' +
            'stationäre Gleichgewichtszustände; Reglerbewegungen sind als ' +
            'quasistatische Parametervariation zu lesen.'
    }));
    box.appendChild(row);
    mount.appendChild(box);
  }

  function buildReadout(mount) {
    readout = MV.ui.readout(mount, {
      items: [
        { key: 'e0',        label: 'E(0) = |E_Z − √(Δ² + μ²)|' },
        { key: 'boundary',  label: '√(Δ² + μ²)' },
        { key: 'criterion', label: 'Kriterium E_Z > √(Δ² + μ²)' },
        { key: 'fermi',     label: 'Fermi-Punkte bei eingestelltem E_Z' }
      ],
      note:
        'Die Fermi-Punkte sind die des ungepaarten Systems, also die Lösungen von ' +
        'ε±(k) = μ nach ' + MV.eq('wire.bands') + '. In Stufe 3 ist das Spektrum gelückt; die Punkte sind ' +
        'dort nur als senkrechte Hilfslinien eingezeichnet. Stufe 1 rechnet mit ' +
        'E_Z = 0 und markiert deshalb im Panel die Fermi-Punkte des feldfreien ' +
        'Systems; die Zahl hier bezieht sich stets auf den eingestellten Wert.'
    });
  }

  /* ====================================================================
     Public interface
     ==================================================================== */

  MV.nanowire = {
    init: function (mount) {
      /* The closed-form spectrum is checked before anything is drawn. */
      var check = selfCheck();
      MV.nanowire.check = check;
      if (!check.passed) {
        var warningEl = MV.ui.el('div', {
          class: 'readout',
          html: '<b>Warnung:</b> Die automatische Prüfung der geschlossenen Form des ' +
                'Quasiteilchenspektrums ist fehlgeschlagen (größte Abweichung bei ' +
                'k = 0: ' + MV.fmtExp(check.maxDeviationAtZero, 2) + ', im Grenzfall ' +
                'Δ → 0: ' + MV.fmtExp(check.maxDeviationDeltaToZero, 2) +
                '). Stufe 3 von Panel A ist nicht belastbar.'
        });
        warningEl.style.borderColor = MV.palette.critical;
        warningEl.style.color = MV.palette.critical;
        mount.appendChild(warningEl);
      }

      buildControls(mount);
      buildReadout(mount);

      var panels = MV.ui.el('div', { class: 'panels' });
      mount.appendChild(panels);

      panelA = MV.ui.panel(panels, {
        title: 'Panel A — Bänder des Drahtes, Aufbau in drei Stufen',
        quote: 'nanodraht.A',
        lead: 'Stufe 1 und 2 zeigen ε(k) mit μ als Linie, Stufe 3 zeigt ±E(k) mit ' +
              'der Nulllinie als Fermi-Niveau.',
        wide: true,
        aspect: 0.42,
        info: INFO_A
      });

      panelB = MV.ui.panel(panels, {
        title: 'Panel B — effektive p-Wellen-Paarung Δ_eff(k)',
        quote: 'nanodraht.B',
        lead: 'Exakter Ausdruck nach Gl. ' + MV.eq('wire.deltaEff') + ' und lineare Näherung für α k ≪ E_Z.',
        aspect: 0.78,
        controls: true,
        info: INFO_B
      });

      panelC = MV.ui.panel(panels, {
        title: 'Panel C — Phasendiagramm',
        quote: 'nanodraht.C',
        lead: 'Topologische Phase oberhalb der Grenze E_Z = √(Δ² + μ²).',
        aspect: 0.78,
        info: INFO_C
      });

      sliders.k = MV.ui.slider(panelB.controls, {
        label: 'k — Auswertungspunkt für die Abweichung',
        min: -K_PANEL_B, max: K_PANEL_B, step: 0.01, value: state.k, digits: 2,
        onInput: function (v) { state.k = v; refresh(); }
      });

      /* Der Regler misst, wie gut die Näherung an einer beliebigen Stelle ist.
         Physikalisch zählt Δ_eff an den Fermi-Punkten — der Knopf springt
         dorthin. Liegt der äußerste außerhalb des Reglerbereichs, wird auf den
         Rand geklemmt und in der Statuszeile steht, wo er wirklich liegt. */
      MV.ui.buttonRow(panelB.controls, {
        buttons: [{
          text: 'k auf k_F setzen',
          onClick: function () {
            var fi = fermiInfo(state);
            if (!fi.rows.length) return;
            var kf = fi.rows[fi.rows.length - 1].k;
            state.k = Math.max(-K_PANEL_B, Math.min(K_PANEL_B, kf));
            sliders.k.set(state.k, true);
            refresh();
          }
        }]
      });

      /* Schnittstelle für die Vorführungen (demo.js). Die Stufe läuft über
         die vorhandene Schaltergruppe, damit sie auch optisch umspringt. */
      MV.demo.provide('nanodraht', MV.demo.handle({
        state: state,
        sliders: sliders,
        refresh: refresh,
        segments: { stage: segStage },
        extra: {
          /* Fermi-Punkte und Δ_eff für gedachte Werte — ohne den Reiter zu
             verstellen; so kann eine Endprüfung Anfang und Ende gegeneinander
             halten, nachdem der Durchlauf vorbei ist. */
          fermiCountAt: function (over) { return fermiPoints(withState(over)).length; },
          deltaEffAt: function (k, over) { return deltaEff(k, withState(over)); },
          boundaryAt: function (over) {
            var s = withState(over);
            return Math.sqrt(s.delta * s.delta + s.mu * s.mu);
          },
          topologicalAt: function (over) { return derive(withState(over)).topological; },
          read: function () { return derive(state); }
        }
      }));
    },

    draw: function () { refresh(); },

    stop: function () { /* nothing runs in the background in this tab */ }
  };

}(MV));
