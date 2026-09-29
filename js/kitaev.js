/* ------------------------------------------------------------------------
   kitaev.js — tab 1: Kitaev chain, section 3.3 of the thesis.

   Panel A: orbit of h(k) in the (h_y, h_z) plane      (figure 3.2c)
   Panel B: bulk spectrum ±E(k)                        (figure 3.2a)

   Conventions follow Alicea (2012) as adopted in the thesis:
     H = −μ Σ c†c − (1/2) Σ [ t c†_x c_{x+1} + Δ e^{iφ} c_x c_{x+1} + h.c. ]   (3.1)
   The factor 1/2 puts the phase transition at |μ| = t, not at 2t.
   The phase is fixed at φ = 0 throughout, hence h_x ≡ 0.
   ------------------------------------------------------------------------ */

var MV = window.MV || (window.MV = {});

(function (MV) {
  'use strict';

  var PI = Math.PI;

  /* Tolerance below which sgn(ε) — and with it ν — counts as undefined.
     The gap closes exactly when ε₀ = 0 or ε_π = 0, so the same tolerance
     decides when panel B reports a gap closure. */
  var TOL = 0.02;

  var N_CURVE = 720;    /* samples for both plotted curves (721 points) */
  var N_WIND = 2000;    /* angle increments for the winding number       */

  /* μ-Durchlauf: quasi-static parameter variation, not time evolution.
     One traverse from −2 to +2 takes SWEEP_SECONDS. */
  var SWEEP_MIN = -2, SWEEP_MAX = 2, SWEEP_SECONDS = 7;

  var state = { mu: 0, t: 1, delta: 1 };

  var sliders = {}, panelA = null, panelB = null, readout = null;
  var sweepButton = null;
  var sweep = { running: false, dir: 1, raf: null, last: 0 };

  /* ====================================================================
     Physics
     ==================================================================== */

  /* Bulk quantities, Gl. (3.7), with φ = 0:
       ε_k  = −t cos k − μ
       Δ̃_k = −i Δ sin k                    →  |Δ̃_k| = Δ |sin k|
     Pseudospin H_k = h(k)·σ:
       h_x = 0,  h_y = −Δ sin k,  h_z = ε_k                            */
  function hVector(k, s) {
    return { y: -s.delta * Math.sin(k), z: -s.t * Math.cos(k) - s.mu };
  }

  /* Excitation spectrum, Gl. (3.8):
       E(k) = sqrt( ε_k² + |Δ̃_k|² ) = sqrt( (t cos k + μ)² + Δ² sin²k ) */
  function energy(k, s) {
    var a = s.t * Math.cos(k) + s.mu;
    var b = s.delta * Math.sin(k);
    return Math.sqrt(a * a + b * b);
  }

  /**
   * Bulk gap Δ_bulk = min_k E(k) — exact, not a grid search.
   *
   * With u = cos k ∈ [−1, 1]:
   *   E² = (t² − Δ²) u² + 2tμ u + (μ² + Δ²) = a u² + b u + c
   * a > 0: convex, minimum at the vertex u* = −b/2a clamped to [−1, 1].
   * a < 0: concave, minimum at one of the two ends.
   * a = 0: linear, minimum at the end selected by the sign of b.
   *
   * Returns { value, us, ks, flat }. `ks` holds every k in [−π, π] that
   * attains the minimum, so that a gap closure at k = ±π marks both edges of
   * the displayed Brillouin zone.
   *
   * For a = b = 0 — that is t = Δ together with μ = 0, the topological limit
   * of Gl. (3.4) — the polynomial is constant and E(k) = Δ for every k. There
   * is then no distinguished minimum; `flat` says so and `ks` stays empty, so
   * that no arbitrary subset of momenta gets marked as "the" minimum.
   */
  function bulkGap(s) {
    var a = s.t * s.t - s.delta * s.delta;
    var b = 2 * s.t * s.mu;
    var c = s.mu * s.mu + s.delta * s.delta;
    var cands;

    if (Math.abs(a) < 1e-12 && Math.abs(b) < 1e-12) {
      return { value: Math.sqrt(Math.max(0, c)), us: [], ks: [], flat: true };
    }

    if (Math.abs(a) < 1e-12) {
      cands = (b > 0) ? [-1] : (b < 0) ? [1] : [-1, 1];
    } else if (a > 0) {
      cands = [Math.min(1, Math.max(-1, -b / (2 * a)))];
    } else {
      cands = [-1, 1];
    }

    var vals = cands.map(function (u) { return a * u * u + b * u + c; });
    var best = Math.min.apply(null, vals);
    var us = cands.filter(function (u, i) { return vals[i] - best < 1e-12; });

    var ks = [];
    us.forEach(function (u) {
      if (u >= 1 - 1e-12) ks.push(0);
      else if (u <= -1 + 1e-12) { ks.push(-PI); ks.push(PI); }
      else { var k = Math.acos(u); ks.push(-k); ks.push(k); }
    });

    return { value: Math.sqrt(Math.max(0, best)), us: us, ks: ks, flat: false };
  }

  /**
   * max_k E(k) for a given μ, exact, from the same polynomial as bulkGap.
   * Used to size the energy axis of panel B independently of the current μ.
   */
  function maxEnergy(s, mu) {
    var a = s.t * s.t - s.delta * s.delta;
    var b = 2 * s.t * mu;
    var c = mu * mu + s.delta * s.delta;
    var cands = [-1, 1];
    if (a < 0) {
      var u = -b / (2 * a);                    /* vertex is a maximum here */
      if (u > -1 && u < 1) cands.push(u);
    }
    var best = -Infinity;
    cands.forEach(function (u) {
      var v = a * u * u + b * u + c;
      if (v > best) best = v;
    });
    return Math.sqrt(Math.max(0, best));
  }

  /**
   * Winding number W of the closed orbit h(k) around the origin.
   *
   * Counted mathematically positive in the (h_y, h_z) plane, i.e. from the
   * angle θ = atan2(h_z, h_y) accumulated while k runs from −π to π. With
   * this choice W = −1 in the topological phase. Only the parity of W enters
   * ν = (−1)^W, so the sign convention has no physical consequence — see the
   * info panel.
   *
   * Returns null while the orbit runs through the origin (W undefined).
   */
  function winding(s) {
    var total = 0, i, k, h, th, prev = null, d;
    for (i = 0; i <= N_WIND; i++) {
      k = -PI + 2 * PI * i / N_WIND;
      h = hVector(k, s);
      if (Math.abs(h.y) < 1e-12 && Math.abs(h.z) < 1e-12) return null;
      th = Math.atan2(h.z, h.y);
      if (prev !== null) {
        d = th - prev;
        while (d > PI) d -= 2 * PI;
        while (d < -PI) d += 2 * PI;
        total += d;
      }
      prev = th;
    }
    return Math.round(total / (2 * PI));
  }

  /** Everything the two panels and the read-out need for one parameter set. */
  function derive(s) {
    var eps0 = -s.t - s.mu;      /* ε₀  = ε_k at k = 0  */
    var epsPi = s.t - s.mu;      /* ε_π = ε_k at k = π  */
    var critical = Math.abs(eps0) < TOL || Math.abs(epsPi) < TOL;

    var s0 = critical && Math.abs(eps0) < TOL ? 0 : (eps0 > 0 ? 1 : -1);
    var sPi = critical && Math.abs(epsPi) < TOL ? 0 : (epsPi > 0 ? 1 : -1);

    return {
      s: s,
      eps0: eps0,
      epsPi: epsPi,
      s0: s0,
      sPi: sPi,
      critical: critical,
      nu: critical ? null : s0 * sPi,          /* Gl. (3.11) */
      gap: bulkGap(s),
      W: critical ? null : winding(s)
    };
  }

  /* ====================================================================
     Formatting helpers
     ==================================================================== */

  function signLabel(v) {
    return v === 0 ? 'undefiniert' : (v > 0 ? '+1' : '−1');
  }

  function phaseName(d) {
    if (d.critical) return 'Übergang, ν undefiniert';
    return d.nu < 0 ? 'topologisch (ν = −1)' : 'trivial (ν = +1)';
  }

  function phaseKind(d) {
    return d.critical ? 'critical' : (d.nu < 0 ? 'topological' : 'trivial');
  }

  function phaseColor(d) {
    return MV.palette[phaseKind(d)];
  }

  /* ====================================================================
     Panel A — orbit of h(k) in the (h_y, h_z) plane   (figure 3.2c)
     ==================================================================== */

  function drawPanelA(d) {
    var p = panelA.plot, s = d.s, i, k, h;

    /* Same number of data units per pixel on both axes, so that the orbit at
       t = Δ shows up as a circle.
       The frame does not depend on μ. μ shifts the orbit rigidly and does not
       deform it; a frame that followed μ would rescale the picture and make
       the orbit appear to shrink and grow, which is exactly the opposite of
       what happens. It is sized to hold both the origin and the whole orbit
       over the range of the μ-Durchlauf, and only widens further if μ is
       dragged beyond that range. */
    var m = p.measure();
    var R = Math.max(
      MV.ceilTo(SWEEP_MAX + s.t, 0.5),
      MV.ceilTo(Math.abs(s.mu) + s.t, 0.5)
    );
    var X = R * m.pw / m.ph;

    p.begin({
      x: [-X, X], y: [-R, R],
      xLabel: 'h_y', yLabel: 'h_z',
      xTickCount: 5, yTickCount: 5
    });

    /* orbit */
    var pts = [];
    for (i = 0; i <= N_CURVE; i++) {
      k = -PI + 2 * PI * i / N_CURVE;
      h = hVector(k, s);
      pts.push([h.y, h.z]);
    }
    p.curve(pts, { color: phaseColor(d), width: 2, legend: 'h(k), Pfeil: wachsendes k' });

    /* direction of increasing k: tangent (dh_y, dh_z)/dk = (−Δ cos k, t sin k),
       drawn at k = ±π/2 where the orbit crosses h_z = −μ */
    [PI / 2, -PI / 2].forEach(function (kk) {
      var hh = hVector(kk, s);
      var dz = s.t * Math.sin(kk);
      p.arrow(hh.y, hh.z, 0, -Math.sign(dz) * 16, { color: phaseColor(d), width: 1.4 });
    });

    /* axes through the origin */
    p.hline(0, { color: MV.palette.line, width: 1 });
    p.vline(0, { color: MV.palette.line, width: 1 });

    /* The two special momenta always sit on the h_z axis, k = π above k = 0
       (their distance is 2t). Their labels normally point away from each
       other; when one of them comes close to the origin cross, that label is
       flipped into the gap between the two so the cross stays legible. */
    var near = 0.25 * R;
    p.marker(0, d.eps0, {
      label: 'k = 0', color: MV.palette.ink, r: 4.5,
      anchor: (d.eps0 >= 0 && d.eps0 < near) ? 'above' : 'below'
    });
    p.marker(0, d.epsPi, {
      label: 'k = π', color: MV.palette.ink, r: 4.5,
      anchor: (d.epsPi <= 0 && -d.epsPi < near) ? 'below' : 'above'
    });
    /* The origin is named in the legend rather than by a label at the marker:
       near the transition it coincides with k = 0 or k = π, and two labels on
       the same spot would cover each other. */
    p.marker(0, 0, {
      color: MV.palette.critical, shape: 'cross', r: 5, width: 2,
      legend: 'Ursprung'
    });

    /* the reading of the figure, in words, under the plot */
    panelA.setStatus(
      d.critical ? 'Bahn läuft durch den Ursprung — Übergang'
                 : (d.nu < 0 ? 'Ursprung innen: ν = −1, topologisch'
                             : 'Ursprung außen: ν = +1, trivial'),
      phaseColor(d)
    );

    p.finish();
  }

  /* ====================================================================
     Panel B — bulk spectrum ±E(k)   (figure 3.2a)
     ==================================================================== */

  function drawPanelB(d) {
    var p = panelB.plot, s = d.s, i, k, e;

    var pos = [], neg = [], maxE = 0;
    for (i = 0; i <= N_CURVE; i++) {
      k = -PI + 2 * PI * i / N_CURVE;
      e = energy(k, s);
      if (e > maxE) maxE = e;
      pos.push([k, e]);
      neg.push([k, -e]);
    }
    /* Energy axis independent of μ, for the same reason as in panel A: an axis
       that follows max E(k) would rescale during the μ-Durchlauf and hide how
       far the bands actually move. Sized to the largest value E(k) takes
       anywhere in the range of the Durchlauf, and widened only if μ is dragged
       beyond that range. */
    var Y = Math.max(
      MV.ceilTo(Math.max(maxEnergy(s, SWEEP_MIN), maxEnergy(s, SWEEP_MAX)), 0.5),
      MV.ceilTo(maxE, 0.5)
    );

    p.begin({
      x: [-PI, PI], y: [-Y, Y],
      xTicks: MV.piTicks(), yTickCount: 5,
      xLabel: 'k', yLabel: '±E(k)'
    });

    var closed = d.gap.value < TOL;

    /* mark the momenta that attain min_k E(k) */
    d.gap.ks.forEach(function (kk) {
      p.vline(kk, { color: closed ? MV.palette.critical : MV.palette.line, dash: [3, 3] });
    });

    p.curve(pos, { color: MV.palette.curve1, width: 2, legend: '+E(k)' });
    p.curve(neg, { color: MV.palette.curve2, width: 2, dash: [6, 3], legend: '−E(k)' });
    p.hline(0, { color: MV.palette.axis, dash: [4, 4], label: 'E = 0', labelAt: 'left' });

    /* Label and dimension line go to whichever minimum sits furthest from the
       edges of the zone, so that neither is clipped at the frame. */
    var iBest = 0, bestDist = -1;
    d.gap.ks.forEach(function (kk, i) {
      var dist = PI - Math.abs(kk);
      if (dist > bestDist) { bestDist = dist; iBest = i; }
    });

    d.gap.ks.forEach(function (kk, idx) {
      p.marker(kk, d.gap.value, {
        color: closed ? MV.palette.critical : MV.palette.ink,
        r: 4,
        /* When the gap is closed the banner in the corner names the location;
           a label at the marker would collide with the E = 0 caption. */
        label: (idx === iBest && !closed) ? 'min E(k)' : '',
        anchor: 'above'
      });
      if (!closed) p.marker(kk, -d.gap.value, { color: MV.palette.ink, r: 3 });
    });

    /* Δ_bulk drawn from E = 0 up to +min E(k) — never as 2·min E(k).
       In the flat case there is no distinguished momentum; the measure is
       then placed at k = 0 purely for legibility. */
    if (!closed) {
      var kDim = d.gap.flat ? 0 : d.gap.ks[iBest];
      p.dimension(kDim, 0, d.gap.value,
        'Δ_bulk = ' + MV.fmt(d.gap.value, 2),
        { color: MV.palette.ink, side: kDim > 0.5 ? 'left' : 'right' });
    }

    panelB.setStatus(
      closed
        ? 'Lückenschluss bei ' + (Math.abs(d.eps0) < TOL ? 'k = 0' : 'k = ±π')
        : d.gap.flat
          ? 'E(k) ist konstant gleich Δ: bei t = Δ und μ = 0 ist die Lücke ' +
            'impulsunabhängig, Δ_bulk = ' + MV.fmt(d.gap.value, 2) +
            ' wird bei jedem k angenommen'
          : 'Lücke Δ_bulk = ' + MV.fmt(d.gap.value, 2) + ', Minimum bei k = ' +
            d.gap.ks.slice().sort(function (a, b) { return a - b; })
              .map(function (kk) { return MV.fmtSigned(kk / PI, 2) + '·π'; }).join(', '),
      closed ? MV.palette.critical : MV.palette.muted
    );

    p.finish();
  }

  /* ====================================================================
     Read-out and update cycle
     ==================================================================== */

  function refresh() {
    var d = derive(state);

    readout.set('eps0', MV.fmtSigned(d.eps0, 2));
    readout.set('epsPi', MV.fmtSigned(d.epsPi, 2));
    readout.set('s0', signLabel(d.s0));
    readout.set('sPi', signLabel(d.sPi));
    readout.set('nu', d.critical ? 'undefiniert' : MV.fmtSigned(d.nu, 0));
    readout.set('gap', MV.fmt(d.gap.value, 2));
    readout.set('W', d.W === null ? 'undefiniert' : MV.fmtSigned(d.W, 0));
    readout.status(phaseKind(d), phaseName(d));

    drawPanelA(d);
    drawPanelB(d);

    panelA.info.update();
    panelB.info.update();
    return d;
  }

  /* ====================================================================
     μ-Durchlauf (quasi-static parameter variation)
     ==================================================================== */

  function setSweep(on) {
    if (on === sweep.running) return;
    sweep.running = on;
    sweepButton.setAttribute('aria-pressed', on ? 'true' : 'false');
    sweepButton.textContent = on ? 'μ-Durchlauf anhalten' : 'μ-Durchlauf (quasistatisch)';
    if (on) {
      if (state.mu >= SWEEP_MAX) { state.mu = SWEEP_MAX; sweep.dir = -1; }
      else if (state.mu <= SWEEP_MIN) { state.mu = SWEEP_MIN; sweep.dir = 1; }
      sliders.mu.set(state.mu, true);
      sweep.last = 0;
      sweep.raf = window.requestAnimationFrame(step);
    } else if (sweep.raf) {
      window.cancelAnimationFrame(sweep.raf);
      sweep.raf = null;
    }
  }

  function step(ts) {
    if (!sweep.running) return;
    var dt = sweep.last ? Math.min(0.1, (ts - sweep.last) / 1000) : 0;
    sweep.last = ts;
    var v = state.mu + sweep.dir * ((SWEEP_MAX - SWEEP_MIN) / SWEEP_SECONDS) * dt;
    if (v >= SWEEP_MAX) { v = SWEEP_MAX; sweep.dir = -1; }
    if (v <= SWEEP_MIN) { v = SWEEP_MIN; sweep.dir = 1; }
    state.mu = v;
    sliders.mu.set(v, true);       /* silent: does not fire onInput */
    refresh();
    sweep.raf = window.requestAnimationFrame(step);
  }

  /* ====================================================================
     Info panels
     ==================================================================== */

  /** Shared live parameter block; `tail` adds the part specific to one panel. */
  function paramLine(tail) {
    var d = derive(state);
    return '<p>μ = <span class="num">' + MV.fmt(state.mu, 2) + '</span>, ' +
      't = <span class="num">' + MV.fmt(state.t, 2) + '</span>, ' +
      'Δ = <span class="num">' + MV.fmt(state.delta, 2) + '</span>, ' +
      'φ = 0 (fest). Alle Größen dimensionslos in Einheiten von t.</p>' +
      '<p>Daraus: ε₀ = <span class="num">' + MV.fmtSigned(d.eps0, 2) + '</span>, ' +
      'ε_π = <span class="num">' + MV.fmtSigned(d.epsPi, 2) + '</span>, ' +
      's₀ = <span class="num">' + signLabel(d.s0) + '</span>, ' +
      's_π = <span class="num">' + signLabel(d.sPi) + '</span>, ' +
      'ν = <span class="num">' + (d.critical ? 'undefiniert' : MV.fmtSigned(d.nu, 0)) + '</span>, ' +
      'Δ_bulk = <span class="num">' + MV.fmt(d.gap.value, 2) + '</span>, ' +
      'W = <span class="num">' + (d.W === null ? 'undefiniert' : MV.fmtSigned(d.W, 0)) + '</span>.</p>' +
      tail(d);
  }

  function paramsA() {
    return paramLine(function () {
      return '<p>Halbachsen der Ellipse: Δ = <span class="num">' + MV.fmt(state.delta, 2) +
        '</span> in h_y-Richtung, t = <span class="num">' + MV.fmt(state.t, 2) +
        '</span> in h_z-Richtung; Zentrum bei (0, ' + MV.fmt(-state.mu, 2) + ').</p>';
    });
  }

  function paramsB() {
    return paramLine(function (d) {
      if (d.gap.flat) {
        return '<p>Wegen t = Δ und μ = 0 ist E(k) impulsunabhängig gleich Δ; das ' +
          'Minimum wird bei jedem k angenommen. Das ist der topologische Grenzfall ' +
          'aus Gl. ' + MV.eq('kitaev.sweetSpot') + ' der Arbeit.</p>';
      }
      var ks = d.gap.ks.slice().sort(function (a, b) { return a - b; })
        .map(function (k) { return MV.fmtSigned(k / PI, 2) + '·π'; }).join(', ');
      return '<p>Das Minimum von E(k) wird angenommen bei k = <span class="num">' + ks +
        '</span>' + (d.gap.value < TOL ? ' (Lückenschluss)' : '') + '.</p>';
    });
  }

  var W_NOTE =
    'Umlaufsinn: h(k) für k von −π nach π, mathematisch positiv gezählt in der ' +
    '(h_y, h_z)-Ebene. Das Vorzeichen von W hängt von Achsenreihenfolge und ' +
    'Umlaufsinn ab; in ν = (−1)^W geht nur die Parität ein, sodass es ' +
    'physikalisch ohne Belang ist. W ist nur bei φ = 0 definiert (Klasse BDI); ' +
    'maßgeblich ist ν (Klasse D).';

  var ALICEA_NOTE =
    'Normierung nach Alicea (2012); der Faktor 1/2 vor Hopping- und Paarungsterm ' +
    'in Gl. ' + MV.eq('kitaev.hamiltonian') + ' legt den Phasenübergang auf |μ| = t, nicht auf 2t.';

  var INFO_A = {
    what:
      '<p>Bahn des Vektors h(k) in der (h_y, h_z)-Ebene, während der Impuls k die ' +
      'Brillouin-Zone von −π bis π durchläuft. Waagerecht ist h_y aufgetragen, ' +
      'senkrecht h_z. Der Ursprung ist als Kreuz markiert, die Punkte k = 0 und ' +
      'k = π als Punkte auf der h_z-Achse. Der Pfeil gibt die Richtung wachsender ' +
      'k an.</p>',
    model:
      '<p>Kitaev-Kette, Gl. ' + MV.eq('kitaev.hamiltonian') + ' der Arbeit:</p>' +
      '<span class="eq">H = −μ Σ_x c†_x c_x − (1/2) Σ_x [ t c†_x c_{x+1} + Δ e^{iφ} c_x c_{x+1} + h.c. ]</span>' +
      '<p>Im Impulsraum lautet der Bloch-Hamiltonian H_k = h(k)·σ mit den ' +
      'Bulk-Größen aus Gl. ' + MV.eq('kitaev.bulk') + '.</p>',
    computed:
      '<span class="eq">ε_k  = −t cos k − μ\n' +
      'Δ̃_k = −i Δ e^{iφ} sin k        (Gl. 3.7)\n\n' +
      'h_x(k) = 0            (wegen φ = 0)\n' +
      'h_y(k) = −Δ sin k\n' +
      'h_z(k) = −t cos k − μ = ε_k</span>' +
      '<p>Die Bahn ist eine Ellipse mit den Halbachsen Δ in h_y-Richtung und t in ' +
      'h_z-Richtung, zentriert bei (0, −μ). Die beiden markierten Punkte liegen ' +
      'stets auf der h_z-Achse, weil h_y bei k = 0 und k = π verschwindet:</p>' +
      '<span class="eq">k = 0:  h = (0, ε₀),   ε₀  = −t − μ\n' +
      'k = π:  h = (0, ε_π),  ε_π =  t − μ</span>' +
      '<p>Invariante nach Gl. ' + MV.eq('kitaev.signs') + ' und ' + MV.eq('kitaev.invariant') + ':</p>' +
      '<span class="eq">s₀ = sgn ε₀,   s_π = sgn ε_π,   ν = s₀ · s_π</span>' +
      '<p>ν = −1 topologisch, ν = +1 trivial. Windungszahl aus derselben Bahn:</p>' +
      '<span class="eq">W = (1/2π) ∮ dθ,   θ(k) = atan2( h_z(k), h_y(k) )</span>',
    params: paramsA,
    numerics:
      '<ul>' +
      '<li>Bahn aus 721 Stützstellen, Δk = 2π/720; die Kurve ist exakt ausgewertet, ' +
      'nicht interpoliert.</li>' +
      '<li>W numerisch als Summe von 2000 Winkelinkrementen, jeweils auf das ' +
      'Intervall (−π, π] reduziert, anschließend auf die nächste ganze Zahl gerundet.</li>' +
      '<li>Übergangstoleranz |ε| &lt; ' + MV.fmt(TOL, 2) + ': darunter gilt sgn ε als ' +
      'undefiniert, und ' +
      'ν sowie W werden als „undefiniert" ausgewiesen statt als 0.</li>' +
      '<li>Beide Achsen haben dieselbe Anzahl Einheiten pro Pixel, damit die Bahn ' +
      'bei t = Δ als Kreis erscheint.</li>' +
      '<li>Der dargestellte Ausschnitt hängt <b>nicht von μ ab</b>. Er ist so ' +
      'bemessen, dass er über den gesamten Bereich des μ-Durchlaufs sowohl den ' +
      'Ursprung als auch die ganze Bahn enthält, und weitet sich nur, wenn μ über ' +
      'diesen Bereich hinaus geschoben wird. Dadurch behält die Bahn ihre Größe ' +
      'auf dem Bildschirm und wird durch μ ausschließlich verschoben — so, wie es ' +
      'die Rechnung sagt. Ein mitwachsender Ausschnitt würde die Bahn scheinbar ' +
      'schrumpfen lassen und damit das Gegenteil suggerieren.</li>' +
      '</ul>',
    convention:
      '<p>' + ALICEA_NOTE + '</p>' +
      '<p>Die Phase ist fest auf φ = 0 gesetzt. Damit ist Δ̃_k rein imaginär, es gilt ' +
      'h_x ≡ 0, und die Bahn liegt vollständig in der (h_y, h_z)-Ebene.</p>' +
      '<p>' + W_NOTE + '</p>' +
      '<p>Abweichung von der gedruckten Abbildung: ' + MV.figShort('kitaev.orbitFig') + ' trägt h_y/t gegen h_z/t ' +
      'auf, weil dort Δ = t fest gewählt ist. Hier ist t ein Regler; die Achsen ' +
      'zeigen deshalb h_y und h_z unnormiert, damit eine Änderung von t sichtbar ' +
      'bleibt.</p>',
    reference:
      '<p>' + MV.sec('3.3') + ' der Arbeit, ' + MV.fig('kitaev.orbitFig') + '. Verwendete Gleichungen: ' +
      '' + MV.eq('kitaev.hamiltonian') + ', ' + MV.eq('kitaev.bulk') + ', ' + MV.eq('kitaev.signs') + ', ' + MV.eq('kitaev.invariant') + '.</p>',
    reading:
      '<p><b>Zieh μ von 0 nach 1,5 und beobachte, wann die Bahn den Ursprung ' +
      'überstreicht.</b> Genau dort springt die Anzeige von ν = −1 auf ν = +1.</p>' +
      '<p>μ verschiebt die Ellipse vertikal, verformt sie nicht: die Halbachsen ' +
      'bleiben Δ und t. Umschließt die Bahn den Ursprung, haben ε₀ und ε_π ' +
      'entgegengesetztes Vorzeichen und es gilt ν = −1; liegt der Ursprung außerhalb, ' +
      'ist ν = +1. Genau am Übergang läuft die Bahn durch den Ursprung — dort ist die ' +
      'Richtung von h undefiniert und die Invariante kann springen.</p>'
  };

  var INFO_B = {
    what:
      '<p>Anregungsspektrum ±E(k) über der Brillouin-Zone k ∈ [−π, π]. Die Linie ' +
      'E = 0 ist eingezeichnet. Senkrechte gestrichelte Linien und Punkte markieren ' +
      'die Impulse, an denen E(k) sein Minimum annimmt; die Maßlinie gibt den ' +
      'Zahlenwert Δ_bulk an. Schließt die Lücke, sind alle betroffenen Stellen ' +
      'markiert.</p>',
    model:
      '<p>Kitaev-Kette, Gl. ' + MV.eq('kitaev.hamiltonian') + ' der Arbeit, mit den Bulk-Größen aus Gl. ' + MV.eq('kitaev.bulk') + ':</p>' +
      '<span class="eq">ε_k  = −t cos k − μ\nΔ̃_k = −i Δ e^{iφ} sin k</span>',
    computed:
      '<p>Anregungsspektrum, Gl. ' + MV.eq('kitaev.spectrum') + ':</p>' +
      '<span class="eq">E(k) = sqrt( ε_k² + |Δ̃_k|² ) = sqrt( (t cos k + μ)² + Δ² sin²k )</span>' +
      '<p>Lücke:</p>' +
      '<span class="eq">Δ_bulk = min_k E(k)</span>' +
      '<p>Kritische Linie, Gl. ' + MV.eq('kitaev.criticalLine') + ': |μ| = t. Der Lückenschluss tritt bei μ = −t ' +
      'an k = 0 auf (dort ε₀ = −t − μ = 0) und bei μ = +t an k = ±π ' +
      '(dort ε_π = t − μ = 0).</p>',
    params: paramsB,
    numerics:
      '<ul>' +
      '<li>Kurve aus 721 Stützstellen, Δk = 2π/720.</li>' +
      '<li>Δ_bulk wird nicht am Gitter gesucht, sondern geschlossen bestimmt. Mit ' +
      'u = cos k ∈ [−1, 1] ist E² = (t² − Δ²)u² + 2tμ u + (μ² + Δ²) ein Polynom ' +
      'zweiten Grades in u: für t &gt; Δ liegt das Minimum am auf [−1, 1] ' +
      'beschnittenen Scheitel u* = −tμ/(t² − Δ²), für t ≤ Δ an einem der beiden ' +
      'Ränder. Die Umkehrung k = ±arccos u liefert alle minimierenden Impulse; bei ' +
      'u = −1 sind das beide Ränder der dargestellten Zone.</li>' +
      '<li>Verschwinden beide Koeffizienten, also t = Δ zusammen mit μ = 0, so ist ' +
      'das Polynom konstant und E(k) = Δ für jedes k. Dann gibt es kein ' +
      'ausgezeichnetes Minimum; es werden keine einzelnen Impulse markiert, und die ' +
      'Maßlinie steht nur der Lesbarkeit halber bei k = 0.</li>' +
      '<li>Als Lückenschluss gilt Δ_bulk &lt; ' + MV.fmt(TOL, 2) + ', dieselbe Toleranz wie für ' +
      'sgn ε in Panel A. Beides ist konsistent, weil E(k) = 0 nur bei sin k = 0 ' +
      'und ε_k = 0 zugleich möglich ist.</li>' +
      '<li>Die Energieachse hängt <b>nicht von μ ab</b>. Sie ist auf den größten ' +
      'Wert bemessen, den E(k) irgendwo im Bereich des μ-Durchlaufs annimmt; dieser ' +
      'wird aus demselben Polynom in u exakt bestimmt wie das Minimum. Sie weitet ' +
      'sich nur, wenn μ über diesen Bereich hinaus geschoben wird. Eine mitwachsende ' +
      'Achse würde die Bänder scheinbar stehen lassen, während sie sich tatsächlich ' +
      'verschieben, und den Lückenschluss dadurch verharmlosen.</li>' +
      '</ul>',
    convention:
      '<p>E(k) ist das Anregungsspektrum: die Energie, ein einzelnes Quasiteilchen zu ' +
      'erzeugen. Die Lücke ist deshalb Δ_bulk = min_k E(k) und ' +
      '<b>nicht</b> 2·min_k E(k). Dargestellt sind ±E(k); der angezeigte Zahlenwert ' +
      'bleibt min_k E(k), und die Maßlinie reicht von E = 0 bis +min_k E(k).</p>' +
      '<p>' + ALICEA_NOTE + '</p>' +
      '<p>φ ist hier die supraleitende Phase des Paarungsterms aus ' + MV.eqn('kitaev.hamiltonian') +
      '; sie ist fest auf null gesetzt und erscheint an keinem Regler. Mit φ = 0 ist ' +
      'Δ̃_k rein imaginär und h_x verschwindet.</p>',
    reference:
      '<p>' + MV.sec('3.3') + ' der Arbeit, ' + MV.fig('kitaev.spectrumFig') + '. Verwendete Gleichungen: ' +
      '' + MV.eq('kitaev.hamiltonian') + ', ' + MV.eq('kitaev.bulk') + ', ' + MV.eq('kitaev.spectrum') + ', ' + MV.eq('kitaev.criticalLine') + '.</p>',
    reading:
      '<p><b>Zieh μ auf −1 und dann auf +1 und beobachte, wo sich die Lücke ' +
      'schließt:</b> bei μ = −1 am Zonenzentrum k = 0, bei μ = +1 am Zonenrand ' +
      'k = ±π. Der Zahlenwert Δ_bulk geht dabei beide Male auf null.</p>' +
      '<p>Die Lücke schließt ausschließlich bei |μ| = t, und zwar bei μ = −t am ' +
      'Zonenzentrum k = 0, bei μ = +t am Zonenrand k = ±π. Die beiden Ränder der ' +
      'dargestellten Zone sind derselbe Punkt der Brillouin-Zone, deshalb sind dort ' +
      'stets beide Stellen markiert.</p>' +
      '<p>Der μ-Durchlauf ist eine Folge statischer Gleichgewichtszustände, keine ' +
      'Zeitentwicklung. Am Lückenschluss verschwindet die Energieskala, die ' +
      'Adiabatizität bricht dort zwangsläufig (Landau-Zener); die Darstellung ' +
      'aufeinanderfolgender Grundzustände ist deshalb nicht als realer Vorgang zu ' +
      'lesen.</p>'
  };

  /* ====================================================================
     Construction
     ==================================================================== */

  function buildControls(mount) {
    var box = MV.ui.el('div', { class: 'controls' });

    function onSlider(key) {
      return function (v, programmatic) {
        state[key] = v;
        if (!programmatic) setSweep(false);
        refresh();
      };
    }

    sliders.mu = MV.ui.slider(box, {
      label: 'μ — chemisches Potential',
      min: -2.5, max: 2.5, step: 0.01, value: state.mu, digits: 2,
      onInput: onSlider('mu')
    });
    sliders.t = MV.ui.slider(box, {
      label: 't — Hopping-Amplitude',
      min: 0.1, max: 2.0, step: 0.01, value: state.t, digits: 2,
      onInput: onSlider('t')
    });
    sliders.delta = MV.ui.slider(box, {
      label: 'Δ — Paarungsamplitude',
      min: 0.1, max: 1.5, step: 0.01, value: state.delta, digits: 2,
      onInput: onSlider('delta')
    });

    var row = MV.ui.el('div', { class: 'controls__wide' });

    function preset(mu, t, delta) {
      return function () {
        setSweep(false);
        if (t !== null) { state.t = t; sliders.t.set(t, true); }
        if (delta !== null) { state.delta = delta; sliders.delta.set(delta, true); }
        var m = typeof mu === 'function' ? mu() : mu;
        m = Math.min(2.5, Math.max(-2.5, m));
        state.mu = m;
        sliders.mu.set(m, true);
        refresh();
      };
    }

    MV.ui.buttonRow(row, {
      label: 'Voreinstellungen:',
      buttons: [
        { text: 'topologisch (μ = 0, t = Δ = 1)', onClick: preset(0, 1, 1) },
        { text: 'kritisch (μ = −t)', onClick: preset(function () { return -state.t; }, null, null) },
        { text: 'kritisch (μ = +t)', onClick: preset(function () { return state.t; }, null, null) },
        { text: 'trivial (μ = ' + MV.fmt(1.8, 1) + ' t)',
          onClick: preset(function () { return 1.8 * state.t; }, null, null) }
      ]
    });

    var sweepRow = MV.ui.buttonRow(row, {
      label: 'μ-Durchlauf:',
      buttons: [{
        text: 'μ-Durchlauf (quasistatisch)',
        pressed: false,
        onClick: function () { setSweep(!sweep.running); }
      }]
    });
    sweepButton = sweepRow.buttons[0];

    row.appendChild(MV.ui.el('p', {
      class: 'note',
      html: 'Der Durchlauf führt μ zwischen −2 und +2 hin und zurück und zeigt eine ' +
            'Folge statischer Gleichgewichtszustände (quasistatische Parametervariation), ' +
            'keine Zeitentwicklung. Jede Reglerbewegung hält ihn an.'
    }));

    box.appendChild(row);
    mount.appendChild(box);
  }

  function buildReadout(mount) {
    readout = MV.ui.readout(mount, {
      items: [
        { key: 'eps0',  label: 'ε₀ = −t − μ' },
        { key: 'epsPi', label: 'ε_π = t − μ' },
        { key: 's0',    label: 's₀ = sgn ε₀' },
        { key: 'sPi',   label: 's_π = sgn ε_π' },
        { key: 'nu',    label: 'ν = s₀ · s_π' },
        { key: 'gap',   label: 'Δ_bulk = min_k E(k)' },
        { key: 'W',     label: 'W (Windungszahl)' }
      ],
      note:
        'W erscheint nur ergänzend: es ist ausschließlich bei φ = 0 definiert ' +
        '(Klasse BDI). Maßgeblich für diese Arbeit ist ν (Klasse D), es gilt ' +
        'ν = (−1)^W. Zur Vorzeichenkonvention siehe „Details zur Abbildung" ' +
        'unter Panel A.'
    });
  }

  /* ====================================================================
     Public interface
     ==================================================================== */

  /* Shared with tab 3, which checks the continuum edge of the finite chain
     against exactly this quantity. One implementation, not two. */
  MV.kitaevBulk = { gap: bulkGap, energy: energy };

  MV.kitaev = {
    init: function (mount) {
      buildControls(mount);
      buildReadout(mount);

      var panels = MV.ui.el('div', { class: 'panels' });
      mount.appendChild(panels);

      panelA = MV.ui.panel(panels, {
        title: 'Panel A — Bahn von h(k) in der (h_y, h_z)-Ebene',
        quote: 'kitaev.A',
        lead: 'Ellipse mit Halbachsen Δ und t, Zentrum bei (0, −μ). ' +
              'Entscheidend ist, ob der Ursprung eingeschlossen wird.',
        aspect: 0.92,
        info: INFO_A
      });

      panelB = MV.ui.panel(panels, {
        title: 'Panel B — Bulk-Spektrum ±E(k)',
        quote: 'kitaev.B',
        lead: 'Anregungsspektrum nach Gl. ' + MV.eq('kitaev.spectrum') + '. Die Lücke ist Δ_bulk = min_k E(k).',
        aspect: 0.92,
        info: INFO_B
      });

      /* Schnittstelle für die Vorführungen an den Leitsätzen (demo.js). Sie
         reicht nur die vorhandenen Regler und den vorhandenen μ-Durchlauf
         heraus; eine Vorführung kann damit nichts tun, was nicht auch von
         Hand ginge. */
      MV.demo.provide('kitaev', MV.demo.handle({
        state: state,
        sliders: sliders,
        refresh: refresh,
        stop: function () { setSweep(false); },
        extra: {
          tol: TOL,
          /* ν, E(k) und die Lücke für ein gedachtes μ — ohne den Reiter zu
             verstellen, damit eine Endprüfung nachrechnen kann, was der
             Durchlauf gezeigt hat. */
          nuAt: function (mu) {
            return derive({ mu: mu, t: state.t, delta: state.delta }).nu;
          },
          energyAt: function (k, mu) {
            return energy(k, { mu: mu, t: state.t, delta: state.delta });
          },
          gapAt: function (mu) {
            return bulkGap({ mu: mu, t: state.t, delta: state.delta }).value;
          },
          read: function () { return derive(state); }
        }
      }));
    },

    draw: function () { refresh(); },

    /* Called when the tab is left: a hidden canvas cannot be measured, so the
       μ-Durchlauf must not keep running in the background. */
    stop: function () { setSweep(false); }
  };

}(MV));
