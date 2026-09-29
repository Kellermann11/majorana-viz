/* ------------------------------------------------------------------------
   demos.js — die fünfzehn Vorführungen, eine je Leitsatz.

   Hier steht, welche Regler eine Vorführung bewegt und woran sie sich am
   Ende messen lässt. Der Ablauf selbst — Knopf, Durchlauf, Anhalten,
   Ergebniszeile — steht in demo.js, die Physik in den Reiter-Modulen. Diese
   Datei rechnet nichts nach, was die Module nicht ohnehin rechnen; jede
   Endprüfung liest ihre Zahlen aus deren Schnittstelle.

   Eine Vorführung ist ein Durchlauf, keine Zeitentwicklung (CLAUDE.md §13).
   Das Tempo ist Anzeigetempo. Echte Zeit kommt in diesen Abläufen nirgends
   vor — auch nicht in #gatter Panel B, wo die vorhandene Schrittfolge läuft.

   Interface language is German, code and comments are English.
   ------------------------------------------------------------------------ */

var MV = window.MV || (window.MV = {});

(function (MV) {
  'use strict';

  var PI = Math.PI;

  /* Anzeigetempo, keine physikalische Zeit. */
  var LEAD = 500;    /* kurzer Vorlauf, damit der Ausgangszustand zu sehen ist */
  var RUN = 4800;    /* ein Durchlauf                                          */
  var HOLD = 1500;   /* ein Haltepunkt                                         */
  var STEP = 700;    /* ein Schritt in einer Folge gleichartiger Schritte      */

  function handle(tab) { return MV.demo.handleOf(tab); }

  /* ====================================================================
     #kitaev
     ==================================================================== */

  function registerKitaev() {
    var k = null;

    /* Wo ν springt — gesucht, nicht behauptet. Zwischen dem letzten μ mit
       ν = −1 und dem ersten mit ν = +1 liegt das Band, in dem der Reiter ν
       als undefiniert führt (|ε| < Toleranz); seine Mitte ist die
       Sprungstelle. */
    function scanJump(lo, hi, n) {
      var last = null, first = null, i, mu, v;
      for (i = 0; i <= n; i++) {
        mu = lo + (hi - lo) * i / n;
        v = k.nuAt(mu);
        if (v === -1) last = mu;
        if (v === 1 && first === null) first = mu;
      }
      return { last: last, first: first };
    }

    /* Erst grob, dann fein im gefundenen Streifen: ν auszuwerten kostet je
       Punkt eine Lücke und eine Windungszahl, und die Vorführung soll am
       Ende nicht merklich stehenbleiben. Die Auflösung ist dieselbe wie bei
       einem durchgehend feinen Raster. */
    function jumpOf(lo, hi) {
      var coarse = scanJump(lo, hi, 60);
      if (coarse.last === null || coarse.first === null) return null;
      var pad = (hi - lo) / 60;
      var fine = scanJump(Math.max(lo, coarse.last - pad),
                          Math.min(hi, coarse.first + pad), 200);
      if (fine.last === null || fine.first === null) return null;
      return (fine.last + fine.first) / 2;
    }

    MV.demo.register('kitaev.A', {
      prepare: function () {
        k = handle('kitaev');
        k.setAll({ t: 1, delta: 1, mu: 0 });
      },
      steps: [
        { wait: LEAD },
        { tween: { from: 0, to: 1.5, ms: RUN,
                   apply: function (v) { k.set('mu', v); } } },
        { wait: HOLD }
      ],
      check: function () {
        var t = k.get('t');
        var a = k.nuAt(0), b = k.nuAt(1.5), j = jumpOf(0, 1.5);
        var ok = a === -1 && b === 1 && j !== null && Math.abs(j - t) < 0.005;
        return {
          ok: ok,
          text: 'ν: ' + MV.fmtSigned(a, 0) + ' → ' + MV.fmtSigned(b, 0) +
                '; der Sprung liegt bei μ = ' + MV.fmt(j, 2) + ', also bei |μ| = t = ' +
                MV.fmt(t, 2) + '. Die Bahn schließt den Ursprung dann nicht mehr ein.'
        };
      }
    });

    /* Die Zielwerte ±t stehen fest, weil prepare t = 1 setzt. */
    MV.demo.register('kitaev.B', {
      prepare: function () {
        k = handle('kitaev');
        k.setAll({ t: 1, delta: 1, mu: 0 });
      },
      steps: [
        { wait: LEAD },
        { tween: { from: 0, to: -1, ms: 2400,
                   apply: function (v) { k.set('mu', v); } } },
        { wait: HOLD },
        { tween: { from: -1, to: 1, ms: 2800,
                   apply: function (v) { k.set('mu', v); } } },
        { wait: HOLD }
      ],
      check: function () {
        var t = k.get('t'), tol = k.tol;
        var e0 = k.energyAt(0, -t);
        var ePi = k.energyAt(PI, t);
        var g1 = k.gapAt(-t), g2 = k.gapAt(t);
        var ok = e0 < tol && ePi < tol && g1 < tol && g2 < tol;
        return {
          ok: ok,
          text: 'Lückenschluss: bei μ = −t ist E(k = 0) = ' + MV.fmt(e0, 3) +
                ', bei μ = +t ist E(k = ±π) = ' + MV.fmt(ePi, 3) +
                ' — beide unter der Toleranz ' + MV.fmt(tol, 2) +
                ', mit der dieser Reiter einen Lückenschluss zählt.'
        };
      }
    });
  }

  /* ====================================================================
     #nanodraht
     ==================================================================== */

  function registerNanodraht() {
    var n = null;

    MV.demo.register('nanodraht.A', {
      prepare: function () {
        n = handle('nanodraht');
        /* Stufe 2 — plus Zeeman: dort sind die nackten Bänder mit der
           Zeeman-Lücke zu sehen, um die es im Leitsatz geht. */
        n.setAll({ alpha: 1.4, delta: 0.6, ez: 1.0, k: 0.8, mu: 1.5, stage: 1 });
      },
      steps: [
        { wait: LEAD },
        { tween: { from: 1.5, to: 0, ms: RUN,
                   apply: function (v) { n.set('mu', v); } } },
        { wait: HOLD }
      ],
      check: function () {
        var before = n.fermiCountAt({ mu: 1.5 });
        var after = n.fermiCountAt({ mu: 0 });
        var now = n.read().fermi.length;
        var ez = n.get('ez');
        return {
          ok: before === 4 && after === 2 && now === 2,
          text: 'Fermi-Punkte: ' + before + ' → ' + after +
                '. Bei μ = 0 liegt μ im Zeeman-Fenster ±' + MV.fmt(ez, 2) +
                ', und nur das untere Band ist besetzt.'
        };
      }
    });

    MV.demo.register('nanodraht.B', {
      prepare: function () {
        n = handle('nanodraht');
        n.setAll({ mu: -0.3, alpha: 1.4, delta: 0.6, k: 0.8, ez: 0.5 });
      },
      steps: [
        { wait: LEAD },
        { tween: { from: 0.5, to: 3, ms: RUN,
                   apply: function (v) { n.set('ez', v); } } },
        { wait: HOLD }
      ],
      check: function () {
        var kk = n.get('k'), N = 60, i, ez, v, prev = Infinity, mono = true;
        for (i = 0; i <= N; i++) {
          ez = 0.5 + (3 - 0.5) * i / N;
          v = n.deltaEffAt(kk, { ez: ez });
          if (v >= prev) mono = false;
          prev = v;
        }
        var a = n.deltaEffAt(kk, { ez: 0.5 }), b = n.deltaEffAt(kk, { ez: 3 });
        return {
          ok: mono && b < a,
          text: 'Δ_eff(k = ' + MV.fmt(kk, 2) + ') fällt von ' + MV.fmt(a, 3) +
                ' auf ' + MV.fmt(b, 3) + ', während E_Z von ' + MV.fmt(0.5, 2) +
                ' auf ' + MV.fmt(3, 2) + ' steigt — monoton über ' + (N + 1) +
                ' Stützstellen.'
        };
      }
    });

    MV.demo.register('nanodraht.C', {
      prepare: function () {
        n = handle('nanodraht');
        n.setAll({ alpha: 1.4, delta: 0.6, ez: 1.0, k: 0.8, mu: 0 });
      },
      steps: [
        { wait: LEAD },
        { tween: { from: 0, to: 1.5, ms: RUN,
                   apply: function (v) { n.set('mu', v); } } },
        { wait: HOLD }
      ],
      check: function () {
        var ez = n.get('ez'), delta = n.get('delta');
        /* Das Kriterium E_Z > √(Δ² + μ²) nach μ aufgelöst; das ist dieselbe
           Ungleichung, nur anders geschrieben. */
        var boundary = ez > delta ? Math.sqrt(ez * ez - delta * delta) : null;
        var start = n.topologicalAt({ mu: 0 });
        var end = n.topologicalAt({ mu: 1.5 });
        var inside = boundary !== null && n.topologicalAt({ mu: boundary - 0.1 });
        var outside = boundary !== null && !n.topologicalAt({ mu: boundary + 0.1 });
        return {
          ok: start === true && end === false && inside && outside,
          text: 'Die Phase wechselt von topologisch zu trivial. Die Grenze liegt bei ' +
                '|μ| = √(E_Z² − Δ²) = ' + MV.fmt(boundary, 2) +
                '; dort ist E_Z = √(Δ² + μ²).'
        };
      }
    });
  }

  /* ====================================================================
     #energieskalen
     ==================================================================== */

  function registerEnergieskalen() {
    var e = null;
    var closedAt = null;     /* Zustand nach dem ersten Durchlauf */
    var k4Before = null;
    var teilA = null;        /* Zwischenstand der Vorführung zu 5.6 */

    /* Über Dekaden wird im Logarithmus gefahren, sonst stünde der Wert die
       halbe Zeit im obersten Zehntel und risse die übrigen sieben Dekaden am
       Ende durch. Der logarithmische Regler bewegt sich damit gleichmäßig. */
    function decades(fromV, toV, ms) {
      return { tween: {
        from: Math.log(fromV) / Math.LN10,
        to: Math.log(toV) / Math.LN10,
        ms: ms,
        apply: function (v) { e.set('TP', Math.pow(10, v)); }
      } };
    }

    MV.demo.register('energieskalen.A', {
      prepare: function () {
        e = handle('energieskalen');
        closedAt = null;
        e.defaults();
      },
      steps: [
        { wait: LEAD },
        { tween: { from: 0.1, to: 0.45, ms: 2600,
                   apply: function (v) { e.set('xi', v); } } },
        { set: function () {
            var d = e.read();
            closedAt = { hasWindow: d.hasWindow, xiClose: d.xiClose, xi: d.s.xi };
          } },
        { wait: HOLD },
        { tween: { from: 0.45, to: 0.1, ms: 1200,
                   apply: function (v) { e.set('xi', v); } } },
        { tween: { from: 100, to: 30, ms: 2600,
                   apply: function (v) { e.set('ecRef', v); } } },
        { wait: HOLD }
      ],
      check: function () {
        var d = e.read();
        var ok = closedAt !== null && closedAt.hasWindow === false &&
                 closedAt.xiClose !== null && closedAt.xi > closedAt.xiClose &&
                 d.hasWindow === true;
        return {
          ok: ok,
          text: 'Bei ξ = ' + MV.fmt(closedAt ? closedAt.xi : 0, 2) +
                ' µm gibt es kein Fenster mehr — es schließt sich bei ξ_zu = ' +
                MV.fmt(closedAt ? closedAt.xiClose : 0, 2) +
                ' µm. Mit E_C^ref = ' + MV.fmt(d.s.ecRef, 0) +
                ' µeV liegt die obere Grenze bei ' + MV.fmt(d.hi, 2) +
                ' µm, das Fenster bei ' + MV.fmtRange(d.lo, d.hi, 2, 'µm') + '.'
        };
      }
    });

    /* §5.6 stellt zwei Sätze nebeneinander: eine große Ladungsenergie
       garantiert keine lange Lebensdauer — und umgekehrt kann die Lebensdauer
       lang sein, ohne dass die Ladungsenergie überhaupt schützt. Der Leitsatz
       ist der zweite Satz; zu sehen ist die Aussage erst, wenn beide Fälle
       nacheinander vorkommen. Teil a zeigt den ersten, Teil b den zweiten.

       Dass T_P dabei bis in den Mikrosekundenbereich läuft, ist eine Eingabe
       und keine Behauptung über einen Mechanismus: welcher Vorgang eine kurze
       Lebensdauer verursacht, sagt dieser Reiter nicht. */
    MV.demo.register('energieskalen.B', {
      prepare: function () {
        e = handle('energieskalen');
        e.defaults();
        e.markCriteria(['K3', 'K4']);
        k4Before = e.ratiosNow().K4;
        teilA = null;
      },
      steps: [
        { wait: LEAD },
        /* ---- a) schwebend: Ladungsschutz da, Parität kurzlebig ---- */
        decades(1, 1e-7, 3000),
        { set: function () {
            teilA = {
              floating: e.get('floating'),
              K3: e.criteriaVerdict('K3'),
              K4: e.criteriaVerdict('K4'),
              TP: e.get('TP'),
              ratio: e.ratiosNow().K4
            };
          } },
        { wait: HOLD },
        /* ---- b) geerdet: kein Ladungsschutz, Parität langlebig ---- */
        decades(1e-7, 1, 800),
        { set: function () { e.holdK3(true); e.showPanel('A'); } },
        { wait: 500 },
        { set: function () { e.set('floating', false); } },
        { wait: 1200 },
        { set: function () { e.showPanel('B'); } },
        { wait: 1000 },
        { set: function () { e.holdK3(false); } },
        { wait: 600 },
        decades(1, 10, 2000),
        { wait: HOLD }
      ],
      check: function () {
        var r = e.ratiosNow();
        var shown = e.criteriaShown();
        var marked = e.markedCriteria();
        var tp = e.get('TP'), tau = e.get('tau');
        /* K4 ist τ_meas/T_P, τ_meas in µs und T_P in s — dieselbe Rechnung
           wie im Reiter, hier nur zur Gegenprobe. */
        var expected = (tau * 1e-6) / tp;

        var aOk = teilA !== null &&
                  teilA.floating === true &&
                  teilA.K3 === 'erfüllt' &&
                  teilA.K4 === 'nicht erfüllt';
        var bOk = e.get('floating') === false &&
                  shown.indexOf('K3') < 0 && r.K3 === null &&
                  shown.indexOf('K4') >= 0 &&
                  e.criteriaVerdict('K4') === 'erfüllt' &&
                  marked.indexOf('K4') >= 0 &&
                  Math.abs(tp - 10) < 1e-9 &&
                  Math.abs(r.K4 - expected) < 1e-18 &&
                  k4Before !== null && r.K4 < k4Before;

        return {
          ok: aOk && bOk,
          text: 'Schwebend mit T_P = ' + (teilA ? e.fmtTP(teilA.TP) : '—') +
                ': Ladungsschutz vorhanden, K4 verletzt (τ_meas/T_P = ' +
                (teilA ? MV.fmt(teilA.ratio, 1) : '—') + '). Geerdet mit T_P = ' +
                e.fmtTP(tp) + ': kein Ladungsschutz, K4 erfüllt (' +
                MV.fmtExp(r.K4, 2) + '). K3 und K4 sind unabhängig.'
        };
      }
    });
  }

  /* ====================================================================
     #auslese
     ==================================================================== */

  function registerAuslese() {
    var a = null;

    /* Die sechs Paare der Reihe nach — derselbe Weg, den ein Klick auf einen
       Chip oder auf zwei Majoranas im Bild nimmt. */
    function pairSteps() {
      var out = [{ wait: LEAD }], i;
      for (i = 0; i < 6; i++) {
        out.push({ set: (function (idx) {
          return function () {
            var p = a.allPairs()[idx];
            a.selectPair(p[0], p[1]);
          };
        }(i)) });
        out.push({ wait: STEP });
      }
      return out;
    }

    MV.demo.register('auslese.A', {
      prepare: function () {
        a = handle('auslese');
        a.forgetVisited();
        a.refresh();
      },
      steps: pairSteps(),
      check: function () {
        var c = a.counts();
        return {
          ok: c.pairs === 6 && c.axes === 3,
          text: c.pairs + ' von 6 Paaren ausprobiert, ' + c.axes +
                ' von 3 Achsen gefunden — vier Majoranas lassen sich auf ' +
                'C(4,2) = 6 Weisen zu Paaren zusammenfassen, und je zwei dieser ' +
                'Paare beschreiben dieselbe Achse.'
        };
      }
    });

    /* Vorgaben von Plugge, Anh. B.3 — dieselben wie der Knopf „Vorgaben". */
    function pluggeDefaults() {
      a = handle('auslese');
      a.setAll({ t0: 5, t1: 1, phi: 0, eps: 0, tau: PI / 4, gamma: 0 });
    }

    MV.demo.register('auslese.B', {
      prepare: pluggeDefaults,
      steps: [
        { wait: LEAD },
        { tween: { from: 5, to: 0, ms: RUN,
                   apply: function (v) { a.set('t0', v); } } },
        { wait: HOLD }
      ],
      check: function () {
        var r = a.read();
        return {
          ok: r.blind === true && Math.abs(r.tplus - r.tminus) < 1e-12,
          text: '|t₊| = ' + MV.fmt(r.tplus, 3) + ' und |t₋| = ' + MV.fmt(r.tminus, 3) +
                ', die Differenz |t₊|² − |t₋|² ist ' + MV.fmtSigned(r.split, 3) +
                ' — nicht unterscheidbar (' + r.why + ').'
        };
      }
    });

    MV.demo.register('auslese.C', {
      prepare: pluggeDefaults,
      steps: [
        { set: function () { a.showPanelC(); } },
        { wait: LEAD },
        { tween: { from: 5, to: 0, ms: RUN,
                   apply: function (v) { a.set('t0', v); } } },
        { wait: HOLD }
      ],
      check: function () {
        var r = a.read();
        var ok = r.deltaGmax < 1e-12 && r.contrastC2 < 1e-12 && r.peakGapC3 < 1e-12;
        return {
          ok: ok,
          text: 'Alle drei Verfahren werden blind: δG_max = ' + MV.fmt(r.deltaGmax, 3) +
                ', Kontrast im Zeitbereich = ' + MV.fmt(r.contrastC2, 3) +
                ', Abstand der beiden Resonanzen = ' + MV.fmt(r.peakGapC3, 3) +
                '. Sie lesen dieselbe Größe, und die trägt hier nichts mehr.'
        };
      }
    });
  }

  /* ====================================================================
     #gatter
     ==================================================================== */

  function registerGatter() {
    var g = null;
    var afterFirst = null;

    MV.demo.register('gatter.A', {
      prepare: function () {
        g = handle('gatter');
        g.setStart(90, 0);          /* x̂ = +1 */
        g.clearSeq();
        afterFirst = null;
      },
      steps: [
        { wait: LEAD },
        { set: function () { g.applyGate('T'); } },
        { set: function () { afterFirst = g.vertexInfo(); } },
        { wait: HOLD },
        { set: function () { g.applyGate('T'); } },
        { wait: HOLD }
      ],
      check: function () {
        var now = g.vertexInfo();
        var ok = afterFirst !== null && afterFirst.onVertex === false &&
                 now.onVertex === true && now.label === 'ŷ = +1';
        return {
          ok: ok,
          text: 'Ecke des Oktaeders: nach dem ersten T-Gatter „nein" (Abstand ' +
                MV.fmt(afterFirst ? afterFirst.d : 0, 3) +
                '), nach dem zweiten „ja" — ' + now.label + ' ' + now.ket +
                '. Eine einzelne Drehung um π/4 führt zwischen die Ecken.'
        };
      }
    });

    MV.demo.register('gatter.B', {
      prepare: function () {
        g = handle('gatter');
        g.set('proto', 0);                        /* Ŝ_z                    */
        g.setStart(90, 0);                        /* x̂ = +1                 */
        g.set('correct', true);
        g.set('outcomes', [null, null, null]);    /* gezogen, nicht gesetzt */
      },
      steps: [
        { wait: LEAD },
        /* Die Wiedergabe von Panel B, so wie der Knopf „Start" sie fährt.
           Die Vorführung wartet nur, bis sie von selbst am Ende steht. */
        { run: function (next, instant, ctx) {
            if (instant) { g.stepToEnd(); next(); return; }
            MV.demo.onCancel(function () { g.stop(); });
            g.startPlayer();
            (function poll() {
              if (ctx.cancelled) return;
              if (!g.isPlaying()) { next(); return; }
              window.requestAnimationFrame(poll);
            }());
          } },
        { wait: HOLD }
      ],
      check: function () {
        var run = g.runNow();
        var ok = g.atEnd() && run.impossible < 0 && run.dev < 1e-12;
        return {
          ok: ok,
          text: 'Die Schrittfolge endet mit einer Abweichung vom Ziel von ' +
                MV.fmtExp(run.dev, 1) + ', also unter ' + MV.fmtExp(1e-12, 0) +
                ' — ' + g.protoName() + ' ist aus zwei Messungen und einer ' +
                'klassisch bedingten Korrektur erzeugt.'
        };
      }
    });
  }

  /* ====================================================================
     #pb-tetron
     ==================================================================== */

  function registerPbTetron() {
    var p = null;

    MV.demo.register('pbtetron.A', {
      prepare: function () {
        p = handle('pbtetron');
        p.set('zoom', 0);
      },
      steps: [
        { wait: 1200 },
        { set: function () { p.set('zoom', 1); } },
        { wait: HOLD },
        { set: function () { p.set('zoom', 2); } },
        { wait: HOLD }
      ],
      check: function () {
        var line = p.readoutStatus();
        var wanted = 'oberer Draht von Tetron BA — geerdet';
        return {
          ok: p.get('zoom') === 2 && line.indexOf(wanted) >= 0,
          text: 'Zoomstufe ' + (p.get('zoom') + 1) + ' von ' + p.zoomCount() + ' — ' +
                p.zoomName() + '. Darunter steht: „' + line + '"'
        };
      }
    });

    MV.demo.register('pbtetron.B', {
      prepare: function () {
        p = handle('pbtetron');
        p.markTestRows(0);
      },
      steps: [
        { wait: 800 },
        /* Der Durchlauf läuft über den Anteil der Zeilen, nicht über eine
           feste Zahl: wie viele es sind, sagt die Tabelle selbst. */
        { tween: { from: 0, to: 1, ms: 2400, apply: function (v) {
            p.markTestRows(Math.round(v * p.testRowCount()));
          } } },
        { wait: HOLD }
      ],
      check: function () {
        var total = p.testRowCount(), marked = p.markedCount();
        var first = p.testRowText(0);
        return {
          ok: total > 0 && marked === total && first.length > 0,
          text: marked + ' von ' + total + ' Zeilen der 3 µm-Teststruktur ' +
                'hervorgehoben — sie gehören zum separaten Draht, nicht zum ' +
                'vermessenen. Nur ein Rahmen: kein Status, keine Wertung.'
        };
      }
    });

    MV.demo.register('pbtetron.C', {
      prepare: function () {
        p = handle('pbtetron');
        p.set('material', 0);
      },
      steps: [
        { wait: 1200 },
        { set: function () { p.set('material', 1); } },
        { wait: HOLD },
        { set: function () { p.set('material', 2); } },
        { wait: HOLD }
      ],
      check: function () {
        var f = p.factors();
        var cols = p.materialColumns();
        var ok = p.get('material') === 2 && cols === 2 &&
                 isFinite(f.deltaT) && f.tauLo > f.deltaT;
        return {
          ok: ok,
          text: 'Beide Plattformen nebeneinander: Δ_T wächst um ' + f.deltaTText +
                ', τ_Z um ' + MV.fmtRange(f.tauLo, f.tauHi, 1, null, MV.fmtExp) +
                '. Bewertet wird das in ' + MV.secOfWork('6.5') + '.'
        };
      }
    });
  }

  /* ====================================================================
     Einbau — von main.js gerufen, nachdem alle Reiter aufgebaut sind
     ==================================================================== */

  MV.demos = {
    install: function () {
      registerKitaev();
      registerNanodraht();
      registerEnergieskalen();
      registerAuslese();
      registerGatter();
      registerPbTetron();
      return MV.demo.keys().length;
    }
  };

}(MV));
