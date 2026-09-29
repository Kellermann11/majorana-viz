/* ------------------------------------------------------------------------
   auslese.js — Reiter #auslese: Pauli-Achsen und die drei Ausleseverfahren.
   Abschnitte 5.2, 5.3 und 5.5 der Arbeit.

   Panel A  Was gemessen wird: eine Pauli-Achse als Majorana-Paar, in der
            Drei-Punkte-Geometrie aus Plugge et al., Fig. 3(a).
   Panel B  Die Primitive t_z = t₀ + t₁z als Zeigerdiagramm.
   Panel C  Drei Ausleseverfahren, die alle dieselbe Größe |t_z| lesen.

   Die Formeln der Verfahren 2 und 3 stammen aus Plugge et al., Gl. (2)–(6)
   und Anhang B; die Arbeit beschreibt sie in §5.5 nur qualitativ. Verfahren 1
   steht als (5.10) auch in der Arbeit.

   Alles dimensionslos, willkürliche Einheiten wie bei Plugge (ℏ = 1);
   G in Einheiten von (e²/h)·ν₁ν₂. Kein Flussquant — φ ist eine reine Phase.
   ------------------------------------------------------------------------ */

var MV = window.MV || (window.MV = {});

(function (MV) {
  'use strict';

  var PI = Math.PI;

  /* ====================================================================
     Complex 4×4 matrices

     Small and explicit: the point of panel C is that the algebra is checked
     numerically rather than asserted, so the arithmetic is written out here
     instead of being taken on trust from anywhere else.
     Layout: re and im are Float64Array(16), row major.
     ==================================================================== */

  function cmat(re, im) {
    return { re: Float64Array.from(re), im: Float64Array.from(im || new Array(16).fill(0)) };
  }

  function czeros() { return { re: new Float64Array(16), im: new Float64Array(16) }; }

  function cid() {
    var m = czeros();
    for (var i = 0; i < 4; i++) m.re[i * 4 + i] = 1;
    return m;
  }

  function cmul(a, b) {
    var out = czeros(), i, j, k, ar, ai, br, bi;
    for (i = 0; i < 4; i++) {
      for (j = 0; j < 4; j++) {
        var sr = 0, si = 0;
        for (k = 0; k < 4; k++) {
          ar = a.re[i * 4 + k]; ai = a.im[i * 4 + k];
          br = b.re[k * 4 + j]; bi = b.im[k * 4 + j];
          sr += ar * br - ai * bi;
          si += ar * bi + ai * br;
        }
        out.re[i * 4 + j] = sr;
        out.im[i * 4 + j] = si;
      }
    }
    return out;
  }

  function cadd(a, b, sb) {
    sb = sb === undefined ? 1 : sb;
    var out = czeros();
    for (var i = 0; i < 16; i++) {
      out.re[i] = a.re[i] + sb * b.re[i];
      out.im[i] = a.im[i] + sb * b.im[i];
    }
    return out;
  }

  /** Multiply by a complex scalar (sr + i·si). */
  function cscale(a, sr, si) {
    si = si || 0;
    var out = czeros();
    for (var i = 0; i < 16; i++) {
      out.re[i] = a.re[i] * sr - a.im[i] * si;
      out.im[i] = a.re[i] * si + a.im[i] * sr;
    }
    return out;
  }

  /** Largest absolute entry of a − b. */
  function cdiff(a, b) {
    var worst = 0, i, d;
    for (i = 0; i < 16; i++) {
      d = Math.sqrt((a.re[i] - b.re[i]) * (a.re[i] - b.re[i]) +
                    (a.im[i] - b.im[i]) * (a.im[i] - b.im[i]));
      if (d > worst) worst = d;
    }
    return worst;
  }

  function cmaxAbs(a) {
    var worst = 0, i, d;
    for (i = 0; i < 16; i++) {
      d = Math.sqrt(a.re[i] * a.re[i] + a.im[i] * a.im[i]);
      if (d > worst) worst = d;
    }
    return worst;
  }

  /** Kronecker product of two 2×2 matrices, given as [re2x2, im2x2]. */
  function kron(a, b) {
    var out = czeros(), i, j, k, l;
    for (i = 0; i < 2; i++) {
      for (j = 0; j < 2; j++) {
        for (k = 0; k < 2; k++) {
          for (l = 0; l < 2; l++) {
            var ar = a.re[i * 2 + j], ai = a.im[i * 2 + j];
            var br = b.re[k * 2 + l], bi = b.im[k * 2 + l];
            var row = i * 2 + k, col = j * 2 + l;
            out.re[row * 4 + col] = ar * br - ai * bi;
            out.im[row * 4 + col] = ar * bi + ai * br;
          }
        }
      }
    }
    return out;
  }

  function m2(re, im) {
    return { re: Float64Array.from(re), im: Float64Array.from(im || [0, 0, 0, 0]) };
  }

  /* Pauli matrices and the 2×2 identity. */
  var S0 = m2([1, 0, 0, 1]);
  var SX = m2([0, 1, 1, 0]);
  var SY = m2([0, 0, 0, 0], [0, -1, 1, 0]);
  var SZ = m2([1, 0, 0, -1]);

  /**
   * Representation of the four Majorana operators.
   *
   *   γ₁ = σ_x ⊗ 1     γ₂ = σ_y ⊗ 1     γ₃ = σ_z ⊗ σ_x     γ₄ = σ_z ⊗ σ_y
   *
   * All four are hermitian and square to the identity; the choice is the
   * standard Jordan–Wigner pair of qubits and is not unique — any other
   * representation of the same Clifford algebra would do. See NOTES.md.
   */
  var GAMMA = [
    kron(SX, S0),
    kron(SY, S0),
    kron(SZ, SX),
    kron(SZ, SY)
  ];

  /* ====================================================================
     The algebra check   (5.5) – (5.8)
     ==================================================================== */

  var TOL = 1e-12;

  /** Bilinear i γ_a γ_b, with 1-based indices as in the thesis. */
  function bilinear(a, b) {
    return cscale(cmul(GAMMA[a - 1], GAMMA[b - 1]), 0, 1);
  }

  /* Axis definitions. `pair` is the form of (5.6), `alt` the second form that
     (5.8) adds; both are equal on the sector selected by (5.5). */
  var AXES = [
    {
      key: 'x',
      symbol: 'x̂',
      pair: [1, 2],
      alt: [3, 4],
      geometry: 'längs — Parität des oberen Drahtes und zugleich die des unteren',
      loopA: 'oberer Draht',
      loopB: 'unterer Draht',
      /* The thesis states this correspondence explicitly in section 5.3;
         für ŷ und ẑ steht dort keine. Die Begründung dafür, was sich daraus
         trotzdem folgern lässt, steht in CONVENTION_NOTE. */
      microsoft: 'Z'
    },
    {
      key: 'y',
      symbol: 'ŷ',
      pair: [3, 1],
      alt: [2, 4],
      geometry: 'diagonal — eine der beiden Diagonalen',
      loopA: 'Diagonale γ₃–γ₁',
      loopB: 'Diagonale γ₂–γ₄',
      microsoft: null
    },
    {
      key: 'z',
      symbol: 'ẑ',
      pair: [2, 3],
      alt: [1, 4],
      geometry: 'quer — Parität des rechten Endenpaares und zugleich die des linken',
      loopA: 'rechtes Endenpaar',
      loopB: 'linkes Endenpaar',
      microsoft: null
    }
  ];

  /* ====================================================================
     Der Sektor P̂ = +1 als Qubit

     (5.5) wählt den Sektor P̂ = −γ₁γ₂γ₃γ₄ = +1 aus; er ist zweidimensional
     und trägt das Qubit. Welche zwei Zustände darin |0⟩ und |1⟩ heißen, legt
     die Arbeit nicht fest — das ist eine Basiswahl, keine Physik.

     Gewählt ist hier die Basis, in der (5.6) wörtlich gilt:

       |0⟩ = ( e₀ − e₃ ) / √2        |1⟩ = − ( e₀ + e₃ ) / √2

     über den Basisvektoren e₀ … e₃ der Darstellung aus GAMMA. In ihr ergeben
     alle sechs Bilinearen exakt die Pauli-Matrizen; nimmt man stattdessen
     einfach die beiden Sektor-Basisvektoren e₀ und e₃ in dieser Reihenfolge,
     kommen −σ_z, −σ_y und −σ_x heraus — dieselbe Algebra, aber eine andere
     Benennung der Achsen. Siehe NOTES.md, „Verifikation der Pauli-Algebra".
     ==================================================================== */

  var QUBIT_BASIS = (function () {
    var s = 1 / Math.SQRT2;
    return [
      { re: [s, 0, 0, -s], im: [0, 0, 0, 0] },
      { re: [-s, 0, 0, -s], im: [0, 0, 0, 0] }
    ];
  }());

  /** ⟨a| M |b⟩ in der Qubit-Basis; liefert [Re, Im]. */
  function sectorElement(M, a, b) {
    var re = 0, im = 0, r, c, ar, ai, bc, bi, mr, mi, xr, xi;
    for (r = 0; r < 4; r++) {
      for (c = 0; c < 4; c++) {
        ar = QUBIT_BASIS[a].re[r]; ai = -QUBIT_BASIS[a].im[r];   /* konjugiert */
        bc = QUBIT_BASIS[b].re[c]; bi = QUBIT_BASIS[b].im[c];
        mr = M.re[r * 4 + c]; mi = M.im[r * 4 + c];
        xr = ar * mr - ai * mi; xi = ar * mi + ai * mr;
        re += xr * bc - xi * bi; im += xr * bi + xi * bc;
      }
    }
    return [re, im];
  }

  /** Das Bilineare i γ_i γ_j als 2×2-Matrix im Sektor, Zeilen zuerst. */
  function sectorMatrix(i, j) {
    var M = bilinear(i, j), out = [], a, b;
    for (a = 0; a < 2; a++) for (b = 0; b < 2; b++) out.push(sectorElement(M, a, b));
    return out;
  }

  /** Das Quadrat eines Bilinearen im Sektor. */
  function sectorSquare(i, j) {
    var M = bilinear(i, j), S = cmul(M, M), out = [], a, b;
    for (a = 0; a < 2; a++) for (b = 0; b < 2; b++) out.push(sectorElement(S, a, b));
    return out;
  }

  var PAULI = {
    x: [[0, 0], [1, 0], [1, 0], [0, 0]],
    y: [[0, 0], [0, -1], [0, 1], [0, 0]],
    z: [[1, 0], [0, 0], [0, 0], [-1, 0]],
    id: [[1, 0], [0, 0], [0, 0], [1, 0]]
  };

  function mat2Diff(m, p) {
    var w = 0, k;
    for (k = 0; k < 4; k++) w = Math.max(w, Math.hypot(m[k][0] - p[k][0], m[k][1] - p[k][1]));
    return w;
  }

  function mat2Mul(a, b) {
    var out = [], r, c, k, re, im, x, y;
    for (r = 0; r < 2; r++) {
      for (c = 0; c < 2; c++) {
        re = 0; im = 0;
        for (k = 0; k < 2; k++) {
          x = a[r * 2 + k]; y = b[k * 2 + c];
          re += x[0] * y[0] - x[1] * y[1];
          im += x[0] * y[1] + x[1] * y[0];
        }
        out.push([re, im]);
      }
    }
    return out;
  }

  function mat2TimesI(m) { return m.map(function (c) { return [-c[1], c[0]]; }); }

  /** Ein Matrixeintrag als kurze Zeichenkette: 0, 1, −1, i, −i, … */
  function cellText(c) {
    var r = Math.abs(c[0]) < 1e-12 ? 0 : c[0];
    var i = Math.abs(c[1]) < 1e-12 ? 0 : c[1];
    function n(v) { return Math.abs(v - Math.round(v)) < 1e-9 ? String(Math.round(v)) : MV.fmt(v, 3); }
    if (i === 0) return n(r).replace('-', '−');
    if (r === 0) return (i === 1 ? 'i' : i === -1 ? '−i' : n(i).replace('-', '−') + 'i');
    return (n(r) + (i > 0 ? '+' : '−') + (Math.abs(i) === 1 ? '' : n(Math.abs(i))) + 'i').replace('-', '−');
  }

  function mat2Html(m, cls) {
    return '<span class="mat2 ' + (cls || '') + '">' +
      '<span class="mat2__grid">' +
      m.map(function (c) { return '<span>' + cellText(c) + '</span>'; }).join('') +
      '</span></span>';
  }

  /* ====================================================================
     Die sechs Paare
     ==================================================================== */

  /** Abweichung als kurze Zahl; exakt null bleibt "0". */
  function fmtDev(v) {
    if (v === 0) return '0';
    return MV.fmtExp(v, 1);
  }

  function pairKey(i, j) { return Math.min(i, j) + '-' + Math.max(i, j); }

  /** Zu welcher Achse gehört ein Paar, und welche der beiden Schleifen ist es? */
  /**
   * z ist der Eigenwert der in Panel A gewählten Achse, nicht irgendeiner.
   * Der Name kommt aus der Auswahl, damit Panel B und Panel A nicht
   * auseinanderlaufen.
   */
  function zName() {
    return 'Eigenwert von ' + AXES[state.axis].symbol;
  }

  function axisOfPair(i, j) {
    var k = pairKey(i, j), a;
    for (a = 0; a < AXES.length; a++) {
      if (pairKey(AXES[a].pair[0], AXES[a].pair[1]) === k) return { idx: a, which: 'pair' };
      if (pairKey(AXES[a].alt[0], AXES[a].alt[1]) === k) return { idx: a, which: 'alt' };
    }
    return null;
  }

  /** Alle sechs Paare, nach Achse gruppiert — die Reihenfolge der Anzeige. */
  function allPairs() {
    var out = [];
    AXES.forEach(function (ax, idx) {
      out.push({ pair: ax.pair, axis: idx, which: 'pair' });
      out.push({ pair: ax.alt, axis: idx, which: 'alt' });
    });
    return out;
  }

  function runAlgebraCheck() {
    var out = {
      anticommutators: [],
      worstAnticommutator: 0,
      squares: [],
      cyclic: [],
      decompositions: [],
      sectorMatrices: [],
      parityEigenvalues: [],
      worst: 0,
      passed: true
    };
    var i, j;

    /* 1 — {γ_i, γ_j} = 2 δ_ij */
    for (i = 1; i <= 4; i++) {
      var row = [];
      for (j = 1; j <= 4; j++) {
        var anti = cadd(cmul(GAMMA[i - 1], GAMMA[j - 1]), cmul(GAMMA[j - 1], GAMMA[i - 1]));
        var want = i === j ? cscale(cid(), 2) : czeros();
        var dev = cdiff(anti, want);
        row.push(dev);
        if (dev > out.worstAnticommutator) out.worstAnticommutator = dev;
      }
      out.anticommutators.push(row);
    }
    out.worst = Math.max(out.worst, out.worstAnticommutator);

    /* 2 — parity operator (5.5) and the projector onto P = +1 */
    var P = cscale(cmul(cmul(GAMMA[0], GAMMA[1]), cmul(GAMMA[2], GAMMA[3])), -1);
    out.paritySquareDeviation = cdiff(cmul(P, P), cid());
    out.worst = Math.max(out.worst, out.paritySquareDeviation);

    var proj = cscale(cadd(cid(), P), 0.5);
    out.projectorIdempotency = cdiff(cmul(proj, proj), proj);
    out.worst = Math.max(out.worst, out.projectorIdempotency);

    /* trace of the projector = dimension of the sector */
    var tr = 0;
    for (i = 0; i < 4; i++) tr += proj.re[i * 4 + i];
    out.sectorDimension = tr;

    /* diagonal of P in this representation — it is already diagonal here */
    for (i = 0; i < 4; i++) out.parityEigenvalues.push(P.re[i * 4 + i]);

    function onSector(m) { return cmul(cmul(proj, m), proj); }

    /* 3 — squares, cyclic relation (5.7), and the two forms of (5.8) */
    AXES.forEach(function (ax) {
      var A = bilinear(ax.pair[0], ax.pair[1]);
      var B = bilinear(ax.alt[0], ax.alt[1]);
      var As = onSector(A);

      var sq = cdiff(onSector(cmul(A, A)), proj);
      out.squares.push({ symbol: ax.symbol, deviation: sq });
      out.worst = Math.max(out.worst, sq);

      var dec = cdiff(As, onSector(B));
      out.decompositions.push({
        symbol: ax.symbol,
        first: 'i γ' + sub(ax.pair[0]) + 'γ' + sub(ax.pair[1]),
        second: 'i γ' + sub(ax.alt[0]) + 'γ' + sub(ax.alt[1]),
        deviation: dec
      });
      out.worst = Math.max(out.worst, dec);
    });

    /* x̂ŷ = iẑ and its cyclic partners */
    [[0, 1, 2], [1, 2, 0], [2, 0, 1]].forEach(function (t) {
      var A = bilinear(AXES[t[0]].pair[0], AXES[t[0]].pair[1]);
      var B = bilinear(AXES[t[1]].pair[0], AXES[t[1]].pair[1]);
      var C = bilinear(AXES[t[2]].pair[0], AXES[t[2]].pair[1]);
      var lhs = onSector(cmul(A, B));
      var rhs = onSector(cscale(C, 0, 1));
      var dev = cdiff(lhs, rhs);
      out.cyclic.push({
        relation: AXES[t[0]].symbol + ' ' + AXES[t[1]].symbol + ' = i ' + AXES[t[2]].symbol,
        deviation: dev
      });
      out.worst = Math.max(out.worst, dev);
    });

    /* 4 — the three operators restricted to the two-dimensional sector,
       in the basis in which P is diagonal with eigenvalue +1 */
    var basis = [];
    for (i = 0; i < 4; i++) if (P.re[i * 4 + i] > 0) basis.push(i);
    AXES.forEach(function (ax) {
      var A = bilinear(ax.pair[0], ax.pair[1]);
      var cells = [];
      basis.forEach(function (r) {
        basis.forEach(function (c) {
          cells.push({ re: A.re[r * 4 + c], im: A.im[r * 4 + c] });
        });
      });
      out.sectorMatrices.push({ symbol: ax.symbol, cells: cells });
    });

    /* 5 — the three axes together have to cover each of the six pairs γ_iγ_j
       exactly once; there are C(4,2) = 6 pairs and 3 axes with 2 pairings. */
    var seen = {}, duplicates = [], count = 0;
    AXES.forEach(function (ax) {
      [ax.pair, ax.alt].forEach(function (pr) {
        var kk = Math.min(pr[0], pr[1]) + '-' + Math.max(pr[0], pr[1]);
        if (seen[kk]) duplicates.push(kk); else seen[kk] = true;
        count++;
      });
    });
    out.pairCoverage = {
      pairsUsed: count,
      distinctPairs: Object.keys(seen).length,
      duplicates: duplicates,
      ok: count === 6 && Object.keys(seen).length === 6 && duplicates.length === 0
    };
    if (!out.pairCoverage.ok) out.passed = false;

    out.passed = out.worst < TOL && out.pairCoverage.ok;
    out.tolerance = TOL;
    return out;
  }

  function sub(n) {
    return ['₀', '₁', '₂', '₃', '₄'][n];
  }

  MV.ausleseAlgebraCheck = runAlgebraCheck;

  /* Die Darstellung wird auch vom Gatter-Reiter gebraucht: dort ist zu
     zeigen, dass der Zopfoperator U₂₃ aus (5.11) und das Gatter Ŝ_z aus (5.13)
     derselbe Operator sind. Beide Reiter müssen dieselben γ verwenden, sonst
     prüft die Rechnung dort nur ihre eigene, zweite Darstellung gegen sich
     selbst. Deshalb hier exportiert und nicht dort neu aufgeschrieben. */
  MV.majoranaRep = {
    GAMMA: GAMMA,
    bilinear: bilinear,
    cmul: cmul, cadd: cadd, cscale: cscale, cid: cid,
    cdiff: cdiff, cmaxAbs: cmaxAbs,
    sectorMatrix: sectorMatrix, sectorElement: sectorElement, QUBIT_BASIS: QUBIT_BASIS
  };

  /* ====================================================================
     Zustand des Reiters

     Vorgaben sind die Illustrationswerte aus Plugge, Anh. B.3: mit t₀ = 5 und
     t₁ = 1 ist ω₊ = 6 und ω₋ = 4. ε ist die Verstimmung der Quantenpunkte
     (Plugge Gl. (4)) und kommt nur in diesem Reiter vor.
     ==================================================================== */

  var state = {
    t0: 5, t1: 1, phi: 0, eps: 0,
    tau: PI / 4,       /* Wartezeit τ_m der Zeitbereichs-Auslese (Plugge).
                          π/4 gibt mit t₀ = 5, t₁ = 1, φ = 0, ε = 0 vollen
                          Kontrast: P₊ = 0, P₋ = 1. Mit π/8 wären es nur 0,5. */
    gamma: 0,          /* Γ_tot der Frequenzbereichs-Auslese                  */
    axis: 2,
    pick: [2, 3], pending: null, visited: { '2-3': true },
    cycA: 0, cycB: 1
  };

  var sliders = {}, panelA = null, panelB = null, panelC = null;
  var readout = null, check = runAlgebraCheck();
  var pairHost = null, sectorHost = null, cyclicHost = null, wiringHost = null;
  var plotC1 = null, plotC2 = null, plotC3 = null;

  var N_CURVE = 480;

  /* Feste Werte der Frequenzbereichs-Auslese nach Plugge, Anh. B.3. */
  var LAMBDA = 2, OMEGA0 = 10, KAPPA_IN = 0.2, KAPPA_OUT = 0.2;

  /* ====================================================================
     Panel A — Geometrie nach Plugge et al., Fig. 3(a)

     Abbildung 5.1 der Arbeit (= Plugge Fig. 1a) zeigt nur die ẑ-Schleife.
     Für alle drei Achsen braucht es die Drei-Punkte-Geometrie aus Plugge
     Fig. 3(a): oberhalb der Box liegt ein zusätzlicher schwebender TS-Draht,
     der als ein einziges, über seine Länge ausgedehntes fermionisches Niveau
     wirkt und die phasenkohärente Verbindung zum fernen Ende der Box (γ₁)
     herstellt. Drei Quantenpunkte sitzen rechts: Punkt 1 am Zusatzdraht,
     Punkt 2 neben γ₂, Punkt 3 neben γ₃.

     γ₄ wird in dieser Geometrie nicht kontaktiert.
     ==================================================================== */

  var GEO = {
    yF: 8.2, yM: 5.6, yB: 3.0,   /* Zeilenmitten: Zusatzdraht, obere, untere */
    half: 0.42,                   /* halbe Höhe eines Drahtes                */
    tsL: 5.0, tsR: 16.0,          /* topologischer Abschnitt                 */
    semiL: 4.0, semiR: 17.0,      /* nicht proximisierte Bereiche            */
    bridgeL: 9.8, bridgeR: 11.0,
    refLx: 3.2, refLw: 0.8,       /* Referenzarm links: Zusatzdraht ↔ Box    */
    dotL: 17.4, dotR: 19.2,
    refRx: 19.8, refRw: 0.8       /* Referenzarm rechts: verbindet die Punkte */
  };

  function rowY(i) { return i === 1 ? GEO.yF : i === 2 ? GEO.yM : GEO.yB; }

  /** Lage der vier Nullmoden. γ₁ und γ₄ links, γ₂ und γ₃ rechts. */
  function siteOf(idx) {
    var xL = GEO.tsL + 0.6, xR = GEO.tsR - 0.6;
    return [null,
      [xL, GEO.yM],   /* γ₁ oben links   */
      [xR, GEO.yM],   /* γ₂ oben rechts  */
      [xR, GEO.yB],   /* γ₃ unten rechts */
      [xL, GEO.yB]    /* γ₄ unten links  */
    ][idx];
  }

  /** Mitte des Quantenpunkts, der diese Nullmode anspricht. */
  function dotOf(i) { return [(GEO.dotL + GEO.dotR) / 2, rowY(i)]; }

  /** Ist das Paar in dieser Geometrie verdrahtet? Alles ohne γ₄ ist es. */
  function isWired(i, j) { return i !== 4 && j !== 4; }

  /**
   * Weg vom Quantenpunkt zu seiner Nullmode. Punkt 2 und 3 sitzen direkt
   * neben γ₂ und γ₃; Punkt 1 erreicht γ₁ nur über den Zusatzdraht und den
   * linken Referenzarm — genau dafür ist der Zusatzdraht da.
   */
  function dotToMode(i) {
    var d = dotOf(i), s = siteOf(i);
    if (i === 1) {
      var rx = GEO.refLx + GEO.refLw / 2;
      return [d, [GEO.semiL + 0.6, GEO.yF], [rx, GEO.yF], [rx, GEO.yM], s];
    }
    return [d, s];
  }

  /** Weg durch die Box: Punkt → Nullmode → Nullmode → Punkt. Amplitude t₁. */
  function boxPath(i, j) {
    var a = dotToMode(i), b = dotToMode(j);
    return a.concat(b.slice().reverse());
  }

  /** Weg über den Referenzarm rechts, der dieselben zwei Punkte verbindet. */
  function refPath(i, j) {
    var rx = GEO.refRx + GEO.refRw / 2;
    return [[GEO.dotR, rowY(i)], [rx, rowY(i)], [rx, rowY(j)], [GEO.dotR, rowY(j)]];
  }

  /* Strichstärke proportional zur Amplitude; eine Amplitude null zeichnet
     keinen Weg. Bezug sind die Reglerobergrenzen 8 und 4. */
  function widthFor(amp, max) { return 0.6 + 3.4 * Math.min(1, Math.abs(amp) / max); }
  function alphaFor(amp) { return amp < 1e-9 ? 0 : 0.25 + 0.75 * Math.min(1, Math.abs(amp) / 4); }

  function drawPanelA() {
    var p = panelA.plot;
    var sel = state.pick;
    var info = axisOfPair(sel[0], sel[1]);
    var ax = AXES[info.idx];
    var other = (info.which === 'pair') ? ax.alt : ax.pair;
    var wired = isWired(sel[0], sel[1]);

    var m = p.measure();
    var X = 20.4;                       /* Inhalt liegt zwischen x = 3,2 und 20,6 */
    var x0 = 1.7;
    var Y = X * m.ph / m.pw;
    p.begin({ x: [x0, x0 + X], y: [(11.2 - Y) / 2, (11.2 + Y) / 2],
              grid: false, frame: false, xTicks: [], yTicks: [] });

    var wireFill = '#cfe2f3', wireEdge = '#5b8fc4';
    var semiFill = '#d6ead2', semiEdge = '#6f9c66';
    var dotFill = '#f3c9cf', dotEdge = '#b5707c';
    var bridgeFill = '#e8873c';
    var floatFill = '#dbe8f5';
    var h = GEO.half;

    /* ---- nicht proximisierte Bereiche und Referenzarme ---- */
    [GEO.yF, GEO.yM, GEO.yB].forEach(function (y) {
      p.rect(GEO.semiL, y - h, GEO.tsL, y + h, { fill: semiFill, stroke: semiEdge });
      p.rect(GEO.tsR, y - h, GEO.semiR, y + h, { fill: semiFill, stroke: semiEdge });
    });
    /* Beide Arme heißen bei Plugge R, tun aber Verschiedenes: der rechte ist
       der Referenzarm mit Amplitude t₀, der linke gehört zum Weg durch die
       Box. In der Legende deshalb getrennt. */
    p.rect(GEO.refLx, GEO.yM - h, GEO.refLx + GEO.refLw, GEO.yF + h,
      { fill: semiFill, stroke: semiEdge,
        legend: 'R links: Verbindung Zusatzdraht → Box-Ende bei γ₁ — Teil des ' +
                'Weges durch die Box' });
    p.rect(GEO.refRx, GEO.yB - h, GEO.refRx + GEO.refRw, GEO.yF + h,
      { fill: semiFill, stroke: semiEdge,
        legend: 'R rechts: Referenzarm, verbindet die drei Punkte, Amplitude t₀' });
    p.text('R', GEO.refLx + GEO.refLw / 2, (GEO.yM + GEO.yF) / 2,
      { color: '#4d7a45', align: 'center', baseline: 'middle', font: 'bold 12px "Segoe UI", Arial' });
    p.text('R', GEO.refRx + GEO.refRw / 2, (GEO.yB + GEO.yF) / 2,
      { color: '#4d7a45', align: 'center', baseline: 'middle', font: 'bold 12px "Segoe UI", Arial' });

    /* ---- Brücke und die drei Drähte ---- */
    p.rect(GEO.bridgeL, GEO.yB + h, GEO.bridgeR, GEO.yM - h,
      { fill: bridgeFill, stroke: bridgeFill, legend: 's-Wellen-Brücke S (trivial)' });
    p.rect(GEO.tsL, GEO.yF - h, GEO.tsR, GEO.yF + h,
      { fill: floatFill, stroke: wireEdge, dash: [5, 3],
        legend: 'schwebender TS-Zusatzdraht' });
    [GEO.yM, GEO.yB].forEach(function (y, k) {
      p.rect(GEO.tsL, y - h, GEO.tsR, y + h, {
        fill: wireFill, stroke: wireEdge,
        legend: k === 0 ? 'topologische Drähte der Box' : undefined
      });
    });
    /* Der Zusatzdraht ist selbst topologisch und trägt an seinen Enden zwei
       eigene Nullmoden. Sie gehören nicht zum Tetron und bekommen deshalb
       keine γ-Nummer; neutral gezeichnet, damit sie nicht mit γ₁ … γ₄
       verwechselt werden. */
    [GEO.tsL, GEO.tsR].forEach(function (x, i) {
      p.marker(x, GEO.yF, {
        r: 4.5, shape: 'cross', color: MV.palette.muted, width: 2,
        legend: i === 0
          ? 'eigene Nullmoden des Zusatzdrahts — er wirkt als ein über die ganze ' +
            'Länge ausgedehntes Niveau (Plugge et al.)'
          : null
      });
    });

    p.text('TS', GEO.tsL + 1.6, GEO.yF, { color: '#2c5f92', align: 'center', baseline: 'middle', font: '11px "Segoe UI", Arial' });
    p.text('TS', GEO.tsL + 1.6, GEO.yM, { color: '#2c5f92', align: 'center', baseline: 'middle', font: '11px "Segoe UI", Arial' });
    p.text('TS', GEO.tsL + 1.6, GEO.yB, { color: '#2c5f92', align: 'center', baseline: 'middle', font: '11px "Segoe UI", Arial' });
    p.text('S', (GEO.bridgeL + GEO.bridgeR) / 2, (GEO.yB + GEO.yM) / 2,
      { color: '#ffffff', align: 'center', baseline: 'middle', font: 'bold 12px "Segoe UI", Arial' });

    /* ---- die drei Quantenpunkte ---- */
    [1, 2, 3].forEach(function (i) {
      var y = rowY(i), open = wired && sel.indexOf(i) >= 0;
      p.rect(GEO.dotL, y - h, GEO.dotR, y + h, {
        fill: open ? MV.palette.topologicalSoft : dotFill,
        stroke: open ? MV.palette.topological : dotEdge,
        legend: i === 1 ? 'Quantenpunkte 1, 2, 3' : undefined
      });
      p.text(String(i), (GEO.dotL + GEO.dotR) / 2, y, {
        color: open ? MV.palette.topological : MV.palette.ink,
        align: 'center', baseline: 'middle', font: 'bold 12px "Segoe UI", Arial'
      });
    });

    /* ---- die zweite Schleife derselben Achse, schwächer ---- */
    var o1 = siteOf(other[0]), o2 = siteOf(other[1]);
    p.curve([o1, o2], {
      color: MV.palette.critical, width: 4, alpha: 0.28, dash: [8, 5],
      legend: 'zweite Schleife derselben Achse nach ' + MV.eq('tetron.decomposition')
    });

    /* ---- die beiden Wege des gewählten Paars ---- */
    if (wired) {
      if (state.t0 > 1e-9) {
        p.curve(refPath(sel[0], sel[1]), {
          color: MV.palette.curve3, width: widthFor(state.t0, 8), alpha: alphaFor(state.t0),
          legend: 'Weg über den Referenzarm, Amplitude t₀'
        });
      }
      if (state.t1 > 1e-9) {
        p.curve(boxPath(sel[0], sel[1]), {
          color: MV.palette.topological, width: widthFor(state.t1, 4), alpha: alphaFor(state.t1),
          legend: 'Weg durch die Box, Amplitude t₁'
        });
      }
      var bp = boxPath(sel[0], sel[1]);
      var rp = refPath(sel[0], sel[1]);
      var cx = (bp[bp.length - 2][0] + rp[1][0]) / 2;
      var cy = (rowY(sel[0]) + rowY(sel[1])) / 2;
      p.text('φ', cx, cy, { color: MV.palette.ink, align: 'center', baseline: 'middle',
        font: 'italic 15px "Segoe UI", Arial', box: true });
    }

    /* ---- die vier Nullmoden ---- */
    [1, 2, 3, 4].forEach(function (i) {
      var s = siteOf(i);
      var on = sel.indexOf(i) >= 0;
      if (on) p.marker(s[0], s[1], { color: MV.palette.topological, r: 10, shape: 'ring', width: 2 });
      p.marker(s[0], s[1], {
        color: on ? MV.palette.topological : MV.palette.muted,
        shape: 'cross', r: 6, width: on ? 3 : 2
      });
      p.text('γ' + sub(i), s[0], s[1] + (i <= 2 ? 0.95 : -0.95), {
        color: on ? MV.palette.topological : MV.palette.ink,
        align: 'center', baseline: 'middle', font: 'bold 13px "Segoe UI", Arial', box: true
      });
    });

    p.text(MV.figShort('tetron.layoutFig') + ' der Arbeit zeigt nur die ẑ-Schleife; die ' +
           'Verdrahtung aller drei Achsen folgt Plugge et al., Fig. 3.',
      x0 + X / 2, (11.2 - Y) / 2 + 0.45, {
        color: MV.palette.muted, align: 'center', baseline: 'middle',
        font: '11px "Segoe UI", Arial'
      });

    p.finish();

    panelA.setStatus(
      wired
        ? 'Gates öffnen die Kontakte an Punkt ' + sel.map(function (i) { return i; }).join(' und ') +
          ' und schließen die Schleife. Beide Wege umschließen den Fluss φ; nach ' +
          MV.eq('tetron.readoutH') + ' addieren sich ihre Amplituden zu t_z = t₀ + t₁z.'
        : 'γ₄ wird in dieser Geometrie nicht kontaktiert. Die Paarung ist nach ' +
          MV.eq('tetron.decomposition') + ' gleichwertig zu ' + ax.symbol +
          ', aber nicht verdrahtet — deshalb ist hier kein Weg gezeichnet.',
      wired ? MV.palette.topological : MV.palette.critical);

    if (wiringHost) {
      wiringHost.hidden = wired;
      wiringHost.innerHTML = wired ? '' :
        '<b>Nicht verdrahtet.</b> Die Drei-Punkte-Geometrie aus Plugge et al., ' +
        'Fig. 3(a), kontaktiert γ₁, γ₂ und γ₃. γ₄ bleibt unkontaktiert, deshalb ' +
        'gibt es für die Paarungen mit γ₄ keinen Interferenzweg. Nach ' +
        MV.eq('tetron.decomposition') + ' beschreiben sie dieselbe Achse wie die jeweils ' +
        'verdrahtete Paarung; gemessen wird die Achse also trotzdem, nur über die ' +
        'andere Schleife.';
    }

    renderPairGrid();
    renderSector();
  }

  /* Both designations side by side, this thesis first and dominant. */
  /* Eine Zeile, keine drei Sätze Metatext. Die Begründung steht einmal für
     das ganze Panel in Punkt 6 des Info-Panels. */
  var schematicHost = null;

  /* ====================================================================
     Auswahl: anklicken, umschalten, zählen
     ==================================================================== */

  /** Ein Paar setzen; die Reihenfolge bleibt die von (5.6) beziehungsweise (5.8). */
  function selectPair(i, j) {
    var info = axisOfPair(i, j);
    if (!info) return;
    var ax = AXES[info.idx];
    state.pick = (info.which === 'pair') ? ax.pair.slice() : ax.alt.slice();
    state.axis = info.idx;
    state.visited[pairKey(i, j)] = true;
    state.pending = null;
    refreshAll();
  }

  /**
   * Klick auf eine Mode. Der erste Klick merkt sich eine Mode, der zweite
   * vervollständigt das Paar. Ein Klick auf die gemerkte Mode nimmt sie
   * zurück; vier Majoranas lassen sich immer paaren, ein Fehlgriff ist also
   * nicht möglich.
   */
  function clickMode(i) {
    if (state.pending === null) {
      if (state.pick.indexOf(i) >= 0 && state.pick.length === 2) {
        /* aus dem aktuellen Paar heraus: die andere Mode bleibt stehen */
        state.pending = state.pick[0] === i ? state.pick[1] : state.pick[0];
        drawPanelA();
        return;
      }
      state.pending = i;
      drawPanelA();
      return;
    }
    if (state.pending === i) { state.pending = null; drawPanelA(); return; }
    selectPair(state.pending, i);
  }

  /** Pixel → Datenkoordinaten; Plot hält x0/x1/px0/px1 nach begin(). */
  function dataAt(p, ev) {
    var r = p.canvas.getBoundingClientRect();
    var px = ev.clientX - r.left, py = ev.clientY - r.top;
    return [
      p.x0 + (px - p.px0) / (p.px1 - p.px0) * (p.x1 - p.x0),
      p.y0 + (p.py1 - py) / (p.py1 - p.py0) * (p.y1 - p.y0)
    ];
  }

  function attachClicks() {
    var cv = panelA.canvas;
    cv.style.cursor = 'pointer';
    cv.addEventListener('click', function (ev) {
      var q = dataAt(panelA.plot, ev);
      var best = null, bestD = Infinity, i, s, d;
      for (i = 1; i <= 4; i++) {
        s = siteOf(i);
        d = Math.hypot(q[0] - s[0], q[1] - s[1]);
        if (d < bestD) { bestD = d; best = i; }
      }
      if (bestD <= 1.1) clickMode(best);
    });
  }

  /* ---------------------------------------------------- die sechs Paare */

  function renderPairGrid() {
    var seenAxes = {}, count = 0, k;
    for (k in state.visited) {
      if (Object.prototype.hasOwnProperty.call(state.visited, k)) {
        count++;
        var pr = k.split('-');
        seenAxes[axisOfPair(+pr[0], +pr[1]).idx] = true;
      }
    }
    var axesSeen = Object.keys(seenAxes).length;

    var html = '<div class="pairs">';
    AXES.forEach(function (ax, idx) {
      html += '<div class="pairs__axis">' +
        '<span class="pairs__name">' + ax.symbol + '<span class="pairs__geo">' +
        ax.geometry.split('—')[0].trim() + '</span></span>';
      [ax.pair, ax.alt].forEach(function (pr) {
        var key = pairKey(pr[0], pr[1]);
        var cur = (pairKey(state.pick[0], state.pick[1]) === key);
        html += '<button type="button" class="pairs__chip' +
          (cur ? ' is-current' : '') + (state.visited[key] ? ' is-visited' : '') +
          '" data-pair="' + pr[0] + ',' + pr[1] + '"' +
          ' aria-pressed="' + (cur ? 'true' : 'false') + '">' +
          'γ' + sub(pr[0]) + 'γ' + sub(pr[1]) + '</button>';
      });
      html += '</div>';
    });
    html += '</div>';

    html += '<p class="pairs__count">' + count + ' von 6 Paaren ausprobiert, ' +
      axesSeen + ' von 3 Achsen gefunden.';
    if (count === 6) {
      html += ' <b>Damit sind alle Möglichkeiten erschöpft:</b> vier Majoranas ' +
        'lassen sich auf C(4,2) = 6 Weisen zu Paaren zusammenfassen, und je zwei ' +
        'dieser Paare beschreiben dieselbe Achse — es bleiben genau drei.';
    }
    html += '</p>';

    pairHost.innerHTML = html;
    Array.prototype.forEach.call(pairHost.querySelectorAll('.pairs__chip'), function (b) {
      b.addEventListener('click', function () {
        var pr = b.getAttribute('data-pair').split(',');
        selectPair(+pr[0], +pr[1]);
      });
    });
  }

  /* ------------------------------------- die Rechnung zum gewählten Paar */

  function renderSector() {
    var sel = state.pick;
    var info = axisOfPair(sel[0], sel[1]);
    var ax = AXES[info.idx];
    var other = (info.which === 'pair') ? ax.alt : ax.pair;

    var M = sectorMatrix(sel[0], sel[1]);
    var Mo = sectorMatrix(other[0], other[1]);
    var S = sectorSquare(sel[0], sel[1]);
    var devP = mat2Diff(M, PAULI[ax.key]);
    var devO = mat2Diff(M, Mo);
    var devS = mat2Diff(S, PAULI.id);

    sectorHost.innerHTML =
      '<div class="sector">' +
      '<div class="sector__row">' +
      '<span class="sector__lab">i γ' + sub(sel[0]) + 'γ' + sub(sel[1]) +
      ' im Sektor P̂ = +1</span>' + mat2Html(M) +
      '<span class="sector__eq">=</span>' +
      mat2Html(PAULI[ax.key], 'mat2--ref') +
      '<span class="sector__lab">' + ax.symbol + ' nach ' + MV.eq('tetron.pauli') +
      '</span>' +
      '<span class="sector__dev">Abweichung ' + fmtDev(devP) + '</span>' +
      '</div>' +
      '<div class="sector__row">' +
      '<span class="sector__lab">zweite Schleife i γ' + sub(other[0]) + 'γ' + sub(other[1]) +
      '</span>' + mat2Html(Mo) +
      '<span class="sector__dev">dieselbe Matrix, Abweichung ' + fmtDev(devO) +
      ' — welche der beiden Schleifen verdrahtet wird, ändert an der gemessenen ' +
      'Größe nichts (' + MV.eq('tetron.decomposition') + ')</span>' +
      '</div>' +
      '<div class="sector__row">' +
      '<span class="sector__lab">Quadrat</span>' + mat2Html(S) +
      '<span class="sector__dev">Einheitsmatrix, Abweichung ' + fmtDev(devS) +
      ' — der Eigenwert ist ±1</span>' +
      '</div>' +
      '</div>';
  }

  /* ------------------------------------------- zyklische Relation (5.7) */

  function renderCyclic() {
    var a = AXES[state.cycA], b = AXES[state.cycB];
    var A = sectorMatrix(a.pair[0], a.pair[1]);
    var B = sectorMatrix(b.pair[0], b.pair[1]);
    var prod = mat2Mul(A, B);

    var html = '<div class="sector sector--cyclic"><div class="sector__row">' +
      '<span class="sector__lab">' + a.symbol + ' · ' + b.symbol + '</span>' +
      mat2Html(A) + '<span class="sector__eq">·</span>' + mat2Html(B) +
      '<span class="sector__eq">=</span>' + mat2Html(prod);

    if (state.cycA === state.cycB) {
      html += '<span class="sector__dev">Gleiche Achse: das Quadrat ist die ' +
        'Einheitsmatrix, Abweichung ' + fmtDev(mat2Diff(prod, PAULI.id)) + '.</span>';
    } else {
      /* die dritte Achse ist die, die weder a noch b ist */
      var third = [0, 1, 2].filter(function (t) { return t !== state.cycA && t !== state.cycB; })[0];
      var C = sectorMatrix(AXES[third].pair[0], AXES[third].pair[1]);
      var cyc = (state.cycB - state.cycA + 3) % 3 === 1;
      var rhs = cyc ? mat2TimesI(C) : mat2TimesI(C).map(function (c) { return [-c[0], -c[1]]; });
      html += '<span class="sector__eq">=</span>' + mat2Html(rhs, 'mat2--ref') +
        '<span class="sector__lab">' + (cyc ? '' : '−') + 'i ' + AXES[third].symbol + '</span>' +
        '<span class="sector__dev">Abweichung ' + fmtDev(mat2Diff(prod, rhs)) +
        (cyc ? ' — zyklisch, wie ' + MV.eq('tetron.cyclic') + ' es verlangt.'
             : ' — antizyklisch, deshalb das Minuszeichen.') + '</span>';
    }
    html += '</div></div>';
    cyclicHost.innerHTML = html;
  }


  /* ====================================================================
     Die Primitive und die drei Verfahren

     Alle Formeln dieses Abschnitts stammen aus Plugge et al., Gl. (2)–(6)
     und Anhang B; die Arbeit beschreibt die Verfahren in §5.5 nur
     qualitativ. Einzige Ausnahme ist (5.10), die Leitfähigkeit, die auch in
     der Arbeit steht.
     ==================================================================== */

  /** t_z = t₀ + t₁ z, mit t₀ reell und t₁ = |t₁| e^{iφ}.  Plugge (2). */
  function tz(z) {
    var re = state.t0 + z * state.t1 * Math.cos(state.phi);
    var im = z * state.t1 * Math.sin(state.phi);
    return { re: re, im: im, abs: Math.hypot(re, im) };
  }

  /** Rabi-Frequenz ω_z = √(ε² + |t_z|²).  Plugge (4). */
  function omegaZ(z) {
    var a = tz(z).abs;
    return Math.sqrt(state.eps * state.eps + a * a);
  }

  /** |t₊|² − |t₋|² = 4|t₀||t₁| cos φ — der einzige Term, der z trägt. */
  function splitAbs2() { return 4 * state.t0 * state.t1 * Math.cos(state.phi); }

  /** Unterscheiden die Verfahren die beiden Paritäten überhaupt? */
  function indistinguishable() { return Math.abs(splitAbs2()) < 1e-12; }

  function whyBlind() {
    if (state.t0 < 1e-9) return 'kein Referenzarm: t₀ = 0';
    if (state.t1 < 1e-9) return 'kein Weg durch die Box: t₁ = 0';
    if (Math.abs(Math.cos(state.phi)) < 1e-9) return 'blinder Arbeitspunkt: φ = π/2';
    return '';
  }

  /* Verfahren 1 — Leitfähigkeit.  Arbeit (5.10), Plugge (3). */
  function conductance(phi, z) {
    return state.t0 * state.t0 + state.t1 * state.t1 +
           2 * state.t0 * state.t1 * z * Math.cos(phi);
  }
  function deltaGmax() { return 4 * state.t0 * state.t1; }

  /* Verfahren 2 — Zeitbereich.  Plugge (B.5). */
  function pReturn(z, tau) {
    var w = omegaZ(z);
    if (w < 1e-12) return 1;                 /* ω_z = 0: der Zustand bleibt */
    var c = Math.cos(w * tau), s = Math.sin(w * tau);
    return c * c + (state.eps * state.eps / (w * w)) * s * s;
  }

  /* Verfahren 3 — Frequenzbereich.  Plugge (5), (6). */
  function amplitude(w, z) {
    var W = omegaZ(z);
    var g = (W < 1e-12) ? 0 : -LAMBDA * tz(z).abs / (2 * W);   /* nackter Resonator */
    /* χ_z = g² / [ i Γ_tot − (2ω_z − ω) ] */
    var dr = -(2 * W - w), di = state.gamma;
    var d2 = dr * dr + di * di;
    var chr = 0, chi = 0;
    if (d2 > 1e-18) { chr = g * g * dr / d2; chi = -g * g * di / d2; }
    else { return { re: 0, im: 0, abs2: 0, phase: 0 }; }        /* ω = 2ω_z, Γ = 0 */
    var br = (OMEGA0 - w) + chr;
    var bi = -(KAPPA_IN + KAPPA_OUT) / 2 + chi;
    var ni = -Math.sqrt(KAPPA_IN * KAPPA_OUT);
    var b2 = br * br + bi * bi;
    if (b2 < 1e-18) return { re: 0, im: 0, abs2: 0, phase: 0 };
    var re = (ni * bi) / b2, im = (ni * br) / b2;
    return { re: re, im: im, abs2: re * re + im * im, phase: -Math.atan2(im, re) };
  }

  /* ====================================================================
     Panel B — Die Primitive als Zeigerdiagramm
     ==================================================================== */

  function drawPanelB() {
    var p = panelB.plot;
    var tp = tz(+1), tm = tz(-1);
    var R = Math.max(state.t0 + state.t1, 1) * 1.25;

    var m = p.measure();
    var xLo = -0.18 * R, xHi = R * 1.08;
    var Yh = (xHi - xLo) * m.ph / m.pw;

    p.begin({ x: [xLo, xHi], y: [-Yh / 2, Yh / 2],
              xLabel: 'Re t', yLabel: 'Im t', xTickCount: 5, yTickCount: 5 });

    p.hline(0, { color: MV.palette.line });
    p.vline(0, { color: MV.palette.line });

    /* t₀ auf der reellen Achse */
    p.curve([[0, 0], [state.t0, 0]], { color: MV.palette.curve3, width: 3, legend: 't₀ (Referenzarm)' });

    /* von der Spitze aus ±t₁ e^{iφ} */
    var dx = state.t1 * Math.cos(state.phi), dy = state.t1 * Math.sin(state.phi);
    p.curve([[state.t0, 0], [state.t0 + dx, dy]], { color: MV.palette.muted, width: 2, dash: [5, 3] });
    p.curve([[state.t0, 0], [state.t0 - dx, -dy]], { color: MV.palette.muted, width: 2, dash: [5, 3],
      legend: '± t₁ e^{iφ} (Box)' });

    /* die beiden Summen */
    p.curve([[0, 0], [tp.re, tp.im]], { color: MV.palette.curve1, width: 2.6, legend: 't₊ = t₀ + t₁e^{iφ}' });
    p.curve([[0, 0], [tm.re, tm.im]], { color: MV.palette.curve2, width: 2.6, dash: [7, 4],
      legend: 't₋ = t₀ − t₁e^{iφ}' });
    p.marker(tp.re, tp.im, { color: MV.palette.curve1, r: 4 });
    p.marker(tm.re, tm.im, { color: MV.palette.curve2, r: 4 });
    p.marker(state.t0, 0, { color: MV.palette.curve3, r: 3.5 });

    /* Bei φ = 0 liegen alle drei Zeiger auf der reellen Achse übereinander.
       Das ist richtig so — die Spitzenbeschriftung hält sie trotzdem
       auseinander. */
    p.text('t₀', state.t0, 0, { color: MV.palette.curve3, align: 'center',
      baseline: 'middle', font: 'bold 11px "Segoe UI", Arial', box: true });
    p.text('t₊', tp.re, tp.im, { color: MV.palette.curve1, align: 'left',
      baseline: 'bottom', font: 'bold 12px "Segoe UI", Arial', box: true });
    p.text('t₋', tm.re, tm.im, { color: MV.palette.curve2, align: 'right',
      baseline: 'top', font: 'bold 12px "Segoe UI", Arial', box: true });

    p.finish();

    var blind = indistinguishable();
    var msg =
      '|t₊| = ' + MV.fmt(tp.abs, 3) + ' · |t₋| = ' + MV.fmt(tm.abs, 3) +
      ' · |t₊|² − |t₋|² = ' + MV.fmtSigned(splitAbs2(), 3) + '   →   ' +
      (blind ? 'nicht unterscheidbar — |t₊| = |t₋| (' + whyBlind() + ')' : 'unterscheidbar') +
      '   ·   z = ' + zName();

    /* Bei φ = 0 liegen t₀, t₊ und t₋ alle auf der reellen Achse; das Dreieck
       der Zeiger ist dann zu einer Strecke entartet und nicht zu sehen.
       φ = 0 bleibt Vorgabe, weil Panel C damit den besten Kontrast zeigt. */
    if (Math.abs(Math.sin(state.phi)) < 0.05) {
      msg += '   ·   Bei φ = 0 liegen alle Zeiger auf der reellen Achse — ' +
             'zieh φ, um das Dreieck zu sehen.';
    }

    panelB.setStatus(msg, blind ? MV.palette.critical : MV.palette.topological);
  }

  /* ====================================================================
     Panel C — Drei Verfahren, eine Größe
     ==================================================================== */

  function drawC1() {
    var p = plotC1, i, phi;
    var xMax = 4 * PI;
    var Y = MV.ceilTo(Math.max(conductance(0, 1), 0.5) * 1.15, 0.5);
    var ticks = [];
    for (i = 0; i <= 4; i++) ticks.push({ v: i * PI, label: i === 0 ? '0' : (i === 1 ? 'π' : i + 'π') });

    p.begin({ x: [0, xMax], y: [0, Y], xTicks: ticks, yTickCount: 5,
              xLabel: 'φ', yLabel: 'G_z / [ (e²/h) ν₁ν₂ ]' });

    var plus = [], minus = [];
    for (i = 0; i <= N_CURVE; i++) {
      phi = xMax * i / N_CURVE;
      plus.push([phi, conductance(phi, +1)]);
      minus.push([phi, conductance(phi, -1)]);
    }
    for (var n = 0; (n + 0.5) * PI <= xMax; n++) {
      p.vline((n + 0.5) * PI, { color: MV.palette.critical, dash: [2, 4] });
    }
    p.curve(plus, { color: MV.palette.curve1, width: 2, legend: 'z = +1' });
    p.curve(minus, { color: MV.palette.curve2, width: 2, dash: [7, 4], legend: 'z = −1' });
    for (n = 0; n * PI <= xMax + 1e-9; n++) {
      p.marker(n * PI, conductance(n * PI, +1), { color: MV.palette.topological, r: 3.5, shape: 'ring', width: 1.6 });
      p.marker(n * PI, conductance(n * PI, -1), { color: MV.palette.topological, r: 3.5, shape: 'ring', width: 1.6 });
    }
    p.vline(state.phi, { color: MV.palette.ink, width: 1.5 });
    p.finish();
  }

  function drawC2() {
    var p = plotC2, i, t;
    var wMax = Math.max(omegaZ(1), omegaZ(-1), 0.5);
    var xMax = Math.max(2 * PI / wMax, 0.2) * 1.6;
    p.begin({ x: [0, xMax], y: [-0.05, 1.08], yTickCount: 5, xTickCount: 5,
              xLabel: 'Wartezeit τ_m', yLabel: 'P_z(2; τ_m)' });

    var plus = [], minus = [];
    for (i = 0; i <= N_CURVE; i++) {
      t = xMax * i / N_CURVE;
      plus.push([t, pReturn(+1, t)]);
      minus.push([t, pReturn(-1, t)]);
    }
    p.curve(plus, { color: MV.palette.curve1, width: 2, legend: 'z = +1' });
    p.curve(minus, { color: MV.palette.curve2, width: 2, dash: [7, 4], legend: 'z = −1' });
    p.vline(state.tau, { color: MV.palette.ink, width: 1.5 });
    p.marker(state.tau, pReturn(+1, state.tau), { color: MV.palette.curve1, r: 4 });
    p.marker(state.tau, pReturn(-1, state.tau), { color: MV.palette.curve2, r: 4 });
    p.finish();
  }

  function drawC3() {
    var i, w;
    var w0 = OMEGA0 - 1.5, w1 = OMEGA0 + 1.5;

    var p = plotC3.abs;
    p.begin({ x: [w0, w1], y: [0, 1.12], xTickCount: 5, yTickCount: 4,
              xLabel: 'ω', yLabel: '|A_ω|²' });
    var ap = [], am = [];
    for (i = 0; i <= N_CURVE; i++) {
      w = w0 + (w1 - w0) * i / N_CURVE;
      ap.push([w, amplitude(w, +1).abs2]);
      am.push([w, amplitude(w, -1).abs2]);
    }
    p.curve(ap, { color: MV.palette.curve1, width: 2, legend: 'z = +1' });
    p.curve(am, { color: MV.palette.curve2, width: 2, dash: [7, 4], legend: 'z = −1' });
    p.vline(OMEGA0, { color: MV.palette.ink, dash: [3, 4] });
    p.finish();

    var q = plotC3.phase;
    q.begin({ x: [w0, w1], y: [-185, 185], xTickCount: 5,
              yTicks: [-180, -90, 0, 90, 180].map(function (v) { return { v: v, label: v + '°' }; }),
              xLabel: 'ω', yLabel: 'Δφ(ω)' });
    var pp = [], pm = [];
    for (i = 0; i <= N_CURVE; i++) {
      w = w0 + (w1 - w0) * i / N_CURVE;
      pp.push([w, amplitude(w, +1).phase * 180 / PI]);
      pm.push([w, amplitude(w, -1).phase * 180 / PI]);
    }
    q.curve(pp, { color: MV.palette.curve1, width: 2, legend: 'z = +1' });
    q.curve(pm, { color: MV.palette.curve2, width: 2, dash: [7, 4], legend: 'z = −1' });
    q.vline(OMEGA0, { color: MV.palette.ink, dash: [3, 4] });
    q.marker(OMEGA0, amplitude(OMEGA0, +1).phase * 180 / PI, { color: MV.palette.curve1, r: 4 });
    q.marker(OMEGA0, amplitude(OMEGA0, -1).phase * 180 / PI, { color: MV.palette.curve2, r: 4 });
    q.finish();
  }

  /** Maximum von |A_ω|² in der Nähe von ω₀, für die Anzeige. */
  function peakOf(z) {
    var lo = OMEGA0 - 1.5, hi = OMEGA0 + 1.5, best = lo, bm = -1, i, w, m;
    for (i = 0; i <= 3000; i++) {
      w = lo + (hi - lo) * i / 3000;
      m = amplitude(w, z).abs2;
      if (m > bm) { bm = m; best = w; }
    }
    return { w: best, abs2: bm };
  }

  function drawPanelC() {
    drawC1(); drawC2(); drawC3();
    var dc = Math.abs(pReturn(+1, state.tau) - pReturn(-1, state.tau));
    var pk = peakOf(+1), mk = peakOf(-1);
    panelC.setStatus(
      'C1: δG_max = ' + MV.fmt(deltaGmax(), 3) +
      ' · C2: Kontrast bei τ_m = ' + MV.fmt(state.tau, 2) + ' ist ' + MV.fmt(dc, 3) +
      ' · C3: Maxima bei ω = ' + MV.fmt(pk.w, 3) + ' und ' + MV.fmt(mk.w, 3) +
      ', Δφ(ω₀) = ' + MV.fmt(amplitude(OMEGA0, +1).phase * 180 / PI, 1) + '° und ' +
      MV.fmt(amplitude(OMEGA0, -1).phase * 180 / PI, 1) + '°',
      indistinguishable() ? MV.palette.critical : MV.palette.muted);
  }

  /* ====================================================================
     Aufbau des Reiters
     ==================================================================== */

  var PLUGGE = 'Plugge et al., New J. Phys. 19, 012001 (2017)';

  var UNITS_NOTE =
    'Willkürliche Einheiten wie bei Plugge (ℏ = 1); G in Einheiten von ' +
    '(e²/h)·ν₁ν₂. ν₁ und ν₂ sind die Zustandsdichten der Zuleitungen, nicht die ' +
    'Z₂-Invariante.';

  var PHI_NOTE =
    'φ ist in diesem Reiter der dimensionslose eingeschlossene Fluss, nicht die ' +
    'supraleitende Phase aus Kapitel 3 und 4. Ein Flussquant wird bewusst nicht ' +
    'angegeben — φ ist hier eine reine Phase.';

  var PARITY_NOTE =
    'P̂ ist hier die Gesamtparität −γ₁γ₂γ₃γ₄ nach ' + MV.eq('tetron.parity') + ', nicht der ' +
    'Operator (−1)^N̂ aus Kapitel 2.';

  function refreshAll() {
    var tp = tz(+1), tm = tz(-1);
    readout.set('tplus', MV.fmt(tp.abs, 3));
    readout.set('tminus', MV.fmt(tm.abs, 3));
    readout.set('omega', MV.fmt(omegaZ(+1), 3) + ' / ' + MV.fmt(omegaZ(-1), 3));
    readout.set('split', MV.fmtSigned(splitAbs2(), 3));
    readout.set('dg', MV.fmt(deltaGmax(), 3));
    readout.status(
      indistinguishable() ? 'critical' : 'topological',
      indistinguishable()
        ? 'Kein Verfahren unterscheidet die Paritäten — ' + whyBlind()
        : 'Alle drei Verfahren lesen dieselbe Größe |t_z|');

    drawPanelA();
    drawPanelB();
    drawPanelC();
    renderCyclic();
    if (panelA.info) panelA.info.update();
    if (panelB.info) panelB.info.update();
    if (panelC.info) panelC.info.update();
  }

  function buildControls(mount) {
    var box = MV.ui.el('div', { class: 'controls' });
    function on(k) { return function (v) { state[k] = v; refreshAll(); }; }

    sliders.t0 = MV.ui.slider(box, {
      label: '|t₀| — Amplitude über den Referenzarm',
      min: 0, max: 8, step: 0.05, value: state.t0, digits: 2, onInput: on('t0')
    });
    sliders.t1 = MV.ui.slider(box, {
      label: '|t₁| — Kotunnelamplitude durch die Box',
      min: 0, max: 4, step: 0.05, value: state.t1, digits: 2, onInput: on('t1')
    });
    sliders.phi = MV.ui.slider(box, {
      label: 'φ — eingeschlossener Fluss',
      min: 0, max: 2 * PI, step: PI / 60, value: state.phi, digits: 2,
      hint: 'Schrittweite π/60, damit π/2 und π genau zu treffen sind.',
      onInput: on('phi')
    });
    sliders.eps = MV.ui.slider(box, {
      label: 'ε — Verstimmung der Quantenpunkte',
      min: -8, max: 8, step: 0.05, value: state.eps, digits: 2,
      hint: 'Nur in den Verfahren 2 und 3; Plugge Gl. (4).',
      onInput: on('eps')
    });
    sliders.tau = MV.ui.slider(box, {
      label: 'τ_m — Wartezeit der Zeitbereichs-Auslese',
      min: 0, max: PI, step: PI / 400, value: state.tau, digits: 3,
      hint: 'Nicht die Messzeit τ_meas aus der Arbeit. Schrittweite π/400, damit ' +
            'π/4 und π/2 genau zu treffen sind.',
      onInput: on('tau')
    });
    sliders.gamma = MV.ui.slider(box, {
      label: 'Γ_tot — Verbreiterung des Doppelpunkts',
      min: 0, max: 1, step: 0.01, value: state.gamma, digits: 2,
      hint: 'Nur Verfahren 3; Plugge Anh. B.3 rechnet mit Γ_tot = 0.',
      onInput: on('gamma')
    });

    var row = MV.ui.el('div', { class: 'controls__wide' });
    MV.ui.buttonRow(row, {
      label: 'Voreinstellungen:',
      buttons: [
        { text: 'Vorgaben (Plugge Fig. 2)', onClick: function () {
            state.t0 = 5; state.t1 = 1; state.phi = 0; state.eps = 0;
            state.tau = PI / 4; state.gamma = 0;
            ['t0', 't1', 'phi', 'eps', 'tau', 'gamma'].forEach(function (k) {
              sliders[k].set(state[k], true);
            });
            refreshAll();
          } },
        { text: 'ohne Referenzarm (|t₀| = 0)', onClick: function () {
            state.t0 = 0; sliders.t0.set(0, true); refreshAll();
          } },
        { text: 'blinder Punkt (φ = π/2)', onClick: function () {
            state.phi = PI / 2; sliders.phi.set(PI / 2, true); refreshAll();
          } },
        { text: 'Zeitbereich τ_m = π/2', onClick: function () {
            state.tau = PI / 2; sliders.tau.set(PI / 2, true); refreshAll();
          } }
      ]
    });
    box.appendChild(row);
    mount.appendChild(box);
  }

  function buildReadout(mount) {
    readout = MV.ui.readout(mount, {
      items: [
        { key: 'tplus',  label: '|t₊|' },
        { key: 'tminus', label: '|t₋|' },
        { key: 'omega',  label: 'ω₊ / ω₋' },
        { key: 'split',  label: '|t₊|² − |t₋|² = 4|t₀||t₁| cos φ' },
        { key: 'dg',     label: 'δG_max = 4|t₀t₁|' }
      ],
      note: 'Alle drei Verfahren lesen dieselbe Größe: |t_z| mit t_z = t₀ + t₁z nach ' +
            MV.eqn('tetron.readoutH') + ' und ' + PLUGGE + ', Gl. (2). Verschwindet ' +
            '|t₊|² − |t₋|², unterscheidet keines der drei die Parität.'
    });
  }

  /* ---- Panel C: drei Teilgrafiken nebeneinander ---- */
  function buildTriptych(body) {
    var row = MV.ui.el('div', { class: 'triptych' });
    function cell(title, sub, two) {
      var c = MV.ui.el('div', { class: 'triptych__cell' });
      c.appendChild(MV.ui.el('h4', { class: 'triptych__title', html: title }));
      c.appendChild(MV.ui.el('p', { class: 'triptych__sub', html: sub }));
      var cv1 = MV.ui.el('canvas', { class: 'panel__canvas' });
      c.appendChild(cv1);
      var cv2 = null;
      if (two) { cv2 = MV.ui.el('canvas', { class: 'panel__canvas' }); c.appendChild(cv2); }
      var lg = MV.ui.el('div', { class: 'legend-inline' });
      c.appendChild(lg);
      row.appendChild(c);
      return { c1: cv1, c2: cv2, legend: lg };
    }
    var a = cell('C1 — Leitfähigkeit',
      'Arbeit ' + MV.eq('tetron.conductance') + ', ' + PLUGGE + ' Gl. (3)', false);
    var b = cell('C2 — Zeitbereich',
      PLUGGE + ', Gl. (4) und Anh. B', false);
    var c = cell('C3 — Frequenzbereich',
      PLUGGE + ', Gl. (5) und (6)', true);
    body.appendChild(row);

    plotC1 = new MV.Plot(a.c1, { aspect: 0.82, margin: { l: 52, r: 12, t: 12, b: 40 } });
    plotC1.setLegendHost(a.legend);
    plotC2 = new MV.Plot(b.c1, { aspect: 0.82, margin: { l: 52, r: 12, t: 12, b: 40 } });
    plotC2.setLegendHost(b.legend);
    plotC3 = {
      abs: new MV.Plot(c.c1, { aspect: 0.40, margin: { l: 52, r: 12, t: 10, b: 34 } }),
      phase: new MV.Plot(c.c2, { aspect: 0.40, margin: { l: 52, r: 12, t: 10, b: 34 } })
    };
    plotC3.abs.setLegendHost(c.legend);
  }

  function comparisonTable() {
    return '<table class="numtable compare"><thead><tr>' +
      '<th scope="col">Verfahren</th><th scope="col">was von z abhängt</th>' +
      '<th scope="col">wodurch begrenzt</th></tr></thead><tbody>' +
      '<tr><th scope="row">C1 Leitfähigkeit</th><td>|t_z|², über den Kreuzterm ' +
      '2|t₀||t₁| z cos φ</td><td class="compare__note">Nach ' + MV.secOfWork('5.5') + ' die ' +
      'langsamste: Kotunnelleitfähigkeiten sind klein, und die Phasenkohärenz im ' +
      'Referenzarm verlangt kleine Spannungen, sodass die Messzeit durch die nötige ' +
      'Datenakkumulation bestimmt wird.</td></tr>' +
      '<tr><th scope="row">C2 Zeitbereich</th><td>ω_z = √(ε² + |t_z|²), die ' +
      'Rabi-Frequenz</td><td class="compare__note">Schneller als C1. Nach ' + MV.sec('5.5') +
      ' durch Ladungsdephasierung begrenzt — hier nicht modelliert, Plugge beschreibt ' +
      'sie nur qualitativ.</td></tr>' +
      '<tr><th scope="row">C3 Frequenzbereich</th><td>die Lage von 2ω_z, also die ' +
      'dispersive Verschiebung</td><td class="compare__note">Nach ' + MV.sec('5.5') +
      ' kann die Dephasierung durch längere Integration kompensiert werden. Das ist ' +
      'das Verfahren, das im Pb-Experiment als Messung der Quantenkapazität zum ' +
      'Einsatz kommt.</td></tr>' +
      '</tbody></table>' +
      '<p class="note"><b>Güte der Primitive</b> nach ' + MV.secOfWork('5.5') + ': die Messzeit ' +
      'muss klein gegen die Paritätslebensdauer bleiben (τ_meas ≪ T_P), und die ' +
      'während einer Operation akkumulierte dynamische Phase klein gegen eins ' +
      '(δE·τ_meas/ℏ ≪ 1). Beide Bedingungen sind Größenvergleiche und im Reiter ' +
      '<a href="#energieskalen">Energieskalen</a> durchgerechnet.</p>';
  }

  /* ====================================================================
     Info-Panels — acht Punkte nach CLAUDE.md §11
     ==================================================================== */

  function paramLine() {
    var tp = tz(+1), tm = tz(-1);
    return '<ul>' +
      '<li>|t₀| = <span class="num">' + MV.fmt(state.t0, 2) + '</span>, ' +
      '|t₁| = <span class="num">' + MV.fmt(state.t1, 2) + '</span>, ' +
      'φ = <span class="num">' + MV.fmt(state.phi, 2) + '</span>, ' +
      'ε = <span class="num">' + MV.fmt(state.eps, 2) + '</span> (willkürliche Einheiten)</li>' +
      '<li>|t₊| = <span class="num">' + MV.fmt(tp.abs, 3) + '</span>, ' +
      '|t₋| = <span class="num">' + MV.fmt(tm.abs, 3) + '</span>, ' +
      'ω₊ = <span class="num">' + MV.fmt(omegaZ(+1), 3) + '</span>, ' +
      'ω₋ = <span class="num">' + MV.fmt(omegaZ(-1), 3) + '</span></li>' +
      '<li>τ_m = <span class="num">' + MV.fmt(state.tau, 2) + '</span>, ' +
      'Γ_tot = <span class="num">' + MV.fmt(state.gamma, 2) + '</span>; fest: λ = 2, ' +
      'ω₀ = 10, κ_in = κ_out = 0,2</li>' +
      '<li>verdrahtetes Paar: γ' + sub(state.pick[0]) + 'γ' + sub(state.pick[1]) +
      ' → ' + AXES[state.axis].symbol + '</li></ul>';
  }

  var ALGEBRA_NOTE =
    '<li><b>Darstellung.</b> Die vier Majorana-Operatoren sind als 4×4-Matrizen ' +
    'ausgeschrieben: γ₁ = σ_x⊗1, γ₂ = σ_y⊗1, γ₃ = σ_z⊗σ_x, γ₄ = σ_z⊗σ_y. Beim Laden ' +
    'werden alle sechzehn Antikommutatoren {γ_i, γ_j} − 2δ_ij gerechnet; größte ' +
    'Abweichung <b>' + fmtDev(check.worstAnticommutator) + '</b> bei Toleranz ' +
    fmtDev(TOL) + '. Ebenso ohne Abweichung: ' + MV.eq('tetron.parity') + ', die Spur des ' +
    'Projektors auf den Rechenraum, die Quadrate der drei Achsen sowie ' +
    MV.eq('tetron.cyclic') + ' und ' + MV.eq('tetron.decomposition') + '.</li>' +
    '<li><b>Qubit-Basis im Sektor.</b> Der Sektor P̂ = +1 ist zweidimensional; welche ' +
    'zwei Zustände darin |0⟩ und |1⟩ heißen, ist eine Basiswahl. Gewählt ist die ' +
    'Basis, in der ' + MV.eq('tetron.pauli') + ' wörtlich gilt. Die angezeigten Matrizen ' +
    'werden bei jedem Klick neu gerechnet.</li>';

  var INFO_A = {
    what:
      '<p>Aufsicht auf das Bauelement und die Frage, <b>was</b> gemessen wird. Zwei ' +
      'der vier Nullmoden werden angeklickt; daraus folgt die Pauli-Achse nach ' +
      MV.eqn('tetron.pauli') + ', das Bilineare in beiden Formen nach ' + MV.eq('tetron.decomposition') +
      ' und die geometrische Lesart.</p>' +
      '<p>Die Geometrie ist die Drei-Punkte-Anordnung aus ' + PLUGGE + ', Fig. 3(a): ' +
      'oberhalb der Box ein schwebender TS-Zusatzdraht, der als ein einziges ' +
      'ausgedehntes fermionisches Niveau die phasenkohärente Verbindung zum fernen ' +
      'Ende der Box herstellt, und drei Quantenpunkte rechts.</p>',
    model:
      '<p>Box nach ' + MV.fig('tetron.layoutFig') + ' der Arbeit: zwei parallele TS-Drähte, eine ' +
      's-Wellen-Brücke, eine schwebende Insel, vier Nullmoden.</p>' +
      '<span class="eq">P̂ = −γ₁γ₂γ₃γ₄ ≐ +1                       ' + MV.eq('tetron.parity') + '\n' +
      'x̂ = iγ₁γ₂ ,  ŷ = iγ₃γ₁ ,  ẑ = iγ₂γ₃       ' + MV.eq('tetron.pauli') + '\n' +
      'x̂ = iγ₁γ₂ = iγ₃γ₄   längs\n' +
      'ŷ = iγ₃γ₁ = iγ₂γ₄   diagonal\n' +
      'ẑ = iγ₂γ₃ = iγ₁γ₄   quer                  ' + MV.eq('tetron.decomposition') + '</span>' +
      '<p>Zuordnung der Punktpaare nach ' + PLUGGE + ', Fig. 3(b): Punkte 1 &amp; 2 → x̂, ' +
      'Punkte 2 &amp; 3 → ẑ, Punkte 3 &amp; 1 → ŷ.</p>',
    computed:
      '<p>Für das gewählte Paar wird iγ_iγ_j in der 4×4-Darstellung gebildet, auf den ' +
      'Sektor P̂ = +1 projiziert und als 2×2-Matrix angezeigt, daneben die zugehörige ' +
      'Pauli-Matrix und die größte Abweichung. Ebenso die zweite Schleife derselben ' +
      'Achse und das Quadrat.</p>',
    params: paramLine,
    numerics: '<ul>' + ALGEBRA_NOTE +
      '<li><b>Schematisch vereinfacht.</b> Die genaue Führung der Wege in Fig. 3(a) ' +
      'ist der Abbildung nicht eindeutig zu entnehmen. Übernommen sind: die drei ' +
      'gestapelten Drähte, die Lage der vier Nullmoden, die drei Quantenpunkte rechts ' +
      'und ihre Zuordnung. Ergänzt sind die Lage der beiden Referenzarme und der ' +
      'Verlauf der Verbindungslinien.</li></ul>',
    convention:
      '<p>' + UNITS_NOTE + '</p><p>' + PHI_NOTE + '</p><p>' + PARITY_NOTE + '</p>' +
      '<p><b>Konventionshinweis.</b> Microsoft Quantum nennt die Parität eines ' +
      'Drahtes — hier x̂ = iγ₁γ₂ — Pauli-Z (' + MV.sec('5.3') + ' der Arbeit). ' +
      'Für ŷ und ẑ gibt die Arbeit keine Übersetzung an.</p>',
    reference:
      '<p>' + MV.secOfWork('5.2', '5.3') + ', ' + MV.fig('tetron.layoutFig') + '. Gleichungen ' +
      MV.eq('tetron.parity') + ', ' + MV.eq('tetron.pauli') + ', ' + MV.eq('tetron.cyclic') + ', ' + MV.eq('tetron.decomposition') + '. ' +
      'Geometrie der drei Achsen: ' + PLUGGE + ', Fig. 3 — nicht in der Arbeit.</p>',
    reading:
      '<p><b>Klick nacheinander verschiedene Paare an.</b> Nach sechs Versuchen hast ' +
      'du nur drei Achsen gefunden. Beachte, welche Paare verdrahtet sind und welche ' +
      'nicht: alles, was γ₄ enthält, ist nach ' + MV.eq('tetron.decomposition') + ' gleichwertig, ' +
      'aber in dieser Geometrie nicht kontaktiert.</p>'
  };

  var INFO_B = {
    what:
      '<p>Die Primitive, auf der alles beruht: die Gesamtamplitude t_z = t₀ + t₁z als ' +
      'Zeiger in der komplexen Ebene. t₀ liegt auf der reellen Achse; von seiner ' +
      'Spitze aus zeigt ±t₁e^{iφ} zu den beiden Summen t₊ und t₋.</p>',
    model:
      '<span class="eq">t_z = t₀ + t₁ z ,   z = ±1 ,   t₀ reell,  t₁ = |t₁| e^{iφ}\n' +
      '|t_z|² = |t₀|² + |t₁|² + 2|t₀||t₁| z cos φ</span>' +
      '<p>' + PLUGGE + ', Gl. (2); in der Arbeit ' + MV.eqn('tetron.readoutH') + '.</p>',
    computed:
      '<p>|t₊|, |t₋| und die Differenz |t₊|² − |t₋|² = 4|t₀||t₁| cos φ. Diese Differenz ' +
      'ist der einzige Term, der den Qubit-Zustand trägt.</p>',
    params: paramLine,
    numerics:
      '<ul><li>Reine Zeigerarithmetik, keine Näherung.</li>' +
      '<li>Für t₀ = 0 oder φ = π/2 ist |t₊| = |t₋| exakt; dann unterscheidet <b>keines</b> ' +
      'der drei Verfahren die Parität.</li></ul>',
    convention: '<p>' + UNITS_NOTE + '</p><p>' + PHI_NOTE + '</p>',
    reference:
      '<p>' + MV.secOfWork('5.3') + ', ' + MV.eqn('tetron.readoutH') + '. Formel aus ' + PLUGGE + ', Gl. (2).</p>',
    reading:
      '<p><b>Zieh |t₀| auf null:</b> Beide Zeiger werden gleich lang. <b>Stell dann φ ' +
      'auf π/2</b> — auch mit Referenzarm sind sie gleich lang. Die Statuszeile nennt ' +
      'jeweils den Grund.</p>'
  };

  var INFO_C = {
    what:
      '<p>Drei Ausleseverfahren nebeneinander. Sie sehen verschieden aus, lesen aber ' +
      'dieselbe Größe: |t_z|. Verschwindet der z-abhängige Term, fällt der Unterschied ' +
      'zwischen den Paritäten in allen dreien zugleich weg.</p>',
    model:
      '<span class="eq">C1   G_z = (e²/h) ν₁ν₂ |t_z|² ,   δG_max ~ 4|t₀t₁|\n' +
      'C2   P_z(2; τ_m) = cos²(ω_z τ_m) + (ε²/ω_z²) sin²(ω_z τ_m)\n' +
      '     ω_z = √(ε² + |t_z|²)\n' +
      'C3   g_z = −λ|t_z|/(2ω_z)\n' +
      '     χ_z = g_z² / [ iΓ_tot − (2ω_z − ω) ]\n' +
      '     A_ω = −i√(κ_in κ_out) / [ −(i/2)(κ_in+κ_out) + (ω₀−ω) + χ_z ]</span>' +
      '<p>C1 steht als ' + MV.eqn('tetron.conductance') + ' auch in der Arbeit; C2 und C3 stammen ' +
      'aus ' + PLUGGE + ', Gl. (4) bis (6) und Anhang B. Die Arbeit beschreibt die ' +
      'Verfahren in ' + MV.sec('5.5') + ' nur qualitativ.</p>',
    computed:
      '<p>C1: G_z(φ) über zwei Perioden, beide Paritäten. C2: Rückkehrwahrscheinlichkeit ' +
      'über der Wartezeit. C3: Transmission |A_ω|² und Phasenverschiebung Δφ(ω) = ' +
      '−arg A_ω über der Resonatorfrequenz.</p>',
    params: paramLine,
    numerics:
      '<ul><li>' + N_CURVE + ' Stützstellen je Kurve.</li>' +
      '<li><b>Keine Dämpfung in C2.</b> Plugge beschreibt die Ladungsdephasierung nur ' +
      'qualitativ; sie ist hier nicht modelliert. Die Kurven laufen deshalb ungedämpft ' +
      'weiter, als sie es im Experiment täten.</li>' +
      '<li>Grenzfälle: ω_z = 0 ergibt in C2 P = 1 und in C3 den nackten Resonator ' +
      '(g_z = 0). Bei ω = 2ω_z und Γ_tot = 0 geht A_ω gegen null; beides wird ' +
      'abgefangen, es entstehen keine NaN.</li>' +
      '<li>Kontrolle gegen ' + PLUGGE + ', Fig. 2(e, f): mit den Vorgaben liegen die ' +
      'Transmissionsmaxima bei ω = 9,586 und 10,414 mit |A|² = 1, und Δφ(ω₀) beträgt ' +
      '−68,2° für z = +1 und +68,2° für z = −1.</li></ul>',
    convention:
      '<p>' + UNITS_NOTE + '</p><p>' + PHI_NOTE + '</p>' +
      '<p>C1 ist 2π-periodisch in φ, die beiden Kurven sind um π verschoben. Ein ' +
      'Flussquant wird hier <b>nicht</b> angegeben — φ ist dimensionslos.</p>',
    reference:
      '<p>' + MV.eqn('tetron.conductance') + ' steht in ' + MV.sec('5.3') + ' der ' +
      'Arbeit; der Vergleich der drei Verfahren in ' + MV.sec('5.5') + '. ' +
      'Verfahren 2 und 3 aus ' + PLUGGE + ', Gl. (4) bis (6), Anhang B.</p>',
    reading:
      '<p><b>Zieh |t₀| auf null und beobachte alle drei Teilgrafiken:</b> Der ' +
      'Unterschied zwischen den Paritäten verschwindet überall zugleich, weil alle ' +
      'dieselbe Größe |t_z| lesen.</p>' +
      '<p><b>Zieh dann τ_m durch den ganzen Bereich und sieh C2 zu.</b> Bei ε = 0 ' +
      'ist P_z = cos²(ω_z τ_m); voller Kontrast herrscht genau dann, wenn ω₊τ_m und ' +
      'ω₋τ_m um ein <b>ungerades Vielfaches von π/2</b> auseinanderliegen — dann ' +
      'steht die eine Kurve auf 0, während die andere auf 1 steht. Mit den ' +
      'Vorgabewerten ist (ω₊ − ω₋) τ_m = ' + MV.fmt(2 * (PI / 4), 3) + ' = π/2, ' +
      'und das ist der erste solche Punkt. Später gäbe es weitere — man misst aber ' +
      'so früh wie möglich, weil jede Wartezeit Gelegenheit zur Dephasierung ist.</p>'
  };

  /* ====================================================================
     Öffentliche Schnittstelle
     ==================================================================== */

  MV.auslese = {
    init: function (mount) {
      MV.auslese.check = check;
      if (!check.passed) {
        var warn = MV.ui.el('div', { class: 'readout',
          html: '<b>Warnung:</b> Die Prüfung der Pauli-Algebra ist fehlgeschlagen ' +
                '(größte Abweichung ' + fmtDev(check.worst) + ').' });
        warn.style.borderColor = MV.palette.critical;
        warn.style.color = MV.palette.critical;
        mount.appendChild(warn);
      }

      buildControls(mount);
      buildReadout(mount);

      var panels = MV.ui.el('div', { class: 'panels' });
      mount.appendChild(panels);

      panelA = MV.ui.panel(panels, {
        title: 'Panel A — Was gemessen wird',
        quote: 'auslese.A',
        lead: '<b>Klick zwei der vier Majoranas an.</b> Die gewählte Schleife wird ' +
              'verdrahtet, die zweite Schleife derselben Achse leuchtet schwächer mit. ' +
              'Die Rechnung darunter folgt der Auswahl.',
        wide: true, aspect: 0.40, controls: true,
        margin: { l: 14, r: 14, t: 10, b: 10 },
        info: INFO_A
      });

      wiringHost = MV.ui.el('div', { class: 'callout callout--soft', hidden: 'hidden' });
      panelA.controls.appendChild(wiringHost);
      pairHost = MV.ui.el('div', { class: 'pairhost' });
      panelA.controls.appendChild(pairHost);
      sectorHost = MV.ui.el('div', { class: 'sectorhost' });
      panelA.controls.appendChild(sectorHost);

      var cycRow = MV.ui.el('div', { class: 'cyclicrow' });
      panelA.controls.appendChild(cycRow);
      MV.ui.segmented(cycRow, {
        label: 'Produkt zweier Achsen nach ' + MV.eq('tetron.cyclic') + ' — erste:',
        options: ['x̂', 'ŷ', 'ẑ'], value: state.cycA,
        onChange: function (i) { state.cycA = i; renderCyclic(); }
      });
      MV.ui.segmented(cycRow, {
        label: 'zweite:', options: ['x̂', 'ŷ', 'ẑ'], value: state.cycB,
        onChange: function (i) { state.cycB = i; renderCyclic(); }
      });
      cyclicHost = MV.ui.el('div', { class: 'cyclichost' });
      panelA.controls.appendChild(cyclicHost);

      attachClicks();

      panelB = MV.ui.panel(panels, {
        title: 'Panel B — Die Primitive',
        quote: 'auslese.B',
        lead: 'Die Gesamtamplitude t_z = t₀ + t₁z als Zeiger; <b>z ist der ' +
              'Eigenwert der in Panel A gewählten Achse</b>. Nur die Differenz ' +
              '|t₊|² − |t₋|² trägt den Qubit-Zustand.',
        aspect: 0.72,
        info: INFO_B
      });

      panelC = MV.ui.panel(panels, {
        title: 'Panel C — Drei Verfahren, eine Größe',
        quote: 'auslese.C',
        lead: 'Alle drei lesen |t_z|. Sie unterscheiden sich darin, <i>wie</i> sie es ' +
              'tun — und damit in ihrer Geschwindigkeit, nicht in der gemessenen Größe.',
        wide: true, canvas: false, controls: true,
        info: INFO_C
      });
      buildTriptych(panelC.body);
      panelC.controls.appendChild(MV.ui.el('div', { html: comparisonTable() }));

      /* Schnittstelle für die Vorführungen (demo.js). Ein Paar wird über
         selectPair gewählt — denselben Weg nimmt ein Klick auf die Chips
         und auf die Majoranas im Bild. */
      MV.demo.provide('auslese', MV.demo.handle({
        state: state,
        sliders: sliders,
        refresh: refreshAll,
        extra: {
          allPairs: function () {
            return allPairs().map(function (e) { return e.pair; });
          },
          selectPair: selectPair,
          forgetVisited: function () { state.visited = {}; },
          counts: function () {
            var seen = {}, n = 0, k, pr;
            for (k in state.visited) {
              if (!Object.prototype.hasOwnProperty.call(state.visited, k)) continue;
              n++;
              pr = k.split('-');
              seen[axisOfPair(+pr[0], +pr[1]).idx] = true;
            }
            return { pairs: n, axes: Object.keys(seen).length };
          },
          read: function () {
            var tp = tz(+1), tm = tz(-1);
            var pk = peakOf(+1), mk = peakOf(-1);
            return {
              tplus: tp.abs,
              tminus: tm.abs,
              split: splitAbs2(),
              blind: indistinguishable(),
              why: whyBlind(),
              deltaGmax: deltaGmax(),
              contrastC2: Math.abs(pReturn(+1, state.tau) - pReturn(-1, state.tau)),
              peakGapC3: Math.abs(pk.w - mk.w)
            };
          },
          showPanelC: function () {
            if (panelC && panelC.root.scrollIntoView) {
              panelC.root.scrollIntoView({ block: 'center' });
            }
          }
        }
      }));
    },

    draw: function () { refreshAll(); },

    stop: function () { /* nichts läuft im Hintergrund */ }
  };

}(MV));
