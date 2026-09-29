/* ------------------------------------------------------------------------
   gatter.js — Reiter #gatter: was mit Clifford-Gattern erreichbar ist und
   was nicht. Abschnitte 5.4 und 5.5 der Arbeit.

   Panel A  Die neun Gatter auf der Bloch-Kugel. Die sechs Eigenzustände von
            x̂, ŷ und ẑ sind die Ecken eines Oktaeders; jedes Clifford-Gatter
            bildet Ecken auf Ecken ab, das T-Gatter nicht. Das ist die
            Anschauung zu Gottesman-Knill, kein Beweis.
   Panel B  Dieselben Gatter aus Messungen: vier Protokolle, gerechnet.

   Drehachse und Drehwinkel jedes Gatters werden aus seiner Matrix bestimmt,
   nicht hingeschrieben — der gezeichnete Bogen ist damit die nachgerechnete
   Wirkung des Operators und nicht eine zweite, unabhängige Behauptung.

   Alles dimensionslos: ein Zustandsraum hat keine Einheiten. Es gibt in
   diesem Reiter keine Zeitentwicklung und keine Rate; ein Gatter wird
   angewandt, nicht gefahren.

   θ und φ sind hier die Winkel des Startzustands auf der Bloch-Kugel —
   weder die supraleitende Phase aus Kapitel 3 noch der eingeschlossene
   Fluss aus Kapitel 5.
   ------------------------------------------------------------------------ */

var MV = window.MV || (window.MV = {});

(function (MV) {
  'use strict';

  var PI = Math.PI;
  var DEG = 180 / PI;
  var TOL = 1e-12;
  var S2 = 1 / Math.SQRT2;

  var NC = 'Nielsen und Chuang, Quantum Computation and Quantum Information, Kap. 10';

  /* ====================================================================
     1 — Komplexe 2×2-Matrizen

     Eine Matrix ist ein Feld aus vier Paaren [Re, Im], zeilenweise:
     [u₀₀, u₀₁, u₁₀, u₁₁].
     ==================================================================== */

  function cm(a, b) { return [a[0] * b[0] - a[1] * b[1], a[0] * b[1] + a[1] * b[0]]; }
  function ca(a, b) { return [a[0] + b[0], a[1] + b[1]]; }
  function cs(a, k) { return [a[0] * k, a[1] * k]; }
  function cc(a) { return [a[0], -a[1]]; }

  function mmul(A, B) {
    var out = [], r, c, k, s;
    for (r = 0; r < 2; r++) {
      for (c = 0; c < 2; c++) {
        s = [0, 0];
        for (k = 0; k < 2; k++) s = ca(s, cm(A[r * 2 + k], B[k * 2 + c]));
        out.push(s);
      }
    }
    return out;
  }

  function mdag(A) { return [cc(A[0]), cc(A[2]), cc(A[1]), cc(A[3])]; }

  function madd(A, B) { return A.map(function (e, i) { return ca(e, B[i]); }); }
  function mscale(A, k) { return A.map(function (e) { return cs(e, k); }); }

  function mmaxAbs(A) {
    return A.reduce(function (m, e) { return Math.max(m, Math.hypot(e[0], e[1])); }, 0);
  }
  function mdiff(A, B) {
    return mmaxAbs(A.map(function (e, i) { return [e[0] - B[i][0], e[1] - B[i][1]]; }));
  }

  var ID = [[1, 0], [0, 0], [0, 0], [1, 0]];
  var SX = [[0, 0], [1, 0], [1, 0], [0, 0]];
  var SY = [[0, 0], [0, -1], [0, 1], [0, 0]];
  var SZ = [[1, 0], [0, 0], [0, 0], [-1, 0]];
  var SIGMA = { x: SX, y: SY, z: SZ };

  /** exp(−i(α/2) n̂·σ̂) = cos(α/2)·1 − i sin(α/2) n̂·σ̂. */
  function rotMatrix(n, alpha) {
    var c = Math.cos(alpha / 2), s = Math.sin(alpha / 2);
    var nx = SIGMA.x, ny = SIGMA.y, nz = SIGMA.z;
    var M = madd(madd(mscale(nx, n[0]), mscale(ny, n[1])), mscale(nz, n[2]));
    /* −i · M */
    var miM = M.map(function (e) { return [e[1], -e[0]]; });
    return madd(mscale(ID, c), mscale(miM, s));
  }

  /**
   * Abweichung zweier Matrizen bis auf eine globale Phase. Ein Operator, der
   * auf Zustände wirkt, ist nur bis auf diese Phase bestimmt; ohne die
   * Anpassung meldete jeder Vergleich eine bedeutungslose Differenz.
   */
  function devUpToPhase(A, B) {
    /* Phase aus dem betragsgrößten Eintrag von B nehmen. */
    var i, best = 0, bm = -1, m;
    for (i = 0; i < 4; i++) {
      m = Math.hypot(B[i][0], B[i][1]);
      if (m > bm) { bm = m; best = i; }
    }
    if (bm < 1e-14) return mmaxAbs(A);
    var q = cm(A[best], cc(B[best]));
    var qm = Math.hypot(q[0], q[1]);
    if (qm < 1e-14) return mdiff(A, B);
    var ph = [q[0] / qm, q[1] / qm];
    return mdiff(A, B.map(function (e) { return cm(e, ph); }));
  }

  /**
   * Drehachse und Drehwinkel aus der Matrix. Die globale Phase wird über die
   * Determinante entfernt; die verbleibende SU(2)-Matrix hat die Form
   * a·1 − i(b σ_x + c σ_y + d σ_z) mit a = cos(α/2), (b,c,d) = sin(α/2) n̂.
   */
  function axisAngle(U) {
    var det = [U[0][0] * U[3][0] - U[0][1] * U[3][1] - (U[1][0] * U[2][0] - U[1][1] * U[2][1]),
               U[0][0] * U[3][1] + U[0][1] * U[3][0] - (U[1][0] * U[2][1] + U[1][1] * U[2][0])];
    var delta = Math.atan2(det[1], det[0]);
    var ph = [Math.cos(-delta / 2), Math.sin(-delta / 2)];
    var V = U.map(function (e) { return cm(e, ph); });

    var a = (V[0][0] + V[3][0]) / 2;
    var b = -(V[1][1] + V[2][1]) / 2;
    var c = (V[2][0] - V[1][0]) / 2;
    var d = (V[3][1] - V[0][1]) / 2;
    if (a < 0) { a = -a; b = -b; c = -c; d = -d; }   /* α in [0, π] wählen */

    var len = Math.hypot(b, c, d);
    var alpha = 2 * Math.atan2(len, a);
    var n = len < 1e-12 ? [0, 0, 1] : [b / len, c / len, d / len];
    return { n: n, alpha: alpha, identity: len < 1e-12 };
  }

  /* ====================================================================
     2 — Die neun Gatter

     Definitionen nach CLAUDE.md, Abschnitt 8: Ŝ_x = e^{−i(π/4)x̂},
     Ŝ_y = Ŝ_z Ŝ_x Ŝ_z†, Ĥ = (x̂ + ẑ)/√2, T-Gatter = diag(1, e^{iπ/4}).
     Ŝ_y und Ĥ werden hier aus den anderen gerechnet, nicht noch einmal
     hingeschrieben.
     ==================================================================== */

  var M_SZ_GATE = [[1, 0], [0, 0], [0, 0], [0, 1]];             /* diag(1, i)  */
  var M_SX_GATE = rotMatrix([1, 0, 0], PI / 2);                 /* e^{−iπx̂/4} */
  var M_SY_GATE = mmul(mmul(M_SZ_GATE, M_SX_GATE), mdag(M_SZ_GATE));
  var M_H = mscale(madd(SX, SZ), S2);
  var M_T = [[1, 0], [0, 0], [0, 0], [Math.cos(PI / 4), Math.sin(PI / 4)]];
  var M_TD = mdag(M_T);

  var GATES = [
    { key: 'x',   label: 'x̂',        U: SX,        clifford: true,
      note: 'Pauli-Gatter. Am Bauelement der bestätigte Elektronentransfer ' +
            'über das zugehörige Kontaktpaar oder reine klassische Buchführung.' },
    { key: 'y',   label: 'ŷ',        U: SY,        clifford: true,
      note: 'Pauli-Gatter, siehe x̂.' },
    { key: 'z',   label: 'ẑ',        U: SZ,        clifford: true,
      note: 'Pauli-Gatter. Für ẑ nennt die Arbeit den Fall ausdrücklich: ' +
            'bestätigter Transfer eines Elektrons durch die Box.' },
    { key: 'Sz',  label: 'Ŝ_z',      U: M_SZ_GATE, clifford: true,
      note: 'Derselbe Operator wie der Zopf U₂₃ aus ' + MV.figShort('ops.braidFig') +
            ' der Arbeit — hier aber aus zwei Messungen und einer klassisch ' +
            'bedingten Pauli-Korrektur erzeugt, ohne dass eine Nullmode bewegt wird.' },
    { key: 'Sx',  label: 'Ŝ_x',      U: M_SX_GATE, clifford: true,
      note: 'Analoges Protokoll mit Ancilla in y_A = −1.' },
    { key: 'Sy',  label: 'Ŝ_y',      U: M_SY_GATE, clifford: true,
      note: 'Nicht eigenständig ausgeführt, sondern als Ŝ_z Ŝ_x Ŝ_z† ' +
            'zusammengesetzt.' },
    { key: 'H',   label: 'Ĥ',        U: M_H,       clifford: true,
      note: 'Hadamard, hier als (x̂ + ẑ)/√2 definiert; bis auf eine Phase ' +
            'gleich Ŝ_z Ŝ_x Ŝ_z.' },
    { key: 'T',   label: 'T-Gatter',  U: M_T,      clifford: false,
      note: 'Kein Clifford-Gatter. Nicht aus Messungen an Majorana-Bilinearen ' +
            'erreichbar; kommt über magische Zustände und ist nicht ' +
            'topologisch geschützt.' },
    { key: 'Td',  label: 'T-Gatter†', U: M_TD,     clifford: false,
      note: 'Das Inverse des T-Gatters, mit derselben Einschränkung.' }
  ];

  var GATE_BY_KEY = {};
  GATES.forEach(function (g) {
    g.rot = axisAngle(g.U);
    GATE_BY_KEY[g.key] = g;
  });

  /* ====================================================================
     3 — Ein Qubit: Zustand und Bloch-Vektor
     ==================================================================== */

  /** |ψ⟩ = cos(θ/2)|0⟩ + e^{iφ} sin(θ/2)|1⟩, θ und φ im Bogenmaß. */
  function ketOf(theta, phi) {
    return [[Math.cos(theta / 2), 0],
            [Math.sin(theta / 2) * Math.cos(phi), Math.sin(theta / 2) * Math.sin(phi)]];
  }

  function applyU(U, psi) {
    return [ca(cm(U[0], psi[0]), cm(U[1], psi[1])),
            ca(cm(U[2], psi[0]), cm(U[3], psi[1]))];
  }

  /** r = (2Re c₀*c₁, 2Im c₀*c₁, |c₀|² − |c₁|²). */
  function blochOf(psi) {
    var p = cm(cc(psi[0]), psi[1]);
    return [2 * p[0], 2 * p[1],
            psi[0][0] * psi[0][0] + psi[0][1] * psi[0][1] -
            psi[1][0] * psi[1][0] - psi[1][1] * psi[1][1]];
  }

  /** Drehung eines Vektors um n̂ um den Winkel α, Rodrigues. */
  function rotate(r, n, alpha) {
    var c = Math.cos(alpha), s = Math.sin(alpha);
    var dot = n[0] * r[0] + n[1] * r[1] + n[2] * r[2];
    var cx = [n[1] * r[2] - n[2] * r[1], n[2] * r[0] - n[0] * r[2], n[0] * r[1] - n[1] * r[0]];
    return [r[0] * c + cx[0] * s + n[0] * dot * (1 - c),
            r[1] * c + cx[1] * s + n[1] * dot * (1 - c),
            r[2] * c + cx[2] * s + n[2] * dot * (1 - c)];
  }

  /* Die sechs Ecken des Oktaeders: die Eigenzustände von x̂, ŷ und ẑ. Die
     Basis ist die der Arbeit, CLAUDE.md Abschnitt 8: |0⟩ ist ẑ = +1. */
  var VERTICES = [
    { r: [0, 0, 1],  label: 'ẑ = +1', ket: '|0⟩',  theta: 0,      phi: 0 },
    { r: [0, 0, -1], label: 'ẑ = −1', ket: '|1⟩',  theta: PI,     phi: 0 },
    { r: [1, 0, 0],  label: 'x̂ = +1', ket: '|+⟩',  theta: PI / 2, phi: 0 },
    { r: [-1, 0, 0], label: 'x̂ = −1', ket: '|−⟩',  theta: PI / 2, phi: PI },
    { r: [0, 1, 0],  label: 'ŷ = +1', ket: '|+i⟩', theta: PI / 2, phi: PI / 2 },
    { r: [0, -1, 0], label: 'ŷ = −1', ket: '|−i⟩', theta: PI / 2, phi: -PI / 2 }
  ];

  /** Nächstgelegene Ecke und ihr Abstand. */
  function nearestVertex(r) {
    var best = null, bd = Infinity;
    VERTICES.forEach(function (v) {
      var d = Math.hypot(r[0] - v.r[0], r[1] - v.r[1], r[2] - v.r[2]);
      if (d < bd) { bd = d; best = v; }
    });
    return { v: best, d: bd };
  }

  /* ====================================================================
     4 — Register aus mehreren Qubits

     Ein Register ist ein komplexer Vektor der Länge 2ⁿ. Qubit 0 ist das
     höchstwertige Bit; das Protokoll für CNOT braucht drei Qubits, das für
     das T-Gatter zwei.
     ==================================================================== */

  function reg(n) { return { n: n, re: new Float64Array(1 << n), im: new Float64Array(1 << n) }; }
  function bitOf(i, k, n) { return (i >> (n - 1 - k)) & 1; }

  function norm2(v) {
    var s = 0, i;
    for (i = 0; i < v.re.length; i++) s += v.re[i] * v.re[i] + v.im[i] * v.im[i];
    return s;
  }

  function normalise(v) {
    var s = Math.sqrt(norm2(v)), i;
    if (s < 1e-14) return v;
    for (i = 0; i < v.re.length; i++) { v.re[i] /= s; v.im[i] /= s; }
    return v;
  }

  /** Ein Ein-Qubit-Zustand als Register. */
  function regOf(psi) {
    var v = reg(1);
    v.re[0] = psi[0][0]; v.im[0] = psi[0][1];
    v.re[1] = psi[1][0]; v.im[1] = psi[1][1];
    return v;
  }

  function tensor(a, b) {
    var out = reg(a.n + b.n), i, j, k;
    for (i = 0; i < (1 << a.n); i++) {
      for (j = 0; j < (1 << b.n); j++) {
        k = (i << b.n) | j;
        out.re[k] = a.re[i] * b.re[j] - a.im[i] * b.im[j];
        out.im[k] = a.re[i] * b.im[j] + a.im[i] * b.re[j];
      }
    }
    return out;
  }

  /** Pauli-Operator auf Qubit k des Registers. */
  function pauli(v, k, ax) {
    var n = v.n, N = 1 << n, out = reg(n), i, b, j, s;
    for (i = 0; i < N; i++) {
      b = bitOf(i, k, n);
      if (ax === 'z') {
        s = b ? -1 : 1;
        out.re[i] += s * v.re[i]; out.im[i] += s * v.im[i];
      } else {
        j = i ^ (1 << (n - 1 - k));
        if (ax === 'x') {
          out.re[j] += v.re[i]; out.im[j] += v.im[i];
        } else {
          /* ŷ: |0⟩ ↦ i|1⟩, |1⟩ ↦ −i|0⟩ */
          s = b ? -1 : 1;
          out.re[j] += -s * v.im[i]; out.im[j] += s * v.re[i];
        }
      }
    }
    return out;
  }

  /** Produkt mehrerer Pauli-Operatoren, terms = [[qubit, achse], …]. */
  function pauliString(v, terms) {
    var w = v, t;
    for (t = 0; t < terms.length; t++) w = pauli(w, terms[t][0], terms[t][1]);
    return w;
  }

  /**
   * Projektive Messung: P_b = (1 + b·Ô)/2. Das Normquadrat des projizierten
   * Vektors vor der Normierung ist die Wahrscheinlichkeit des Ergebnisses.
   */
  function project(v, terms, b) {
    var Ov = pauliString(v, terms), out = reg(v.n), i;
    for (i = 0; i < out.re.length; i++) {
      out.re[i] = 0.5 * (v.re[i] + b * Ov.re[i]);
      out.im[i] = 0.5 * (v.im[i] + b * Ov.im[i]);
    }
    return { vec: out, prob: norm2(out) };
  }

  /** Ein 2×2-Gatter auf Qubit k des Registers. */
  function gate1(v, k, U) {
    var n = v.n, N = 1 << n, out = reg(n), m = 1 << (n - 1 - k), i, j, a0, a1, r0, r1;
    for (i = 0; i < N; i++) {
      if (i & m) continue;
      j = i | m;
      a0 = [v.re[i], v.im[i]];
      a1 = [v.re[j], v.im[j]];
      r0 = ca(cm(U[0], a0), cm(U[1], a1));
      r1 = ca(cm(U[2], a0), cm(U[3], a1));
      out.re[i] = r0[0]; out.im[i] = r0[1];
      out.re[j] = r1[0]; out.im[j] = r1[1];
    }
    return out;
  }

  /**
   * Reduzierte Dichtematrix über eine Teilmenge der Qubits. Für ein mit dem
   * Rest verschränktes Qubit ist sie gemischt — genau das macht den
   * Bloch-Vektor kürzer als 1.
   */
  function reduced(v, qubits) {
    var n = v.n, N = 1 << n, d = 1 << qubits.length;
    var re = new Float64Array(d * d), im = new Float64Array(d * d);
    var sIdx = new Int32Array(N), rIdx = new Int32Array(N);
    var rest = [], q, i, j, si, ri, a, b;
    for (q = 0; q < n; q++) if (qubits.indexOf(q) < 0) rest.push(q);
    for (i = 0; i < N; i++) {
      si = 0; qubits.forEach(function (k) { si = (si << 1) | bitOf(i, k, n); });
      ri = 0; rest.forEach(function (k) { ri = (ri << 1) | bitOf(i, k, n); });
      sIdx[i] = si; rIdx[i] = ri;
    }
    for (i = 0; i < N; i++) {
      for (j = 0; j < N; j++) {
        if (rIdx[i] !== rIdx[j]) continue;
        a = sIdx[i]; b = sIdx[j];
        /* ρ_ab = Σ v_a conj(v_b) */
        re[a * d + b] += v.re[i] * v.re[j] + v.im[i] * v.im[j];
        im[a * d + b] += v.im[i] * v.re[j] - v.re[i] * v.im[j];
      }
    }
    return { d: d, re: re, im: im };
  }

  /** Bloch-Vektor aus einer 2×2-Dichtematrix. */
  function blochOfRho(rho) {
    return [2 * rho.re[1], -2 * rho.im[1], rho.re[0] - rho.re[3]];
  }

  /** Dichtematrix eines reinen Registerzustands. */
  function rhoOfPure(v) {
    var d = v.re.length, re = new Float64Array(d * d), im = new Float64Array(d * d), a, b;
    for (a = 0; a < d; a++) {
      for (b = 0; b < d; b++) {
        re[a * d + b] = v.re[a] * v.re[b] + v.im[a] * v.im[b];
        im[a * d + b] = v.im[a] * v.re[b] - v.re[a] * v.im[b];
      }
    }
    return { d: d, re: re, im: im };
  }

  /**
   * Abstand zweier Dichtematrizen, größter Eintragsbetrag. Der Vergleich
   * läuft bewusst über ρ und nicht über den Zustandsvektor: ρ ist gegen die
   * globale Phase unempfindlich, sodass keine Phasenanpassung nötig ist,
   * die einen echten Unterschied verdecken könnte.
   */
  function rhoDev(A, B) {
    var m = 0, i;
    for (i = 0; i < A.re.length; i++) {
      m = Math.max(m, Math.hypot(A.re[i] - B.re[i], A.im[i] - B.im[i]));
    }
    return m;
  }

  /* ====================================================================
     5 — Die vier Protokolle

     Wortlaut aus CLAUDE.md, Abschnitt 8; die ersten drei stehen so auch in
     Abschnitt 5.5 der Arbeit. Das vierte ist die Standardkonstruktion über
     einen magischen Zustand und wird in der Arbeit nicht ausgeführt.
     ==================================================================== */

  var KET_YP = regOf([[S2, 0], [0, S2]]);                        /* |+i⟩        */
  var KET_YM = regOf([[S2, 0], [0, -S2]]);                       /* |−i⟩        */
  var KET_0 = regOf([[1, 0], [0, 0]]);
  var KET_1 = regOf([[0, 0], [1, 0]]);
  var KET_M = regOf([[S2, 0], [S2 * Math.cos(PI / 4), S2 * Math.sin(PI / 4)]]);

  /** Ĉ aus CLAUDE.md §8: ½(1+ẑ)_C ⊗ 1_T + ½(1−ẑ)_C ⊗ x̂_T. */
  function cnotOn(v, ctrl, targ) {
    var pPlus = project(v, [[ctrl, 'z']], +1).vec;
    var pMinus = project(v, [[ctrl, 'z']], -1).vec;
    var flipped = pauli(pMinus, targ, 'x');
    var out = reg(v.n), i;
    for (i = 0; i < out.re.length; i++) {
      out.re[i] = pPlus.re[i] + flipped.re[i];
      out.im[i] = pPlus.im[i] + flipped.im[i];
    }
    return out;
  }

  /**
   * Wendet eine Liste von Korrekturoperatoren nacheinander an und schreibt
   * nach jedem den Registerzustand mit. Der Endzustand ist derselbe wie
   * vorher — die Zwischenstufen sind neu und nur für die Wiedergabe da.
   */
  function applyOps(v, list) {
    var ops = [], w = v;
    list.forEach(function (o) {
      w = gate1(w, o.qubit, o.U);
      ops.push({ name: o.name, qubit: o.qubit, U: o.U, after: w });
    });
    return { v: w, ops: ops };
  }

  var PROTOCOLS = [
    {
      key: 'Sz',
      name: 'Ŝ_z',
      title: 'Ŝ_z — Drehung um π/2 um die z-Achse',
      roles: ['C', 'A'],
      out: [0],
      outName: 'C',
      inQubits: 1,
      source: MV.secOfWork('5.5') + ', ' + MV.eq('ops.sz') + '; Plugge et al., Fig. 5(b)',
      prepText: 'Präpariere die Ancilla A im Eigenzustand y_A = +1.',
      prepShort: '|+i⟩_A',
      prep: function (psi) { return tensor(regOf(psi), KET_YP); },
      steps: [
        { text: 'Miss die gemeinsame Parität ẑ_C ẑ_A.', terms: [[0, 'z'], [1, 'z']], sym: 'b₁' },
        { text: 'Lies x̂_A aus.', terms: [[1, 'x']], sym: 'b₂' }
      ],
      corrText: 'Wende ẑ_C an, falls b₁b₂ = −1.',
      noCorrCond: 'b₁b₂ = +1',
      correct: function (v, r) {
        return applyOps(v, (r[0] * r[1] === -1) ? [{ name: 'ẑ_C', qubit: 0, U: SZ }] : []);
      },
      corrJoin: function (names) { return names[0]; },
      target: function (psi) { return regOf(applyU(M_SZ_GATE, psi)); }
    },
    {
      key: 'Sx',
      name: 'Ŝ_x',
      title: 'Ŝ_x — Drehung um π/2 um die x-Achse',
      roles: ['C', 'A'],
      out: [0],
      outName: 'C',
      inQubits: 1,
      source: MV.secOfWork('5.5') + '; Plugge et al., Fig. 5(c)',
      prepText: 'Präpariere die Ancilla A im Eigenzustand y_A = −1.',
      prepShort: '|−i⟩_A',
      prep: function (psi) { return tensor(regOf(psi), KET_YM); },
      steps: [
        { text: 'Miss die gemeinsame Parität x̂_C x̂_A.', terms: [[0, 'x'], [1, 'x']], sym: 'c₁' },
        { text: 'Lies ẑ_A aus.', terms: [[1, 'z']], sym: 'c₂' }
      ],
      corrText: 'Wende x̂_C an, falls c₁c₂ = −1.',
      noCorrCond: 'c₁c₂ = +1',
      correct: function (v, r) {
        return applyOps(v, (r[0] * r[1] === -1) ? [{ name: 'x̂_C', qubit: 0, U: SX }] : []);
      },
      corrJoin: function (names) { return names[0]; },
      target: function (psi) { return regOf(applyU(M_SX_GATE, psi)); }
    },
    {
      key: 'CNOT',
      name: 'CNOT',
      title: 'CNOT — Kontrolle C, Ziel T',
      roles: ['C', 'A', 'T'],
      out: [0, 2],
      outName: 'C und T',
      inQubits: 2,
      source: MV.secOfWork('5.5') + '; Plugge et al., Fig. 5(a)',
      prepText: 'Präpariere die Ancilla A in |0⟩.',
      prepShort: '|0⟩_A',
      prep: function (psi, t0) {
        return tensor(tensor(regOf(psi), KET_0), t0 ? KET_1 : KET_0);
      },
      steps: [
        { text: 'Miss die gemeinsame Parität x̂_A x̂_T.', terms: [[1, 'x'], [2, 'x']], sym: 'a₁' },
        { text: 'Miss die gemeinsame Parität ẑ_C ẑ_A.', terms: [[0, 'z'], [1, 'z']], sym: 'a₂' },
        { text: 'Lies x̂_A aus.', terms: [[1, 'x']], sym: 'a₃' }
      ],
      corrText: 'Wende ẑ_C an, falls a₁a₃ = −1, und x̂_T, falls a₂ = −1.',
      noCorrCond: 'a₁a₃ = +1 und a₂ = +1',
      correct: function (v, r) {
        var list = [];
        if (r[0] * r[2] === -1) list.push({ name: 'ẑ_C', qubit: 0, U: SZ });
        if (r[1] === -1) list.push({ name: 'x̂_T', qubit: 2, U: SX });
        return applyOps(v, list);
      },
      corrJoin: function (names) { return names.join(' und '); },
      target: function (psi, t0) {
        return cnotOn(tensor(regOf(psi), t0 ? KET_1 : KET_0), 0, 1);
      }
    },
    {
      key: 'T',
      name: 'T-Gatter',
      title: 'T-Gatter über einen magischen Zustand',
      roles: ['D', 'M'],
      out: [1],
      outName: 'M',
      inQubits: 1,
      magic: true,
      source: 'Standardkonstruktion, ' + NC + '. In der Arbeit nicht ausgeführt; ' +
              MV.sec('5.5') + ' nennt den Weg über magische Zustände und hält fest, ' +
              'dass er nicht topologisch geschützt ist.',
      prepText: 'Halte den magischen Zustand |M⟩ = (|0⟩ + e^{iπ/4}|1⟩)/√2 auf Qubit M bereit.',
      prepShort: '|M⟩_M',
      prep: function (psi) { return tensor(regOf(psi), KET_M); },
      steps: [
        { text: 'Miss die gemeinsame Parität ẑ_D ẑ_M.', terms: [[0, 'z'], [1, 'z']], sym: 's' },
        { text: 'Lies x̂_D aus.', terms: [[0, 'x']], sym: 'x' }
      ],
      corrText: 'Auf M: x̂ falls s = −1, danach ẑ falls x = −1, danach Ŝ_z falls s = −1.',
      noCorrCond: 's = +1 und x = +1',
      correct: function (v, r) {
        var list = [];
        if (r[0] === -1) list.push({ name: 'x̂', qubit: 1, U: SX });
        if (r[1] === -1) list.push({ name: 'ẑ', qubit: 1, U: SZ });
        if (r[0] === -1) list.push({ name: 'Ŝ_z', qubit: 1, U: M_SZ_GATE });
        return applyOps(v, list);
      },
      corrJoin: function (names) { return names.join(', ') + ' auf M'; },
      target: function (psi) { return regOf(applyU(M_T, psi)); }
    }
  ];

  var PROTO_BY_KEY = {};
  PROTOCOLS.forEach(function (p) { PROTO_BY_KEY[p.key] = p; });

  /**
   * Ein Durchlauf. `outcomes` ist ein Feld aus ±1 oder null; null heißt
   * "nach Born ziehen". `correct` schaltet den Korrekturschritt ab.
   */
  function runProtocol(p, psi, t0, outcomes, withCorrection, rand) {
    var v = p.prep(psi, t0);
    var stages = [v], results = [], probs = [], probPlus = [], drawn = [];
    var s, pr, pPlus, b, res;

    for (s = 0; s < p.steps.length; s++) {
      pPlus = project(v, p.steps[s].terms, +1).prob;
      b = outcomes[s];
      if (b === null || b === undefined) {
        b = (rand() < pPlus) ? +1 : -1;
        drawn.push(true);
      } else {
        drawn.push(false);
      }
      pr = project(v, p.steps[s].terms, b);
      probs.push(b === +1 ? pPlus : 1 - pPlus);
      probPlus.push(pPlus);
      if (pr.prob < 1e-12) {
        /* Ein unmögliches Ergebnis — der Zustand bliebe undefiniert. */
        results.push(b);
        stages.push(v);
        return {
          impossible: s, stages: stages, results: results,
          probs: probs, probPlus: probPlus, drawn: drawn, corrOps: []
        };
      }
      v = normalise(pr.vec);
      results.push(b);
      stages.push(v);
    }

    var applied = null, corrOps = [];
    if (withCorrection) {
      res = p.correct(v, results);
      corrOps = res.ops;
      v = res.v;
      applied = corrOps.length
        ? p.corrJoin(corrOps.map(function (o) { return o.name; }))
        : null;
    }
    stages.push(v);

    var rhoOut = reduced(v, p.out);
    var rhoGoal = rhoOfPure(p.target(psi, t0));
    return {
      impossible: -1,
      stages: stages,
      results: results,
      probs: probs,
      probPlus: probPlus,
      drawn: drawn,
      applied: applied,
      corrOps: corrOps,
      dev: rhoDev(rhoOut, rhoGoal)
    };
  }

  /* Deterministischer Zufall, damit ein gezeigtes Bild reproduzierbar ist. */
  function mulberry32(a) {
    return function () {
      a |= 0; a = (a + 0x6D2B79F5) | 0;
      var t = Math.imul(a ^ (a >>> 15), 1 | a);
      t = (t + Math.imul(t ^ (t >>> 7), 61 | t)) ^ t;
      return ((t ^ (t >>> 14)) >>> 0) / 4294967296;
    };
  }

  /* ====================================================================
     6 — Selbstprüfung

     CLAUDE.md, Abschnitt 13: was behauptet wird, wird gerechnet. Die
     Ergebnisse stehen im Info-Panel unter "Numerik"; schlägt etwas fehl,
     sagt der Reiter das sichtbar, statt es zu verschweigen.
     ==================================================================== */

  function runSelfCheck() {
    var parts = [], worst = 0;

    function part(name, ok, detail, dev) {
      parts.push({ name: name, ok: ok, detail: detail });
      if (dev !== undefined && isFinite(dev)) worst = Math.max(worst, dev);
    }

    /* (a) Unitarität aller zehn Matrizen */
    var uni = 0;
    GATES.forEach(function (g) { uni = Math.max(uni, mdiff(mmul(mdag(g.U), g.U), ID)); });
    part('Unitarität aller ' + GATES.length + ' Gatter', uni < TOL,
         'größte Abweichung von U†U = 1: ' + MV.fmtExp(uni, 1), uni);

    /* (b) Die Identitäten aus CLAUDE.md, Abschnitt 8 */
    var ids = [
      ['Ŝ_z ≐ e^{−i(π/4)ẑ}', M_SZ_GATE, rotMatrix([0, 0, 1], PI / 2)],
      ['Ŝ_x ≐ e^{−i(π/4)x̂}', M_SX_GATE, rotMatrix([1, 0, 0], PI / 2)],
      ['Ŝ_y = Ŝ_z Ŝ_x Ŝ_z† ≐ e^{−i(π/4)ŷ}', M_SY_GATE, rotMatrix([0, 1, 0], PI / 2)],
      ['Ĥ = (x̂ + ẑ)/√2 ≐ Ŝ_z Ŝ_x Ŝ_z', M_H, mmul(mmul(M_SZ_GATE, M_SX_GATE), M_SZ_GATE)],
      ['T-Gatter ≐ e^{−i(π/8)ẑ}', M_T, rotMatrix([0, 0, 1], PI / 4)],
      ['T-Gatter · T-Gatter ≐ Ŝ_z', mmul(M_T, M_T), M_SZ_GATE],
      ['Ŝ_z · Ŝ_z ≐ ẑ', mmul(M_SZ_GATE, M_SZ_GATE), SZ],
      ['Ĥ · Ĥ = 1', mmul(M_H, M_H), ID]
    ];
    var idWorst = 0, idBad = [];
    ids.forEach(function (t) {
      var d = devUpToPhase(t[1], t[2]);
      idWorst = Math.max(idWorst, d);
      if (d >= TOL) idBad.push(t[0]);
    });
    part('Die acht Gatter-Identitäten aus CLAUDE.md, Abschnitt 8', idBad.length === 0,
         idBad.length ? 'verletzt: ' + idBad.join('; ')
                      : 'alle erfüllt, größte Abweichung ' + MV.fmtExp(idWorst, 1), idWorst);

    /* (c) Der Zopfoperator U₂₃ und Ŝ_z sind derselbe Operator. Gerechnet in
           der γ-Darstellung aus dem Reiter Auslese, nicht in einer zweiten
           eigenen — sonst prüfte die Rechnung nur sich selbst. */
    var braid = null;
    if (MV.majoranaRep && MV.majoranaRep.sectorMatrix) {
      /* U₂₃ = (1/√2)(1 + γ₂γ₃) = (1/√2)(1 − i·(iγ₂γ₃)) = (1/√2)(1 − i ẑ) */
      var zMat = MV.majoranaRep.sectorMatrix(2, 3);
      var miZ = zMat.map(function (e) { return [e[1], -e[0]]; });
      var U23 = mscale(madd(ID, miZ), S2);
      braid = devUpToPhase(U23, M_SZ_GATE);
      part('Zopf ' + MV.eq('ops.braid') + ' und Gatter ' + MV.eq('ops.sz') +
           ' sind derselbe Operator',
           braid < TOL, 'Abweichung bis auf globale Phase: ' + MV.fmtExp(braid, 1), braid);
    } else {
      part('Zopf ' + MV.eq('ops.braid') + ' und Gatter ' + MV.eq('ops.sz'), null,
           'nicht geprüft — die γ-Darstellung aus dem Reiter Auslese ist nicht geladen');
    }

    /* (d) Drehachse und Drehwinkel aus der Matrix beschreiben wirklich die
           Wirkung auf den Bloch-Vektor. Ohne diese Prüfung wäre der
           gezeichnete Bogen nur eine zweite Behauptung. */
    var rotWorst = 0;
    GATES.forEach(function (g) {
      [[0.7, 0.4], [2.1, -1.2], [PI / 2, 0], [0, 0], [PI, 0], [1.3, 2.9]].forEach(function (a) {
        var psi = ketOf(a[0], a[1]);
        var viaState = blochOf(applyU(g.U, psi));
        var viaRot = rotate(blochOf(psi), g.rot.n, g.rot.alpha);
        rotWorst = Math.max(rotWorst,
          Math.hypot(viaState[0] - viaRot[0], viaState[1] - viaRot[1], viaState[2] - viaRot[2]));
      });
    });
    part('Gezeichneter Bogen = Wirkung des Operators', rotWorst < TOL,
         'größte Abweichung zwischen Drehung und Zustandsrechnung über sechs ' +
         'Startzustände und alle Gatter: ' + MV.fmtExp(rotWorst, 1), rotWorst);

    /* (e) Clifford-Gatter bilden Ecken auf Ecken ab, das T-Gatter nicht. */
    var cliffOk = true, tMoves = false, tImage = null, offWorst = 0;
    GATES.forEach(function (g) {
      VERTICES.forEach(function (v) {
        var img = rotate(v.r, g.rot.n, g.rot.alpha);
        var near = nearestVertex(img);
        if (g.clifford) {
          if (near.d >= 1e-9) cliffOk = false;
          offWorst = Math.max(offWorst, near.d);
        } else if (g.key === 'T') {
          if (near.d >= 1e-9) tMoves = true;
          if (v.label === 'x̂ = +1') tImage = img;
        }
      });
    });
    part('Clifford-Gatter bilden die sechs Ecken auf Ecken ab', cliffOk,
         'größte Abweichung von einer Ecke über alle sieben Clifford-Gatter und ' +
         'alle sechs Ecken: ' + MV.fmtExp(offWorst, 1), offWorst);
    part('Das T-Gatter tut das nicht', tMoves,
         tImage ? 'x̂ = +1 landet auf (' + MV.fmt(tImage[0], 4) + '; ' +
                  MV.fmt(tImage[1], 4) + '; ' + MV.fmt(tImage[2], 4) +
                  ') — zwischen den Ecken, wie in CLAUDE.md, Abschnitt 8, angegeben'
                : 'nicht bestimmt');

    /* (f) Auf die Pole von ẑ wirkt das T-Gatter nur als Phase. */
    var poleMove = 0;
    [VERTICES[0], VERTICES[1]].forEach(function (v) {
      var img = rotate(v.r, GATE_BY_KEY.T.rot.n, GATE_BY_KEY.T.rot.alpha);
      poleMove = Math.max(poleMove, Math.hypot(img[0] - v.r[0], img[1] - v.r[1], img[2] - v.r[2]));
    });
    part('Auf die Pole von ẑ wirkt das T-Gatter nur als Phase', poleMove < TOL,
         'Verschiebung beider Pole: ' + MV.fmtExp(poleMove, 1), poleMove);

    /* (g) Die vier Protokolle über alle Kombinationen von Messergebnissen
           und mehrere Startzustände. */
    var rnd = mulberry32(7);
    var protoWorst = 0, protoBad = [], noCorrCount = {};
    var startStates = [[0.9, 0.6], [2.4, -1.1], [PI / 2, 0], [0, 0], [PI, 0], [0.3, 2.2]];

    var generic = ketOf(0.9, 0.6);

    PROTOCOLS.forEach(function (p) {
      var m = p.steps.length, total = 1 << m, wrong = 0, k, j, outs, t0s, ti;
      t0s = p.inQubits === 2 ? [false, true] : [false];
      for (ti = 0; ti < t0s.length; ti++) {
        for (k = 0; k < total; k++) {
          outs = [];
          for (j = 0; j < m; j++) outs.push((k >> j) & 1 ? -1 : +1);
          startStates.forEach(function (a) {
            var psi = ketOf(a[0], a[1]);
            var r = runProtocol(p, psi, t0s[ti], outs, true, rnd);
            if (r.impossible >= 0) return;
            protoWorst = Math.max(protoWorst, r.dev);
            if (r.dev >= TOL) protoBad.push(p.name + ' bei ' + outs.join(','));
          });
          /* Nur der generische Startzustand zählt für die Aussage, wie viele
             Zweige die Korrektur wirklich braucht. */
          if (ti === 0) {
            var r2 = runProtocol(p, generic, false, outs, false, rnd);
            if (r2.impossible < 0 && r2.dev > 1e-6) wrong++;
          }
        }
      }
      noCorrCount[p.key] = { wrong: wrong, total: total };
    });

    part('Alle vier Protokolle, jede Kombination von Messergebnissen, ' +
         startStates.length + ' Startzustände',
         protoBad.length === 0,
         protoBad.length ? 'fehlgeschlagen: ' + protoBad.slice(0, 6).join('; ')
                         : 'größte Abweichung vom Ziel: ' + MV.fmtExp(protoWorst, 1),
         protoWorst);

    var noCorrTxt = PROTOCOLS.map(function (p) {
      var c = noCorrCount[p.key];
      return p.name + ' ' + c.wrong + ' von ' + c.total;
    }).join(', ');
    part('Ohne die Korrektur ist das Ergebnis in einem Teil der Zweige falsch',
         true, noCorrTxt + ', gezählt an einem Startzustand, der Eigenzustand ' +
         'keiner der drei Achsen ist — die Korrektur trägt also wirklich etwas bei');

    /* (h) Die Wahrscheinlichkeiten der beiden Ausgänge summieren sich zu 1. */
    var pSum = 0;
    PROTOCOLS.forEach(function (p) {
      var psi = ketOf(1.1, 0.5), v = p.prep(psi, false), s;
      for (s = 0; s < p.steps.length; s++) {
        var pp = project(v, p.steps[s].terms, +1);
        var pm = project(v, p.steps[s].terms, -1);
        pSum = Math.max(pSum, Math.abs(pp.prob + pm.prob - 1));
        v = normalise(pp.prob > pm.prob ? pp.vec : pm.vec);
      }
    });
    part('P(+1) + P(−1) = 1 in jedem Messschritt', pSum < TOL,
         'größte Abweichung: ' + MV.fmtExp(pSum, 1), pSum);

    /* (i) Der Bogen der Korrektur ist ihre Wirkung.

       Ein Korrekturoperator wirkt nur auf ein Qubit; seine reduzierte
       Dichtematrix geht damit in U ρ U† über, und das ist genau eine Drehung
       des Bloch-Vektors um Achse und Winkel des Operators. Gezeichnet wird
       diese Drehung — geprüft wird, dass sie mit der Zustandsrechnung
       übereinstimmt, sonst zeigte der Bogen etwas anderes als die Rechnung. */
    var arcWorst = 0, arcBad = [], arcCount = 0;
    PROTOCOLS.forEach(function (p) {
      var m = p.steps.length, total = 1 << m, k, j, outs;
      var t0s = p.inQubits === 2 ? [false, true] : [false];
      t0s.forEach(function (t0) {
        for (k = 0; k < total; k++) {
          outs = [];
          for (j = 0; j < m; j++) outs.push((k >> j) & 1 ? -1 : +1);
          startStates.forEach(function (a) {
            var r = runProtocol(p, ketOf(a[0], a[1]), t0, outs, true, rnd);
            if (r.impossible >= 0) return;
            var before = r.stages[m];
            r.corrOps.forEach(function (o) {
              var rot = axisAngle(o.U);
              var rb = blochOfRho(reduced(before, [o.qubit]));
              var ra = blochOfRho(reduced(o.after, [o.qubit]));
              var got = rotate(rb, rot.n, rot.alpha);
              var d = Math.hypot(got[0] - ra[0], got[1] - ra[1], got[2] - ra[2]);
              arcWorst = Math.max(arcWorst, d);
              arcCount++;
              if (d >= TOL) arcBad.push(p.name + ' ' + o.name + ' bei ' + outs.join(','));
              before = o.after;
            });
          });
        }
      });
    });
    part('Bogen = Wirkung der Korrektur', arcBad.length === 0,
         arcBad.length ? 'verletzt: ' + arcBad.slice(0, 5).join('; ')
           : arcCount + ' Korrekturoperatoren über alle vier Protokolle, alle ' +
             'Ausgänge und ' + startStates.length + ' Startzustände; größte ' +
             'Abweichung ' + MV.fmtExp(arcWorst, 1),
         arcWorst);

    /* (ii) Kontrollwerte des Ŝ_z-Protokolls aus CLAUDE.md, Abschnitt 8. */
    var pSz = PROTO_BY_KEY.Sz, ctlBad = [], ctlWorst = 0;
    function near(a, b) { var d = Math.abs(a - b); ctlWorst = Math.max(ctlWorst, d); return d < 1e-9; }
    function len(r) { return Math.hypot(r[0], r[1], r[2]); }

    [[1, 1], [1, -1], [-1, 1], [-1, -1]].forEach(function (o) {
      var r = runProtocol(pSz, ketOf(PI / 2, 0), false, o, true, rnd);
      var c1 = blochOfRho(reduced(r.stages[1], [0]));
      var a1 = blochOfRho(reduced(r.stages[1], [1]));
      if (!near(len(c1), 0) || !near(len(a1), 0)) ctlBad.push('Länge nach Messung 1 bei ' + o);
      var c2 = blochOfRho(reduced(r.stages[2], [0]));
      if (!near(c2[0], 0) || !near(Math.abs(c2[1]), 1) || !near(c2[2], 0)) {
        ctlBad.push('C nach Messung 2 bei ' + o);
      }
      var c3 = blochOfRho(reduced(r.stages[3], [0]));
      if (!near(c3[0], 0) || !near(c3[1], 1) || !near(c3[2], 0)) ctlBad.push('Ende bei ' + o);
    });

    var th60 = 60 / DEG, cos60 = Math.cos(th60);
    [1, -1].forEach(function (b1) {
      var r = runProtocol(pSz, ketOf(th60, 0.8), false, [b1, 1], true, rnd);
      var c = blochOfRho(reduced(r.stages[1], [0]));
      var a = blochOfRho(reduced(r.stages[1], [1]));
      if (!near(c[0], 0) || !near(c[1], 0) || !near(c[2], cos60)) ctlBad.push('C bei θ = 60°');
      if (!near(a[0], 0) || !near(a[1], 0) || !near(a[2], b1 * cos60)) ctlBad.push('A bei θ = 60°');
      if (!near(r.probPlus[0], 0.5)) ctlBad.push('P(+1) bei θ = 60°');
    });

    part('Kontrollwerte des Ŝ_z-Protokolls', ctlBad.length === 0,
         ctlBad.length ? 'verletzt: ' + ctlBad.slice(0, 4).join('; ')
           : 'Eingang x̂ = +1: beide Längen 0 nach Messung 1, C auf (0, ±1, 0) nach ' +
             'Messung 2, in allen vier Zweigen (0, 1, 0) am Ende. Eingang θ = 60°: ' +
             'C auf (0, 0, ' + MV.fmt(cos60, 3) + '), A auf (0, 0, b₁ · ' +
             MV.fmt(cos60, 3) + '), P(+1) = 0,500. Größte Abweichung ' +
             MV.fmtExp(ctlWorst, 1),
         ctlWorst);

    /* (iii) Die Wiedergabe rechnet nichts zweites: ihr Endbild ist die letzte
       Stufe aus runProtocol. */
    var fWorst = 0, fBad = [];
    PROTOCOLS.forEach(function (p) {
      var outs = [], j;
      for (j = 0; j < p.steps.length; j++) outs.push(j % 2 ? -1 : 1);
      var r = runProtocol(p, ketOf(0.9, 0.6), false, outs, true, rnd);
      if (r.impossible >= 0) return;
      var lastS = p.steps.length + 1;
      var f = frameOf(r, p, lastS, 1);
      var vEnd = r.stages[r.stages.length - 1], q;
      for (q = 0; q < p.roles.length; q++) {
        var want = blochOfRho(reduced(vEnd, [q]));
        var d = Math.hypot(f.vecs[q][0] - want[0], f.vecs[q][1] - want[1],
                           f.vecs[q][2] - want[2]);
        fWorst = Math.max(fWorst, d);
        if (d >= TOL) fBad.push(p.name + ', Qubit ' + p.roles[q]);
      }
    });
    part('Wiedergabe = Rechnung', fBad.length === 0,
         fBad.length ? 'weicht ab: ' + fBad.join('; ')
           : 'das Endbild der Wiedergabe ist in allen vier Protokollen die letzte ' +
             'Stufe aus runProtocol, Abweichung ' + MV.fmtExp(fWorst, 1),
         fWorst);

    /* (iv) Eine Messung ist ein Sprung — nie eine Bewegung.

       CLAUDE.md §8: "Es gibt keinen Weg zwischen dem Zustand vorher und
       nachher, also werden keine Zwischenzustände interpoliert." Geprüft wird
       das an der Wiedergabe selbst: In einem Messschritt muss jedes Einzelbild
       entweder genau den Zustand vor der Messung oder genau den danach zeigen
       — nichts dazwischen. Eine Interpolation fiele hier sofort auf. */
    var jumpBad = [], jumpFrames = 0;
    PROTOCOLS.forEach(function (p) {
      var outs = [], j;
      for (j = 0; j < p.steps.length; j++) outs.push(j % 2 ? -1 : 1);
      var r = runProtocol(p, ketOf(0.9, 0.6), false, outs, true, rnd);
      if (r.impossible >= 0) return;
      var i2, u2, q2, f2, before, after, k2;
      for (i2 = 1; i2 <= p.steps.length; i2++) {
        before = r.stages[i2 - 1];
        after = r.stages[i2];
        for (k2 = 0; k2 <= 40; k2++) {
          u2 = k2 / 40;
          f2 = frameOf(r, p, i2, u2);
          jumpFrames++;
          for (q2 = 0; q2 < p.roles.length; q2++) {
            var rb = blochOfRho(reduced(before, [q2]));
            var ra = blochOfRho(reduced(after, [q2]));
            var got = f2.vecs[q2];
            var db = Math.hypot(got[0] - rb[0], got[1] - rb[1], got[2] - rb[2]);
            var da = Math.hypot(got[0] - ra[0], got[1] - ra[1], got[2] - ra[2]);
            if (Math.min(db, da) >= TOL) {
              jumpBad.push(p.name + ', Schritt ' + i2 + ', u = ' + MV.fmt(u2, 2));
            }
          }
        }
      }
    });
    part('Eine Messung ist ein Sprung, keine Bewegung', jumpBad.length === 0,
         jumpBad.length
           ? 'interpolierter Zwischenzustand bei: ' + jumpBad.slice(0, 4).join('; ')
           : jumpFrames + ' Einzelbilder über alle Messschritte aller vier ' +
             'Protokolle zeigen ausnahmslos den Zustand vor oder den nach der ' +
             'Messung — nie etwas dazwischen');

    var passed = parts.every(function (p) { return p.ok !== false; });
    return { parts: parts, passed: passed, worst: worst, noCorr: noCorrCount };
  }

  /* ====================================================================
     7 — Die Bloch-Kugel

     Eigene orthografische Projektion, keine Bibliothek. Die beiden
     Bildschirmachsen sind orthonormal, es ist also eine echte
     Parallelprojektion: die Einheitskugel wird genau in den gezeichneten
     Umrisskreis abgebildet, und ein Bloch-Vektor der Länge 1 kann nie
     darüber hinausragen. Preis dafür ist die Verkürzung — ein Vektor, der
     auf den Betrachter zu zeigt, erscheint kurz. Weil die Kugel gedreht
     werden kann, ist das hier kein Problem mehr: man dreht sie weg.
     ==================================================================== */

  var cam = { az: 55 * PI / 180, el: 20 * PI / 180 };

  function basis() {
    var az = cam.az, el = cam.el;
    return {
      u: [-Math.sin(az), Math.cos(az), 0],
      w: [-Math.cos(az) * Math.sin(el), -Math.sin(az) * Math.sin(el), Math.cos(el)],
      v: [Math.cos(az) * Math.cos(el), Math.sin(az) * Math.cos(el), Math.sin(el)]
    };
  }

  function proj(B, r) {
    return {
      u: B.u[0] * r[0] + B.u[1] * r[1] + B.u[2] * r[2],
      w: B.w[0] * r[0] + B.w[1] * r[1] + B.w[2] * r[2],
      v: B.v[0] * r[0] + B.v[1] * r[1] + B.v[2] * r[2]      /* > 0: vorn      */
    };
  }

  /** Ein Zeichner für eine Kugel an fester Stelle im Bild. */
  function sphere(ctx, cx, cy, R) {
    var B = basis();
    var api = {
      px: function (r) {
        var s = proj(B, r);
        return [cx + R * s.u, cy - R * s.w];
      },
      front: function (r) { return proj(B, r).v >= 0; }
    };

    /* Umriss */
    ctx.save();
    ctx.strokeStyle = MV.palette.line;
    ctx.lineWidth = 1.4;
    ctx.beginPath();
    ctx.arc(cx, cy, R, 0, 2 * PI);
    ctx.stroke();

    /* Drei Großkreise: Äquator und die beiden Meridiane durch die Achsen.
       Die hinteren Hälften blasser, sonst ist nicht zu sehen, wie die Kugel
       im Raum liegt. */
    [[[1, 0, 0], [0, 1, 0]], [[0, 0, 1], [1, 0, 0]], [[0, 0, 1], [0, 1, 0]]].forEach(function (pair) {
      var i, t, r, p, wasFront = null;
      ctx.lineWidth = 1;
      ctx.beginPath();
      for (i = 0; i <= 144; i++) {
        t = i / 144 * 2 * PI;
        r = [pair[0][0] * Math.cos(t) + pair[1][0] * Math.sin(t),
             pair[0][1] * Math.cos(t) + pair[1][1] * Math.sin(t),
             pair[0][2] * Math.cos(t) + pair[1][2] * Math.sin(t)];
        p = api.px(r);
        var f = api.front(r);
        if (f !== wasFront) {
          if (wasFront !== null) ctx.stroke();
          ctx.beginPath();
          ctx.globalAlpha = f ? 0.9 : 0.35;
          ctx.strokeStyle = MV.palette.line;
          ctx.setLineDash(f ? [] : [3, 4]);
          ctx.moveTo(p[0], p[1]);
          wasFront = f;
        } else {
          ctx.lineTo(p[0], p[1]);
        }
      }
      ctx.stroke();
      ctx.globalAlpha = 1;
      ctx.setLineDash([]);
    });
    ctx.restore();
    return api;
  }

  /** Einen Vektor als Pfeil vom Mittelpunkt aus. */
  function drawVector(ctx, sp, cx, cy, r, color, width) {
    var e = sp.px(r);
    ctx.save();
    ctx.strokeStyle = color;
    ctx.fillStyle = color;
    ctx.lineWidth = width || 2.4;
    ctx.lineCap = 'round';
    ctx.beginPath();
    ctx.moveTo(cx, cy);
    ctx.lineTo(e[0], e[1]);
    ctx.stroke();
    ctx.beginPath();
    ctx.arc(e[0], e[1], 4, 0, 2 * PI);
    ctx.fill();
    ctx.restore();
    return e;
  }

  /* ====================================================================
     8 — Zustand des Reiters
     ==================================================================== */

  var state = {
    theta: 90,          /* Startzustand, Grad                                */
    phi: 0,
    seq: [],            /* angewandte Gatter, als Schlüssel                  */

    proto: 0,           /* Index in PROTOCOLS                                */
    target1: false,     /* Ziel-Qubit des CNOT startet in |1⟩                */
    outcomes: [null, null, null],
    correct: true,
    step: 0,            /* angezeigter Schritt in Panel B, 0 = Start         */
    seed: 7
  };

  var sliders = {}, panelA = null, panelB = null;
  var readoutA = null, gateTableHost = null, seqHost = null;
  var ladderHost = null, readoutB = null, outcomeRow = null, stepSeg = null;
  var counterHost = null, playButton = null;
  var segProto = null, corrButtons = null, targetButtons = null;

  /* ---------------------------------------------------- abgeleitete Größen */

  /** Der Weg des Bloch-Vektors: Startpunkt und je ein Bogen pro Gatter. */
  function path() {
    var psi = ketOf(state.theta / DEG, state.phi / DEG);
    var pts = [blochOf(psi)], arcs = [], i, g;
    for (i = 0; i < state.seq.length; i++) {
      g = GATE_BY_KEY[state.seq[i]];
      var from = pts[pts.length - 1];
      psi = applyU(g.U, psi);
      var to = blochOf(psi);
      arcs.push({ gate: g, from: from, to: to });
      pts.push(to);
    }
    return { psi: psi, pts: pts, arcs: arcs, r: pts[pts.length - 1] };
  }

  /* ====================================================================
     9 — Panel A: die Gatter auf der Kugel
     ==================================================================== */

  function drawA() {
    var p = panelA.plot;
    if (!p) return;
    var ctx = p.ctx;

    p.begin({ x: [0, 10], y: [0, 10], grid: false, frame: false, xTicks: [], yTicks: [] });

    var cx = (p.px0 + p.px1) / 2;
    var cy = (p.py0 + p.py1) / 2;
    var R = Math.min((p.px1 - p.px0) * 0.30, (p.py1 - p.py0) * 0.40);

    var sp = sphere(ctx, cx, cy, R);
    var w = path();

    /* Die zwölf Kanten des Oktaeders. Sie machen sichtbar, dass die sechs
       Eigenzustände einen Körper aufspannen und ein Clifford-Gatter ihn auf
       sich selbst abbildet. */
    ctx.save();
    ctx.lineWidth = 1;
    VERTICES.forEach(function (a, i) {
      VERTICES.forEach(function (b, j) {
        if (j <= i) return;
        /* Gegenüberliegende Ecken sind keine Kante. */
        if (Math.abs(a.r[0] + b.r[0]) < 1e-9 &&
            Math.abs(a.r[1] + b.r[1]) < 1e-9 &&
            Math.abs(a.r[2] + b.r[2]) < 1e-9) return;
        var mid = [(a.r[0] + b.r[0]) / 2, (a.r[1] + b.r[1]) / 2, (a.r[2] + b.r[2]) / 2];
        var f = sp.front(mid);
        var pa = sp.px(a.r), pb = sp.px(b.r);
        ctx.globalAlpha = f ? 0.55 : 0.18;
        ctx.strokeStyle = MV.palette.curve1;
        ctx.setLineDash(f ? [] : [3, 4]);
        ctx.beginPath();
        ctx.moveTo(pa[0], pa[1]);
        ctx.lineTo(pb[0], pb[1]);
        ctx.stroke();
      });
    });
    ctx.globalAlpha = 1;
    ctx.setLineDash([]);
    ctx.restore();

    /* Die Drehachse des zuletzt angewandten Gatters */
    if (w.arcs.length) {
      var g = w.arcs[w.arcs.length - 1].gate;
      var a1 = sp.px([g.rot.n[0] * 1.24, g.rot.n[1] * 1.24, g.rot.n[2] * 1.24]);
      var a2 = sp.px([-g.rot.n[0] * 1.24, -g.rot.n[1] * 1.24, -g.rot.n[2] * 1.24]);
      ctx.save();
      ctx.strokeStyle = MV.palette.critical;
      ctx.lineWidth = 1.2;
      ctx.setLineDash([5, 4]);
      ctx.globalAlpha = 0.8;
      ctx.beginPath();
      ctx.moveTo(a1[0], a1[1]);
      ctx.lineTo(a2[0], a2[1]);
      ctx.stroke();
      ctx.restore();
      p.label('Drehachse ' + g.label, a1[0], a1[1] - 7, {
        color: MV.palette.critical, align: 'center', baseline: 'bottom', box: true
      });
    }

    /* Die Bögen. Der letzte kräftig, die früheren blass — so bleibt die
       Reihenfolge lesbar, ohne dass der Weg zum Knäuel wird. */
    w.arcs.forEach(function (arc, idx) {
      var last = idx === w.arcs.length - 1;
      var n = arc.gate.rot.n, al = arc.gate.rot.alpha;
      var steps = 64, i, r, q, f, wasFront = null;
      ctx.save();
      ctx.lineWidth = last ? 2.6 : 1.6;
      ctx.lineCap = 'round';
      for (i = 0; i <= steps; i++) {
        r = rotate(arc.from, n, al * i / steps);
        q = sp.px(r);
        f = sp.front(r);
        if (f !== wasFront) {
          if (wasFront !== null) ctx.stroke();
          ctx.beginPath();
          ctx.strokeStyle = last ? MV.palette.topological : MV.palette.muted;
          ctx.globalAlpha = (last ? 1 : 0.45) * (f ? 1 : 0.4);
          ctx.setLineDash(f ? [] : [3, 3]);
          ctx.moveTo(q[0], q[1]);
          wasFront = f;
        } else {
          ctx.lineTo(q[0], q[1]);
        }
      }
      ctx.stroke();
      ctx.restore();
    });

    /* Die Zwischenpunkte des Weges */
    ctx.save();
    w.pts.forEach(function (r, i) {
      if (i === w.pts.length - 1) return;
      var q = sp.px(r);
      ctx.globalAlpha = sp.front(r) ? 0.75 : 0.3;
      ctx.fillStyle = MV.palette.bg;
      ctx.strokeStyle = i === 0 ? MV.palette.trivial : MV.palette.muted;
      ctx.lineWidth = 1.4;
      ctx.beginPath();
      ctx.arc(q[0], q[1], i === 0 ? 4 : 3, 0, 2 * PI);
      ctx.fill();
      ctx.stroke();
    });
    ctx.restore();

    /* Die sechs Ecken */
    VERTICES.forEach(function (v) {
      var q = sp.px(v.r);
      var f = sp.front(v.r);
      ctx.save();
      ctx.globalAlpha = f ? 1 : 0.4;
      ctx.fillStyle = MV.palette.curve1;
      ctx.beginPath();
      ctx.arc(q[0], q[1], 3.4, 0, 2 * PI);
      ctx.fill();
      ctx.restore();

      var dx = q[0] - cx, dy = q[1] - cy;
      var len = Math.hypot(dx, dy);
      var ox = len < 1 ? 0 : dx / len * 15;
      var oy = len < 1 ? -15 : dy / len * 15;
      p.label(v.label + '  ' + v.ket, q[0] + ox, q[1] + oy, {
        color: f ? MV.palette.curve1 : MV.palette.muted,
        align: Math.abs(dx) < 6 ? 'center' : (dx > 0 ? 'left' : 'right'),
        baseline: Math.abs(dy) < 6 ? 'middle' : (dy > 0 ? 'top' : 'bottom'),
        box: true,
        font: '11px system-ui'
      });
    });

    /* Der aktuelle Zustand */
    var end = drawVector(ctx, sp, cx, cy, w.r, MV.palette.topological, 2.6);
    var nv = nearestVertex(w.r);
    if (nv.d >= 1e-9) {
      p.label('zwischen den Ecken', end[0], end[1] - 10, {
        color: MV.palette.topological, align: 'center', baseline: 'bottom', box: true
      });
    }

    p.finish();

    panelA.setStatus(
      state.seq.length === 0
        ? 'Noch kein Gatter angewandt — der Zeiger steht auf dem Startzustand. ' +
          'Die Kugel lässt sich mit der Maus drehen.'
        : 'Folge: ' + state.seq.map(function (k) { return GATE_BY_KEY[k].label; }).join(' → ') +
          '. ' + (nv.d < 1e-9
            ? 'Der Zustand steht auf der Ecke ' + nv.v.label + '.'
            : 'Der Zustand steht zwischen den Ecken — das kann kein Clifford-Gatter erzeugen.'),
      state.seq.length && nv.d >= 1e-9 ? MV.palette.critical : MV.palette.muted);
  }

  function updateReadoutA() {
    var w = path();
    var r = w.r;
    var th = Math.acos(Math.max(-1, Math.min(1, r[2]))) * DEG;
    var ph = Math.atan2(r[1], r[0]) * DEG;
    var nv = nearestVertex(r);

    readoutA.set('rx', MV.fmt(r[0], 4));
    readoutA.set('ry', MV.fmt(r[1], 4));
    readoutA.set('rz', MV.fmt(r[2], 4));
    readoutA.set('theta', MV.fmt(th, 1) + '°');
    readoutA.set('phi', MV.fmt(ph, 1) + '°');
    readoutA.set('len', MV.fmt(Math.hypot(r[0], r[1], r[2]), 6));
    readoutA.set('vertex', nv.d < 1e-9 ? 'ja — ' + nv.v.label + ' ' + nv.v.ket : 'nein');
    readoutA.set('count', String(state.seq.length));

    if (state.seq.length === 0) {
      readoutA.status('trivial', 'Startzustand, noch kein Gatter angewandt');
    } else if (nv.d < 1e-9) {
      readoutA.status('topological', 'auf einer Ecke des Oktaeders');
    } else {
      readoutA.status('critical', 'zwischen den Ecken — nur mit einem T-Gatter erreichbar');
    }

    seqHost.innerHTML = state.seq.length
      ? '<b>Angewandt:</b> ' + state.seq.map(function (k, i) {
          return '<span class="gateseq__item">' + (i + 1) + '. ' + GATE_BY_KEY[k].label + '</span>';
        }).join('')
      : '<span class="gateseq__empty">Noch kein Gatter angewandt.</span>';
  }

  /** Ein Eintrag der komplexen 2×2-Matrix als kurze Zeichenkette. */
  function fmtC(e) {
    var re = e[0], im = e[1];
    var near = function (a, b) { return Math.abs(a - b) < 1e-12; };
    if (near(re, 0) && near(im, 0)) return '0';
    if (near(im, 0)) return near(re, 1) ? '1' : near(re, -1) ? '−1' : MV.fmt(re, 3);
    if (near(re, 0)) return near(im, 1) ? 'i' : near(im, -1) ? '−i' : MV.fmt(im, 3) + 'i';
    return MV.fmt(re, 3) + (im < 0 ? '−' : '+') + MV.fmt(Math.abs(im), 3) + 'i';
  }

  function fmtMatrix(U) {
    return '<span class="mat2"><span class="mat2__grid">' + U.map(function (e) {
      return '<span>' + fmtC(e) + '</span>';
    }).join('') + '</span></span>';
  }

  /** Wo eine Ecke landet — als Eckenname, sonst als Koordinaten. */
  function imageOf(g, v) {
    var img = rotate(v.r, g.rot.n, g.rot.alpha);
    var near = nearestVertex(img);
    if (near.d < 1e-9) return near.v.label;
    return '(' + MV.fmt(img[0], 3) + '; ' + MV.fmt(img[1], 3) + '; ' + MV.fmt(img[2], 3) + ')';
  }

  function axisName(n) {
    var names = [[[1, 0, 0], 'x̂'], [[0, 1, 0], 'ŷ'], [[0, 0, 1], 'ẑ'],
                 [[S2, 0, S2], '(x̂+ẑ)/√2']];
    var hit = null;
    names.forEach(function (t) {
      if (Math.hypot(n[0] - t[0][0], n[1] - t[0][1], n[2] - t[0][2]) < 1e-9) hit = t[1];
      if (Math.hypot(n[0] + t[0][0], n[1] + t[0][1], n[2] + t[0][2]) < 1e-9) hit = '−' + t[1];
    });
    return hit || '(' + MV.fmt(n[0], 3) + '; ' + MV.fmt(n[1], 3) + '; ' + MV.fmt(n[2], 3) + ')';
  }

  /**
   * Die Gattertabelle. Jede Zelle ist gerechnet — Drehachse und Drehwinkel
   * aus der Matrix, die beiden Bildspalten durch Anwenden der Drehung.
   */
  function renderGateTable() {
    var html = '<table class="numtable gatetable"><caption>Alle neun Gatter, ' +
      'gerechnet aus der jeweiligen Matrix</caption><thead><tr>' +
      '<th scope="col">Gatter</th><th scope="col">Matrix</th>' +
      '<th scope="col">Drehachse</th><th scope="col" class="num">Winkel</th>' +
      '<th scope="col">x̂ = +1 landet auf</th>' +
      '<th scope="col">Ecken auf Ecken</th></tr></thead><tbody>';

    GATES.forEach(function (g) {
      var allVertices = VERTICES.every(function (v) {
        return nearestVertex(rotate(v.r, g.rot.n, g.rot.alpha)).d < 1e-9;
      });
      html += '<tr' + (g.clifford ? '' : ' class="numtable__binding"') + '>' +
        '<th scope="row">' + g.label + '</th>' +
        '<td>' + fmtMatrix(g.U) + '</td>' +
        '<td>' + axisName(g.rot.n) + '</td>' +
        '<td class="num">' + MV.fmt(g.rot.alpha * DEG, 0) + '°</td>' +
        '<td>' + imageOf(g, VERTICES[2]) + '</td>' +
        '<td>' + (allVertices ? 'ja' : '<b>nein</b>') + '</td></tr>';
    });

    html += '</tbody></table>' +
      '<p class="note">Hervorgehoben sind die beiden Zeilen, die keine ' +
      'Clifford-Gatter sind. Sie sind die einzigen, in denen die letzte Spalte ' +
      '„nein" trägt — und die einzigen, die am Bauelement nicht aus Messungen ' +
      'an Majorana-Bilinearen erzeugt werden können.</p>';
    gateTableHost.innerHTML = html;
  }

  /* ------------------------------------------------------- Drehen per Maus */

  function attachDrag(canvas) {
    var dragging = false, lastX = 0, lastY = 0;
    canvas.style.cursor = 'grab';

    canvas.addEventListener('pointerdown', function (ev) {
      dragging = true;
      lastX = ev.clientX; lastY = ev.clientY;
      canvas.style.cursor = 'grabbing';
      if (canvas.setPointerCapture) canvas.setPointerCapture(ev.pointerId);
      ev.preventDefault();
    });

    canvas.addEventListener('pointermove', function (ev) {
      if (!dragging) return;
      cam.az -= (ev.clientX - lastX) * 0.010;
      cam.el += (ev.clientY - lastY) * 0.010;
      /* Über den Pol hinaus zu drehen würde das Bild spiegeln, ohne einen
         neuen Blickwinkel zu erschließen. */
      cam.el = Math.max(-1.45, Math.min(1.45, cam.el));
      lastX = ev.clientX; lastY = ev.clientY;
      drawA();
    });

    function stop(ev) {
      if (!dragging) return;
      dragging = false;
      canvas.style.cursor = 'grab';
      if (canvas.releasePointerCapture && ev.pointerId !== undefined) {
        try { canvas.releasePointerCapture(ev.pointerId); } catch (e) { /* egal */ }
      }
    }
    canvas.addEventListener('pointerup', stop);
    canvas.addEventListener('pointercancel', stop);
    canvas.addEventListener('pointerleave', stop);
  }

  /* ====================================================================
     10 — Panel B: dieselben Gatter aus Messungen

     Die Wiedergabe ist eine Schrittfolge, kein Film (CLAUDE.md §8). Jede
     Messung ist ein Sprung — es gibt im Modell keinen Weg zwischen dem
     Zustand davor und danach, also wird auch keiner gezeichnet. Nur die
     Korrektur ist ein unitäres Gatter und darf als Bogen laufen. Das Tempo
     ist ein Anzeigetempo und keine physikalische Zeit.

     Gerechnet wird nichts zweites: die Wiedergabe liest die Stufen aus
     runProtocol und stellt sie nacheinander dar.
     ==================================================================== */

  var AX_SYM = { x: 'x̂', y: 'ŷ', z: 'ẑ' };

  /**
   * Das Zielqubit des CNOT heißt in der Arbeit (§5.5) und bei Plugge „T" — und
   * kollidiert damit mit dem T-Gatter im selben Abschnitt. Umbenannt wird es
   * nicht, die Seite folgt der Arbeit; überall, wo es auftaucht, steht aber
   * „T (Ziel)". Das Gatter heißt stets „T-Gatter".
   */
  function roleLabel(role) { return role === 'T' ? 'T (Ziel)' : role; }

  /* Zeitaufteilung innerhalb eines Messschritts (CLAUDE.md §8, Richtwert) */
  var U_RESULT = 1 / 3;     /* ab hier steht das Ergebnis da                 */
  var U_JUMP = 1 / 2;       /* ab hier ist gesprungen                        */
  var FADE_S = 0.3;         /* Überblendung am Ort, höchstens                */

  var play = {
    running: false,
    u: 1,                   /* Fortschritt im angezeigten Schritt, 0 … 1     */
    tempo: 1.5,             /* Sekunden je Schritt                           */
    raf: 0,
    last: 0,
    runs: 0,
    hits: 0,
    counted: false
  };

  var runCache = { key: null, run: null };

  function reducedMotion() {
    return window.matchMedia &&
           window.matchMedia('(prefers-reduced-motion: reduce)').matches;
  }

  /** Der gezeichnete Fortschritt. Bei reduzierter Bewegung immer der Endstand. */
  function drawU() { return reducedMotion() ? 1 : play.u; }

  function currentRun() {
    var p = PROTOCOLS[state.proto];
    var key = [state.proto, state.theta, state.phi, state.target1,
               state.outcomes.slice(0, p.steps.length).join(','),
               state.correct, state.seed].join('|');
    if (runCache.key === key) return runCache.run;
    var psi = ketOf(state.theta / DEG, state.phi / DEG);
    runCache.key = key;
    runCache.run = runProtocol(p, psi, state.target1,
      state.outcomes.slice(0, p.steps.length), state.correct, mulberry32(state.seed));
    return runCache.run;
  }

  /** Der letzte anzeigbare Schritt — bei unmöglichem Ausgang früher. */
  function lastStep(run, proto) {
    return run.impossible >= 0 ? run.impossible + 1 : proto.steps.length + 1;
  }

  function blochAt(v, k) { return blochOfRho(reduced(v, [k])); }

  /**
   * Die Abweichung vom Ziel. Unterhalb der Toleranz ist sie Rundung und keine
   * Physik — dann steht „0" da, und der wahre Wert daneben in Klammern, damit
   * niemand eine Null für ein Weglassen hält.
   */
  function devText(dev) {
    return dev < TOL
      ? '0 (Rundung, ' + MV.fmtExp(dev, 1) + ')'
      : MV.fmtExp(dev, 2);
  }

  /** Der Bloch-Vektor, auf dem das Ausgabequbit landen soll. */
  function targetBloch(proto, k) {
    var psi = ketOf(state.theta / DEG, state.phi / DEG);
    var tgt = proto.target(psi, state.target1);
    return blochOfRho(reduced(tgt, [proto.out.indexOf(k)]));
  }

  function termsText(proto, terms) {
    return terms.map(function (t) {
      return AX_SYM[t[1]] + '_' + proto.roles[t[0]];
    }).join(' ');
  }

  /**
   * Alles, was ein Einzelbild braucht — aus Schritt und Fortschritt.
   * Hier wird nur gelesen und interpoliert, nie gerechnet.
   */
  function frameOf(run, proto, step, u) {
    var n = proto.roles.length;
    var f = {
      step: step, u: u, vecs: [], shadow: null, arcs: [],
      title: '', caption: '', highlight: null,
      showResult: false, jumped: false, fade: 1
    };
    var q;

    if (step <= 0) {
      for (q = 0; q < n; q++) f.vecs.push(blochAt(run.stages[0], q));
      f.title = 'Präparation';
      f.caption = proto.prepText;
      return f;
    }

    var m = proto.steps.length;

    if (step <= m) {
      var i = step - 1;
      var before = run.stages[i];
      var after = run.stages[Math.min(i + 1, run.stages.length - 1)];
      var st = proto.steps[i];
      var impossible = run.impossible === i;

      f.title = 'Messung ' + termsText(proto, st.terms);
      f.highlight = st.terms;
      f.showResult = u >= U_RESULT;
      f.jumped = !impossible && u >= U_JUMP;

      for (q = 0; q < n; q++) f.vecs.push(blochAt(f.jumped ? after : before, q));
      if (f.jumped) {
        f.shadow = [];
        for (q = 0; q < n; q++) f.shadow.push(blochAt(before, q));
        var span = Math.min(FADE_S / Math.max(play.tempo, 1e-6), 1 - U_JUMP);
        f.fade = reducedMotion() || span <= 0
          ? 1 : Math.min(1, (u - U_JUMP) / span);
      }

      var pPlus = run.probPlus[i];
      f.caption = 'P(+1) = ' + MV.fmt(pPlus, 3) + ' · P(−1) = ' + MV.fmt(1 - pPlus, 3);
      if (f.showResult) {
        f.caption += '  →  ' + st.sym + ' = ' + MV.fmtSigned(run.results[i], 0) +
          ' (' + (run.drawn[i] ? 'gezogen' : 'von Hand') + ')';
      }
      if (impossible && f.showResult) {
        f.caption += '  —  Wahrscheinlichkeit null, dieser Ausgang kann nicht eintreten.';
      } else if (f.jumped) {
        f.caption += '  —  Sprung durch Messung: Projektion (1 ± Ô)/2, dann Normierung';
      }
      return f;
    }

    /* Korrekturschritt */
    var base = run.stages[m];
    for (q = 0; q < n; q++) f.vecs.push(blochAt(base, q));
    f.title = 'Korrektur';
    f.showResult = true;

    var ops = run.corrOps || [];
    if (!ops.length) {
      f.caption = !state.correct
        ? 'Korrektur weggelassen'
        : 'keine Korrektur nötig (' + proto.noCorrCond + ')';
      return f;
    }

    var slice = 1 / ops.length;
    var j, rot, from, w, cur = Math.min(ops.length - 1, Math.floor(u / slice));
    if (u >= 1) cur = ops.length;
    for (j = 0; j < ops.length; j++) {
      rot = axisAngle(ops[j].U);
      from = f.vecs[ops[j].qubit];
      if (j < cur) {
        f.vecs[ops[j].qubit] = blochAt(ops[j].after, ops[j].qubit);
        f.arcs.push({ qubit: ops[j].qubit, from: from, n: rot.n, alpha: rot.alpha, upto: 1 });
      } else if (j === cur) {
        w = Math.max(0, Math.min(1, (u - j * slice) / slice));
        f.vecs[ops[j].qubit] = rotate(from, rot.n, rot.alpha * w);
        f.arcs.push({ qubit: ops[j].qubit, from: from, n: rot.n, alpha: rot.alpha, upto: w });
      }
    }
    f.caption = ops.map(function (o) { return o.name; }).join(' → ') +
      '  —  der Bogen zeigt, was die Korrektur am Zustand bewirkt, keine ' +
      'Zeitentwicklung; am Bauelement ist sie ein bestätigter Elektronentransfer ' +
      'oder reine Buchführung.';
    return f;
  }

  /* ------------------------------------------------------------ Zeichnen */

  function drawB() {
    var p = panelB.plot;
    if (!p) return;
    var ctx = p.ctx;
    var proto = PROTOCOLS[state.proto];
    var run = currentRun();
    var step = Math.min(state.step, lastStep(run, proto));
    var f = frameOf(run, proto, step, drawU());

    p.begin({ x: [0, 10], y: [0, 10], grid: false, frame: false, xTicks: [], yTicks: [] });

    var n = proto.roles.length;
    var R = Math.min((p.px1 - p.px0) / (n * 2.7), (p.py1 - p.py0) * 0.33);
    var gap = (p.px1 - p.px0) / n;
    var cy = (p.py0 + p.py1) / 2 + 4;

    p.label(f.title, (p.px0 + p.px1) / 2, p.py0 - 8, {
      color: MV.palette.ink, align: 'center', baseline: 'bottom',
      font: '600 14px system-ui'
    });

    proto.roles.forEach(function (role, k) {
      var cx = p.px0 + gap * (k + 0.5);
      var sp = sphere(ctx, cx, cy, R);
      var isOut = proto.out.indexOf(k) >= 0;

      /* Achsen; die gemessene Achse dieses Schritts hervorgehoben. */
      [[[1, 0, 0], 'x̂', 'x'], [[0, 1, 0], 'ŷ', 'y'], [[0, 0, 1], 'ẑ', 'z']].forEach(function (ax) {
        var lit = f.highlight && f.highlight.some(function (t) {
          return t[0] === k && t[1] === ax[2];
        });
        var e = sp.px(ax[0]), e2 = sp.px([-ax[0][0], -ax[0][1], -ax[0][2]]);
        ctx.save();
        ctx.strokeStyle = lit ? MV.palette.critical : MV.palette.line;
        ctx.lineWidth = lit ? 2.2 : 1;
        ctx.beginPath();
        ctx.moveTo(e2[0], e2[1]);
        ctx.lineTo(e[0], e[1]);
        ctx.stroke();
        ctx.restore();
        var dx = e[0] - cx, dy = e[1] - cy, len = Math.hypot(dx, dy) || 1;
        p.label(ax[1], e[0] + dx / len * 10, e[1] + dy / len * 10, {
          color: lit ? MV.palette.critical : MV.palette.muted,
          align: 'center', baseline: 'middle',
          font: (lit ? '600 11px ' : '10px ') + 'system-ui'
        });
      });

      /* Das Ziel, von Anfang an sichtbar. */
      if (isOut) {
        var tr = targetBloch(proto, k);
        var tl = Math.hypot(tr[0], tr[1], tr[2]);
        var tp = tl < 1e-9 ? [cx, cy] : sp.px(tr);
        ctx.save();
        ctx.strokeStyle = MV.palette.curve3;
        ctx.lineWidth = 1.8;
        ctx.beginPath();
        ctx.arc(tp[0], tp[1], 6.5, 0, 2 * PI);
        ctx.stroke();
        ctx.restore();
        p.label('Ziel', tp[0], tp[1] - 10, {
          color: MV.palette.curve3, align: 'center', baseline: 'bottom',
          font: '10px system-ui', box: true
        });
      }

      /* Der Schatten des Zustands vor dem Sprung. */
      if (f.shadow) {
        var sr = f.shadow[k];
        if (Math.hypot(sr[0], sr[1], sr[2]) > 1e-9) {
          var se = sp.px(sr);
          ctx.save();
          ctx.strokeStyle = MV.palette.muted;
          ctx.globalAlpha = 0.45;
          ctx.lineWidth = 1.8;
          ctx.setLineDash([4, 4]);
          ctx.beginPath();
          ctx.moveTo(cx, cy);
          ctx.lineTo(se[0], se[1]);
          ctx.stroke();
          ctx.restore();
        }
      }

      /* Die Bögen der Korrektur. */
      f.arcs.filter(function (a) { return a.qubit === k; }).forEach(function (a) {
        var steps = 48, i2, r2, q2, wasFront = null, ff;
        ctx.save();
        ctx.lineWidth = 2.2;
        ctx.lineCap = 'round';
        for (i2 = 0; i2 <= steps; i2++) {
          r2 = rotate(a.from, a.n, a.alpha * a.upto * i2 / steps);
          q2 = sp.px(r2);
          ff = sp.front(r2);
          if (ff !== wasFront) {
            if (wasFront !== null) ctx.stroke();
            ctx.beginPath();
            ctx.strokeStyle = MV.palette.critical;
            ctx.globalAlpha = ff ? 0.9 : 0.35;
            ctx.setLineDash(ff ? [] : [3, 3]);
            ctx.moveTo(q2[0], q2[1]);
            wasFront = ff;
          } else {
            ctx.lineTo(q2[0], q2[1]);
          }
        }
        ctx.stroke();
        ctx.restore();
      });

      /* Der Zustandsvektor. */
      var r = f.vecs[k];
      var L = Math.hypot(r[0], r[1], r[2]);
      if (L > 1e-9) {
        ctx.save();
        ctx.globalAlpha = f.fade;
        drawVector(ctx, sp, cx, cy, r, isOut ? MV.palette.topological : MV.palette.curve2,
                   isOut ? 2.6 : 2.0);
        ctx.restore();
      } else {
        ctx.save();
        ctx.fillStyle = MV.palette.muted;
        ctx.beginPath();
        ctx.arc(cx, cy, 3, 0, 2 * PI);
        ctx.fill();
        ctx.restore();
      }

      p.label('Qubit ' + roleLabel(role) + (isOut ? ' — Ausgabe' : ''), cx, cy - R - 16, {
        color: MV.palette.ink, align: 'center', baseline: 'bottom', font: '13px system-ui'
      });
      p.label('Länge ' + MV.fmt(L, 3) +
              (L < 1 - 1e-9 ? ' — verschränkt' : ' — reiner Zustand'),
              cx, cy + R + 16, {
        color: MV.palette.muted, align: 'center', baseline: 'top', font: '11px system-ui'
      });
    });

    p.label(f.caption, (p.px0 + p.px1) / 2, p.py1 + 34, {
      color: MV.palette.muted, align: 'center', baseline: 'bottom', font: '11px system-ui'
    });

    p.finish();

    var lastS = lastStep(run, proto);
    var atEnd = step >= lastS && play.u >= 1;
    panelB.setStatus(
      (step === 0 ? 'Präparation'
        : step > proto.steps.length ? 'Korrekturschritt'
        : 'Messschritt ' + step + ' von ' + proto.steps.length) +
      (play.running ? ' — läuft' : '') +
      (atEnd
        ? (run.impossible >= 0
            ? ' — Wahrscheinlichkeit null, die Wiedergabe hält hier an.'
            : ' — Abweichung vom Ziel: ' + devText(run.dev) +
              (run.dev < TOL ? ', das Gatter ist erzeugt.' : ', das Ziel ist nicht erreicht.'))
        : ''),
      atEnd && run.impossible < 0
        ? (run.dev < TOL ? MV.palette.topological : MV.palette.critical)
        : MV.palette.muted);
  }

  /* ------------------------------------------------------------- Leiter */

  function renderLadder() {
    var proto = PROTOCOLS[state.proto];
    var run = currentRun();
    var step = Math.min(state.step, lastStep(run, proto));
    var u = drawU();
    var items = [];

    items.push({ t: proto.prepText, r: proto.prepShort, cur: step === 0 });

    proto.steps.forEach(function (st, i) {
      var res = '—';
      var reached = step > i + 1 || (step === i + 1 && u >= U_RESULT);
      if (reached && i < run.results.length) {
        res = st.sym + ' = ' + MV.fmtSigned(run.results[i], 0) +
              ', P = ' + MV.fmt(run.probs[i], 3) +
              (run.drawn[i] ? '' : ' (von Hand)');
      } else if (step === i + 1) {
        res = 'P(+1) = ' + MV.fmt(run.probPlus[i], 3);
      }
      items.push({ t: st.text, r: res, cur: step === i + 1 });
    });

    var corrDone = step > proto.steps.length && run.impossible < 0;
    items.push({
      t: proto.corrText,
      r: !corrDone ? '—'
        : !state.correct ? 'weggelassen'
        : run.applied ? run.applied + ' angewandt' : 'nichts zu tun',
      cur: step === proto.steps.length + 1
    });

    ladderHost.innerHTML = '<ol class="ladder">' + items.map(function (it, i) {
      return '<li class="ladder__item' + (it.cur ? ' is-current' : '') + '">' +
        '<span class="ladder__n">' + (i + 1) + '</span>' +
        '<span class="ladder__t">' + it.t + '</span>' +
        '<span class="ladder__r">' + it.r + '</span></li>';
    }).join('') + '</ol>';
  }

  function updateReadoutB() {
    var proto = PROTOCOLS[state.proto];
    var run = currentRun();
    var step = Math.min(state.step, lastStep(run, proto));
    var u = drawU();
    var atEnd = step >= lastStep(run, proto) && play.u >= 1;

    readoutB.set('proto', proto.name);
    readoutB.set('out', 'Qubit ' + (proto.outName === 'C und T'
      ? 'C und T (Ziel)' : roleLabel(proto.outName)));

    var shown = run.results.filter(function (b, i) {
      return step > i + 1 || (step === i + 1 && u >= U_RESULT);
    }).map(function (b, i) {
      return proto.steps[i].sym + ' = ' + MV.fmtSigned(b, 0);
    });
    readoutB.set('results', shown.join(', ') || '—');

    readoutB.set('dev', !atEnd ? '—'
      : run.impossible >= 0 ? 'unmöglicher Ausgang' : devText(run.dev));
    /* Die Schalterstellung allein sagt nicht, ob in diesem Zweig überhaupt
       etwas zu tun war. Die Liste liefert correct()/applyOps ohnehin. */
    readoutB.set('corr', !state.correct ? 'weggelassen'
      : !atEnd ? 'an'
      : 'an — angewandt: ' + (run.applied || 'keine'));

    if (!atEnd) {
      readoutB.status('trivial', 'Lauf noch nicht zu Ende — die Abweichung steht am Schluss');
    } else if (run.impossible >= 0) {
      readoutB.status('trivial', 'Dieser Ausgang hat Wahrscheinlichkeit null');
    } else if (run.dev < TOL) {
      readoutB.status('topological', 'Das Gatter ' + proto.name + ' ist erzeugt');
    } else {
      readoutB.status('critical', 'Das Ziel ist nicht erreicht — ' +
        (state.correct ? 'das sollte nicht vorkommen' : 'die Korrektur fehlt'));
    }
  }

  /* --------------------------------------------------------- Wiedergabe */

  function schedule() {
    if (play.raf) return;
    play.raf = window.requestAnimationFrame(tick);
  }

  function tick(ts) {
    play.raf = 0;
    if (!play.running) return;
    if (!play.last) play.last = ts;
    var dt = Math.min(0.25, (ts - play.last) / 1000);
    play.last = ts;

    var proto = PROTOCOLS[state.proto];
    var run = currentRun();
    var lastS = lastStep(run, proto);

    play.u += dt / Math.max(play.tempo, 1e-6);
    if (play.u >= 1) {
      if (state.step >= lastS) {
        play.u = 1;
        play.running = false;
      } else {
        state.step += 1;
        play.u = 0;
      }
    }
    refreshB();
    if (play.running) schedule();
  }

  function pause() {
    play.running = false;
    if (play.raf) { window.cancelAnimationFrame(play.raf); play.raf = 0; }
    play.last = 0;
  }

  /** Neuer Lauf: Anzeige auf Start, für gezogene Messwerte neue Ergebnisse. */
  function newRun() {
    pause();
    state.step = 0;
    play.u = 1;
    play.counted = false;
    if (state.outcomes.some(function (o) { return o === null; })) {
      state.seed = (state.seed * 1103515245 + 12345) & 0x7fffffff;
    }
  }

  function resetCounter() {
    var proto = PROTOCOLS[state.proto];
    var run = currentRun();
    play.runs = 0;
    play.hits = 0;
    play.counted = state.step >= lastStep(run, proto) && play.u >= 1;
  }

  function maybeCount() {
    if (play.counted) return;
    var proto = PROTOCOLS[state.proto];
    var run = currentRun();
    if (state.step < lastStep(run, proto) || play.u < 1) return;
    play.counted = true;
    play.runs += 1;
    if (run.impossible < 0 && run.dev < TOL) play.hits += 1;
  }

  function renderCounter() {
    if (!counterHost) return;
    counterHost.textContent = 'Läufe: ' + play.runs + ' · Ziel erreicht: ' + play.hits;
  }

  function onStartPause() {
    var proto = PROTOCOLS[state.proto];
    var run = currentRun();
    if (play.running) { pause(); refreshB(); return; }
    if (state.step >= lastStep(run, proto) && play.u >= 1) newRun();
    play.running = true;
    play.last = 0;
    schedule();
    refreshB();
  }

  function stepForward() {
    var proto = PROTOCOLS[state.proto];
    var run = currentRun();
    pause();
    if (play.u < 1) play.u = 1;
    else if (state.step < lastStep(run, proto)) { state.step += 1; play.u = 1; }
    refreshB();
  }

  function stepBack() {
    pause();
    if (play.u < 1) play.u = 1;
    else if (state.step > 0) state.step -= 1;
    play.u = 1;
    refreshB();
  }

  function toStart() {
    pause();
    state.step = 0;
    play.u = 1;
    refreshB();
  }

  /** Die Wiedergabeleiste — direkt unter der Grafik, die sie abspielt. */
  function buildPlayBar(mount) {
    var bar = MV.ui.el('div', { class: 'playbar' });

    var row = MV.ui.el('div', { class: 'playbar__row' });
    var g = MV.ui.buttonRow(row, {
      buttons: [
        { text: 'Start', onClick: onStartPause },
        { text: 'Schritt zurück', onClick: stepBack },
        { text: 'Schritt vor', onClick: stepForward },
        { text: 'Anfang', onClick: toStart }
      ]
    });
    playButton = g.buttons[0];
    bar.appendChild(row);

    MV.ui.slider(bar, {
      label: 'Tempo — Sekunden je Schritt',
      min: 0.5, max: 4, step: 0.1, value: play.tempo, digits: 1,
      hint: 'Anzeigetempo — keine physikalische Zeit.',
      onInput: function (v) { play.tempo = v; }
    });

    counterHost = MV.ui.el('div', { class: 'counter' });
    bar.appendChild(counterHost);

    mount.appendChild(bar);
  }

  /* Die Selbstprüfung steht hier, weil sie die Wiedergabe mitprüft und dafür
     den Zustand aus Abschnitt 10 braucht. Ihr Ergebnis geht in Punkt 5 der
     Info-Panels, die weiter unten gebaut werden. */
  var check = runSelfCheck();

  /* ====================================================================
     11 — Aufbau des Reiters
     ==================================================================== */

  function refreshA() {
    drawA();
    updateReadoutA();
    if (panelA.info) panelA.info.update();
  }

  function refreshB() {
    var proto = PROTOCOLS[state.proto];
    var run = currentRun();
    if (state.step > lastStep(run, proto)) state.step = lastStep(run, proto);
    maybeCount();
    drawB();
    renderLadder();
    updateReadoutB();
    renderCounter();
    if (playButton) playButton.textContent = play.running ? 'Pause' : 'Start';
    if (panelB.info) panelB.info.update();
  }

  function refreshAll() { refreshA(); refreshB(); }

  /**
   * Eine Änderung am Eingang hält eine laufende Wiedergabe an — der
   * angezeigte Schritt bleibt stehen — und setzt den Zähler zurück: die
   * bisherigen Läufe gehörten zu einem anderen Eingangszustand.
   */
  function onInputChanged() {
    pause();
    resetCounter();
    refreshAll();
  }

  /* --------------------------------------------------------- Einstellungen

     Knopf und Vorführung (demo.js) gehen durch dieselbe Funktion; sonst
     stünde die Beschriftung irgendwann anders als der Zustand. */

  function setStart(theta, phi) {
    state.theta = theta;
    state.phi = phi;
    sliders.theta.set(theta, true);
    sliders.phi.set(phi, true);
    onInputChanged();
  }

  function setCorrect(on) {
    pause();
    state.correct = !!on;
    if (corrButtons) {
      corrButtons[0].setAttribute('aria-pressed', state.correct ? 'true' : 'false');
      corrButtons[1].setAttribute('aria-pressed', state.correct ? 'false' : 'true');
    }
    resetCounter();
    refreshB();
  }

  function setTarget(one) {
    state.target1 = !!one;
    if (targetButtons) {
      targetButtons[0].setAttribute('aria-pressed', state.target1 ? 'false' : 'true');
      targetButtons[1].setAttribute('aria-pressed', state.target1 ? 'true' : 'false');
    }
    onInputChanged();
  }

  function setOutcomes(list) {
    state.outcomes = list.slice();
    newRun();
    buildOutcomeRow();
    refreshB();
  }

  /**
   * Die Knöpfe für die Messergebnisse. Einmal je Protokoll gebaut, nicht bei
   * jeder Auffrischung: sonst würde ein Knopf mitten in der Behandlung seines
   * eigenen Klicks aus dem Dokument entfernt, und jedes Ziehen eines Reglers
   * risse die halbe Leiste neu auf.
   */
  function buildOutcomeRow() {
    var proto = PROTOCOLS[state.proto];
    outcomeRow.innerHTML = '';
    proto.steps.forEach(function (st, i) {
      MV.ui.segmented(outcomeRow, {
        label: st.sym + ':',
        options: ['gezogen', '+1', '−1'],
        value: state.outcomes[i] === null ? 0 : (state.outcomes[i] === 1 ? 1 : 2),
        onChange: function (k) {
          pause();
          state.outcomes[i] = k === 0 ? null : (k === 1 ? 1 : -1);
          play.counted = false;
          refreshB();
        }
      });
    });
  }

  /* ------------------------------------------------------------ Info-Panels */

  function numericsHtml() {
    return '<ul>' + check.parts.map(function (p) {
      var mark = p.ok === true ? '✓' : p.ok === false ? '✗' : '!';
      return '<li><b>' + mark + '</b> ' + p.name + ' — ' + p.detail + '</li>';
    }).join('') + '</ul>' +
      '<p><b>Zur Größenordnung:</b> Werte um 10⁻¹⁶ sind Rundung der ' +
      'Gleitkommarechnung und kein physikalischer Unterschied — die Anzeige ' +
      'schreibt dort „0" und nennt den wahren Wert in Klammern. Ein Zweig, in dem ' +
      'das Protokoll das Ziel wirklich verfehlt, liegt dagegen in der ' +
      'Größenordnung <b>1</b>. Dazwischen gibt es nichts: entweder das Gatter ' +
      'entsteht exakt, oder es entsteht gar nicht.</p>' +
      '<p>Toleranz ' + MV.fmtExp(TOL, 0) + '. Die Prüfung läuft beim Laden der ' +
      'Seite, nicht als Vorberechnung: schlägt sie fehl, steht das über den ' +
      'Abbildungen. Ohne die Korrektur ist das Ergebnis in ' + noCorrText() +
      ' Zweigen falsch — gezählt an einem Startzustand, der Eigenzustand keiner ' +
      'der drei Achsen ist.</p>';
  }

  /** "Ŝ_z 2 von 4, …" — aus der Selbstprüfung, nicht hingeschrieben. */
  function noCorrText() {
    return PROTOCOLS.map(function (p) {
      var c = check.noCorr[p.key];
      return p.name + ' ' + c.wrong + ' von ' + c.total;
    }).join(', ');
  }

  var PHI_NOTE =
    'θ und φ sind in diesem Reiter die Winkel des Startzustands auf der ' +
    'Bloch-Kugel — weder die supraleitende Phase aus Kapitel 3 noch der ' +
    'eingeschlossene Fluss aus Kapitel 5. Diese Bedeutung stammt aus der ' +
    'Darstellung, nicht aus der Arbeit.';

  var NO_TIME_NOTE =
    'Ein Gatter wird hier angewandt, nicht gefahren: es gibt keine Rate, ' +
    'keine Pulsform und kein Zeitintervall. Der gezeichnete Bogen ist der Weg, ' +
    'den die Drehung nimmt, und keine Zeitentwicklung. Eine Messung ist im ' +
    'Modell eine ideale Projektion; ihre Dauer τ_meas kommt hier nicht vor. ' +
    'Das Tempo der Wiedergabe ist frei gewählt und steht für nichts.';

  var INFO_A = {
    what: 'Die Bloch-Kugel eines Qubits im Sektor P̂ = +1. Die sechs ' +
          'Eigenzustände von x̂, ŷ und ẑ nach ' + MV.eq('tetron.pauli') + ' sind ' +
          'die Ecken eines Oktaeders; die Kanten sind eingezeichnet. Jeder ' +
          'Knopf wendet ein Gatter an und zeichnet den Bogen, den der Zustand ' +
          'dabei beschreibt.',
    model: 'Zu jedem Gatter werden Drehachse n̂ und Drehwinkel α aus der ' +
           'Matrix bestimmt: die globale Phase wird über die Determinante ' +
           'entfernt, die verbleibende SU(2)-Matrix hat die Form ' +
           'cos(α/2)·1 − i sin(α/2) n̂·σ̂. Der Bloch-Vektor dreht sich dann um ' +
           'n̂ um α. Ŝ_y und Ĥ werden aus Ŝ_z und Ŝ_x zusammengesetzt, nicht ' +
           'eigens hingeschrieben.',
    computed: 'Bloch-Vektor und seine Länge, θ und φ des aktuellen Zustands, ' +
              'ob er auf einer Ecke des Oktaeders steht, sowie für jedes ' +
              'Gatter die Tabelle aus Matrix, Drehachse, Winkel und Bild der ' +
              'sechs Ecken.',
    params: function () {
      var w = path();
      var nv = nearestVertex(w.r);
      return '<ul><li>Startzustand: θ = ' + MV.fmt(state.theta, 0) + '°, φ = ' +
        MV.fmt(state.phi, 0) + '°</li>' +
        '<li>angewandte Gatter: ' + (state.seq.length
          ? state.seq.map(function (k) { return GATE_BY_KEY[k].label; }).join(' → ')
          : 'keines') + '</li>' +
        '<li>Bloch-Vektor: (' + MV.fmt(w.r[0], 4) + '; ' + MV.fmt(w.r[1], 4) + '; ' +
        MV.fmt(w.r[2], 4) + '), Länge ' + MV.fmt(Math.hypot(w.r[0], w.r[1], w.r[2]), 6) + '</li>' +
        '<li>auf einer Ecke: ' + (nv.d < 1e-9 ? 'ja, ' + nv.v.label : 'nein, Abstand ' +
          MV.fmt(nv.d, 4) + ' zur nächsten') + '</li>' +
        '<li>Blickwinkel: Azimut ' + MV.fmt(cam.az * DEG, 0) + '°, Höhe ' +
        MV.fmt(cam.el * DEG, 0) + '°</li></ul>';
    },
    numerics: numericsHtml(),
    convention: 'Basis |0⟩ ist ẑ = +1, |1⟩ ist ẑ = −1, wie in ' +
                MV.eq('tetron.pauli') + '. ' + PHI_NOTE + ' ' + NO_TIME_NOTE +
                ' Die Projektion ist orthografisch und orthonormal, der Umriss ' +
                'ist also genau der Ort aller reinen Zustände; ein Vektor kann ' +
                'nie darüber hinausragen, erscheint aber verkürzt, wenn er auf ' +
                'den Betrachter zu zeigt. Dreh die Kugel, dann ist er wieder lang.',
    reference: MV.secOfWork('5.4') + ' und ' + MV.sec('5.5') + '; ' +
               MV.eq('ops.braid') + ', ' + MV.eq('ops.braidZ') + ' und ' +
               MV.eq('ops.sz') + '. Die Vertauschung selbst wird hier nicht ' +
               'gezeigt — sie steht als ' + MV.fig('ops.braidFig') + ' in der ' +
               'Arbeit; dieser Reiter setzt hinter ihrem Ergebnis an. ' +
               'Die Zusammenstellung der neun Gatter und ' +
               'das Oktaeder-Bild sind Standard der Quanteninformation (' + NC +
               ') und stehen so nicht in der Arbeit; sie sind die Anschauung zu ' +
               'dem Satz über Gottesman-Knill in ' + MV.sec('5.4') + ', kein Beweis.',
    reading: '<p><b>Setze den Startzustand auf x̂ = +1 und drücke ' +
             '„T-Gatter“:</b> die Zeile <b>Ecke des Oktaeders</b> springt von ' +
             '<i>ja</i> auf <i>nein</i>, und der Bloch-Vektor steht auf ' +
             '(0,7071; 0,7071; 0) — mitten zwischen zwei Ecken. Drück ' +
             '„T-Gatter“ ein zweites Mal: jetzt steht dort wieder <i>ja</i>, ' +
             'nämlich ŷ = +1, denn zwei T-Gatter ergeben Ŝ_z. Mit jedem der ' +
             'sieben Clifford-Knöpfe bleibt die Zeile dagegen immer auf ' +
             '<i>ja</i>, ganz gleich wie oft man drückt.</p>'
  };

  var INFO_B = {
    what: 'Ein Durchlauf des gewählten Protokolls, Schritt für Schritt ' +
          'wiedergegeben. Für jedes Qubit des Registers steht eine Bloch-Kugel; ' +
          'gezeichnet ist der Vektor der reduzierten Dichtematrix. Ist ein Qubit ' +
          'mit den übrigen verschränkt, wird sein Vektor kürzer als 1 — es hat ' +
          'dann für sich genommen keinen Zustand. Auf jeder Ausgabe-Kugel steht ' +
          'von Anfang an ein Ring an der Spitze des Ziel-Vektors; hat das Ziel ' +
          'die Länge null, steht der Ring im Mittelpunkt. Die Wiedergabe läuft ' +
          'nur auf „Start" und hält am Ende des Laufs an.',
    model: 'Jede Messung ist eine Projektion P_b = (1 + b·Ô)/2 mit ' +
           'anschließender Normierung; das Normquadrat davor ist die ' +
           'Wahrscheinlichkeit des Ergebnisses. <b>Eine Messung ist deshalb ein ' +
           'Sprung:</b> zwischen dem Zustand davor und danach gibt es im Modell ' +
           'keinen Weg, also wird auch keiner gezeichnet — kein Gleiten, kein ' +
           'Bogen, keine Linie zwischen alter und neuer Spitze. Der alte Vektor ' +
           'bleibt als blasser, gestrichelter Schatten stehen, bis der nächste ' +
           'Schritt beginnt. <b>Die Korrektur dagegen ist eine Drehung:</b> ein ' +
           'unitäres Gatter, das als Bogen um Achse und Winkel des Operators ' +
           'läuft, beim T-Protokoll nacheinander für x̂, ẑ und Ŝ_z. Die drei ' +
           'ersten Protokolle stehen in ' + MV.sec('5.5') + '; das vierte ist ' +
           'die Standardkonstruktion über einen magischen Zustand und wird in ' +
           'der Arbeit nicht ausgeführt.' +
           '<p><b>Wie eine gemeinsame Parität physikalisch gemessen wird</b> ' +
           '(Plugge et al., Fig. 3(c) und Anhang C — nicht in der Arbeit, als ' +
           'Quelle gekennzeichnet): Zwei Quantenpunkte koppeln an beide Tetrons ' +
           'zugleich, und die Referenzarme werden abgeschaltet. Der Übergang ' +
           'zwischen den beiden Punkten hat dann die Amplitude ' +
           't̂ = t_a ẑ_a + t_b ẑ_b — einen Weg durch die eine Box, einen durch ' +
           'die andere. Die Rabi-Frequenz hängt davon nur über das Produkt ' +
           'z_a z_b ab, also über die <b>gemeinsame</b> Parität und nicht über ' +
           'die einzelnen. Jede Box dient dabei als Referenzweg für die andere: ' +
           'der Referenzarm aus dem Einzelqubit-Fall wird nicht gebraucht, weil ' +
           'die zweite Box seine Rolle übernimmt.</p>',
    computed: 'Wahrscheinlichkeit jedes Messergebnisses, der Registerzustand ' +
              'nach jedem Schritt, die Länge jedes Bloch-Vektors und die ' +
              'Abweichung des Ausgabe-Teilsystems vom Ziel. Verglichen werden ' +
              'reduzierte Dichtematrizen und nicht Zustandsvektoren: ρ ist ' +
              'gegen die globale Phase unempfindlich, sodass keine ' +
              'Phasenanpassung nötig ist, die einen echten Unterschied ' +
              'verdecken könnte.',
    params: function () {
      var proto = PROTOCOLS[state.proto];
      var run = currentRun();
      var atEnd = state.step >= lastStep(run, proto) && play.u >= 1;
      return '<ul><li>Protokoll: ' + proto.name + ', Register ' +
        proto.roles.map(roleLabel).join(', ') + '</li>' +
        '<li>Eingangszustand: θ = ' + MV.fmt(state.theta, 0) + '°, φ = ' +
        MV.fmt(state.phi, 0) + '°' +
        (proto.inQubits === 2 ? ', Qubit T (Ziel) in |' + (state.target1 ? 1 : 0) + '⟩' : '') +
        '</li>' +
        '<li>Messergebnisse: ' + (run.results.filter(function (b, i) {
          return state.step > i + 1 || (state.step === i + 1 && drawU() >= U_RESULT);
        }).map(function (b, i) {
          return proto.steps[i].sym + ' = ' + MV.fmtSigned(b, 0) +
            ' (' + (run.drawn[i] ? 'gezogen' : 'von Hand') + ')';
        }).join(', ') || '—') + '</li>' +
        '<li>Korrektur: ' + (state.correct ? 'an' : 'weggelassen') + '</li>' +
        '<li>angezeigter Schritt: ' + state.step + ' von ' +
        lastStep(run, proto) + (play.running ? ', läuft' : ', angehalten') + '</li>' +
        '<li>Anzeigetempo: ' + MV.fmt(play.tempo, 1) + ' s je Schritt — ' +
        'keine physikalische Zeit</li>' +
        '<li>Zähler: ' + play.runs + ' Läufe, davon ' + play.hits +
        ' am Ziel</li>' +
        '<li>Abweichung vom Ziel: ' + (atEnd
          ? (run.impossible >= 0 ? 'unmöglicher Ausgang' : devText(run.dev))
          : 'steht erst am Ende des Laufs') +
        '</li></ul>';
    },
    numerics: numericsHtml(),
    convention: PHI_NOTE + ' Das Protokoll für das T-Gatter ist keine Aussage ' +
                'der Arbeit, sondern die Standardkonstruktion nach ' + NC + '; ' +
                MV.sec('5.5') + ' nennt den Weg über magische Zustände und hält ' +
                'fest, dass er nicht topologisch geschützt ist. ' + NO_TIME_NOTE +
                ' <b>Das Tempo ist keine Zeit:</b> es bestimmt nur, wie lange ' +
                'ein Schritt stehen bleibt. Auch der Bogen der Korrektur ist ' +
                'keine Zeitentwicklung — er zeigt, was der Operator am Zustand ' +
                'bewirkt; am Bauelement ist die Korrektur ein bestätigter ' +
                'Elektronentransfer oder reine Buchführung. Bei eingeschalteter ' +
                'Bewegungsreduzierung des Systems läuft kein Bogen und blendet ' +
                'nichts über: jeder Schritt zeigt sofort seinen Endzustand.',
    reference: MV.secOfWork('5.5') + ', ' + MV.eq('ops.sz') + '; die Protokolle ' +
               'nach Plugge et al., Fig. 5(a) bis (c). Der Gedanke geht auf ' +
               'Bonderson, Freedman und Nayak zurück.',
    reading: '<p><b>Wähle Ŝ_z, setz den Startzustand auf x̂ = +1 und drück ' +
             'Start:</b> Nach der ersten Messung schrumpfen beide Vektoren auf ' +
             'die Länge 0, nach der zweiten steht C auf ŷ = +1 oder ŷ = −1 — nur ' +
             'im zweiten Fall folgt die Korrektur, als einziger Bogen. Drück ' +
             'Start noch einmal: andere Messwerte, dasselbe Ziel. Schalt dann die ' +
             'Korrektur ab und lass zehn Läufe laufen: Der Zähler zeigt, dass nur ' +
             'etwa jeder zweite das Ziel erreicht.</p>'
  };

  /* --------------------------------------------------------------- Aufbau */

  function buildControls(mount) {
    var box = MV.ui.el('div', { class: 'controls' });

    box.appendChild(MV.ui.el('h4', { class: 'controls__head', text: 'Startzustand' }));

    sliders.theta = MV.ui.slider(box, {
      label: 'θ — Polarwinkel',
      min: 0, max: 180, step: 1, value: state.theta, digits: 0, prefix: '',
      onInput: function (v) { state.theta = v; onInputChanged(); }
    });
    sliders.phi = MV.ui.slider(box, {
      label: 'φ — Azimutwinkel',
      min: -180, max: 180, step: 1, value: state.phi, digits: 0,
      hint: 'Winkel auf der Bloch-Kugel, nicht der Fluss aus Kapitel 5.',
      onInput: function (v) { state.phi = v; onInputChanged(); }
    });

    var vrow = MV.ui.el('div', { class: 'controls__wide' });
    MV.ui.buttonRow(vrow, {
      label: 'auf eine Ecke setzen:',
      buttons: VERTICES.map(function (v) {
        return {
          text: v.label,
          onClick: function () {
            state.theta = Math.round(v.theta * DEG);
            state.phi = Math.round(v.phi * DEG);
            sliders.theta.set(state.theta, true);
            sliders.phi.set(state.phi, true);
            onInputChanged();
          }
        };
      })
    });
    box.appendChild(vrow);

    box.appendChild(MV.ui.el('h4', { class: 'controls__head', text: 'Gatter anwenden' }));

    var grow = MV.ui.el('div', { class: 'controls__wide' });
    MV.ui.buttonRow(grow, {
      buttons: GATES.map(function (g) {
        return {
          text: g.label,
          onClick: function () { state.seq.push(g.key); refreshA(); }
        };
      }).concat([
        { text: 'letztes zurück', onClick: function () { state.seq.pop(); refreshA(); } },
        { text: 'Zurücksetzen', onClick: function () { state.seq = []; refreshA(); } }
      ])
    });
    box.appendChild(grow);

    seqHost = MV.ui.el('div', { class: 'gateseq' });
    box.appendChild(seqHost);

    box.appendChild(MV.ui.el('p', {
      class: 'note',
      html: 'Die Kugel lässt sich mit der Maus drehen. Der Startzustand wirkt ' +
            'auf den ganzen Weg zurück: zieht man θ oder φ, wandert die ' +
            'angewandte Folge mit.'
    }));

    box.appendChild(MV.ui.el('h4', { class: 'controls__head', text: 'Protokoll in Panel B' }));

    segProto = MV.ui.segmented(box, {
      label: 'Protokoll:',
      options: PROTOCOLS.map(function (p) { return p.name; }),
      value: state.proto,
      onChange: function (i) {
        state.proto = i;
        state.outcomes = [null, null, null];
        resetCounter();
        newRun();
        buildOutcomeRow();
        buildStepSeg();
        refreshB();
      }
    });

    outcomeRow = MV.ui.el('div', { class: 'controls__wide' });
    box.appendChild(outcomeRow);

    var srow = MV.ui.el('div', { class: 'controls__wide' });
    corrButtons = MV.ui.buttonRow(srow, {
      label: 'Korrektur:',
      buttons: [
        { text: 'anwenden', pressed: true, onClick: function () { setCorrect(true); } },
        { text: 'weglassen', pressed: false, onClick: function () { setCorrect(false); } }
      ]
    }).buttons;
    box.appendChild(srow);

    var trow = MV.ui.el('div', { class: 'controls__wide' });
    MV.ui.buttonRow(trow, {
      label: 'neu ziehen:',
      buttons: [{
        text: 'andere Messergebnisse',
        onClick: function () {
          state.seed = (state.seed * 1103515245 + 12345) & 0x7fffffff;
          newRun();
          refreshB();
        }
      }]
    });
    box.appendChild(trow);

    mount.appendChild(box);
  }

  function buildStepSeg() {
    var proto = PROTOCOLS[state.proto];
    var opts = ['Start'].concat(proto.steps.map(function (s, i) {
      return (i + 1) + ' ' + s.sym;
    })).concat(['Korrektur']);
    stepSeg.innerHTML = '';
    MV.ui.segmented(stepSeg, {
      label: 'angezeigter Schritt:',
      options: opts,
      value: Math.min(state.step, opts.length - 1),
      onChange: function (i) { pause(); state.step = i; play.u = 1; refreshB(); }
    });
  }

  /* ====================================================================
     Öffentliche Schnittstelle
     ==================================================================== */

  MV.gatter = {
    check: check,

    init: function (mount) {
      if (!check.passed) {
        var bad = check.parts.filter(function (p) { return p.ok === false; });
        var warn = MV.ui.el('div', {
          class: 'readout',
          html: '<b>Warnung:</b> ' + bad.length + ' Prüfung(en) dieses Reiters sind ' +
                'fehlgeschlagen: ' + bad.map(function (p) { return p.name; }).join('; ') +
                '. Die Abbildungen sind nicht belastbar.'
        });
        warn.style.borderColor = MV.palette.critical;
        warn.style.color = MV.palette.critical;
        mount.appendChild(warn);
      }

      buildControls(mount);

      var panels = MV.ui.el('div', { class: 'panels' });
      mount.appendChild(panels);

      panelA = MV.ui.panel(panels, {
        title: 'Panel A — Die Gatter auf der Bloch-Kugel',
        quote: 'gatter.A',
        lead: 'Die sechs Eigenzustände von x̂, ŷ und ẑ spannen ein Oktaeder auf. ' +
              '<b>Drücke die Gatter-Knöpfe und sieh zu, wo der Zeiger landet:</b> ' +
              'die sieben Clifford-Gatter bilden Ecken auf Ecken ab, die beiden ' +
              'T-Gatter nicht. Die Kugel lässt sich mit der Maus drehen.',
        wide: true, aspect: 0.56, controls: true,
        margin: { l: 16, r: 16, t: 14, b: 14 },
        info: INFO_A
      });

      readoutA = MV.ui.readout(panelA.controls, {
        items: [
          { key: 'rx', label: 'Bloch x' },
          { key: 'ry', label: 'Bloch y' },
          { key: 'rz', label: 'Bloch z' },
          { key: 'len', label: 'Länge' },
          { key: 'theta', label: 'θ' },
          { key: 'phi', label: 'φ' },
          { key: 'vertex', label: 'Ecke des Oktaeders' },
          { key: 'count', label: 'angewandte Gatter' }
        ]
      });

      gateTableHost = MV.ui.el('div', { class: 'gatetable-host' });
      panelA.controls.appendChild(gateTableHost);

      attachDrag(panelA.canvas);

      panelB = MV.ui.panel(panels, {
        title: 'Panel B — Dieselben Gatter aus Messungen',
        quote: 'gatter.B',
        lead: 'Kein Film, sondern eine Schrittfolge: Jede Messung ist ein ' +
              'Sprung, nur die Korrektur ist eine Drehung. <b>Drück Start</b> — ' +
              'das Tempo regelt die Anzeige, nicht die Physik.',
        wide: true, aspect: 0.46, controls: true,
        margin: { l: 16, r: 16, t: 34, b: 48 },
        info: INFO_B
      });

      /* Die Wiedergabeleiste sitzt direkt unter der Grafik, die sie
         abspielt — die einzige Ausnahme vom Links-Layout (CLAUDE.md §12). */
      buildPlayBar(panelB.controls);

      stepSeg = MV.ui.el('div');
      panelB.controls.appendChild(stepSeg);

      ladderHost = MV.ui.el('div', { class: 'ladder-host' });
      panelB.controls.appendChild(ladderHost);

      readoutB = MV.ui.readout(panelB.controls, {
        items: [
          { key: 'proto', label: 'Protokoll' },
          { key: 'out', label: 'Ausgabe auf' },
          { key: 'results', label: 'Messergebnisse' },
          { key: 'dev', label: 'Abweichung vom Ziel' },
          { key: 'corr', label: 'Korrektur' }
        ]
      });

      panelB.controls.appendChild(MV.ui.el('p', {
        class: 'note',
        html: 'Beim CNOT ist der Eingang zweifach: Qubit C trägt den Zustand aus ' +
              'θ und φ, Qubit T startet in |0⟩ oder |1⟩. Mit C auf x̂ = +1 und ' +
              'T in |0⟩ entsteht ein Bell-Zustand — beide Bloch-Vektoren ' +
              'schrumpfen dann auf null.'
      }));

      targetButtons = MV.ui.buttonRow(panelB.controls, {
        label: 'Qubit T (Ziel) — nur CNOT:',
        buttons: [
          { text: '|0⟩', pressed: true, onClick: function () { setTarget(false); } },
          { text: '|1⟩', pressed: false, onClick: function () { setTarget(true); } }
        ]
      }).buttons;

      /* Schnittstelle für die Vorführungen (demo.js). Panel B fährt die
         Wiedergabe, die dort schon eingebaut ist — die Vorführung drückt nur
         auf Start und wartet, bis sie von selbst am Ende steht. */
      MV.demo.provide('gatter', MV.demo.handle({
        state: state,
        sliders: sliders,
        refresh: refreshAll,
        stop: function () { pause(); if (panelB) refreshB(); },
        segments: { proto: segProto },
        setters: {
          theta: function (v) { setStart(v, state.phi); },
          phi: function (v) { setStart(state.theta, v); },
          correct: setCorrect,
          target1: setTarget,
          outcomes: setOutcomes
        },
        extra: {
          setStart: setStart,
          applyGate: function (key) { state.seq.push(key); refreshA(); },
          clearSeq: function () { state.seq = []; refreshA(); },
          vertexInfo: function () {
            var nv = nearestVertex(path().r);
            return { onVertex: nv.d < 1e-9, label: nv.v.label, ket: nv.v.ket, d: nv.d };
          },
          protoName: function () { return PROTOCOLS[state.proto].name; },
          startPlayer: function () { if (!play.running) onStartPause(); },
          isPlaying: function () { return play.running; },
          atEnd: function () {
            var proto = PROTOCOLS[state.proto];
            return state.step >= lastStep(currentRun(), proto) && play.u >= 1;
          },
          /* Ohne Durchlauf: Schritt für Schritt bis ans Ende — derselbe Weg,
             den der Knopf „Schritt vor" nimmt. */
          stepToEnd: function () {
            var proto = PROTOCOLS[state.proto], guard = 0;
            while (guard++ < 20 &&
                   !(state.step >= lastStep(currentRun(), proto) && play.u >= 1)) {
              stepForward();
            }
          },
          runNow: function () { return currentRun(); },
          tolerance: TOL
        }
      }));

      state.step = 0;
      buildOutcomeRow();
      buildStepSeg();
      renderGateTable();
    },

    draw: function () {
      if (!panelA) return;
      refreshAll();
    },

    /* Vom Router beim Verlassen des Reiters gerufen. */
    stop: function () {
      pause();
      if (panelB) refreshB();
    }
  };

}(MV));
