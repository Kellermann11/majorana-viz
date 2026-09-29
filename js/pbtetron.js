/* ------------------------------------------------------------------------
   pbtetron.js — Reiter #pb-tetron: das InAs–Pb-Tetron aus Kapitel 6.

   Panel A  Vom Feld zur Messschleife, drei Zoomstufen.
   Panel B  Zwei Objekte: die 3 µm-Teststruktur und der 3,5 µm-Draht in BA.
   Panel C  Die Lückenkaskade, Al gegen Pb.

   Dieser Reiter urteilt nicht. Er zeigt das Bauelement so, wie §6.2 der
   Arbeit es beschreibt, und macht die beiden Tatsachen sichtbar, auf denen
   die Befunde in §6.5 beruhen: gemessen wurde an einem **geerdeten**
   Tetron, und der **charakterisierte** Draht ist nicht der **vermessene**.
   Bewertet wird das alles in §6.5 der Arbeit, nicht hier.

   Alle Zahlen stammen aus §6.1, §6.2 und Tab. 6.1 der Arbeit und stehen
   genau so, wie sie dort stehen. Die Pb-Werte gehen auf [9] zurück, die
   Al-Werte nach der Unterschrift von Tab. 6.1 auf [38]. [9] liegt nur als
   Preprint vor; das steht in der Kopfzeile.

   Die Zeichnungen sind eigene Schemata nach der Textbeschreibung, keine
   nachgezeichneten Abbildungen aus [9] und keine nachgebauten Messkurven.
   Sie sind als „schematisch" gekennzeichnet.

   Einheiten: Energien in µeV, Längen in µm bzw. nm, Zeiten in ms und s,
   Magnetfeld in mT, Temperatur in mK.
   ------------------------------------------------------------------------ */

var MV = window.MV || (window.MV = {});

(function (MV) {
  'use strict';

  var PI = Math.PI;
  var TOL = 1e-9;

  var SRC = '[9] = Microsoft Quantum, arXiv:2606.03884 (Preprint)';
  var SRC_AL = '[38], wie in der Unterschrift von ' + MV.tab('materials') + ' angegeben';

  /* Hinweis, der überall dort stehen muss, wo „Z" vorkommt: CLAUDE.md §9. */
  var Z_NOTE =
    'Die Bezeichnung „Pauli-Z" und der Index von τ_Z folgen der Konvention ' +
    'von [9]: gemeint ist die Parität eines Drahtes, nach ' +
    MV.eq('tetron.decomposition') + ' also die hier x̂ genannte Achse.';

  /* ====================================================================
     1 — Die Zahlen aus Kapitel 6

     Eine Stelle, an der jede Zahl aus der Arbeit steht. Faktoren und
     abgeleitete Größen werden daraus gerechnet, nicht danebengeschrieben.
     ==================================================================== */

  var K_B = 86.17;          /* µeV/K, CLAUDE.md §2                          */
  var T_MK = 50;            /* mK, Betriebstemperatur in §6.2               */

  /* Tab. 6.1. `al` und `pb` in µeV; null heißt: die Arbeit gibt dort „–" an. */
  var SCALES = [
    { key: 'gapMother', name: 'Δ_Mutter', al: 300, pb: 1300, approx: true,
      meaning: 'Obergrenze der Kaskade' },
    { key: 'gapInd', name: 'Δ_ind (Draht)', al: null, pb: 570, approx: true,
      meaning: 'Schutz gegen Anregungen' },
    { key: 'gapT', name: 'Δ_T (ob. Quintil)', al: 30, pb: 70, approx: true,
      meaning: 'topologische Lücke, K1' },
    { key: 'eM', name: 'E_M', al: null, pb: 1, bound: true,
      meaning: 'Restaufspaltung, K2' },
    { key: 'kT', name: 'k_BT bei ' + T_MK + ' mK', al: 4.3, pb: 4.3, both: true, approx: true,
      meaning: 'Störskala' }
  ];

  /* τ_Z steht in anderen Einheiten und bekommt deshalb eine eigene Achse. */
  var TAU = {
    alMin: 1, alMax: 12,        /* ms                                        */
    pb: 22000, pbErr: 1000,     /* ms, entspricht 22 ± 1 s, Gl. (6.1)        */
    meaning: 'Paritätslebensdauer, K4'
  };

  /* §6.2, Bauelement. Längen in µm, Breiten in nm. */
  var GEO = {
    wireLen: 3.5, wireWidth: 35,
    backLen: 1.0, backWidth: 20,
    xiCl: 0.1,                  /* µm, „ξ_cl ∼ 100 nm"                       */
    locLen: 1.0,                /* µm, „übersteigt 1 µm"                     */
    testLen: 3.0                /* µm, separate Teststruktur                 */
  };

  /* §6.2, Auslese. */
  var READOUT = {
    inject: 0.5,                /* kHz                                        */
    periodMeas: 1.3, periodMeasErr: 0.3,   /* mT                              */
    periodExp: 1.0, periodExpErr: 0.3,     /* mT                              */
    dwells: 324,
    integration: 200,          /* ms, Tab. 6.2 und §6.5                      */
    resolution: 1               /* µeV, „Auflösung liegt bei etwa 1 µeV"      */
  };

  function fac(a, b) { return b / a; }

  function factorText(a, b, digits) {
    return '≈ × ' + MV.fmt(fac(a, b), digits === undefined ? 1 : digits);
  }

  /* ====================================================================
     2 — Selbstprüfung

     Was der Reiter als Faktor oder abgeleitete Zahl zeigt, wird gerechnet
     und gegen die Angabe der Arbeit gehalten.
     ==================================================================== */

  function runCheck() {
    var parts = [], worst = 0;

    function part(name, ok, detail, dev) {
      parts.push({ name: name, ok: ok, detail: detail });
      if (dev !== undefined && isFinite(dev)) worst = Math.max(worst, dev);
    }

    var kT = K_B * T_MK / 1000;
    var dKT = Math.abs(kT - 4.3);
    part('k_BT bei ' + T_MK + ' mK aus k_B = ' + MV.fmt(K_B, 2) + ' µeV/K',
         dKT < 0.05, MV.fmt(kT, 4) + ' µeV gegen ' + MV.tab('materials') + ': ≈ 4,3 µeV', dKT);

    var dc = 1.76 * kT;
    var dDC = Math.abs(dc - 7.6);
    part('Grenze des Gleichstrom-Transports, 1,76 k_BT',
         dDC < 0.05, MV.fmt(dc, 3) + ' µeV gegen ' + MV.secOfWork('6.2') + ': ≈ 7,6 µeV', dDC);

    var fM = fac(300, 1300), fT = fac(30, 70);
    part('Faktor in Δ_Mutter', Math.abs(fM - 4.3) < 0.05,
         MV.fmt(fM, 4) + ' — angezeigt als ' + factorText(300, 1300), Math.abs(fM - 4.3));
    part('Faktor in Δ_T', Math.abs(fT - 2.3) < 0.05,
         MV.fmt(fT, 4) + ' — angezeigt als ' + factorText(30, 70), Math.abs(fT - 2.3));
    part('Der Gewinn in Δ_T ist kleiner als der in Δ_Mutter', fT < fM,
         MV.fmt(fT, 2) + ' gegen ' + MV.fmt(fM, 2) + ' — wie die Unterschrift von ' +
         MV.tab('materials') + ' es sagt');

    var tLo = fac(TAU.alMax, TAU.pb), tHi = fac(TAU.alMin, TAU.pb);
    part('Faktor in τ_Z als Bereich', tLo > 1e3 && tHi < 1e5,
         MV.fmtRange(tLo, tHi, 1, null, MV.fmtExp) +
         ' — „mehr als drei Größenordnungen" nach ' + MV.secOfWork('6.2'));
    part('Der Gewinn in τ_Z ist größer als beide Lückengewinne', tLo > fM && tLo > fT,
         'kleinster τ-Faktor ' + MV.fmtExp(tLo) + ' gegen ' + MV.fmt(fM, 1) +
         ' und ' + MV.fmt(fT, 1));

    /* Die Längen in Panel A, Stufe 2 sind zueinander maßstäblich; das ist
       eine Aussage über die Zeichnung und wird an ihr geprüft. */
    var ratio = GEO.wireLen / GEO.backLen;
    part('Drahtlänge zu Rückgratlänge in der Zeichnung', Math.abs(ratio - 3.5) < TOL,
         MV.fmt(ratio, 4) + ' — die Längen sind zueinander maßstäblich, die ' +
         'Breiten ausdrücklich nicht');

    var ratioXi = GEO.wireLen / GEO.xiCl;
    part('Drahtlänge zu ξ_cl', ratioXi > 10,
         'Faktor ' + MV.fmt(ratioXi, 0) + ' — „Die Drahtlänge übersteigt die ' +
         'erwartete Kohärenzlänge deutlich" nach ' + MV.secOfWork('6.2'));

    /* Die beiden Flussperioden überlappen innerhalb ihrer Unsicherheiten.
       Das ist eine Rechnung, kein Urteil — bewertet wird in §6.5 der Arbeit. */
    var lo1 = READOUT.periodMeas - READOUT.periodMeasErr;
    var hi2 = READOUT.periodExp + READOUT.periodExpErr;
    part('Gemessene und erwartete Flussperiode überlappen', lo1 <= hi2,
         'gemessen ' + MV.fmtRange(lo1, READOUT.periodMeas + READOUT.periodMeasErr, 1, 'mT') +
         ', erwartet ' + MV.fmtRange(READOUT.periodExp - READOUT.periodExpErr, hi2, 1, 'mT') +
         ' — die Arbeit stellt beide Zahlen nebeneinander, ' +
         'ohne daraus einen Schluss zu ziehen');

    var passed = parts.every(function (p) { return p.ok !== false; });
    return { parts: parts, passed: passed, worst: worst };
  }

  var check = runCheck();

  /* ====================================================================
     3 — Zustand des Reiters
     ==================================================================== */

  var segZoom = null, segMaterial = null, testCell = null;

  var state = {
    zoom: 0,           /* 0 Elementarzelle · 1 Tetron BA · 2 Messschleife   */
    material: 2        /* 0 nur Al · 1 nur Pb · 2 beide                     */
  };

  var panelA = null, panelB = null, panelC = null;
  var plotEnergy = null, plotTime = null;
  var readoutA = null, groundLine = null;

  /* ====================================================================
     4 — Panel A, Zeichenbausteine
     ==================================================================== */

  /** Erdungssymbol: kurzer Stiel, darunter drei kürzer werdende Striche. */
  function ground(p, x, y, scale) {
    var ctx = p.ctx;
    var px = p.X(x), py = p.Y(y);
    var s = scale || 1;
    ctx.save();
    ctx.strokeStyle = MV.palette.axis;
    ctx.lineWidth = 1.4;
    ctx.lineCap = 'round';
    ctx.beginPath();
    ctx.moveTo(px, py);
    ctx.lineTo(px, py + 7 * s);
    ctx.stroke();
    [[9, 0], [6, 4], [3, 8]].forEach(function (d) {
      ctx.beginPath();
      ctx.moveTo(px - d[0] * s, py + (7 + d[1]) * s);
      ctx.lineTo(px + d[0] * s, py + (7 + d[1]) * s);
      ctx.stroke();
    });
    ctx.restore();
  }

  /**
   * Ein Tetron in H-Form: zwei waagerechte Drähte, dazwischen ein
   * senkrechtes Rückgrat in der Mitte. Alle Maße in Datenkoordinaten.
   */
  function tetron(p, cx, cy, halfW, halfH, o) {
    var ctx = p.ctx;
    o = o || {};
    var col = o.color || MV.palette.trivial;
    var wireW = o.wireW || 5;
    var backW = o.backW || 3;

    ctx.save();
    ctx.strokeStyle = col;
    ctx.lineCap = 'round';
    if (o.alpha !== undefined) ctx.globalAlpha = o.alpha;

    /* Rückgrat zuerst, damit die Drähte darüber liegen */
    ctx.lineWidth = backW;
    ctx.beginPath();
    ctx.moveTo(p.X(cx), p.Y(cy - halfH));
    ctx.lineTo(p.X(cx), p.Y(cy + halfH));
    ctx.stroke();

    ctx.lineWidth = wireW;
    [cy + halfH, cy - halfH].forEach(function (y) {
      ctx.beginPath();
      ctx.moveTo(p.X(cx - halfW), p.Y(y));
      ctx.lineTo(p.X(cx + halfW), p.Y(y));
      ctx.stroke();
    });
    ctx.restore();

    if (o.modes) {
      [[cx - halfW, cy + halfH], [cx + halfW, cy + halfH],
       [cx + halfW, cy - halfH], [cx - halfW, cy - halfH]].forEach(function (q, i) {
        p.marker(q[0], q[1], {
          r: o.modeR || 4, color: MV.palette.topological, shape: 'ring',
          label: o.modeLabels ? o.modeLabels[i] : null,
          anchor: i === 0 || i === 3 ? 'left' : 'right',
          labelBox: false
        });
      });
    }
  }

  /* ---- Stufe 1: die Elementarzelle ---- */

  function drawCell(p) {
    /* Reichlich Luft zwischen den beiden Reihen: unter jedem Tetron stehen
       Erdungssymbol und Zustandswort, über dem nächsten sein Name. Im ersten
       Entwurf lagen beide aufeinander. */
    p.begin({ x: [0, 20], y: [0, 13], grid: false, frame: false, xTicks: [], yTicks: [] });

    /* Lage der vier Tetronen nach [9], Fig. 1a, wie in CLAUDE.md §9
       angegeben. Die Arbeit selbst nennt nur, welche geerdet sind. */
    var HW = 2.6, HH = 1.0;
    var cells = [
      { name: 'AA', x: 5.4, y: 9.4, grounded: false },
      { name: 'AB', x: 14.6, y: 9.4, grounded: true },
      { name: 'BA', x: 5.4, y: 3.4, grounded: true, focus: true },
      { name: 'BB', x: 14.6, y: 3.4, grounded: false }
    ];

    cells.forEach(function (c) {
      var col = c.focus ? MV.palette.topological : MV.palette.trivial;
      var top = p.Y(c.y + HH), bot = p.Y(c.y - HH);

      tetron(p, c.x, c.y, HW, HH, {
        color: col, alpha: c.focus ? 1 : 0.5, wireW: c.focus ? 7 : 5
      });

      p.label(c.name, p.X(c.x), top - 12, {
        color: col, align: 'center', baseline: 'bottom',
        font: (c.focus ? '600 16px ' : '14px ') + 'system-ui'
      });

      /* Erdungssymbol direkt unter dem unteren Draht, Zustandswort darunter —
         bei allen vier auf derselben Höhe, damit die Reihe lesbar bleibt. */
      if (c.grounded) ground(p, c.x, c.y - HH, 1);
      p.label(c.grounded ? 'geerdet' : 'schwebend', p.X(c.x), bot + 30, {
        color: c.grounded ? MV.palette.axis : MV.palette.muted,
        align: 'center', baseline: 'top', font: '12px system-ui'
      });

      if (c.focus) {
        p.label('hier wird gemessen', p.X(c.x), bot + 48, {
          color: MV.palette.topological, align: 'center', baseline: 'top',
          font: '600 12px system-ui'
        });
      }
    });

    /* Rahmen um die Elementarzelle */
    p.ctx.save();
    p.ctx.strokeStyle = MV.palette.line;
    p.ctx.setLineDash([5, 5]);
    p.ctx.lineWidth = 1;
    p.ctx.strokeRect(p.X(1.2), p.Y(12.0), p.X(18.8) - p.X(1.2), p.Y(0.6) - p.Y(12.0));
    p.ctx.restore();
    p.label('Elementarzelle, vier Tetronen', p.X(1.2), p.Y(12.0) - 5, {
      color: MV.palette.muted, align: 'left', baseline: 'bottom', font: '12px system-ui'
    });

    p.finish();
  }

  /* ---- Stufe 2: Tetron BA mit Maßen ---- */

  function drawDevice(p) {
    /* Datenkoordinaten in µm. Die Längen sind zueinander maßstäblich. */
    var L = GEO.wireLen, H = GEO.backLen;
    p.begin({
      x: [-1.5, L + 1.6], y: [-1.15, H + 1.35],
      grid: false, frame: false, xTicks: [], yTicks: []
    });

    tetron(p, L / 2, H / 2, L / 2, H / 2, {
      color: MV.palette.topological, wireW: 7, backW: 4,
      modes: true, modeR: 5, modeLabels: ['γ₁', 'γ₂', 'γ₃', 'γ₄']
    });

    /* Maße */
    var ctx = p.ctx;
    function hDim(y, x0, x1, text) {
      var py = p.Y(y);
      ctx.save();
      ctx.strokeStyle = MV.palette.axis;
      ctx.fillStyle = MV.palette.axis;
      ctx.lineWidth = 1;
      ctx.beginPath();
      ctx.moveTo(p.X(x0), py); ctx.lineTo(p.X(x1), py);
      ctx.stroke();
      [[p.X(x0), 1], [p.X(x1), -1]].forEach(function (e) {
        ctx.beginPath();
        ctx.moveTo(e[0], py);
        ctx.lineTo(e[0] + 5 * e[1], py - 3);
        ctx.lineTo(e[0] + 5 * e[1], py + 3);
        ctx.closePath();
        ctx.fill();
      });
      ctx.restore();
      p.label(text, (p.X(x0) + p.X(x1)) / 2, py - 5, {
        color: MV.palette.axis, align: 'center', baseline: 'bottom',
        font: '12px system-ui', box: true
      });
    }

    hDim(H + 0.55, 0, L, 'Draht ' + MV.fmt(L, 1) + ' µm lang, ' + GEO.wireWidth + ' nm breit');

    /* Senkrechtes Maß für das Rückgrat */
    var px = p.X(L / 2) + 26;
    ctx.save();
    ctx.strokeStyle = MV.palette.axis;
    ctx.fillStyle = MV.palette.axis;
    ctx.lineWidth = 1;
    ctx.beginPath();
    ctx.moveTo(px, p.Y(0)); ctx.lineTo(px, p.Y(H));
    ctx.stroke();
    [[p.Y(0), -1], [p.Y(H), 1]].forEach(function (e) {
      ctx.beginPath();
      ctx.moveTo(px, e[0]);
      ctx.lineTo(px - 3, e[0] + 5 * e[1]);
      ctx.lineTo(px + 3, e[0] + 5 * e[1]);
      ctx.closePath();
      ctx.fill();
    });
    ctx.restore();
    p.label('Rückgrat ' + MV.fmt(GEO.backLen, 1) + ' µm, ' + GEO.backWidth + ' nm breit',
            px + 7, (p.Y(0) + p.Y(H)) / 2, {
      color: MV.palette.axis, align: 'left', baseline: 'middle', font: '12px system-ui'
    });

    /* Die beiden Längenskalen, maßstäblich neben dem Draht. Genau darum
       geht es in §6.2: die Drahtlänge übersteigt die Kohärenzlänge
       deutlich, und nur deshalb ist E_M exponentiell unterdrückt. */
    function scaleBar(y, len, text, color) {
      var py = p.Y(y);
      ctx.save();
      ctx.strokeStyle = color;
      ctx.lineWidth = 4;
      ctx.lineCap = 'butt';
      ctx.beginPath();
      ctx.moveTo(p.X(0), py);
      ctx.lineTo(p.X(len), py);
      ctx.stroke();
      ctx.restore();
      p.label(text, p.X(len) + 8, py, {
        color: color, align: 'left', baseline: 'middle', font: '12px system-ui'
      });
    }
    scaleBar(-0.42, GEO.xiCl, 'ξ_cl ∼ ' + Math.round(GEO.xiCl * 1000) + ' nm — ' +
             'saubere Abschätzung der Kohärenzlänge', MV.palette.curve2);
    scaleBar(-0.82, GEO.locLen, 'gemessene Lokalisierungslänge > ' +
             MV.fmt(GEO.locLen, 0) + ' µm', MV.palette.curve1);

    p.label('beide maßstäblich zur Drahtlänge', p.X(0), p.Y(-1.06), {
      color: MV.palette.muted, align: 'left', baseline: 'middle', font: '11px system-ui'
    });

    p.finish();
  }

  /* ---- Stufe 3: die Messschleife am oberen Draht ----

     Die meisten Bezeichnungen stehen in der Legende unter der Abbildung und
     nicht im Bild. Im ersten Entwurf lagen sie im Bild und überschrieben
     einander: eine flache Zeichnung hat keine Ecke, die bei jeder Breite
     frei bleibt.                                                            */

  function drawLoop(p) {
    var L = GEO.wireLen;
    p.begin({
      x: [-1.9, L + 1.9], y: [-2.0, 2.5],
      grid: false, frame: false, xTicks: [], yTicks: []
    });
    var ctx = p.ctx;
    var yLow = -1.5;

    /* Der obere Draht, hervorgehoben; der untere nur angedeutet. */
    ctx.save();
    ctx.strokeStyle = MV.palette.topological;
    ctx.lineWidth = 8;
    ctx.lineCap = 'round';
    ctx.beginPath();
    ctx.moveTo(p.X(0), p.Y(0)); ctx.lineTo(p.X(L), p.Y(0));
    ctx.stroke();
    ctx.globalAlpha = 0.25;
    ctx.lineWidth = 6;
    ctx.beginPath();
    ctx.moveTo(p.X(0), p.Y(yLow)); ctx.lineTo(p.X(L), p.Y(yLow));
    ctx.stroke();
    ctx.lineWidth = 4;
    ctx.beginPath();
    ctx.moveTo(p.X(L / 2), p.Y(yLow)); ctx.lineTo(p.X(L / 2), p.Y(0));
    ctx.stroke();
    ctx.restore();

    p.label('oberer Draht', p.X(0) - 10, p.Y(0), {
      color: MV.palette.topological, align: 'right', baseline: 'middle',
      font: '600 12px system-ui'
    });
    p.label('unterer Draht', p.X(0) - 10, p.Y(yLow), {
      color: MV.palette.muted, align: 'right', baseline: 'middle', font: '12px system-ui'
    });

    /* Die beiden Nullmoden an den Enden des oberen Drahtes */
    p.marker(0, 0, { r: 5, color: MV.palette.topological, shape: 'ring',
      legend: 'Nullmoden an den Enden des oberen Drahtes' });
    p.marker(L, 0, { r: 5, color: MV.palette.topological, shape: 'ring' });

    /* Der lange Quantenpunkt QDL, parallel über dem Draht */
    var qy = 1.75, qx0 = 0.7, qx1 = L - 0.7;
    ctx.save();
    ctx.strokeStyle = MV.palette.curve1;
    ctx.fillStyle = 'rgba(43, 93, 158, 0.10)';
    ctx.lineWidth = 1.8;
    var rx = p.X(qx0), rw = p.X(qx1) - p.X(qx0);
    var ry = p.Y(qy) - 9, rh = 18;
    ctx.beginPath();
    ctx.rect(rx, ry, rw, rh);
    ctx.fill();
    ctx.stroke();
    ctx.restore();
    p.label('QDL', rx + rw / 2, p.Y(qy), {
      color: MV.palette.curve1, align: 'center', baseline: 'middle',
      font: '600 12px system-ui'
    });

    /* QD1 und QD2 — ein Punkt je Drahtende, und zwar **dieselben** Punkte
       für beide Aufgaben: bei der Z-Auslese koppeln sie QD_L an die
       Nullmoden, bei der Charakterisierung wird einer dispersiv ausgelesen
       und der andere injiziert, dann werden die Rollen getauscht ([9]).

       Links QD1, rechts QD2: [9] schreibt „Fig. 5(a) presents the
       aggregation of measurements on QD1 and Fig. 5(b) on QD2", und die
       Unterschrift zu Fig. 5 sagt „(a,b) … the signal measured on the left
       (right) QD". Links und rechts folgen damit der Orientierung der
       Abbildungen in [9]; die Zeichnung hier ist schematisch. */
    [[0, qx0, 'QD1'], [L, qx1, 'QD2']].forEach(function (pair, i) {
      var mx = (pair[0] + pair[1]) / 2;
      ctx.save();
      ctx.strokeStyle = MV.palette.curve1;
      ctx.lineWidth = 1.4;
      ctx.beginPath();
      ctx.moveTo(p.X(pair[0]), p.Y(0));
      ctx.lineTo(p.X(mx), p.Y(0.85));
      ctx.lineTo(p.X(pair[1]), p.Y(qy));
      ctx.stroke();
      ctx.restore();
      p.marker(mx, 0.85, {
        r: 5, color: MV.palette.curve1, shape: 'square',
        label: pair[2], anchor: i === 0 ? 'left' : 'right', labelBox: false,
        legend: i === 0
          ? 'QD1, QD2 — ein Punkt je Drahtende. Bei der Z-Auslese koppeln ' +
            'beide QD_L an die Nullmoden; bei der Charakterisierung wird ' +
            'einer dispersiv ausgelesen, der andere injiziert, dann werden ' +
            'die Rollen getauscht — daraus E_M, Auflösung ≈ ' +
            MV.fmt(READOUT.resolution, 0) + ' µeV'
          : null
      });
    });

    /* Der eingeschlossene Fluss, mitten in der Schleife */
    p.label('Φ', p.X(L / 2), p.Y(1.05), {
      color: MV.palette.critical, align: 'center', baseline: 'middle',
      font: '600 17px system-ui'
    });
    p.legendItems.push({
      label: 'Φ — eingeschlossener Fluss; die Sichtbarkeit der Bimodalität ' +
             'oszilliert mit ihm, Periode ' + MV.fmt(READOUT.periodMeas, 1) + ' ± ' +
             MV.fmt(READOUT.periodMeasErr, 1) + ' mT',
      color: MV.palette.critical, swatch: 'dot'
    });

    /* QD4, außerhalb der Schleife */
    ctx.save();
    ctx.strokeStyle = MV.palette.curve2;
    ctx.lineWidth = 1.4;
    ctx.setLineDash([4, 3]);
    ctx.beginPath();
    ctx.moveTo(p.X(L + 1.2), p.Y(0));
    ctx.lineTo(p.X(L + 0.25), p.Y(0));
    ctx.stroke();
    ctx.restore();
    p.arrow(L + 0.55, 0, -14, 0, { color: MV.palette.curve2, head: 5 });
    p.marker(L + 1.2, 0, { r: 6, color: MV.palette.curve2, shape: 'square',
      legend: 'QD4 — Injektor, ' + MV.fmt(READOUT.inject, 1) +
              ' kHz, nicht Teil der Schleife' });
    p.label('QD4', p.X(L + 1.2), p.Y(0) - 12, {
      color: MV.palette.curve2, align: 'center', baseline: 'bottom', font: '12px system-ui'
    });

    p.finish();
  }

  var ZOOM = [
    { name: 'Elementarzelle', draw: drawCell, aspect: 0.56,
      margin: { l: 24, r: 24, t: 26, b: 26 },
      note: 'Vier Tetronen. Zwei sind geerdet (AB und BA), zwei schwebend ' +
            '(AA und BB). Die schwebenden werden gebraucht, weil bei einer ' +
            'Zweiqubit-Messung mindestens eine Insel eine Ladungsenergie ' +
            'tragen muss.' },
    { name: 'Tetron BA', draw: drawDevice, aspect: 0.44,
      margin: { l: 44, r: 190, t: 36, b: 40 },
      note: 'Die H-Form aus ' + MV.sec('5.2') + '. Die Längen sind zueinander ' +
            'maßstäblich, die Breiten ausdrücklich <b>nicht</b> — 35 nm wären ' +
            'bei dieser Länge kaum ein Haar breit. Weil die Drahtlänge die ' +
            'Kohärenzlänge weit übersteigt, ist E_M nach ' +
            MV.eq('kitaev.splitting') + ' exponentiell unterdrückt.' },
    { name: 'Messschleife', draw: drawLoop, aspect: 0.42,
      margin: { l: 106, r: 60, t: 28, b: 30 },
      note: 'Die Paritätsinformation erscheint als Bimodalität der ' +
            'Quantenkapazität, deren Sichtbarkeit vom eingeschlossenen Fluss Φ ' +
            'abhängt. <b>Stück für Stück derselbe Aufbau wie im Reiter ' +
            'Auslese:</b> QD_L entspricht dem Referenzweg im Grenzfall ' +
            't₀ ≫ t₁, in dem die beiden Punkte zu einem hybridisieren; QD1 und ' +
            'QD2 koppeln an die beiden Nullmoden des oberen Drahtes — dieselbe ' +
            'Rolle wie die Quantenpunkte dort; und die Auslese über die ' +
            'Quantenkapazität ist das dispersive Verfahren 3, so auch ' +
            MV.secOfWork('5.5') + '. <b>Dieselben beiden Punkte tragen auch die ' +
            'Charakterisierung:</b> dort wird einer dispersiv ausgelesen und der ' +
            'andere injiziert, dann werden die Rollen getauscht — daraus E_M.' }
  ];

  function drawA() {
    var p = panelA.plot;
    if (!p) return;
    var z = ZOOM[state.zoom];
    p.aspect = z.aspect;
    p.margin = z.margin;
    z.draw(p);
    panelA.setStatus(
      'Stufe ' + (state.zoom + 1) + ' von 3 — ' + z.name +
      '. Schematisch, nach der Textbeschreibung in ' + MV.secOfWork('6.2') +
      '; keine nachgezeichnete Abbildung aus [9].',
      MV.palette.muted);
  }

  /* ====================================================================
     5 — Panel B: zwei Objekte

     Der charakterisierte Draht ist nicht der vermessene. Die Tabelle stellt
     beide nebeneinander und **urteilt nicht** — kein Status, keine Wertung.
     ==================================================================== */

  var OBJECTS = {
    test: {
      title: '3 µm-Teststruktur',
      sub: 'separate Nanodraht-Teststruktur',
      len: GEO.testLen,
      rows: [
        'Phasendiagramm per Topological-Gap-Protokoll',
        'Δ_T ≈ ' + SCALES[2].pb + ' µeV (oberstes Quintil)',
        'Ausdehnung der topologischen Phase über 1,1 mV·T',
        'Nullspannungsmaxima an beiden Drahtenden',
        'Schließen und Wiederöffnen der nichtlokalen Leitfähigkeit'
      ]
    },
    device: {
      title: '3,5 µm-Draht in Tetron BA',
      sub: 'der vermessene Draht — geerdetes Bauelement',
      len: GEO.wireLen,
      rows: [
        'E_M < ' + SCALES[3].pb + ' µeV, auflösungsbegrenzt (Auflösung ≈ ' +
          MV.fmt(READOUT.resolution, 0) + ' µeV)',
        'Flussperiode ' + MV.fmt(READOUT.periodMeas, 1) + ' ± ' +
          MV.fmt(READOUT.periodMeasErr, 1) + ' mT; aus der Schleifenfläche erwartet ' +
          MV.fmt(READOUT.periodExp, 1) + ' ± ' + MV.fmt(READOUT.periodExpErr, 1) +
          ' mT für h/2e',
        'τ_Z = ' + MV.fmt(TAU.pb / 1000, 0) + ' ± ' + MV.fmt(TAU.pbErr / 1000, 0) +
          ' s nach ' + MV.eq('pb.tauZ') + ', aus ' + READOUT.dwells +
          ' Verweildauern (Gaußsches Mischmodell mit doppelter Schwelle, ' +
          'Einzelexponentialanpassung), bei ' + READOUT.integration +
          ' ms Integrationszeit'
      ]
    }
  };

  function drawWire(p, len, maxLen, color, label) {
    p.begin({
      x: [-0.25, maxLen + 0.25], y: [-1, 1],
      grid: false, frame: false, xTicks: [], yTicks: []
    });
    var ctx = p.ctx;
    ctx.save();
    ctx.strokeStyle = color;
    ctx.lineWidth = 9;
    ctx.lineCap = 'round';
    ctx.beginPath();
    ctx.moveTo(p.X(0), p.Y(0));
    ctx.lineTo(p.X(len), p.Y(0));
    ctx.stroke();
    ctx.restore();
    p.marker(0, 0, { r: 4.5, color: color, shape: 'ring' });
    p.marker(len, 0, { r: 4.5, color: color, shape: 'ring' });
    p.label(label, p.X(len / 2), p.Y(0) - 14, {
      color: color, align: 'center', baseline: 'bottom', font: '12px system-ui'
    });
    p.label('maßstäblich zueinander', p.X(maxLen + 0.25), p.Y(0) + 15, {
      color: MV.palette.muted, align: 'right', baseline: 'top', font: '10px system-ui'
    });
    p.finish();
  }

  function objectHtml(o) {
    return '<ul class="objlist">' + o.rows.map(function (r) {
      return '<li>' + r + '</li>';
    }).join('') + '</ul>';
  }

  /* ====================================================================
     6 — Panel C: die Lückenkaskade
     ==================================================================== */

  function showAl() { return state.material === 0 || state.material === 2; }
  function showPb() { return state.material === 1 || state.material === 2; }

  var COL_AL = MV.palette.curve1, COL_PB = MV.palette.topological;

  function logTicks(lo, hi, unit) {
    var out = [], e;
    for (e = Math.ceil(Math.log10(lo)); e <= Math.floor(Math.log10(hi)); e++) {
      out.push({
        v: e,
        label: (e < 0 ? MV.fmt(Math.pow(10, e), -e) : String(Math.pow(10, e))) +
               (unit ? ' ' + unit : '')
      });
    }
    return out;
  }

  function drawEnergy() {
    var p = plotEnergy;
    if (!p) return;
    var rows = SCALES;
    var n = rows.length;

    p.begin({
      x: [Math.log10(0.4), Math.log10(3000)],
      y: [-0.6, n - 0.4],
      grid: true,
      xTicks: logTicks(1, 1000, 'µeV'),
      yTicks: [],
      frame: true
    });
    var ctx = p.ctx;

    /* k_BT als senkrechte Störskala durch das ganze Bild */
    var kT = SCALES[4].pb;
    p.vline(Math.log10(kT), {
      color: MV.palette.critical, dash: [4, 4], width: 1.2,
      label: 'k_BT', legend: 'k_BT bei ' + T_MK + ' mK ≈ ' + MV.fmt(kT, 1) + ' µeV'
    });

    rows.forEach(function (r, i) {
      var y = n - 1 - i;

      p.label(r.name, p.px0 - 8, p.Y(y), {
        color: MV.palette.ink, align: 'right', baseline: 'middle', font: '12px system-ui'
      });

      if (r.key === 'kT') return;   /* schon als senkrechte Linie gezeichnet */

      var hasAl = r.al !== null && showAl();
      var hasPb = r.pb !== null && showPb();

      if (hasAl && hasPb) {
        ctx.save();
        ctx.strokeStyle = MV.palette.muted;
        ctx.lineWidth = 1.2;
        ctx.beginPath();
        ctx.moveTo(p.X(Math.log10(r.al)), p.Y(y));
        ctx.lineTo(p.X(Math.log10(r.pb)), p.Y(y));
        ctx.stroke();
        ctx.restore();
        p.label(factorText(r.al, r.pb),
                (p.X(Math.log10(r.al)) + p.X(Math.log10(r.pb))) / 2, p.Y(y) - 7, {
          color: MV.palette.muted, align: 'center', baseline: 'bottom',
          font: '11px system-ui', box: true
        });
      }

      if (hasAl) {
        p.marker(Math.log10(r.al), y, {
          r: 5, color: COL_AL, shape: 'square',
          legend: i === 0 ? 'InAs–Al' : null
        });
      }
      if (hasPb) {
        if (r.bound) {
          /* E_M ist eine obere Schranke, kein Wert. */
          p.marker(Math.log10(r.pb), y, { r: 5, color: COL_PB, shape: 'ring',
            legend: 'obere Schranke' });
          p.arrow(Math.log10(r.pb), y, -26, 0, { color: COL_PB, head: 5, width: 1.4 });
          p.label('< ' + MV.fmt(r.pb, 0) + ' µeV', p.X(Math.log10(r.pb)) + 9, p.Y(y), {
            color: COL_PB, align: 'left', baseline: 'middle', font: '11px system-ui'
          });
        } else {
          p.marker(Math.log10(r.pb), y, {
            r: 5, color: COL_PB, legend: i === 0 ? 'InAs–Pb' : null
          });
        }
      }
    });

    p.finish();
  }

  function drawTime() {
    var p = plotTime;
    if (!p) return;

    p.begin({
      x: [Math.log10(0.4), Math.log10(120000)],
      y: [-0.6, 0.6],
      grid: true,
      xTicks: [
        { v: 0, label: '1 ms' }, { v: 1, label: '10 ms' }, { v: 2, label: '100 ms' },
        { v: 3, label: '1 s' }, { v: 4, label: '10 s' }
      ],
      yTicks: [],
      frame: true
    });
    var ctx = p.ctx;

    p.label('τ_Z', p.px0 - 8, p.Y(0), {
      color: MV.palette.ink, align: 'right', baseline: 'middle', font: '12px system-ui'
    });

    if (showAl()) {
      ctx.save();
      ctx.strokeStyle = COL_AL;
      ctx.lineWidth = 7;
      ctx.lineCap = 'butt';
      ctx.beginPath();
      ctx.moveTo(p.X(Math.log10(TAU.alMin)), p.Y(0));
      ctx.lineTo(p.X(Math.log10(TAU.alMax)), p.Y(0));
      ctx.stroke();
      ctx.restore();
      p.label(TAU.alMin + '–' + TAU.alMax + ' ms',
              (p.X(Math.log10(TAU.alMin)) + p.X(Math.log10(TAU.alMax))) / 2, p.Y(0) + 11, {
        color: COL_AL, align: 'center', baseline: 'top', font: '11px system-ui'
      });
    }

    if (showPb()) {
      p.marker(Math.log10(TAU.pb), 0, { r: 5, color: COL_PB });
      ctx.save();
      ctx.strokeStyle = COL_PB;
      ctx.lineWidth = 2;
      ctx.beginPath();
      ctx.moveTo(p.X(Math.log10(TAU.pb - TAU.pbErr)), p.Y(0));
      ctx.lineTo(p.X(Math.log10(TAU.pb + TAU.pbErr)), p.Y(0));
      ctx.stroke();
      ctx.restore();
      p.label(MV.fmt(TAU.pb / 1000, 0) + ' ± ' + MV.fmt(TAU.pbErr / 1000, 0) + ' s',
              p.X(Math.log10(TAU.pb)), p.Y(0) + 11, {
        color: COL_PB, align: 'center', baseline: 'top', font: '11px system-ui'
      });
    }

    if (showAl() && showPb()) {
      ctx.save();
      ctx.strokeStyle = MV.palette.muted;
      ctx.lineWidth = 1.2;
      ctx.beginPath();
      ctx.moveTo(p.X(Math.log10(TAU.alMax)), p.Y(0.22));
      ctx.lineTo(p.X(Math.log10(TAU.pb)), p.Y(0.22));
      ctx.stroke();
      ctx.restore();
      p.label('≈ × ' + MV.fmtRange(fac(TAU.alMax, TAU.pb), fac(TAU.alMin, TAU.pb),
                              1, null, MV.fmtExp),
              (p.X(Math.log10(TAU.alMax)) + p.X(Math.log10(TAU.pb))) / 2, p.Y(0.22) - 5, {
        color: MV.palette.muted, align: 'center', baseline: 'bottom',
        font: '11px system-ui', box: true
      });
    }

    p.finish();
  }

  function materialTable() {
    var cols = [];
    if (showAl()) cols.push('al');
    if (showPb()) cols.push('pb');

    function cell(r, k) {
      if (r.both) return (k === 'al' ? '≈ ' + MV.fmt(r.al, 1) : '—');
      var v = r[k];
      if (v === null) return '–';
      return (r.bound ? '< ' : r.approx ? '∼ ' : '') + MV.fmt(v, 0);
    }

    var html = '<table class="numtable"><caption>' + MV.tab('materials') +
      ' der Arbeit — Energieskalen in µeV, Lebensdauern in eigenen Einheiten. ' +
      'Die Pb-Werte stammen aus [9], die Al-Werte aus ' + SRC_AL +
      '.</caption><thead><tr><th scope="col">Größe</th>' +
      cols.map(function (k) {
        return '<th scope="col" class="num">InAs–' + (k === 'al' ? 'Al' : 'Pb') + '</th>';
      }).join('') +
      '<th scope="col">Bedeutung</th></tr></thead><tbody>';

    SCALES.forEach(function (r) {
      html += '<tr><th scope="row">' + r.name + '</th>' +
        cols.map(function (k) {
          if (r.both) {
            return k === cols[0]
              ? '<td class="num" colspan="' + cols.length + '">≈ ' + MV.fmt(r.al, 1) + '</td>'
              : '';
          }
          return '<td class="num">' + cell(r, k) + '</td>';
        }).join('') +
        '<td>' + r.meaning + '</td></tr>';
    });

    html += '<tr class="numtable__sub"><th scope="row">τ_Z</th>' +
      cols.map(function (k) {
        return '<td class="num">' + (k === 'al'
          ? TAU.alMin + '–' + TAU.alMax + ' ms'
          : MV.fmt(TAU.pb / 1000, 0) + ' ± ' + MV.fmt(TAU.pbErr / 1000, 0) + ' s') + '</td>';
      }).join('') +
      '<td>' + TAU.meaning + '</td></tr>';

    html += '</tbody></table>';
    return html;
  }

  /* ====================================================================
     7 — Auffrischen
     ==================================================================== */

  function refreshA() {
    drawA();
    readoutA.set('zoom', (state.zoom + 1) + ' von 3 — ' + ZOOM[state.zoom].name);
    readoutA.set('device', 'Tetron BA');
    readoutA.set('ground', 'geerdet');
    readoutA.set('wire', 'oberer Nanodraht');
    readoutA.status('trivial', 'Alle Messungen: oberer Draht von Tetron BA — geerdet');
    if (groundLine) groundLine.innerHTML = ZOOM[state.zoom].note;
    if (panelA.info) panelA.info.update();
  }

  function refreshC() {
    drawEnergy();
    drawTime();
    var host = panelC.body.querySelector('.materials-host');
    if (host) host.innerHTML = materialTable();
    panelC.setStatus(
      state.material === 0 ? 'Nur InAs–Al.'
        : state.material === 1 ? 'Nur InAs–Pb.'
        : 'Der Gewinn in Δ_T (' + factorText(30, 70) + ') ist deutlich kleiner als der ' +
          'in Δ_Mutter (' + factorText(300, 1300) + '); der Gewinn in τ_Z ist ungleich ' +
          'größer als beide.',
      MV.palette.muted);
    if (panelC.info) panelC.info.update();
  }

  function refreshAll() {
    refreshA();
    drawWire(MV.pbTetron._wireTest, GEO.testLen, GEO.wireLen, MV.palette.curve1,
             MV.fmt(GEO.testLen, 0) + ' µm');
    drawWire(MV.pbTetron._wireDev, GEO.wireLen, GEO.wireLen, MV.palette.topological,
             MV.fmt(GEO.wireLen, 1) + ' µm');
    refreshC();
    if (panelB.info) panelB.info.update();
  }

  /* ====================================================================
     8 — Info-Panels, acht Punkte nach CLAUDE.md §11
     ==================================================================== */

  function numericsHtml() {
    return '<ul>' + check.parts.map(function (p) {
      var mark = p.ok === true ? '✓' : p.ok === false ? '✗' : '!';
      return '<li><b>' + mark + '</b> ' + p.name + ' — ' + p.detail + '</li>';
    }).join('') + '</ul>' +
      '<p>Gerechnet wird aus den Zahlen der Arbeit; die Faktoren stehen ' +
      'nirgends doppelt. Die Zeichnungen enthalten keine Messkurve: es gibt ' +
      'in diesem Reiter nichts zu simulieren, nur etwas darzustellen.</p>';
  }

  var SCHEMATIC_NOTE =
    'Alle Zeichnungen dieses Reiters sind <b>eigene Schemata</b> nach der ' +
    'Textbeschreibung in ' + MV.secOfWork('6.2') + ' — keine nachgezeichneten ' +
    'Abbildungen aus [9] und keine nachgebauten Messkurven.';

  var INFO_A = {
    what: 'Das vermessene Bauelement in drei Zoomstufen: die Elementarzelle aus ' +
          'vier Tetronen, das Tetron BA in H-Form mit seinen Maßen und die ' +
          'Messschleife am oberen Draht. In Stufe 3 sitzt an jedem Drahtende ' +
          '<b>ein</b> Quantenpunkt — QD1 links, QD2 rechts. Es sind dieselben ' +
          'beiden Punkte für die Auslese und für die Charakterisierung; ein ' +
          'zweites Punktpaar gibt es nicht.',
    model: 'Kein Modell, sondern eine Beschreibung. Die Geometrie stammt aus ' +
           MV.secOfWork('6.2') + ' und aus [9], Fig. 1a für die Anordnung der ' +
           'vier Tetronen in der Zelle. Physik steckt hier nur in einer Aussage: ' +
           'weil die Drahtlänge die Kohärenzlänge weit übersteigt, ist E_M nach ' +
           MV.eqOfWork('kitaev.splitting') + ' exponentiell unterdrückt.',
    computed: 'Nichts wird simuliert. Gerechnet werden nur die Verhältnisse, die ' +
              'die Zeichnung maßstäblich halten: Drahtlänge zu Rückgratlänge und ' +
              'Drahtlänge zu ξ_cl.',
    params: function () {
      return '<ul><li>Zoomstufe: ' + ZOOM[state.zoom].name + '</li>' +
        '<li>Draht ' + MV.fmt(GEO.wireLen, 1) + ' µm × ' + GEO.wireWidth +
        ' nm, Rückgrat ' + MV.fmt(GEO.backLen, 1) + ' µm × ' + GEO.backWidth + ' nm</li>' +
        '<li>ξ_cl ∼ ' + Math.round(GEO.xiCl * 1000) + ' nm, gemessene ' +
        'Lokalisierungslänge > ' + MV.fmt(GEO.locLen, 0) + ' µm</li>' +
        '<li>Verhältnis Drahtlänge zu ξ_cl: ' + MV.fmt(GEO.wireLen / GEO.xiCl, 0) + '</li>' +
        '<li>Injektionsrate QD4: ' + MV.fmt(READOUT.inject, 1) + ' kHz</li></ul>';
    },
    numerics: numericsHtml(),
    convention: '<b>Abweichung von der Arbeit.</b> Die Arbeit beschreibt ' +
                'Charakterisierung und Auslese in getrennten Absätzen; nach [9] ' +
                'sind es dieselben Punkte QD1 und QD2. Die Seite folgt hier ' +
                '[9] — auf ausdrücklichen Auftrag und nur an dieser einen ' +
                'Stelle; sonst gilt weiterhin die Arbeit. ' +
                SCHEMATIC_NOTE + ' Die Anordnung der vier Tetronen in der Zelle ' +
                'folgt [9], Fig. 1a; die Arbeit selbst nennt nur, welche geerdet ' +
                'und welche schwebend sind. Längen sind in Stufe 2 zueinander ' +
                'maßstäblich, Breiten nicht. ' + Z_NOTE,
    reference: MV.secOfWork('6.2') + '. Quelle aller Angaben: ' + SRC + '. Das ' +
               'Ausleseprinzip ist eine Variante des interferometrischen ' +
               'Verfahrens aus dem Reiter Auslese — der Fall t₀ ≫ t₁, in dem die ' +
               'beiden Quantenpunkte zu einem hybridisieren, dessen Energie von ' +
               'der Parität abhängt und über die Quantenkapazität gelesen wird. ' +
               'Dass QD1 und QD2 beide Aufgaben tragen, steht in [9]: „A Pauli-Z ' +
               'measurement is performed by coupling the top wire in the qubit ' +
               'to quantum dot QDL via QD1 and QD2" und, zur Charakterisierung, ' +
               '„… for the case where QD2 is the readout QD and QD1 the injector ' +
               'QD. … Repeating the measurement with the role of QD1 and QD2 ' +
               'reversed …".',
    reading: '<p><b>Zoome von der Zelle bis zur Schleife</b> und achte darauf, an ' +
             'welchem Tetron gemessen wird und ob es geerdet ist: In Stufe 1 ' +
             'tragen AB und BA ein Erdungssymbol, AA und BB nicht. In Stufe 2 und ' +
             '3 bleibt die Zeile unter der Abbildung auf <b>oberer Draht von ' +
             'Tetron BA — geerdet</b> stehen. <b>Sieh dir in Stufe 3 die Punkte ' +
             'an den Drahtenden an:</b> es sind genau zwei, QD1 und QD2, und sie ' +
             'tragen die Auslese und die Charakterisierung gemeinsam. Auf welchem ' +
             'Draht die topologische Lücke gemessen wurde, zeigt Panel B.</p>'
  };

  var INFO_B = {
    what: 'Die beiden Objekte, aus denen die Befunde stammen: links eine separate ' +
          '3 µm-Nanodraht-Teststruktur, rechts der 3,5 µm-Draht im Tetron BA. ' +
          'Die Drahtlängen sind zueinander maßstäblich gezeichnet.',
    model: 'Kein Modell. Die Zeilen sind die Angaben aus ' + MV.secOfWork('6.2') +
           ', jede auf der Seite, auf der sie gewonnen wurde.',
    computed: 'Nichts. Die einzige Rechnung ist der Vergleich der beiden ' +
              'Flussperioden, und der steht unter „Numerik“.',
    params: function () {
      return '<ul><li>Teststruktur: ' + MV.fmt(GEO.testLen, 0) + ' µm</li>' +
        '<li>vermessener Draht: ' + MV.fmt(GEO.wireLen, 1) + ' µm, in Tetron BA, geerdet</li>' +
        '<li>Flussperiode gemessen ' + MV.fmt(READOUT.periodMeas, 1) + ' ± ' +
        MV.fmt(READOUT.periodMeasErr, 1) + ' mT, erwartet ' +
        MV.fmt(READOUT.periodExp, 1) + ' ± ' + MV.fmt(READOUT.periodExpErr, 1) + ' mT</li>' +
        '<li>τ_Z = ' + MV.fmt(TAU.pb / 1000, 0) + ' ± ' + MV.fmt(TAU.pbErr / 1000, 0) +
        ' s aus ' + READOUT.dwells + ' Verweildauern, ' + READOUT.integration +
        ' ms Integrationszeit</li></ul>';
    },
    numerics: numericsHtml(),
    convention: 'Diese Seite <b>bewertet nicht</b>: es gibt keine Statusspalte und ' +
                'kein Urteil. Sie hält nur auseinander, welche Größe von welchem ' +
                'Draht stammt. Was daraus folgt, wird in ' + MV.secOfWork('6.5') +
                ' bewertet, nicht hier. ' + Z_NOTE +
                ' Die Integrationszeit von ' + READOUT.integration + ' ms steht nicht ' +
                'in ' + MV.sec('6.2') + ', sondern in ' + MV.tab('criteriaCheck') +
                ' und in ' + MV.secOfWork('6.5') + '; dort ist sie der Nenner des ' +
                'Verhältnisses von Lebensdauer zu Messzeit.',
    reference: MV.secOfWork('6.2') + ', ' + MV.eq('pb.tauZ') + '. Quelle: ' + SRC + '. ' +
               'Das Kriterium für die topologische Phase ist ' +
               MV.eq('wire.criterion') + ' aus dem Reiter Rashba-Nanodraht.',
    reading: '<p><b>Vergleiche die beiden Spalten: Welche Größe stammt von welchem ' +
             'Draht?</b> Die topologische Lücke Δ_T ≈ 70 µeV, die Ausdehnung über ' +
             '1,1 mV·T, die Nullspannungsmaxima und das Schließen und ' +
             'Wiederöffnen der nichtlokalen Leitfähigkeit stehen <b>alle links</b>, ' +
             'an der Teststruktur. Rechts, am vermessenen Draht, stehen nur E_M, ' +
             'die Flussperiode und τ_Z. Der charakterisierte Draht ist nicht der ' +
             'vermessene.</p>'
  };

  var INFO_C = {
    what: 'Die Energieskalen und die Paritätslebensdauer beider ' +
          'Materialplattformen auf logarithmischen Achsen. Links Energien in ' +
          'µeV, rechts Zeiten — zwei Achsen, weil es zwei Einheiten sind.',
    model: 'Die Lückenkaskade aus ' + MV.secOfWork('6.1') + ': die Lücke des ' +
           'Muttersupraleiters setzt die Obergrenze für die induzierte Lücke im ' +
           'Halbleiter, diese wiederum für die topologische Lücke Δ_T. Die ' +
           'Kaskade ist nicht proportional — zwischen Δ und Δ_T stehen die ' +
           'Spin-Bahn-Kopplung α, der g-Faktor und die Unordnung, und ein ' +
           'dickerer Supraleiter renormiert α und g nach unten. Das ist der ' +
           'Zielkonflikt aus ' + MV.sec('4.3') + '.',
    computed: 'Die Faktoren zwischen den beiden Plattformen, aus den Zahlen von ' +
              MV.tab('materials') + ' gerechnet: Δ_Mutter ' + factorText(300, 1300) +
              ', Δ_T ' + factorText(30, 70) + ', τ_Z ≈ × ' +
              MV.fmtRange(fac(TAU.alMax, TAU.pb), fac(TAU.alMin, TAU.pb),
                          1, null, MV.fmtExp) + '.',
    params: function () {
      return '<ul><li>angezeigt: ' + (state.material === 0 ? 'nur InAs–Al'
        : state.material === 1 ? 'nur InAs–Pb' : 'beide Plattformen') + '</li>' +
        '<li>Δ_Mutter ' + SCALES[0].al + ' → ' + SCALES[0].pb + ' µeV, ' +
        factorText(300, 1300) + '</li>' +
        '<li>Δ_T ' + SCALES[2].al + ' → ' + SCALES[2].pb + ' µeV, ' +
        factorText(30, 70) + '</li>' +
        '<li>τ_Z ' + TAU.alMin + '–' + TAU.alMax + ' ms → ' +
        MV.fmt(TAU.pb / 1000, 0) + ' ± ' + MV.fmt(TAU.pbErr / 1000, 0) + ' s</li>' +
        '<li>k_BT bei ' + T_MK + ' mK ≈ ' + MV.fmt(K_B * T_MK / 1000, 1) + ' µeV</li></ul>';
    },
    numerics: numericsHtml(),
    convention: 'Beide Achsen sind logarithmisch — anders wären vier ' +
                'Größenordnungen nicht in ein Bild zu bringen. E_M ist eine ' +
                '<b>obere Schranke</b> und kein Wert; das zeigt der Pfeil nach ' +
                'links. k_BT gilt für beide Plattformen und steht deshalb als ' +
                'senkrechte Linie. ' + SCHEMATIC_NOTE,
    reference: MV.secOfWork('6.1') + ' und ' + MV.tab('materials') + '. Die ' +
               'Pb-Werte stammen aus [9], die Al-Werte aus ' + SRC_AL +
               '. Quelle: ' + SRC + '. Gl. (1) von [9] ' +
               'lautet H = E_M iγ₁γ₂ + …; die Aufspaltung der beiden ' +
               'Paritätszustände ist dort 2E_M, während die Arbeit δE = E_M ' +
               'setzt. Diese Seite zeigt die Zahlen so, wie die Arbeit sie angibt.',
    reading: '<p><b>Schalte zwischen Al und Pb um:</b> Die Lücken wachsen um ' +
             'Faktoren zwischen zwei und vier — Δ_Mutter ' + factorText(300, 1300) +
             ', Δ_T nur ' + factorText(30, 70) + ' —, die Lebensdauer dagegen um ' +
             'mehr als drei Größenordnungen. Der Materialwechsel wirkt also weit ' +
             'stärker auf den Vergiftungskanal als auf den Anregungskanal. Die ' +
             'Autoren führen das ausdrücklich auf die geringere Dichte von ' +
             'Nichtgleichgewichts-Quasiteilchen zurück und nicht auf die größere ' +
             'topologische Lücke.</p>'
  };

  /* ====================================================================
     9 — Aufbau
     ==================================================================== */

  function buildControls(mount) {
    var box = MV.ui.el('div', { class: 'controls' });

    segZoom = MV.ui.segmented(box, {
      label: 'Zoomstufe:',
      options: ZOOM.map(function (z, i) { return (i + 1) + ' ' + z.name; }),
      value: state.zoom,
      onChange: function (i) { state.zoom = i; refreshA(); }
    });

    groundLine = MV.ui.el('div', { class: 'callout callout--soft' });
    box.appendChild(groundLine);

    segMaterial = MV.ui.segmented(box, {
      label: 'Material in Panel C:',
      options: ['nur Al', 'nur Pb', 'beide'],
      value: state.material,
      onChange: function (i) { state.material = i; refreshC(); }
    });

    box.appendChild(MV.ui.el('p', {
      class: 'note',
      html: 'Quelle aller Angaben zum Bauelement: ' + SRC + '. ' + Z_NOTE
    }));

    mount.appendChild(box);
  }

  /* Hervorhebung der Zeilen zur Teststruktur in Panel B. Ein Rahmen, sonst
     nichts: kein Status, keine Ampel, keine Wertung (CLAUDE.md §9, §13).
     Sie zeigt, welche Zeilen zu welchem Draht gehören — mehr sagt sie nicht. */
  function testRows() {
    return testCell ? testCell.querySelectorAll('.objlist li') : [];
  }

  function markTestRows(n) {
    Array.prototype.forEach.call(testRows(), function (li, i) {
      if (i < n) li.classList.add('is-marked');
      else li.classList.remove('is-marked');
    });
  }

  function markedCount() {
    return testCell ? testCell.querySelectorAll('.objlist li.is-marked').length : 0;
  }

  /* Schnittstelle für die Vorführungen (demo.js) — nur vorhandene Schalter. */
  function buildDemoHandle() {
    MV.demo.provide('pbtetron', MV.demo.handle({
      state: state,
      refresh: refreshAll,
      stop: function () { markTestRows(0); },
      segments: { zoom: segZoom, material: segMaterial },
      extra: {
        zoomName: function () { return ZOOM[state.zoom].name; },
        zoomCount: function () { return ZOOM.length; },
        /* Die Zeilen, wie sie wirklich unter Panel A stehen — gelesen, nicht
           noch einmal zusammengesetzt. */
        readoutStatus: function () {
          var e = readoutA && readoutA.root.querySelector('.status');
          return e ? e.textContent : '';
        },
        panelCStatus: function () {
          var e = panelC && panelC.root.querySelector('.panel__status');
          return e ? e.textContent : '';
        },
        markTestRows: markTestRows,
        markedCount: markedCount,
        testRowCount: function () { return testRows().length; },
        testRowText: function (i) {
          var r = testRows();
          return r[i] ? r[i].textContent : '';
        },
        /* Die Faktoren aus denselben Zahlen wie die Selbstprüfung. */
        factors: function () {
          return {
            deltaT: fac(SCALES[2].al, SCALES[2].pb),
            deltaTText: factorText(SCALES[2].al, SCALES[2].pb),
            tauLo: fac(TAU.alMax, TAU.pb),
            tauHi: fac(TAU.alMin, TAU.pb)
          };
        },
        materialColumns: function () {
          var host = panelC && panelC.body.querySelector('.materials-host');
          return host ? host.querySelectorAll('thead th.num').length : 0;
        }
      }
    }));
  }

  MV.pbTetron = {
    check: check,

    init: function (mount) {
      if (!check.passed) {
        var bad = check.parts.filter(function (p) { return p.ok === false; });
        var warn = MV.ui.el('div', {
          class: 'readout',
          html: '<b>Warnung:</b> ' + bad.length + ' Prüfung(en) dieses Reiters sind ' +
                'fehlgeschlagen: ' + bad.map(function (p) { return p.name; }).join('; ')
        });
        warn.style.borderColor = MV.palette.critical;
        warn.style.color = MV.palette.critical;
        mount.appendChild(warn);
      }

      buildControls(mount);

      var panels = MV.ui.el('div', { class: 'panels' });
      mount.appendChild(panels);

      /* -------------------------------------------------------- Panel A */
      panelA = MV.ui.panel(panels, {
        title: 'Panel A — Vom Feld zur Messschleife',
        quote: 'pbtetron.A',
        lead: 'Drei Zoomstufen: die Elementarzelle, das Tetron BA und die ' +
              'Messschleife an seinem oberen Draht. <b>Zoome durch und achte ' +
              'darauf, an welchem Tetron gemessen wird und ob es geerdet ist.</b>',
        wide: true, aspect: 0.50, controls: true,
        margin: { l: 70, r: 70, t: 30, b: 40 },
        info: INFO_A
      });

      readoutA = MV.ui.readout(panelA.controls, {
        items: [
          { key: 'zoom', label: 'Zoomstufe' },
          { key: 'device', label: 'Bauelement' },
          { key: 'ground', label: 'Ladungszustand' },
          { key: 'wire', label: 'vermessener Draht' }
        ]
      });

      /* -------------------------------------------------------- Panel B */
      panelB = MV.ui.panel(panels, {
        title: 'Panel B — Zwei Objekte',
        quote: 'pbtetron.B',
        lead: 'Links die separate Teststruktur, rechts der vermessene Draht im ' +
              'Tetron BA. <b>Vergleiche die beiden Spalten: Welche Größe stammt ' +
              'von welchem Draht?</b> Diese Seite bewertet nicht — sie hält nur ' +
              'auseinander.',
        wide: true, canvas: false, controls: true,
        info: INFO_B
      });

      var row = MV.ui.el('div', { class: 'triptych triptych--two' });
      var wires = ['test', 'device'].map(function (k) {
        var o = OBJECTS[k];
        var c = MV.ui.el('div', { class: 'triptych__cell' });
        c.appendChild(MV.ui.el('h4', { class: 'triptych__title', html: o.title }));
        c.appendChild(MV.ui.el('p', { class: 'triptych__sub', html: o.sub }));
        var cv = MV.ui.el('canvas', { class: 'panel__canvas' });
        c.appendChild(cv);
        c.appendChild(MV.ui.el('div', { html: objectHtml(o) }));
        row.appendChild(c);
        if (k === 'test') testCell = c;
        return cv;
      });
      panelB.body.appendChild(row);

      MV.pbTetron._wireTest = new MV.Plot(wires[0], {
        aspect: 0.16, margin: { l: 10, r: 10, t: 22, b: 22 }
      });
      MV.pbTetron._wireDev = new MV.Plot(wires[1], {
        aspect: 0.16, margin: { l: 10, r: 10, t: 22, b: 22 }
      });

      panelB.controls.appendChild(MV.ui.el('p', {
        class: 'note',
        html: 'Keine Statusspalte, kein Urteil. Was aus dieser Aufteilung folgt, ' +
              'wird in ' + MV.secOfWork('6.5') + ' bewertet.'
      }));

      /* -------------------------------------------------------- Panel C */
      panelC = MV.ui.panel(panels, {
        title: 'Panel C — Die Lückenkaskade',
        quote: 'pbtetron.C',
        lead: 'Energien links, Lebensdauer rechts — zwei Achsen, weil es zwei ' +
              'Einheiten sind, beide logarithmisch. <b>Schalte zwischen Al und ' +
              'Pb um</b> und sieh zu, wo der Gewinn liegt.',
        wide: true, canvas: false, controls: true,
        info: INFO_C
      });

      var crow = MV.ui.el('div', { class: 'triptych triptych--two' });
      function ccell(title, sub) {
        var c = MV.ui.el('div', { class: 'triptych__cell' });
        c.appendChild(MV.ui.el('h4', { class: 'triptych__title', html: title }));
        c.appendChild(MV.ui.el('p', { class: 'triptych__sub', html: sub }));
        var cv = MV.ui.el('canvas', { class: 'panel__canvas' });
        c.appendChild(cv);
        var lg = MV.ui.el('div', { class: 'legend-inline' });
        c.appendChild(lg);
        crow.appendChild(c);
        return { cv: cv, lg: lg };
      }
      var ce = ccell('Energieskalen', 'in µeV, logarithmisch');
      var ct = ccell('Paritätslebensdauer τ_Z', 'logarithmisch, K4');
      panelC.body.appendChild(crow);
      panelC.body.appendChild(MV.ui.el('div', { class: 'materials-host' }));

      plotEnergy = new MV.Plot(ce.cv, { aspect: 0.62, margin: { l: 96, r: 16, t: 12, b: 34 } });
      plotEnergy.setLegendHost(ce.lg);
      plotTime = new MV.Plot(ct.cv, { aspect: 0.62, margin: { l: 40, r: 16, t: 12, b: 34 } });
      plotTime.setLegendHost(ct.lg);

      buildDemoHandle();
      MV.demo.onHandInput(function () { markTestRows(0); });

      panelC.controls.appendChild(MV.ui.el('p', {
        class: 'note',
        html: 'Die Autoren führen den Zuwachs der Lebensdauer auf die geringere ' +
              'Dichte von Nichtgleichgewichts-Quasiteilchen zurück — größere ' +
              'Mutterlücke, schnellere Rekombination in Blei — und ausdrücklich ' +
              '<b>nicht</b> auf die größere topologische Lücke. ' +
              MV.secOfWork('6.1') + ' hält fest, dass die Vergiftungsrate damit ' +
              'eine eigenständige, unabhängig zu messende Kenngröße ist.'
      }));
    },

    draw: function () {
      if (!panelA) return;
      refreshAll();
    },

    /* Vom Router beim Verlassen des Reiters gerufen: die Hervorhebung aus
       Panel B ist Teil einer Vorführung und bleibt nicht stehen. */
    stop: function () { markTestRows(0); }
  };

}(MV));
