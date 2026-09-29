/* ------------------------------------------------------------------------
   plot.js — shared drawing layer for all tabs.

   No tab brings its own axis routine; everything goes through MV.Plot.
   Plain <script> (no ES modules) so that index.html also works when it is
   opened directly from the file system.
   ------------------------------------------------------------------------ */

var MV = window.MV || (window.MV = {});

(function (MV) {
  'use strict';

  var FONT = '"Segoe UI", "Helvetica Neue", Arial, sans-serif';

  /* Colour roles. Duplicated in style.css as CSS custom properties; both must
     stay in sync. The role colours (topological / trivial / critical) are the
     same in every tab and are never the only carrier of information — every
     place that uses them also prints the state as text. */
  MV.palette = {
    ink:        '#1b1b1b',
    muted:      '#6a7178',
    axis:       '#3a3f45',
    line:       '#d6d9dc',
    grid:       '#e9ebed',
    bg:         '#ffffff',
    bgSoft:     '#f7f8f9',

    topological:     '#9d1b3c',
    topologicalSoft: '#f6e7eb',
    trivial:         '#55606a',
    trivialSoft:     '#eef0f1',
    critical:        '#b5670c',
    criticalSoft:    '#fbf0e2',

    curve1: '#2b5d9e',
    curve2: '#d97a26',
    curve3: '#3f7d54',

    /* Neutral highlight for a region that is neither topological nor trivial
       by itself — the Zeeman window in tab 2. Deliberately not one of the
       three state colours. */
    neutralSoft: '#edeff2'
  };

  /* ------------------------------------------------------------ formatting */

  /**
   * Eine Zahl für die Oberfläche. **Dezimalkomma**, wie in der Arbeit, und
   * der Minusstrich U+2212 statt des Bindestrichs.
   *
   * Gleichungs-, Abschnitts-, Abbildungs- und Tabellennummern laufen nicht
   * hier durch — sie kommen aus `equations.js` und behalten ihren Punkt.
   */
  MV.fmt = function (v, digits) {
    if (v === null || v === undefined || !isFinite(v)) return '—';
    var s = v.toFixed(digits === undefined ? 2 : digits);
    if (/^-0(\.0*)?$/.test(s)) s = s.slice(1);
    return s.replace('-', '−').replace('.', ',');
  };

  /* Signed value. A value that rounds to zero gets no sign — "+0.00" would
     suggest a sign where the quantity has none, which matters at the critical
     point where ε₀ or ε_π vanishes exactly. */
  /* Die Null bekommt kein Vorzeichen. Geprüft wird gegen die bereits
     formatierte Zeichenkette, die ein Dezimalkomma trägt — mit einem Punkt im
     Muster liefe die Prüfung ins Leere und aus 0 würde „+0,000". */
  MV.fmtSigned = function (v, digits) {
    var s = MV.fmt(v, digits);
    if (s === '—' || s.charAt(0) === '−') return s;
    if (/^0(,0*)?$/.test(s)) return s;
    return '+' + s;
  };

  MV.piTicks = function () {
    var P = Math.PI;
    return [
      { v: -P,     label: '−π'   },
      { v: -P / 2, label: '−π/2' },
      { v: 0,      label: '0'              },
      { v: P / 2,  label: 'π/2'       },
      { v: P,      label: 'π'         }
    ];
  };

  /* Round up to the next multiple of `step`. Used to keep an auto-scaled
     frame from jittering while a slider moves: the view then changes in
     visible discrete steps instead of continuously. */
  MV.ceilTo = function (v, step) {
    return Math.ceil(v / step - 1e-9) * step;
  };

  /* Superscript digits, for decade labels like 10⁻¹². */
  MV.sup = function (n) {
    var map = { '-': '⁻', '0': '⁰', '1': '¹', '2': '²', '3': '³', '4': '⁴',
                '5': '⁵', '6': '⁶', '7': '⁷', '8': '⁸', '9': '⁹' };
    return String(n).split('').map(function (ch) { return map[ch] || ch; }).join('');
  };

  /**
   * Eine Zehnerpotenz als „1,6 · 10⁻¹⁶" statt „1.6e-16". Für Größen, die über
   * viele Dekaden laufen: Abweichungen, Faktoren, Energien weit unter einem
   * µeV. Landet die gerundete Mantisse auf 10, rückt die Dekade eine weiter,
   * damit nie „10,0 · 10³" dasteht.
   */
  MV.fmtExp = function (v, digits) {
    if (v === null || v === undefined || !isFinite(v)) return '—';
    if (v === 0) return '0';
    var d = digits === undefined ? 1 : digits;
    var e = Math.floor(Math.log(Math.abs(v)) / Math.LN10);
    var m = v / Math.pow(10, e);
    var p = Math.pow(10, d);
    if (Math.round(Math.abs(m) * p) / p >= 10) { m /= 10; e += 1; }
    return MV.fmt(m, d) + ' · 10' + MV.sup(e);
  };

  /**
   * Ein Bereich: „1,42–5,80 µm", Einheit einmal am Ende.
   *
   * Statt des Halbgeviertstrichs trennt „bis", wenn er falsch gelesen würde:
   * bei einer negativen Grenze („−3–3" sähe aus wie eine Rechnung) und
   * zwischen Zehnerpotenzen, wo der Strich direkt hinter einer Hochzahl die
   * Fortsetzung des Exponenten zu sein scheint.
   *
   * `fmt` ist der Formatierer für die Grenzen, normalerweise MV.fmt.
   */
  MV.fmtRange = function (lo, hi, digits, unit, fmt) {
    var f = fmt || MV.fmt;
    var sep = (lo < 0 || hi < 0 || f === MV.fmtExp) ? ' bis ' : '–';
    return f(lo, digits) + sep + f(hi, digits) + (unit ? ' ' + unit : '');
  };

  /**
   * Ticks for an axis on which log₁₀(value) is plotted linearly, labelled as
   * powers of ten. lo and hi are the exponents.
   */
  MV.decadeTicks = function (lo, hi, maxTicks) {
    lo = Math.floor(lo); hi = Math.ceil(hi);
    var stepDec = Math.max(1, Math.ceil((hi - lo) / (maxTicks || 6)));
    var out = [], e;
    for (e = Math.ceil(lo / stepDec) * stepDec; e <= hi; e += stepDec) {
      out.push({ v: e, label: '10' + MV.sup(e) });
    }
    return out;
  };

  /* "Nice" tick positions covering [lo, hi] with roughly `target` steps. */
  MV.niceTicks = function (lo, hi, target) {
    target = target || 6;
    var span = hi - lo;
    if (!(span > 0)) return [{ v: lo, label: String(lo) }];
    var raw = span / target;
    var mag = Math.pow(10, Math.floor(Math.log(raw) / Math.LN10));
    var norm = raw / mag;
    var step = (norm < 1.5 ? 1 : norm < 3 ? 2 : norm < 7 ? 5 : 10) * mag;
    var decimals = Math.max(0, -Math.floor(Math.log(step) / Math.LN10 + 1e-9));
    var out = [];
    var start = Math.ceil(lo / step - 1e-9) * step;
    for (var v = start; v <= hi + step * 1e-9; v += step) {
      var val = Math.abs(v) < step * 1e-9 ? 0 : v;
      out.push({ v: val, label: MV.fmt(val, decimals) });
    }
    return out;
  };

  /* ----------------------------------------------------------------- Plot */

  /**
   * One cartesian panel on one <canvas>.
   *
   * Per frame:
   *   plot.begin({ x:[x0,x1], y:[y0,y1], xLabel:'k', yLabel:'E(k)' });
   *   ... spans, curves, markers ...
   *   plot.finish();
   *
   * begin()  clears, sizes the backing store for the current devicePixelRatio
   *          and draws the grid.
   * finish() draws the frame, ticks, axis labels and the legend on top.
   */
  function Plot(canvas, opts) {
    opts = opts || {};
    this.canvas = canvas;
    this.ctx = canvas.getContext('2d');
    this.aspect = opts.aspect || 0.68;
    this.margin = { l: 58, r: 16, t: 14, b: 44 };
    if (opts.margin) {
      for (var k in opts.margin) {
        if (Object.prototype.hasOwnProperty.call(opts.margin, k)) this.margin[k] = opts.margin[k];
      }
    }
    this.legendItems = [];
    this.cfg = {};
  }

  Plot.prototype._resize = function () {
    var dpr = window.devicePixelRatio || 1;
    var cssW = this.canvas.clientWidth;
    if (!cssW && this.canvas.parentNode) cssW = this.canvas.parentNode.clientWidth;
    cssW = Math.max(240, cssW || 480);
    var cssH = Math.round(cssW * this.aspect);
    var needW = Math.round(cssW * dpr);
    var needH = Math.round(cssH * dpr);
    if (this.canvas.width !== needW || this.canvas.height !== needH) {
      this.canvas.width = needW;
      this.canvas.height = needH;
    }
    this.canvas.style.height = cssH + 'px';
    this.ctx.setTransform(dpr, 0, 0, dpr, 0, 0);
    this.w = cssW;
    this.h = cssH;
  };

  /**
   * Size of the data area in pixels, available before begin(). Needed when a
   * panel wants the same number of data units per pixel on both axes (an
   * orbit that is a circle has to look like a circle).
   */
  Plot.prototype.measure = function () {
    this._resize();
    return {
      pw: (this.w - this.margin.l - this.margin.r),
      ph: (this.h - this.margin.t - this.margin.b)
    };
  };

  Plot.prototype.begin = function (cfg) {
    this.cfg = cfg = cfg || {};
    this._resize();
    var ctx = this.ctx;
    var m = this.margin;

    ctx.clearRect(0, 0, this.w, this.h);
    ctx.fillStyle = MV.palette.bg;
    ctx.fillRect(0, 0, this.w, this.h);

    this.x0 = cfg.x[0]; this.x1 = cfg.x[1];
    this.y0 = cfg.y[0]; this.y1 = cfg.y[1];

    this.px0 = m.l;
    this.px1 = this.w - m.r;
    this.py0 = m.t;
    this.py1 = this.h - m.b;

    this.xTicks = cfg.xTicks || MV.niceTicks(this.x0, this.x1, cfg.xTickCount || 6);
    this.yTicks = cfg.yTicks || MV.niceTicks(this.y0, this.y1, cfg.yTickCount || 5);
    this.legendItems = [];

    if (cfg.grid !== false) this._grid();
    return this;
  };

  Plot.prototype.X = function (x) {
    return this.px0 + (x - this.x0) / (this.x1 - this.x0) * (this.px1 - this.px0);
  };

  Plot.prototype.Y = function (y) {
    return this.py1 - (y - this.y0) / (this.y1 - this.y0) * (this.py1 - this.py0);
  };

  Plot.prototype._grid = function () {
    var ctx = this.ctx, i, p;
    ctx.save();
    ctx.strokeStyle = MV.palette.grid;
    ctx.lineWidth = 1;
    ctx.beginPath();
    for (i = 0; i < this.xTicks.length; i++) {
      p = Math.round(this.X(this.xTicks[i].v)) + 0.5;
      ctx.moveTo(p, this.py0); ctx.lineTo(p, this.py1);
    }
    for (i = 0; i < this.yTicks.length; i++) {
      p = Math.round(this.Y(this.yTicks[i].v)) + 0.5;
      ctx.moveTo(this.px0, p); ctx.lineTo(this.px1, p);
    }
    ctx.stroke();
    ctx.restore();
  };

  Plot.prototype.clip = function () {
    var ctx = this.ctx;
    ctx.save();
    ctx.beginPath();
    ctx.rect(this.px0, this.py0, this.px1 - this.px0, this.py1 - this.py0);
    ctx.clip();
    return this;
  };

  Plot.prototype.unclip = function () { this.ctx.restore(); return this; };

  /* -------------------------------------------------------------- drawing */

  /** pts: array of [x, y]; a null entry breaks the line. */
  Plot.prototype.curve = function (pts, o) {
    o = o || {};
    var ctx = this.ctx, started = false, i, p;
    this.clip();
    ctx.strokeStyle = o.color || MV.palette.ink;
    ctx.lineWidth = o.width || 1.8;
    ctx.lineJoin = 'round';
    ctx.lineCap = 'round';
    ctx.setLineDash(o.dash || []);
    if (o.alpha !== undefined) ctx.globalAlpha = o.alpha;
    ctx.beginPath();
    for (i = 0; i < pts.length; i++) {
      p = pts[i];
      if (!p || !isFinite(p[0]) || !isFinite(p[1])) { started = false; continue; }
      if (!started) { ctx.moveTo(this.X(p[0]), this.Y(p[1])); started = true; }
      else ctx.lineTo(this.X(p[0]), this.Y(p[1]));
    }
    ctx.stroke();
    ctx.setLineDash([]);
    ctx.globalAlpha = 1;
    this.unclip();
    if (o.legend) {
      this.legendItems.push({ label: o.legend, color: o.color || MV.palette.ink, dash: o.dash });
    }
    return this;
  };

  /** Sample f over [a, b] with n+1 points and draw it. */
  Plot.prototype.fn = function (f, a, b, n, o) {
    var pts = [], i, x;
    n = n || 400;
    for (i = 0; i <= n; i++) {
      x = a + (b - a) * i / n;
      pts.push([x, f(x)]);
    }
    return this.curve(pts, o);
  };

  /**
   * Stem plot: a vertical line from a baseline up to each point, all in one
   * path, with an optional dot on top. Used for lattice-site weights, where a
   * connecting curve would suggest values between the sites.
   */
  Plot.prototype.stems = function (pts, o) {
    o = o || {};
    var ctx = this.ctx, i, p, base = this.Y(o.base === undefined ? 0 : o.base);
    this.clip();
    ctx.strokeStyle = o.color || MV.palette.ink;
    ctx.lineWidth = o.width || 2;
    ctx.lineCap = 'butt';
    ctx.beginPath();
    for (i = 0; i < pts.length; i++) {
      p = pts[i];
      if (!p || !isFinite(p[0]) || !isFinite(p[1])) continue;
      ctx.moveTo(this.X(p[0]), base);
      ctx.lineTo(this.X(p[0]), this.Y(p[1]));
    }
    ctx.stroke();
    if (o.dot) {
      ctx.fillStyle = o.color || MV.palette.ink;
      for (i = 0; i < pts.length; i++) {
        p = pts[i];
        if (!p || !isFinite(p[0]) || !isFinite(p[1])) continue;
        ctx.beginPath();
        ctx.arc(this.X(p[0]), this.Y(p[1]), o.dot, 0, 2 * Math.PI);
        ctx.fill();
      }
    }
    this.unclip();
    if (o.legend) {
      this.legendItems.push({ label: o.legend, color: o.color || MV.palette.ink });
    }
    return this;
  };

  /** Filled and/or outlined rectangle in data coordinates. */
  Plot.prototype.rect = function (x0, y0, x1, y1, o) {
    o = o || {};
    var ctx = this.ctx;
    var ax = this.X(Math.min(x0, x1)), bx = this.X(Math.max(x0, x1));
    var ay = this.Y(Math.max(y0, y1)), by = this.Y(Math.min(y0, y1));
    this.clip();
    if (o.radius) {
      var r = Math.min(o.radius, (bx - ax) / 2, (by - ay) / 2);
      ctx.beginPath();
      ctx.moveTo(ax + r, ay);
      ctx.lineTo(bx - r, ay); ctx.quadraticCurveTo(bx, ay, bx, ay + r);
      ctx.lineTo(bx, by - r); ctx.quadraticCurveTo(bx, by, bx - r, by);
      ctx.lineTo(ax + r, by); ctx.quadraticCurveTo(ax, by, ax, by - r);
      ctx.lineTo(ax, ay + r); ctx.quadraticCurveTo(ax, ay, ax + r, ay);
      ctx.closePath();
    } else {
      ctx.beginPath();
      ctx.rect(ax, ay, bx - ax, by - ay);
    }
    if (o.fill) { ctx.fillStyle = o.fill; ctx.fill(); }
    if (o.stroke) {
      ctx.strokeStyle = o.stroke;
      ctx.lineWidth = o.width || 1;
      ctx.setLineDash(o.dash || []);
      ctx.stroke();
      ctx.setLineDash([]);
    }
    this.unclip();
    if (o.legend) {
      this.legendItems.push({
        label: o.legend, color: o.fill || o.stroke, swatch: 'box',
        fill: o.fill || o.stroke
      });
    }
    return this;
  };

  Plot.prototype.hline = function (y, o) {
    o = o || {};
    var ctx = this.ctx, py = Math.round(this.Y(y)) + 0.5;
    if (py < this.py0 - 1 || py > this.py1 + 1) return this;
    ctx.save();
    ctx.strokeStyle = o.color || MV.palette.axis;
    ctx.lineWidth = o.width || 1;
    ctx.setLineDash(o.dash || []);
    ctx.beginPath();
    ctx.moveTo(this.px0, py); ctx.lineTo(this.px1, py);
    ctx.stroke();
    ctx.restore();
    if (o.label) {
      this.label(o.label,
        o.labelAt === 'left' ? this.px0 + 6 : this.px1 - 6, py - 5,
        {
          color: o.labelColor || o.color || MV.palette.muted,
          align: o.labelAt === 'left' ? 'left' : 'right',
          baseline: 'bottom',
          box: o.labelBox !== false
        });
    }
    if (o.legend) {
      this.legendItems.push({ label: o.legend, color: o.color || MV.palette.axis, dash: o.dash });
    }
    return this;
  };

  Plot.prototype.vline = function (x, o) {
    o = o || {};
    var ctx = this.ctx, px = Math.round(this.X(x)) + 0.5;
    if (px < this.px0 - 1 || px > this.px1 + 1) return this;
    ctx.save();
    ctx.strokeStyle = o.color || MV.palette.axis;
    ctx.lineWidth = o.width || 1;
    ctx.setLineDash(o.dash || []);
    ctx.beginPath();
    ctx.moveTo(px, this.py0); ctx.lineTo(px, this.py1);
    ctx.stroke();
    ctx.restore();
    if (o.legend) {
      this.legendItems.push({ label: o.legend, color: o.color || MV.palette.axis, dash: o.dash });
    }
    return this;
  };

  /** Horizontal band between two y values, spanning the full width. */
  Plot.prototype.hspan = function (ya, yb, o) {
    o = o || {};
    var ctx = this.ctx;
    var pa = this.Y(Math.max(ya, yb));
    var pb = this.Y(Math.min(ya, yb));
    pa = Math.max(pa, this.py0);
    pb = Math.min(pb, this.py1);
    if (pb <= pa) return this;
    ctx.save();
    ctx.fillStyle = o.fill || MV.palette.criticalSoft;
    ctx.fillRect(this.px0, pa, this.px1 - this.px0, pb - pa);
    if (o.edge) {
      ctx.strokeStyle = o.edge;
      ctx.lineWidth = 1;
      ctx.setLineDash(o.dash || [3, 3]);
      ctx.beginPath();
      ctx.moveTo(this.px0, Math.round(pa) + 0.5); ctx.lineTo(this.px1, Math.round(pa) + 0.5);
      ctx.moveTo(this.px0, Math.round(pb) + 0.5); ctx.lineTo(this.px1, Math.round(pb) + 0.5);
      ctx.stroke();
    }
    ctx.restore();
    if (o.legend) {
      this.legendItems.push({ label: o.legend, color: o.edge || o.fill, swatch: 'box', fill: o.fill });
    }
    return this;
  };

  /** Vertical band between two x values, spanning the full height. */
  Plot.prototype.vspan = function (xa, xb, o) {
    o = o || {};
    var ctx = this.ctx;
    var pa = this.X(Math.min(xa, xb));
    var pb = this.X(Math.max(xa, xb));
    pa = Math.max(pa, this.px0);
    pb = Math.min(pb, this.px1);
    if (pb <= pa) return this;
    ctx.save();
    ctx.fillStyle = o.fill || MV.palette.neutralSoft;
    ctx.fillRect(pa, this.py0, pb - pa, this.py1 - this.py0);
    if (o.edge) {
      ctx.strokeStyle = o.edge;
      ctx.lineWidth = 1;
      ctx.setLineDash(o.dash || [3, 3]);
      ctx.beginPath();
      ctx.moveTo(Math.round(pa) + 0.5, this.py0); ctx.lineTo(Math.round(pa) + 0.5, this.py1);
      ctx.moveTo(Math.round(pb) + 0.5, this.py0); ctx.lineTo(Math.round(pb) + 0.5, this.py1);
      ctx.stroke();
    }
    ctx.restore();
    if (o.legend) {
      this.legendItems.push({ label: o.legend, color: o.edge || o.fill, swatch: 'box', fill: o.fill });
    }
    return this;
  };

  /** Fill between a curve and the top edge of the frame. */
  Plot.prototype.areaAbove = function (pts, o) {
    return this._area(pts, o, this.py0);
  };

  /** Fill between a curve and the bottom edge of the frame. */
  Plot.prototype.areaBelow = function (pts, o) {
    return this._area(pts, o, this.py1);
  };

  Plot.prototype._area = function (pts, o, edgePy) {
    o = o || {};
    var ctx = this.ctx, i;
    if (!pts.length) return this;
    this.clip();
    ctx.fillStyle = o.fill || MV.palette.topologicalSoft;
    ctx.beginPath();
    ctx.moveTo(this.X(pts[0][0]), this.Y(pts[0][1]));
    for (i = 1; i < pts.length; i++) ctx.lineTo(this.X(pts[i][0]), this.Y(pts[i][1]));
    ctx.lineTo(this.X(pts[pts.length - 1][0]), edgePy);
    ctx.lineTo(this.X(pts[0][0]), edgePy);
    ctx.closePath();
    ctx.fill();
    this.unclip();
    if (o.legend) {
      this.legendItems.push({ label: o.legend, color: o.fill, swatch: 'box', fill: o.fill });
    }
    return this;
  };

  /**
   * Point marker with optional label.
   * o.anchor: 'above' | 'below' | 'left' | 'right'   (default 'above')
   * o.shape:  'dot' | 'ring' | 'cross' | 'square'
   */
  Plot.prototype.marker = function (x, y, o) {
    o = o || {};
    var ctx = this.ctx;
    var px = this.X(x), py = this.Y(y);
    var r = o.r || 4;
    var col = o.color || MV.palette.ink;
    if (px < this.px0 - 14 || px > this.px1 + 14 || py < this.py0 - 14 || py > this.py1 + 14) return this;

    ctx.save();
    ctx.lineWidth = o.width || 1.8;
    ctx.strokeStyle = col;
    ctx.fillStyle = o.shape === 'ring' ? MV.palette.bg : col;
    if (o.shape === 'cross') {
      ctx.beginPath();
      ctx.moveTo(px - r, py - r); ctx.lineTo(px + r, py + r);
      ctx.moveTo(px - r, py + r); ctx.lineTo(px + r, py - r);
      ctx.stroke();
    } else if (o.shape === 'square') {
      ctx.fillRect(px - r, py - r, 2 * r, 2 * r);
      ctx.strokeRect(px - r, py - r, 2 * r, 2 * r);
    } else {
      ctx.beginPath();
      ctx.arc(px, py, r, 0, 2 * Math.PI);
      ctx.fill();
      if (o.shape === 'ring') ctx.stroke();
    }
    ctx.restore();

    if (o.label) {
      var a = o.anchor || 'above';
      var dx = a === 'left' ? -(r + 5) : a === 'right' ? (r + 5) : 0;
      var dy = a === 'above' ? -(r + 4) : a === 'below' ? (r + 4) : 0;
      this.label(o.label, px + dx, py + dy, {
        color: o.labelColor || col,
        align: a === 'left' ? 'right' : a === 'right' ? 'left' : 'center',
        baseline: a === 'above' ? 'bottom' : a === 'below' ? 'top' : 'middle',
        box: o.labelBox !== false,
        clampTo: true
      });
    }
    if (o.legend) {
      this.legendItems.push({ label: o.legend, color: col, swatch: 'dot', shape: o.shape });
    }
    return this;
  };

  /**
   * Vertical dimension line with arrow heads at both ends (e.g. the 2*E_Z
   * Zeeman gap). x, ya, yb in data coordinates.
   */
  Plot.prototype.dimension = function (x, ya, yb, text, o) {
    o = o || {};
    var ctx = this.ctx;
    var px = this.X(x), pa = this.Y(ya), pb = this.Y(yb);
    var col = o.color || MV.palette.ink;
    var head = 5;
    this.clip();
    ctx.strokeStyle = col;
    ctx.fillStyle = col;
    ctx.lineWidth = 1.2;
    ctx.beginPath();
    ctx.moveTo(px, pa); ctx.lineTo(px, pb);
    ctx.stroke();
    [[pa, pa < pb ? 1 : -1], [pb, pb < pa ? 1 : -1]].forEach(function (e) {
      ctx.beginPath();
      ctx.moveTo(px, e[0]);
      ctx.lineTo(px - head * 0.7, e[0] + head * e[1]);
      ctx.lineTo(px + head * 0.7, e[0] + head * e[1]);
      ctx.closePath();
      ctx.fill();
    });
    if (o.ticks !== false) {
      ctx.setLineDash([2, 3]);
      ctx.beginPath();
      ctx.moveTo(px - 14, pa); ctx.lineTo(px + 14, pa);
      ctx.moveTo(px - 14, pb); ctx.lineTo(px + 14, pb);
      ctx.stroke();
      ctx.setLineDash([]);
    }
    this.unclip();
    if (text) {
      /* o.offset lifts the caption off the midpoint, which is where a zero
         line or a μ line often runs. */
      this.label(text, px + (o.side === 'left' ? -8 : 8), (pa + pb) / 2 + (o.offset || 0), {
        color: col,
        align: o.side === 'left' ? 'right' : 'left',
        baseline: 'middle',
        box: true
      });
    }
    return this;
  };

  /** Short arrow, direction given in pixels (used for spinor hints). */
  Plot.prototype.arrow = function (x, y, dxPx, dyPx, o) {
    o = o || {};
    var ctx = this.ctx;
    var px = this.X(x), py = this.Y(y);
    var ex = px + dxPx, ey = py + dyPx;
    var len = Math.sqrt(dxPx * dxPx + dyPx * dyPx) || 1;
    var ux = dxPx / len, uy = dyPx / len, head = o.head || 5;
    this.clip();
    ctx.strokeStyle = o.color || MV.palette.ink;
    ctx.fillStyle = o.color || MV.palette.ink;
    ctx.lineWidth = o.width || 1.4;
    ctx.beginPath();
    ctx.moveTo(px, py); ctx.lineTo(ex, ey);
    ctx.stroke();
    ctx.beginPath();
    ctx.moveTo(ex, ey);
    ctx.lineTo(ex - head * (ux + uy * 0.6), ey - head * (uy - ux * 0.6));
    ctx.lineTo(ex - head * (ux - uy * 0.6), ey - head * (uy + ux * 0.6));
    ctx.closePath();
    ctx.fill();
    this.unclip();
    return this;
  };

  /** Text in pixel coordinates. o.box draws a translucent plate underneath. */
  Plot.prototype.label = function (text, px, py, o) {
    o = o || {};
    var ctx = this.ctx;
    ctx.save();
    ctx.font = o.font || '12px ' + FONT;
    ctx.textAlign = o.align || 'left';
    ctx.textBaseline = o.baseline || 'alphabetic';

    if (o.clampTo) {
      var wt = ctx.measureText(text).width;
      var half = o.align === 'center' ? wt / 2 : o.align === 'right' ? wt : 0;
      var lo = this.px0 + (o.align === 'left' ? 2 : half + 2);
      var hi = this.px1 - (o.align === 'right' ? 2 : (o.align === 'center' ? wt / 2 : wt) + 2);
      if (px < lo) px = lo;
      if (px > hi) px = hi;
      py = Math.min(Math.max(py, this.py0 + 12), this.py1 - 2);
    }

    if (o.box) {
      var w = ctx.measureText(text).width;
      var h = 13;
      var bx = o.align === 'right' ? px - w - 3 : o.align === 'center' ? px - w / 2 - 3 : px - 3;
      var by = o.baseline === 'bottom' ? py - h - 1 : o.baseline === 'middle' ? py - h / 2 - 1 : py - 1;
      ctx.fillStyle = 'rgba(255,255,255,0.84)';
      ctx.fillRect(bx, by, w + 6, h + 3);
    }
    ctx.fillStyle = o.color || MV.palette.ink;
    ctx.fillText(text, px, py);
    ctx.restore();
    return this;
  };

  /** Text in data coordinates. */
  Plot.prototype.text = function (text, x, y, o) {
    return this.label(text, this.X(x), this.Y(y), o);
  };

  /* ---------------------------------------------------------------- frame */

  Plot.prototype.finish = function () {
    var ctx = this.ctx, cfg = this.cfg, i, t, p;

    ctx.save();
    ctx.strokeStyle = MV.palette.axis;
    ctx.lineWidth = 1;
    /* cfg.frame === false: a schematic drawing wants no axis box around it */
    if (cfg.frame !== false) {
      ctx.strokeRect(
        Math.round(this.px0) + 0.5, Math.round(this.py0) + 0.5,
        Math.round(this.px1 - this.px0), Math.round(this.py1 - this.py0)
      );
    }

    ctx.beginPath();
    for (i = 0; i < this.xTicks.length; i++) {
      p = Math.round(this.X(this.xTicks[i].v)) + 0.5;
      ctx.moveTo(p, this.py1); ctx.lineTo(p, this.py1 + 4);
    }
    for (i = 0; i < this.yTicks.length; i++) {
      p = Math.round(this.Y(this.yTicks[i].v)) + 0.5;
      ctx.moveTo(this.px0 - 4, p); ctx.lineTo(this.px0, p);
    }
    ctx.stroke();
    ctx.restore();

    ctx.save();
    ctx.font = '11px ' + FONT;
    ctx.fillStyle = MV.palette.muted;
    ctx.textAlign = 'center';
    ctx.textBaseline = 'top';
    for (i = 0; i < this.xTicks.length; i++) {
      t = this.xTicks[i];
      ctx.fillText(t.label, this.X(t.v), this.py1 + 7);
    }
    ctx.textAlign = 'right';
    ctx.textBaseline = 'middle';
    for (i = 0; i < this.yTicks.length; i++) {
      t = this.yTicks[i];
      ctx.fillText(t.label, this.px0 - 7, this.Y(t.v));
    }

    /* Optional second row of x tick labels with its own axis caption, for an
       axis that carries two readings of the same quantity — a length in µm
       and the same length in units of ξ, say. */
    if (cfg.xTicks2) {
      ctx.font = '11px ' + FONT;
      ctx.fillStyle = MV.palette.muted;
      ctx.textAlign = 'center';
      ctx.textBaseline = 'top';
      for (i = 0; i < cfg.xTicks2.length; i++) {
        t = cfg.xTicks2[i];
        ctx.fillText(t.label, this.X(t.v), this.py1 + 22);
      }
    }

    ctx.fillStyle = MV.palette.ink;
    ctx.font = '12px ' + FONT;
    if (cfg.xLabel) {
      ctx.textAlign = 'center';
      ctx.textBaseline = 'bottom';
      ctx.fillText(cfg.xLabel, (this.px0 + this.px1) / 2, this.h - 5);
    }
    if (cfg.xLabel2) {
      ctx.font = '11px ' + FONT;
      ctx.fillStyle = MV.palette.muted;
      ctx.textAlign = 'center';
      ctx.textBaseline = 'bottom';
      ctx.fillText(cfg.xLabel2, (this.px0 + this.px1) / 2, this.h - 20);
    }
    if (cfg.yLabel) {
      ctx.save();
      ctx.translate(13, (this.py0 + this.py1) / 2);
      ctx.rotate(-Math.PI / 2);
      ctx.textAlign = 'center';
      ctx.textBaseline = 'top';
      ctx.fillText(cfg.yLabel, 0, 0);
      ctx.restore();
    }
    ctx.restore();

    this._renderLegend();
    return this;
  };

  /**
   * The legend is written into a DOM node below the canvas, not painted into
   * the plot: a wide, flat panel has no corner that stays free for every
   * parameter setting, and an overlaid box hides data. The host node is
   * supplied by MV.ui.panel.
   */
  Plot.prototype.setLegendHost = function (el) { this.legendHost = el; return this; };

  Plot.prototype._renderLegend = function () {
    if (!this.legendHost) return;
    var html = this.legendItems.map(function (it) {
      var kind = it.swatch === 'box' ? 'box'
        : it.swatch === 'dot' ? (it.shape === 'cross' ? 'cross' : it.shape === 'ring' ? 'ring' : 'dot')
        : (it.dash && it.dash.length) ? 'dash' : 'line';
      return '<span class="lg lg--' + kind + '" style="--c:' +
        (it.fill || it.color) + '">' + it.label + '</span>';
    }).join('');
    this.legendHost.innerHTML = html;
  };

  MV.Plot = Plot;

  /* Redraw hooks: window resizes, layout shifts and DPI changes go through
     here. A layout shift after the first draw is the common case — adding the
     panels can bring up the page scrollbar, which narrows every canvas after
     it was already sized. */
  MV.redrawHooks = [];
  MV.onResize = function (fn) { MV.redrawHooks.push(fn); };

  var redrawing = false;
  MV.redrawAll = function () {
    if (redrawing) return;            /* a redraw may itself change the layout */
    redrawing = true;
    try {
      MV.redrawHooks.forEach(function (fn) { fn(); });
    } finally {
      redrawing = false;
    }
  };

  var resizeTimer = null;
  function scheduleRedraw() {
    if (redrawing) return;
    if (resizeTimer) clearTimeout(resizeTimer);
    resizeTimer = setTimeout(MV.redrawAll, 80);
  }

  window.addEventListener('resize', scheduleRedraw);

  if (typeof window.ResizeObserver === 'function') {
    window.addEventListener('DOMContentLoaded', function () {
      new window.ResizeObserver(scheduleRedraw).observe(document.body);
    });
  }

}(MV));
