/* ------------------------------------------------------------------------
   demo.js — die Vorführungen an den Leitsätzen.

   Die Arbeit verweist mit einer Klammer der Form "{Reiter: Handlung —
   Beobachtung}" auf die Seite. Jede dieser Stellen ist ein Satz, der hier als
   Leitsatz über einem Panel steht (leitsaetze.js). Neben dem Leitsatz sitzt
   ein Knopf, der genau diese Handlung vorführt: Startwerte setzen, die
   Handlung als Durchlauf abspielen, im Zielzustand anhalten, das Ergebnis in
   einer Zeile nennen.

   Was eine Vorführung darf und was nicht:

   * Sie bewegt nur Regler, Schalter und Knöpfe, die es schon gibt, über die
     Schnittstellen der Module (`MV.demo.handle`). Keine neuen Formeln, keine
     neuen Regler, keine neue Physik.
   * Sie ist ein Durchlauf, keine Zeitentwicklung. Wo ein Reiter dafür schon
     eine Sprache hat — "μ-Durchlauf, quasistatisch" in #kitaev —, gilt sie
     weiter. Einzige Ausnahme ist #gatter Panel B: dort läuft die vorhandene
     Schrittfolge mit Anzeigetempo (CLAUDE.md §8).
   * Sie prüft am Ende selbst nach, ob der Zielzustand erreicht ist. Schlägt
     die Prüfung fehl, wird der Knopf gesperrt statt die Aussage zu
     beschönigen.

   Bewegung läuft über requestAnimationFrame, auch das Warten zwischen zwei
   Schritten; setInterval kommt nicht vor (CLAUDE.md §12). Bei
   prefers-reduced-motion: reduce gibt es keine Durchläufe — dann wird der
   Endzustand direkt gesetzt und anschließend geprüft.

   Interface language is German, code and comments are English.
   ------------------------------------------------------------------------ */

var MV = window.MV || (window.MV = {});

(function (MV) {
  'use strict';

  var demo = MV.demo = {};

  /* Reiter-Schlüssel des Leitsatzes → Hash-Name des Reiters. Beide sind nicht
     überall gleich: der Leitsatz heißt `pbtetron.A`, der Reiter `#pb-tetron`. */
  var TAB_HASH = {
    kitaev: 'kitaev',
    nanodraht: 'nanodraht',
    energieskalen: 'energieskalen',
    auslese: 'auslese',
    gatter: 'gatter',
    pbtetron: 'pb-tetron'
  };

  /* Reihenfolge der Reiter — die der Arbeit, wie in der Navigation. */
  var TAB_ORDER = ['kitaev', 'nanodraht', 'energieskalen', 'auslese',
                   'gatter', 'pbtetron'];

  var specs = {};      /* Schlüssel → spec                                   */
  var hosts = {};      /* Schlüssel → Liste der Knopfplätze im Markup        */
  var handles = {};    /* Reiter → Schnittstelle des Moduls                  */
  var locked = {};     /* Schlüssel → Grund der Sperre                       */
  var current = null;  /* die eine laufende Vorführung, sonst null           */
  var showTab = null;  /* vom Router gesetzt: Reiter öffnen                  */
  var handInput = [];  /* was beim ersten Handeingriff aufzuräumen ist       */
  var listening = false;

  function tabOf(key) { return key.split('.')[0]; }
  function letterOf(key) { return key.split('.')[1]; }
  function hashOf(key) { return '#' + TAB_HASH[tabOf(key)]; }

  /** "Abschn. 4.1" — die Nummer läuft durch die Prüfung in equations.js. */
  function secShort(n) { MV.sec(n); return 'Abschn. ' + n; }

  function reducedMotion() {
    return !!(window.matchMedia &&
              window.matchMedia('(prefers-reduced-motion: reduce)').matches);
  }

  /* ====================================================================
     Schnittstelle eines Moduls

     Ein Reiter reicht seinen Zustand, seine Regler und sein refresh() herein
     und bekommt dafür set/get/snapshot/restore. Mehr braucht eine Vorführung
     nicht — und mehr soll sie auch nicht können.
     ==================================================================== */

  function clone(v) {
    if (Array.isArray(v)) return v.slice();
    if (v && typeof v === 'object') {
      var o = {}, k;
      for (k in v) if (Object.prototype.hasOwnProperty.call(v, k)) o[k] = v[k];
      return o;
    }
    return v;
  }

  /**
   * cfg: { state, refresh, sliders, segments, setters, stop, extra }
   *
   * `segments` und `setters` gehen vor `sliders`: ein Schalter wird über
   * seinen eigenen Setter gestellt, damit Beschriftung und aria-pressed
   * mitwandern und der Reiter nichts von der Vorführung merkt.
   */
  /**
   * Ein Regler kennt nur seine Schrittweite. Ein Durchlauf, der dazwischen
   * landet, würde den Zustand feiner stellen, als der Regler ihn anzeigt —
   * Zeiger und Zahl liefen auseinander. Deshalb rastet jeder Zwischenwert auf
   * das Raster des Reglers ein; eine Vorführung stellt damit genau das ein,
   * was auch eine Hand einstellen könnte.
   */
  function snapToSlider(sl, v) {
    if (typeof v !== 'number' || !sl) return v;
    /* Der Regler weiß selbst, welche Werte er annehmen kann — auch ein
       logarithmischer, dessen Raster in Dekaden zählt. */
    if (sl.snap) return sl.snap(v);
    if (!sl.input) return v;
    var min = parseFloat(sl.input.min), step = parseFloat(sl.input.step);
    if (!isFinite(min) || !isFinite(step) || step <= 0) return v;
    var snapped = min + Math.round((v - min) / step) * step;
    return Math.round(snapped * 1e9) / 1e9;
  }

  demo.handle = function (cfg) {
    function put(k, v) {
      if (cfg.setters && cfg.setters[k]) { cfg.setters[k](v); return true; }
      if (cfg.segments && cfg.segments[k]) { cfg.segments[k].set(v); return true; }
      var sl = cfg.sliders && cfg.sliders[k];
      if (sl) v = snapToSlider(sl, v);
      cfg.state[k] = v;
      if (sl) sl.set(v, true);
      return false;
    }

    var h = {
      state: cfg.state,
      get: function (k) { return cfg.state[k]; },
      set: function (k, v) { if (!put(k, v)) cfg.refresh(); },
      /** Mehrere Größen auf einmal — ein einziges refresh() am Schluss. */
      setAll: function (obj) {
        Object.keys(obj).forEach(function (k) { put(k, obj[k]); });
        cfg.refresh();
      },
      refresh: function () { cfg.refresh(); },
      stop: cfg.stop || function () {},
      snapshot: function () {
        var s = {}, k;
        for (k in cfg.state) {
          if (Object.prototype.hasOwnProperty.call(cfg.state, k)) s[k] = clone(cfg.state[k]);
        }
        return s;
      },
      restore: function (snap) {
        h.stop();
        Object.keys(snap).forEach(function (k) { put(k, snap[k]); });
        cfg.refresh();
      }
    };
    if (cfg.extra) {
      Object.keys(cfg.extra).forEach(function (k) { h[k] = cfg.extra[k]; });
    }
    return h;
  };

  demo.provide = function (tab, handle) { handles[tab] = handle; };
  demo.handleOf = function (tab) { return handles[tab]; };

  /* ====================================================================
     Registrierung und Knopf
     ==================================================================== */

  /**
   * spec: { section?, prepare?, steps, check }
   *
   * `section` kommt aus dem Leitsatz, wenn sie nicht angegeben ist; ein
   * Schlüssel ohne Leitsatz ist ein Fehler und kein stiller Sonderfall.
   */
  demo.register = function (key, spec) {
    var e = MV.leitsaetze[key];
    if (!e) throw new Error('Vorführung ohne Leitsatz: ' + key);
    spec.key = key;
    if (!spec.section) spec.section = e.section;
    specs[key] = spec;
    renderHosts(key);
    listen();
  };

  demo.has = function (key) { return !!specs[key]; };

  /** Alle registrierten Schlüssel, in der Reihenfolge der Reiter. */
  demo.keys = function () {
    return Object.keys(specs).sort(function (a, b) {
      var d = TAB_ORDER.indexOf(tabOf(a)) - TAB_ORDER.indexOf(tabOf(b));
      return d !== 0 ? d : (a < b ? -1 : 1);
    });
  };

  /** Link-Fragment für die Arbeit: `nanodraht/A`, `pb-tetron/A`. */
  demo.fragmentOf = function (key) {
    return TAB_HASH[tabOf(key)] + '/' + letterOf(key);
  };

  /** Umgekehrt: aus `pb-tetron` und `A` wird `pbtetron.A`. */
  demo.keyOf = function (hashName, letter) {
    var tab = null;
    Object.keys(TAB_HASH).forEach(function (t) {
      if (TAB_HASH[t] === hashName) tab = t;
    });
    return tab && letter ? tab + '.' + letter.toUpperCase() : null;
  };

  /**
   * Der Knopfplatz in der Leitsatz-Zeile. ui.panel ruft das für jedes Panel
   * mit Zitat auf — auch für Panels ohne Vorführung, damit die Reihenfolge
   * von Aufbau und Registrierung keine Rolle spielt. Solange nichts
   * registriert ist, bleibt der Platz leer und unsichtbar.
   */
  demo.attach = function (key, quoteEl) {
    var wrap = document.createElement('span');
    wrap.className = 'panel__demo';
    wrap.hidden = true;

    var btn = document.createElement('button');
    btn.type = 'button';
    btn.className = 'demo-btn';
    btn.addEventListener('click', function () { onButton(key); });

    var res = document.createElement('span');
    res.className = 'demo-result';
    res.hidden = true;

    wrap.appendChild(btn);
    wrap.appendChild(res);
    quoteEl.appendChild(wrap);

    (hosts[key] || (hosts[key] = [])).push({ wrap: wrap, button: btn, result: res });
    renderHosts(key);
  };

  function eachHost(key, fn) { (hosts[key] || []).forEach(fn); }

  function renderHosts(key) {
    var spec = specs[key];
    var running = !!(current && current.key === key);
    eachHost(key, function (h) {
      if (!spec) { h.wrap.hidden = true; return; }
      h.wrap.hidden = false;
      if (locked[key]) {
        h.button.textContent = 'Vorführung gesperrt — Prüfung fehlgeschlagen';
        h.button.disabled = true;
        h.button.classList.add('demo-btn--locked');
        h.button.setAttribute('aria-label',
          'Vorführung zu ' + secShort(spec.section) +
          ' gesperrt, die Endprüfung ist fehlgeschlagen');
        return;
      }
      h.button.disabled = false;
      h.button.classList.remove('demo-btn--locked');
      h.button.textContent = running
        ? 'Vorführung anhalten'
        : 'Vorführen – ' + secShort(spec.section);
      h.button.setAttribute('aria-pressed', running ? 'true' : 'false');
      h.button.setAttribute('aria-label', running
        ? 'laufende Vorführung anhalten'
        : 'Vorführung zum Leitsatz aus ' + secShort(spec.section) + ' starten');
    });
  }

  /** Die Ergebniszeile unter dem Leitsatz; `r` = null löscht sie. */
  function setResult(key, r) {
    eachHost(key, function (h) {
      if (!r) { h.result.hidden = true; h.result.textContent = ''; return; }
      h.result.hidden = false;
      h.result.textContent = r.text;
      h.result.className = 'demo-result' + (r.ok ? '' : ' demo-result--bad');
    });
  }

  function onButton(key) {
    if (current && current.key === key) { demo.stop(); return; }
    demo.run(key);
  }

  /* ====================================================================
     Ablauf

     Schritte: { set: fn } sofort, { tween: { from, to, ms, apply } } als
     gleichmäßiger Durchlauf, { wait: ms } als Pause. { run: fn(next,
     instant) } bleibt Schritten vorbehalten, die ein vorhandener Player des
     Moduls selbst zu Ende führt — in #gatter Panel B ist das die
     Schrittfolge, die dort schon eingebaut ist.
     ==================================================================== */

  function frameWait(ctx, ms, next) {
    if (ctx.instant || ms <= 0) { next(); return; }
    var t0 = null;
    function frame(ts) {
      if (ctx.cancelled) return;
      if (t0 === null) t0 = ts;
      if (ts - t0 >= ms) { next(); return; }
      ctx.raf = window.requestAnimationFrame(frame);
    }
    ctx.raf = window.requestAnimationFrame(frame);
  }

  function frameTween(ctx, tw, next) {
    if (ctx.instant) { tw.apply(tw.to); next(); return; }
    var t0 = null, span = Math.max(1, tw.ms);
    function frame(ts) {
      if (ctx.cancelled) return;
      if (t0 === null) t0 = ts;
      var u = Math.min(1, (ts - t0) / span);
      tw.apply(tw.from + (tw.to - tw.from) * u);
      if (u >= 1) { next(); return; }
      ctx.raf = window.requestAnimationFrame(frame);
    }
    ctx.raf = window.requestAnimationFrame(frame);
  }

  function exec(list, i, ctx, done) {
    if (ctx.cancelled) return;
    if (i >= list.length) { done(); return; }
    var st = list[i];
    function next() { exec(list, i + 1, ctx, done); }
    if (st.set) { st.set(); next(); return; }
    if (st.tween) { frameTween(ctx, st.tween, next); return; }
    if (st.run) { st.run(next, ctx.instant, ctx); return; }
    frameWait(ctx, st.wait || 0, next);
  }

  /**
   * opts: { instant } erzwingt den Endzustand ohne Durchlauf,
   *       { silent }  lässt die Ergebniszeile unberührt (für verifyAll).
   */
  demo.run = function (key, opts) {
    var spec = specs[key];
    if (!spec || locked[key]) return null;
    opts = opts || {};

    demo.stop();

    /* Ein eigener Durchlauf des Moduls — der μ-Durchlauf in #kitaev, der
       Player in #gatter — würde der Vorführung in die Regler greifen. */
    var h = handles[tabOf(key)];
    if (h) h.stop();

    if (!opts.silent) setResult(key, null);

    var ctx = {
      key: key,
      cancelled: false,
      raf: 0,
      instant: opts.instant || reducedMotion(),
      silent: !!opts.silent
    };
    current = ctx;
    renderHosts(key);

    if (spec.prepare) spec.prepare();
    exec(spec.steps || [], 0, ctx, function () {
      if (ctx.cancelled) return;
      current = null;
      finish(key, spec, ctx);
      renderHosts(key);
    });
    return ctx;
  };

  function finish(key, spec, ctx) {
    var r;
    try {
      r = spec.check();
    } catch (err) {
      r = { ok: false, text: 'Prüfung nicht ausführbar: ' + err.message };
    }
    if (!r.ok) locked[key] = r.text;
    if (!ctx.silent) setResult(key, r);
    ctx.result = r;
  }

  /** Hält die laufende Vorführung an — Knopf, Handeingriff, Reiterwechsel. */
  demo.stop = function () {
    if (!current) return;
    var key = current.key;
    current.cancelled = true;
    if (current.raf) window.cancelAnimationFrame(current.raf);
    if (current.onCancel) current.onCancel();
    current = null;
    renderHosts(key);
  };

  demo.stopAll = function () { demo.stop(); };
  demo.running = function () { return current ? current.key : null; };

  /** Ein laufender Schritt meldet hier, wie er sich abbrechen lässt. */
  demo.onCancel = function (fn) { if (current) current.onCancel = fn; };

  /**
   * Ein Modul meldet hier, was es beim ersten Handeingriff aufzuräumen hat —
   * eine Hervorhebung, die eine Vorführung stehengelassen hat. Das gilt auch,
   * wenn gerade keine Vorführung läuft: die Hervorhebung gehört zur
   * Vorführung, und der Reiter gehört wieder dem Nutzer, sobald er ihn
   * anfasst.
   */
  demo.onHandInput = function (fn) {
    handInput.push(fn);
    listen();
  };

  /* Jeder Handeingriff hält an: ein Klick auf einen Regler, einen Schalter
     oder eine Zeichenfläche, jede Taste. Nur die Knöpfe der Vorführung selbst
     sind ausgenommen — sonst hielte der Knopf sich selbst an. */
  function listen() {
    if (listening) return;
    listening = true;
    ['pointerdown', 'keydown'].forEach(function (type) {
      document.addEventListener(type, function (ev) {
        var t = ev.target;
        if (t && t.closest && t.closest('.panel__demo')) return;
        demo.stop();
        handInput.forEach(function (fn) { fn(); });
      }, true);
    });
  }

  /* ====================================================================
     Reiter öffnen — vom Router gesetzt, damit demo.js den Router nicht kennt
     ==================================================================== */

  /** fn(hash, silent) — `silent` zeigt den Reiter, ohne den Verlauf zu ändern. */
  demo.setTabOpener = function (fn) { showTab = fn; };

  demo.open = function (key, andRun) {
    if (!specs[key]) return false;
    if (showTab) showTab(hashOf(key));
    return demo.showAt(key, andRun);
  };

  /**
   * Zum Panel scrollen und dort vorführen, ohne den Reiter zu wechseln —
   * für die Direktlinks `#<reiter>/<Panelbuchstabe>`, bei denen der Router
   * den Reiter schon geöffnet hat.
   */
  demo.showAt = function (key, andRun) {
    if (!specs[key]) return false;
    var h = (hosts[key] || [])[0];
    if (h && h.wrap.scrollIntoView) h.wrap.scrollIntoView({ block: 'center' });
    if (andRun) demo.run(key);
    return true;
  };

  /* ====================================================================
     Gegenprobe für die Konsole

     Führt jede Vorführung ohne Durchlauf aus, sammelt die Endprüfungen und
     stellt danach den vorherigen Zustand jedes Reiters wieder her.
     ==================================================================== */

  demo.verifyAll = function () {
    var keys = demo.keys();
    var back = window.location.hash;
    var snaps = {};
    Object.keys(handles).forEach(function (t) { snaps[t] = handles[t].snapshot(); });

    var rows = keys.map(function (k) {
      var ok, text;
      try {
        if (showTab) showTab(hashOf(k), true);
        var ctx = demo.run(k, { instant: true, silent: true });
        if (!ctx) { ok = false; text = locked[k] || 'nicht ausführbar'; }
        else { ok = ctx.result.ok; text = ctx.result.text; }
      } catch (err) {
        ok = false;
        text = 'Fehler: ' + err.message;
      }
      return {
        'Schlüssel': k,
        'Abschnitt': specs[k].section,
        'Link': demo.fragmentOf(k),
        'ok': ok ? 'ja' : 'NEIN',
        'Ergebnis': text
      };
    });

    Object.keys(handles).forEach(function (t) {
      if (showTab) showTab('#' + TAB_HASH[t], true);
      handles[t].restore(snaps[t]);
    });
    keys.forEach(function (k) { setResult(k, null); });
    if (showTab) showTab(back || '#kitaev', true);

    var good = rows.filter(function (r) { return r.ok === 'ja'; }).length;
    if (console.table) console.table(rows);
    else rows.forEach(function (r) { console.log(r); });
    console.log(good + ' von ' + rows.length +
                ' Vorführungen erreichen ihren Zielzustand.');
    return { rows: rows, ok: good, total: rows.length, passed: good === rows.length };
  };

}(MV));
