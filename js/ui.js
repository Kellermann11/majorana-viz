/* ------------------------------------------------------------------------
   ui.js — shared UI building blocks: sliders, buttons, read-outs, panels,
   the eight-part info panel and the hash router.

   Interface language is German, code and comments are English.
   ------------------------------------------------------------------------ */

var MV = window.MV || (window.MV = {});

(function (MV) {
  'use strict';

  var ui = MV.ui = {};

  /* --------------------------------------------------------------- helpers */

  function el(tag, attrs, children) {
    var node = document.createElement(tag), k;
    if (attrs) {
      for (k in attrs) {
        if (!Object.prototype.hasOwnProperty.call(attrs, k)) continue;
        if (k === 'class') node.className = attrs[k];
        else if (k === 'html') node.innerHTML = attrs[k];
        else if (k === 'text') node.textContent = attrs[k];
        else node.setAttribute(k, attrs[k]);
      }
    }
    (children || []).forEach(function (c) {
      if (c) node.appendChild(typeof c === 'string' ? document.createTextNode(c) : c);
    });
    return node;
  }
  ui.el = el;

  /* ---------------------------------------------------------------- slider */

  /**
   * spec: { label, symbol, hint, min, max, step, value, digits, prefix,
   *         scale, format, onInput }
   *
   * `scale: 'log'` legt den Regler über den Zehnerlogarithmus des Werts.
   * `min`, `max` und `value` stehen trotzdem als Werte in der Spezifikation,
   * und `step` zählt dann in Dekaden — ein Reiter hantiert nie mit
   * Exponenten. `format(v)` ersetzt die Zahlenausgabe, etwa um die Einheit
   * mitwandern zu lassen.
   *
   * Angezeigt wird immer der zuletzt gesetzte Wert, nicht die Stellung des
   * Reglers. Beide stimmen normalerweise überein; sie fallen auseinander,
   * wenn ein Wert außerhalb des Bereichs dieses Reglers gesetzt wird — dann
   * steht der Griff am Anschlag, die Zahl bleibt aber die richtige.
   * `inRange(v)` sagt, ob das der Fall ist.
   *
   * Returns { value(), set(v, silent), snap(v), inRange(v), setHint(html), root }.
   */
  ui.slider = function (mount, spec) {
    var digits = spec.digits === undefined ? 2 : spec.digits;
    var log = spec.scale === 'log';

    function toRaw(v) { return log ? Math.log(v) / Math.LN10 : v; }
    function toVal(r) { return log ? Math.pow(10, r) : r; }

    var rawMin = toRaw(spec.min), rawMax = toRaw(spec.max);
    var current = spec.value;

    var valueEl = el('span', { class: 'slider__value' });
    var input = el('input', {
      type: 'range',
      min: String(rawMin),
      max: String(rawMax),
      step: String(spec.step),
      value: String(toRaw(spec.value)),
      'aria-label': spec.label
    });
    var head = el('div', { class: 'slider__head' }, [
      el('span', { class: 'slider__name', html: spec.label }),
      valueEl
    ]);
    var root = el('div', { class: 'slider' }, [head, input]);
    var hintEl = el('div', { class: 'slider__hint', html: spec.hint || '' });
    if (!spec.hint) hintEl.hidden = true;
    root.appendChild(hintEl);
    mount.appendChild(root);

    /* Der Wert, den dieser Regler wirklich annehmen kann: in den Bereich
       geklemmt und auf sein Raster gerundet. Eine Vorführung stellt damit
       genau das ein, was auch eine Hand einstellen könnte (demo.js). */
    function snap(v) {
      var r = Math.min(rawMax, Math.max(rawMin, toRaw(v)));
      var n = rawMin + Math.round((r - rawMin) / spec.step) * spec.step;
      n = Math.round(n * 1e9) / 1e9;
      return toVal(Math.min(rawMax, Math.max(rawMin, n)));
    }

    var api = {
      root: root,
      input: input,
      value: function () { return current; },
      set: function (v, silent) {
        current = v;
        /* Der Griff wird geklemmt, die Zahl nicht. */
        input.value = String(Math.min(rawMax, Math.max(rawMin, toRaw(v))));
        render();
        if (!silent && spec.onInput) spec.onInput(current, true);
      },
      snap: snap,
      inRange: function (v) {
        var w = v === undefined ? current : v;
        return w >= spec.min - 1e-12 && w <= spec.max + 1e-12;
      },
      setHint: function (html) {
        hintEl.innerHTML = html || '';
        hintEl.hidden = !html;
      }
    };

    function render() {
      valueEl.textContent = spec.format
        ? spec.format(current)
        : (spec.prefix || '') + MV.fmt(current, digits);
    }

    input.addEventListener('input', function () {
      current = toVal(parseFloat(input.value));
      render();
      if (spec.onInput) spec.onInput(current, false);
    });

    render();
    return api;
  };

  /* --------------------------------------------------------------- buttons */

  /** spec: { label, buttons: [{ text, onClick, pressed }] } */
  ui.buttonRow = function (mount, spec) {
    var row = el('div', { class: 'buttonrow' });
    if (spec.label) row.appendChild(el('span', { class: 'buttonrow__label', html: spec.label }));
    var made = spec.buttons.map(function (b) {
      var btn = el('button', { type: 'button', html: b.text });
      if (b.pressed !== undefined) btn.setAttribute('aria-pressed', b.pressed ? 'true' : 'false');
      btn.addEventListener('click', function () { b.onClick(btn); });
      row.appendChild(btn);
      return btn;
    });
    mount.appendChild(row);
    return { root: row, buttons: made };
  };

  /** Mutually exclusive button group; onChange receives the index. */
  ui.segmented = function (mount, spec) {
    var current = spec.value || 0;
    var group = ui.buttonRow(mount, {
      label: spec.label,
      buttons: spec.options.map(function (text, i) {
        return {
          text: text,
          pressed: i === current,
          onClick: function () { api.set(i); }
        };
      })
    });
    var api = {
      root: group.root,
      value: function () { return current; },
      set: function (i) {
        current = i;
        group.buttons.forEach(function (b, j) {
          b.setAttribute('aria-pressed', j === i ? 'true' : 'false');
        });
        if (spec.onChange) spec.onChange(i);
      }
    };
    return api;
  };

  /* --------------------------------------------------------------- readout */

  /**
   * items: [{ key, label }]. Returns { set(key, text), status(kind, text) }.
   * `kind` is 'topological' | 'trivial' | 'critical'; the text is always
   * printed as well, so colour is never the only carrier of the state.
   */
  ui.readout = function (mount, spec) {
    var cells = {};
    var grid = el('div', { class: 'readout__grid' });
    spec.items.forEach(function (it) {
      var val = el('div', { class: 'readout__val', text: '—' });
      grid.appendChild(el('div', { class: 'readout__item' }, [
        el('div', { class: 'readout__key', html: it.label }),
        val
      ]));
      cells[it.key] = val;
    });
    var statusEl = el('div', { class: 'status status--trivial', text: '—' });
    var root = el('div', { class: 'readout' }, [grid]);
    if (spec.status !== false) root.appendChild(statusEl);
    if (spec.note) root.appendChild(el('p', { class: 'note', html: spec.note }));
    mount.appendChild(root);

    return {
      root: root,
      set: function (key, text) {
        if (cells[key]) cells[key].textContent = text;
      },
      status: function (kind, text) {
        statusEl.className = 'status status--' + kind;
        statusEl.textContent = text;
      }
    };
  };

  /* ------------------------------------------------------------ info panel */

  var INFO_SECTIONS = [
    ['what',       'Was dargestellt ist'],
    ['model',      'Zugrundeliegendes Modell'],
    ['computed',   'Berechnete Größen'],
    ['params',     'Aktuelle Parameter'],
    ['numerics',   'Numerik'],
    ['convention', 'Konventionshinweis'],
    ['reference',  'Bezug zur Arbeit'],
    ['reading',    'Worauf zu achten ist']
  ];

  /**
   * Eight-part figure caption, collapsed by default.
   * spec holds one HTML string per section key; `params` is a function that
   * returns HTML and is re-evaluated whenever the panel is open.
   */
  ui.info = function (mount, spec) {
    var body = el('div', { class: 'info__body' });
    var paramsHost = null;

    INFO_SECTIONS.forEach(function (s) {
      var key = s[0], heading = s[1];
      var content = spec[key];
      if (content === undefined) return;
      body.appendChild(el('h4', { text: heading }));
      if (key === 'params') {
        paramsHost = el('div');
        body.appendChild(paramsHost);
      } else {
        body.appendChild(el('div', { html: content }));
      }
    });

    var details = el('details', { class: 'info' }, [
      el('summary', { text: 'Details zur Abbildung' }),
      body
    ]);
    mount.appendChild(details);

    var api = {
      root: details,
      update: function (force) {
        if (!paramsHost || (!details.open && !force)) return;
        paramsHost.innerHTML = spec.params();
      }
    };

    details.addEventListener('toggle', function () { api.update(); });
    return api;
  };

  /* ----------------------------------------------------------------- panel */

  /**
   * spec: { title, lead, wide, aspect, margin, controls (bool), info }
   * Returns { root, canvas, plot, controls, info }.
   */
  ui.panel = function (mount, spec) {
    var root = el('section', { class: 'panel' + (spec.wide ? ' panel--wide' : '') });
    root.appendChild(el('h3', { class: 'panel__title', html: spec.title }));

    /* CLAUDE.md, Abschnitt 11: der Satz aus der Arbeit, den dieses Panel
       sichtbar macht — über dem Bild, ohne Klick, nie im Info-Panel
       versteckt. Ohne Zitat steht dort sichtbar "Leitsatz offen". */
    if (spec.quote) {
      var quoteEl = el('p', { class: 'panel__quote', html: MV.leitsatzHtml(spec.quote) });
      root.appendChild(quoteEl);
      /* Der Knopfplatz der Vorführung sitzt in derselben Leiste wie der Satz,
         den er vorführt. Er bleibt leer, solange für diesen Leitsatz keine
         Vorführung registriert ist (demo.js). */
      if (MV.demo) MV.demo.attach(spec.quote, quoteEl);
    }

    if (spec.lead) root.appendChild(el('p', { class: 'panel__lead', html: spec.lead }));

    /* A panel normally holds a plot. `canvas: false` gives an empty body
       instead, for content that is a table rather than a drawing. */
    var canvas = null, body = null;
    if (spec.canvas === false) {
      body = el('div', { class: 'panel__body' });
      root.appendChild(body);
    } else {
      canvas = el('canvas', { class: 'panel__canvas' });
      root.appendChild(canvas);
    }

    /* One-line reading aid for the current parameters, and the legend. Both
       sit below the plot rather than inside it: which corner of a panel is
       free depends on the parameters, so anything painted into the frame
       eventually covers data. */
    var status = el('p', { class: 'panel__status' });
    root.appendChild(status);

    var legend = el('div', { class: 'legend-inline' });
    root.appendChild(legend);

    var controls = null;
    if (spec.controls) {
      controls = el('div', { class: 'panel__controls' });
      root.appendChild(controls);
    }

    mount.appendChild(root);

    var plot = null;
    if (canvas) {
      plot = new MV.Plot(canvas, { aspect: spec.aspect, margin: spec.margin });
      plot.setLegendHost(legend);
    }
    var info = spec.info ? ui.info(root, spec.info) : null;

    return {
      root: root, canvas: canvas, plot: plot, body: body,
      legend: legend, controls: controls, info: info,
      /** Reading aid under the plot; `color` may be a palette value. */
      setStatus: function (text, color) {
        status.textContent = text || '';
        status.style.color = color || MV.palette.muted;
      }
    };
  };

  /* ------------------------------------------------------------- layout */

  /**
   * Legt die Bedienleiste eines Reiters in eine schmale Spalte links neben
   * die Abbildungen. Die Spalte bleibt beim Scrollen stehen, sodass die
   * Grafik sichtbar ist, während ein Regler gezogen wird; auf schmalen
   * Fenstern klappt sie über die Abbildungen.
   *
   * Wird einmal je Reiter aufgerufen, nachdem dessen Modul seinen Inhalt
   * gebaut hat. Die sechs Reiter-Module hängen ihre Teile deshalb weiterhin
   * einfach der Reihe nach an den Reiter und wissen von der Spalte nichts.
   */
  ui.twoColumn = function (tabEl) {
    var layout = el('div', { class: 'tab-layout' });
    var aside = el('div', { class: 'tab-layout__controls' });
    var main = el('div', { class: 'tab-layout__main' });
    layout.appendChild(aside);
    layout.appendChild(main);

    /* Erst die Liste einfrieren: das Verschieben verändert childNodes. */
    var kids = Array.prototype.slice.call(tabEl.childNodes);
    kids.forEach(function (node) {
      if (node.nodeType === 1 && node.classList.contains('controls')) aside.appendChild(node);
      else main.appendChild(node);
    });

    tabEl.appendChild(layout);
    return layout;
  };

  /* ---------------------------------------------------------------- router */

  /**
   * routes: [{ hash, tabEl, linkEl, chapter, convention, onEnter }]
   * Direct calls of #kitaev / #nanodraht open the right tab; the browser back
   * button works because we only listen to hashchange.
   */
  ui.router = function (routes, contextEl) {
    var previous = null;

    /**
     * `explicit` zeigt einen Reiter, ohne die Adresszeile zu lesen — das
     * braucht die Gegenprobe MV.demo.verifyAll(), die alle Reiter der Reihe
     * nach durchgeht und dabei weder Verlauf noch Zurück-Button anfassen soll.
     */
    function apply(explicit) {
      var raw = explicit || window.location.hash;

      /* Direktlinks aus der Arbeit tragen hinter dem Reiter den Buchstaben
         des Panels: #nanodraht/A öffnet den Reiter, springt zum Panel und
         startet dessen Vorführung. Ein unbekannter Buchstabe öffnet nur den
         Reiter — erfunden wird nichts. */
      var cut = raw.indexOf('/');
      var hash = cut >= 0 ? raw.slice(0, cut) : raw;
      var letter = cut >= 0 ? raw.slice(cut + 1) : '';

      var route = null, i;
      for (i = 0; i < routes.length; i++) {
        if (routes[i].hash === hash) { route = routes[i]; break; }
      }
      /* Ein unbekannter Reiter führt zum ersten. Ein Panelbuchstabe, der an
         einem unbekannten Reiter hing, gilt dann nicht mehr — sonst startete
         #bewertung/A eine Vorführung, nach der niemand gefragt hat. */
      if (!route) { route = routes[0]; letter = ''; }

      if (previous !== route) {
        /* Eine laufende Vorführung gehört zum Reiter, in dem sie läuft. */
        if (MV.demo) MV.demo.stopAll();
        if (previous && previous.onLeave) previous.onLeave();
      }
      previous = route;

      routes.forEach(function (r) {
        var active = r === route;
        r.tabEl.hidden = !active;
        r.linkEl.setAttribute('aria-selected', active ? 'true' : 'false');
      });

      contextEl.innerHTML =
        '<b>' + route.chapter + '</b><span class="sep">|</span>' + route.convention;

      if (route.onEnter) route.onEnter();

      if (letter && MV.demo) {
        var key = MV.demo.keyOf(route.hash.slice(1), letter);
        if (key) MV.demo.showAt(key, true);
      }
    }

    window.addEventListener('hashchange', function () { apply(); });
    apply();
    return { apply: apply, show: function (hash) { apply(hash); } };
  };

}(MV));
