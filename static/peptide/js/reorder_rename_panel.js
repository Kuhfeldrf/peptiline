/*
 * reorder_rename_panel.js — shared "arrange & rename" companion for a
 * <select multiple> on the plotting dashboards.
 *
 * The <select> stays the picker. Whatever is selected shows up in a side panel
 * as a drag-to-reorder, click-to-rename list. That ordered/renamed list is what
 * the figure uses:
 *
 *   panel.getOrder()          -> ["Late Torpor", "Early Torpor", ...]   (real values only)
 *   panel.getLabelOverrides() -> { "P02666": "β-casein" }               (renamed only)
 *
 * A rename is display-only: getOrder() still returns the original value, so the
 * server keeps colour / lookup / stats keyed by the real name and only the
 * printed label changes.
 *
 * No dependencies. Instantiate once per <select>, after its <option>s exist.
 */
(function (global) {
'use strict';

function ReorderRenamePanel(rootEl, opts) {
    opts = opts || {};
    this.root = rootEl;
    this.select = opts.selectEl;
    this.isRealValue = opts.isRealValue || function () { return true; };
    this.onChange = opts.onChange || function () {};
    this.noun = opts.noun || 'item';
    // Optional: given a value, return the label it should show *before* any
    // rename — used when an outside control (e.g. "Strip protein name") changes
    // what the canonical name is. Return null/undefined to fall back to the
    // <option> text. Call refresh() after the outside control changes.
    this.originalLabelFor = opts.originalLabelFor || null;

    this.order = [];      // array of option values, panel-controlled
    this.labels = {};     // { value: customLabel } — only when != original

    this._dragValue = null;
    this._buildScaffold();
    this._wireSelect();
    this.sync();
}

ReorderRenamePanel.prototype._buildScaffold = function () {
    this.root.classList.add('rr-panel');
    this.root.innerHTML =
        '<div class="rr-head">' +
            '<span class="rr-head-label">On the figure &middot; drag to reorder, click a name to rename</span>' +
            '<span class="rr-head-links">' +
                '<a href="#" data-rr="clear">Clear</a>' +
                '<a href="#" data-rr="reset-order">Reset order</a>' +
                '<a href="#" data-rr="reset-names">Reset names</a>' +
            '</span>' +
        '</div>' +
        '<ul class="rr-rows" data-rr="rows"></ul>' +
        '<div class="rr-empty" data-rr="empty"></div>';

    this.rowsEl = this.root.querySelector('[data-rr="rows"]');
    this.emptyEl = this.root.querySelector('[data-rr="empty"]');
    this.emptyEl.textContent = 'Select one or more ' + this.noun + 's to arrange them.';

    var self = this;
    this.root.querySelector('[data-rr="clear"]').addEventListener('click', function (e) {
        e.preventDefault();
        Array.from(self.select.options).forEach(function (o) {
            if (self.isRealValue(o.value)) o.selected = false;
        });
        self.select.dispatchEvent(new Event('change'));
        self.sync();
    });
    this.root.querySelector('[data-rr="reset-order"]').addEventListener('click', function (e) {
        e.preventDefault();
        var sel = new Set(self._selectedRealValues());
        self.order = self._optionValues().filter(function (v) { return sel.has(v); });
        self._render();
    });
    this.root.querySelector('[data-rr="reset-names"]').addEventListener('click', function (e) {
        e.preventDefault();
        self.labels = {};
        self._render();
    });

    this._wireDrag();
};

ReorderRenamePanel.prototype._wireSelect = function () {
    var self = this;
    this.select.addEventListener('change', function () { self.sync(); });
};

ReorderRenamePanel.prototype._optionValues = function () {
    return Array.from(this.select.options).map(function (o) { return o.value; });
};

ReorderRenamePanel.prototype._selectedRealValues = function () {
    var self = this;
    return Array.from(this.select.selectedOptions)
        .map(function (o) { return o.value; })
        .filter(function (v) { return self.isRealValue(v); });
};

ReorderRenamePanel.prototype._optionByValue = function (v) {
    return Array.from(this.select.options).find(function (o) { return o.value === v; });
};

ReorderRenamePanel.prototype._originalLabel = function (v) {
    if (this.originalLabelFor) {
        var forced = this.originalLabelFor(v);
        if (forced != null && forced !== '') return forced;
    }
    var opt = this._optionByValue(v);
    // The <option> text sometimes carries an accession suffix ("Name (P02666)")
    // or helper punctuation; for renaming we want what the user recognises, so
    // fall back to the option's text but strip a trailing " (VALUE)".
    if (!opt) return v;
    var txt = opt.textContent.trim();
    var suffix = ' (' + v + ')';
    if (txt.endsWith(suffix)) txt = txt.slice(0, -suffix.length).trim();
    return txt || v;
};

ReorderRenamePanel.prototype.displayName = function (v) {
    return Object.prototype.hasOwnProperty.call(this.labels, v)
        ? this.labels[v]
        : this._originalLabel(v);
};

// Reconcile panel state with the current <select> state, then re-render.
ReorderRenamePanel.prototype.sync = function () {
    var selected = this._selectedRealValues();
    var selSet = new Set(selected);

    this.order = this.order.filter(function (v) { return selSet.has(v); });
    selected.forEach(function (v) {
        if (this.order.indexOf(v) === -1) this.order.push(v);
    }, this);

    Object.keys(this.labels).forEach(function (v) {
        if (!selSet.has(v)) delete this.labels[v];
    }, this);

    this._render();
};

ReorderRenamePanel.prototype.reset = function () {
    this.order = [];
    this.labels = {};
    this.sync();
};

ReorderRenamePanel.prototype.refresh = function () { this._render(); };

ReorderRenamePanel.prototype.getOrder = function () {
    return this.order.slice();
};

ReorderRenamePanel.prototype.getLabelOverrides = function () {
    var out = {};
    Object.keys(this.labels).forEach(function (v) {
        var name = String(this.labels[v]).trim();
        if (name && name !== this._originalLabel(v)) out[v] = name;
    }, this);
    return out;
};

ReorderRenamePanel.prototype._render = function () {
    var n = this.order.length;
    this.emptyEl.style.display = n ? 'none' : 'block';
    this.rowsEl.style.display = n ? 'flex' : 'none';
    this.root.querySelector('[data-rr="clear"]').classList.toggle('rr-link-off', n === 0);
    this.root.querySelector('[data-rr="reset-order"]').classList.toggle('rr-link-off', n < 2);
    this.root.querySelector('[data-rr="reset-names"]')
        .classList.toggle('rr-link-off', Object.keys(this.labels).length === 0);

    var self = this;
    this.rowsEl.innerHTML = '';
    this.order.forEach(function (val, i) {
        var orig = self._originalLabel(val);
        var renamed = Object.prototype.hasOwnProperty.call(self.labels, val)
            && self.labels[val] !== orig;

        var li = document.createElement('li');
        li.className = 'rr-row';
        li.draggable = true;
        li.dataset.value = val;

        var idx = document.createElement('span');
        idx.className = 'rr-idx';
        idx.textContent = i + 1;

        var grip = document.createElement('i');
        grip.className = 'fas fa-grip-vertical rr-grip';

        var names = document.createElement('span');
        names.className = 'rr-names';
        var disp = document.createElement('span');
        disp.className = 'rr-disp';
        disp.textContent = self.displayName(val);
        disp.title = 'Click to rename';
        names.appendChild(disp);
        if (renamed) {
            var was = document.createElement('span');
            was.className = 'rr-was';
            was.textContent = 'was: ' + orig;
            names.appendChild(was);
        }

        var drop = document.createElement('i');
        drop.className = 'fas fa-xmark rr-drop';
        drop.title = 'Remove from selection';

        li.appendChild(grip);
        li.appendChild(idx);
        li.appendChild(names);
        if (renamed) {
            var badge = document.createElement('span');
            badge.className = 'rr-badge';
            badge.textContent = 'renamed';
            li.appendChild(badge);
        }
        li.appendChild(drop);

        disp.addEventListener('click', function () { self._startRename(li, val, orig); });
        drop.addEventListener('click', function () {
            var opt = self._optionByValue(val);
            if (opt) opt.selected = false;
            self.select.dispatchEvent(new Event('change'));
            self.sync();
        });

        self.rowsEl.appendChild(li);
    });

    this.onChange();
};

ReorderRenamePanel.prototype._startRename = function (li, val, orig) {
    var self = this;
    var names = li.querySelector('.rr-names');
    var current = this.displayName(val);
    names.innerHTML = '';
    var input = document.createElement('input');
    input.className = 'rr-rename';
    input.value = current;
    names.appendChild(input);
    input.focus();
    input.select();

    var commit = function () {
        var next = input.value.trim();
        if (!next || next === orig) delete self.labels[val];
        else self.labels[val] = next;
        self._render();
    };
    input.addEventListener('blur', commit);
    input.addEventListener('keydown', function (e) {
        if (e.key === 'Enter') { e.preventDefault(); input.blur(); }
        else if (e.key === 'Escape') { input.value = current; input.blur(); }
    });
};

ReorderRenamePanel.prototype._wireDrag = function () {
    var self = this;
    this.rowsEl.addEventListener('dragstart', function (e) {
        var li = e.target.closest('.rr-row');
        if (!li) return;
        self._dragValue = li.dataset.value;
        li.classList.add('rr-dragging');
    });
    this.rowsEl.addEventListener('dragend', function () {
        self.rowsEl.querySelectorAll('.rr-dragging, .rr-drop-target')
            .forEach(function (el) { el.classList.remove('rr-dragging', 'rr-drop-target'); });
    });
    this.rowsEl.addEventListener('dragover', function (e) {
        e.preventDefault();
        var li = e.target.closest('.rr-row');
        self.rowsEl.querySelectorAll('.rr-drop-target')
            .forEach(function (el) { el.classList.remove('rr-drop-target'); });
        if (li && li.dataset.value !== self._dragValue) li.classList.add('rr-drop-target');
    });
    this.rowsEl.addEventListener('drop', function (e) {
        e.preventDefault();
        var li = e.target.closest('.rr-row');
        var dragValue = self._dragValue;
        if (!li || li.dataset.value === dragValue) return;
        self.order = self.order.filter(function (v) { return v !== dragValue; });
        var targetIdx = self.order.indexOf(li.dataset.value);
        var rect = li.getBoundingClientRect();
        var after = e.clientY > rect.top + rect.height / 2;
        self.order.splice(targetIdx + (after ? 1 : 0), 0, dragValue);
        self._render();
    });
};

// Wire a plain text input to filter a <select>'s options by substring
// (matches option text or value). Selected options are never hidden.
ReorderRenamePanel.filterSelect = function (inputEl, selectEl) {
    inputEl.addEventListener('input', function () {
        var q = inputEl.value.trim().toLowerCase();
        Array.from(selectEl.options).forEach(function (o) {
            if (!q || o.selected) { o.hidden = false; return; }
            o.hidden = o.textContent.toLowerCase().indexOf(q) === -1
                    && o.value.toLowerCase().indexOf(q) === -1;
        });
    });
};

global.ReorderRenamePanel = ReorderRenamePanel;

})(window);
