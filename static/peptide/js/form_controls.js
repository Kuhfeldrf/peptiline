/*
 * form_controls.js — shared UI helpers for the plotting dashboards.
 *
 * One normalized way to enable/disable an option based on other selections,
 * used by both the Heatmap and Data Analysis dashboards so every gated control
 * looks and behaves the same (no more per-call inline `style.opacity` tweaks).
 */

/**
 * Enable or disable a form control and dim its "row" consistently.
 *
 * Sets `el.disabled` and toggles the `.control-disabled` class (opacity + a
 * not-allowed cursor, defined once in base.html) on the control's visual row.
 *
 * @param {Element|null} el         the input/select being gated (no-op if null)
 * @param {boolean}      enabled    true → interactive; false → disabled + greyed
 * @param {Object}       [opts]
 * @param {Element}      [opts.row] element to dim; defaults to the control's
 *                                  wrapping <label>, or the control itself
 *                                  (pass the control for bare <select>s).
 * @param {boolean}      [opts.clearOnDisable] when disabling, also uncheck a
 *                                  checkbox/radio (e.g. a mode that no longer applies).
 */
function setControlEnabled(el, enabled, opts) {
    if (!el) return;
    opts = opts || {};
    el.disabled = !enabled;
    const row = opts.row || el.closest('label') || el;
    if (row) row.classList.toggle('control-disabled', !enabled);
    if (!enabled && opts.clearOnDisable && (el.type === 'checkbox' || el.type === 'radio')) {
        el.checked = false;
    }
}
