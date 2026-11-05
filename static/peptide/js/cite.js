/* Copy buttons in the nav "Cite" menu.
   Each button names, via data-cite-source, the id of the off-screen element
   that holds the exact text to copy (formatted citation or BibTeX record). */
(function () {
    'use strict';

    function copyFrom(source) {
        var text = source.value !== undefined ? source.value : source.textContent;

        if (navigator.clipboard && navigator.clipboard.writeText) {
            return navigator.clipboard.writeText(text);
        }

        // Fallback for older browsers, and for pages served over plain HTTP
        // where the async clipboard API is unavailable.
        return new Promise(function (resolve, reject) {
            try {
                source.select();
                if (source.setSelectionRange) {
                    source.setSelectionRange(0, text.length);
                }
                if (document.execCommand('copy')) {
                    resolve();
                } else {
                    reject(new Error('execCommand("copy") returned false'));
                }
            } catch (err) {
                reject(err);
            }
        });
    }

    function flash(button, message) {
        if (!button.getAttribute('data-cite-label')) {
            button.setAttribute('data-cite-label', button.textContent);
        }
        button.textContent = message;
        button.classList.add('pl-cite-copied');
        window.clearTimeout(button.plCiteTimer);
        button.plCiteTimer = window.setTimeout(function () {
            button.textContent = button.getAttribute('data-cite-label');
            button.classList.remove('pl-cite-copied');
        }, 2000);
    }

    function init() {
        var buttons = document.querySelectorAll('.pl-cite-copy');

        Array.prototype.forEach.call(buttons, function (button) {
            button.addEventListener('click', function (event) {
                event.preventDefault();
                // Keep the click from reaching Bootstrap's document handler,
                // which would close the dropdown out from under the button.
                event.stopPropagation();

                var source = document.getElementById(button.getAttribute('data-cite-source'));
                if (!source) {
                    return;
                }

                copyFrom(source).then(
                    function () { flash(button, 'Copied'); },
                    function () { flash(button, 'Copy failed'); }
                );
            });
        });
    }

    if (document.readyState === 'loading') {
        document.addEventListener('DOMContentLoaded', init);
    } else {
        init();
    }
}());
