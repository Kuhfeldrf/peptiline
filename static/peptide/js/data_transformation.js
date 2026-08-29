/**
 * Data Transformation Wizard - AJAX navigation and progress polling.
 */
(function() {
    'use strict';

    var currentStep = 1;
    var definedGroups = [];
    var availableColumns = [];
    var selectedColumns = [];  // columns selected for the current group being built
    var csrfToken = document.querySelector('[name=csrfmiddlewaretoken]').value;

    // -----------------------------------------------------------------------
    // Drag & drop file inputs
    // -----------------------------------------------------------------------
    document.querySelectorAll('.dt-dropzone').forEach(function(zone) {
        var input = zone.querySelector('input[type="file"]');
        var dropText = zone.querySelector('.drop-text');
        var fileName = zone.querySelector('.file-name');
        var mergeFileList = zone.querySelector('.merge-file-list');
        var isMulti = input && input.hasAttribute('multiple');

        function showFile(name) {
            dropText.classList.add('hidden');
            fileName.textContent = name;
            fileName.classList.remove('hidden');
            if (mergeFileList) mergeFileList.classList.add('hidden');
        }

        function showMultipleFiles(files) {
            dropText.classList.add('hidden');
            fileName.classList.add('hidden');
            if (!mergeFileList) return;
            mergeFileList.innerHTML = '';
            for (var i = 0; i < files.length; i++) {
                var item = document.createElement('div');
                item.className = 'merge-file-item';
                item.textContent = files[i].name;
                mergeFileList.appendChild(item);
            }
            mergeFileList.classList.remove('hidden');
        }

        input.addEventListener('change', function() {
            if (!this.files.length) return;
            if (isMulti && this.files.length > 1) {
                showMultipleFiles(this.files);
            } else {
                showFile(this.files[0].name);
            }
        });

        zone.addEventListener('dragover', function(e) {
            e.preventDefault();
            zone.classList.add('dragover');
        });
        zone.addEventListener('dragleave', function() {
            zone.classList.remove('dragover');
        });
        zone.addEventListener('drop', function(e) {
            e.preventDefault();
            zone.classList.remove('dragover');
            if (e.dataTransfer.files.length) {
                input.files = e.dataTransfer.files;
                if (isMulti && e.dataTransfer.files.length > 1) {
                    showMultipleFiles(e.dataTransfer.files);
                } else {
                    showFile(e.dataTransfer.files[0].name);
                }
                input.dispatchEvent(new Event('change'));
            }
        });
    });

    // -----------------------------------------------------------------------
    // "Load example" links — fetch static example file and inject into input
    // -----------------------------------------------------------------------
    document.querySelectorAll('.dt-load-example').forEach(function(link) {
        link.addEventListener('click', function(e) {
            e.preventDefault();
            var url = link.getAttribute('data-url');
            var inputId = link.getAttribute('data-input-id');
            var input = document.getElementById(inputId);
            var fileName = url.split('/').pop().replace(/\.[0-9a-f]{8,12}(\.[^.]+)$/, '$1');
            var origHtml = link.innerHTML;

            link.innerHTML = '<i class="fas fa-spinner fa-spin" style="font-size:0.75rem;"></i> Loading...';
            link.classList.add('loading');

            fetch(url)
                .then(function(r) {
                    if (!r.ok) throw new Error('HTTP ' + r.status);
                    return r.blob();
                })
                .then(function(blob) {
                    var file = new File([blob], fileName, {type: blob.type || 'application/octet-stream'});
                    var dt = new DataTransfer();
                    dt.items.add(file);
                    input.files = dt.files;
                    input.dispatchEvent(new Event('change'));
                    link.innerHTML = '<i class="fas fa-check" style="font-size:0.75rem;"></i> Loaded';
                    link.classList.remove('loading');
                    link.classList.add('loaded');
                    // Reset after 3 s so it can be clicked again
                    setTimeout(function() {
                        link.innerHTML = origHtml;
                        link.classList.remove('loaded');
                    }, 3000);
                })
                .catch(function() {
                    link.innerHTML = '<i class="fas fa-times" style="font-size:0.75rem;"></i> Failed';
                    link.classList.remove('loading');
                    link.style.color = '#dc3545';
                    setTimeout(function() {
                        link.innerHTML = origHtml;
                        link.style.color = '';
                    }, 3000);
                });
        });
    });

    // Utility: AJAX helper
    function ajax(method, url, data, callback, errorCallback) {
        var xhr = new XMLHttpRequest();
        xhr.open(method, url);
        xhr.setRequestHeader('X-CSRFToken', csrfToken);
        xhr.setRequestHeader('X-Requested-With', 'XMLHttpRequest');

        xhr.onload = function() {
            var resp;
            try {
                resp = JSON.parse(xhr.responseText);
            } catch(e) {
                // Server returned non-JSON — surface the raw text so user can diagnose
                var raw = xhr.responseText ? xhr.responseText.replace(/<[^>]+>/g, ' ').replace(/\s+/g, ' ').trim() : '';
                var preview = raw.substring(0, 300);
                resp = {error: 'Server returned an unexpected response (HTTP ' + xhr.status + '). ' +
                    (preview ? 'Details: ' + preview : 'Check server logs.')};
            }
            if (xhr.status >= 200 && xhr.status < 300) {
                callback(resp);
            } else {
                (errorCallback || showError)(resp.error || 'Request failed', resp);
            }
        };
        xhr.onerror = function() {
            (errorCallback || showError)('Network error');
        };

        if (data instanceof FormData) {
            xhr.send(data);
        } else if (data) {
            xhr.setRequestHeader('Content-Type', 'application/json');
            xhr.send(JSON.stringify(data));
        } else {
            xhr.send();
        }
    }

    function showError(msg) {
        var el = document.getElementById('dt-error');
        el.textContent = msg;
        el.classList.remove('hidden');
        // Long, actionable messages (e.g. lists of invalid peptides) stay visible
        // until the user acts; short errors auto-hide so they don't linger.
        if (!msg || msg.length < 200) {
            setTimeout(function() { el.classList.add('hidden'); }, 10000);
        }
    }

    // Step 2 errors are shown inline within the relevant sub-section so they
    // stay visible no matter how far the user has scrolled. `section` is one of
    // 'rename' | 'techrep' | 'groups' | 'actions'; omitted falls back to the
    // section-wide banner near the top of the card.
    var STEP2_ERROR_IDS = {
        rename: 'dt-error-rename',
        techrep: 'dt-error-techrep',
        groups: 'dt-error-groups',
        actions: 'dt-error-step2-actions'
    };

    function clearStep2Errors() {
        ['dt-error-step2'].concat(Object.keys(STEP2_ERROR_IDS).map(function(k) {
            return STEP2_ERROR_IDS[k];
        })).forEach(function(id) {
            var el = document.getElementById(id);
            if (el) el.classList.add('hidden');
        });
    }

    function showStep2Error(msg, section) {
        clearStep2Errors();
        var el = document.getElementById(STEP2_ERROR_IDS[section] || 'dt-error-step2');
        if (!el) el = document.getElementById('dt-error-step2');
        el.textContent = msg;
        el.classList.remove('hidden');
        el.scrollIntoView({behavior: 'smooth', block: 'nearest'});
        setTimeout(function() { el.classList.add('hidden'); }, 10000);
    }

    function goToStep(step) {
        currentStep = step;
        var indicators = document.querySelectorAll('.dt-step-indicator');
        var contents = document.querySelectorAll('.dt-step-content');
        indicators.forEach(function(ind) {
            var s = parseInt(ind.getAttribute('data-step'));
            ind.classList.remove('active', 'completed');
            if (s < step) ind.classList.add('completed');
            else if (s === step) ind.classList.add('active');
        });
        contents.forEach(function(c) {
            var s = parseInt(c.getAttribute('data-step'));
            c.classList.toggle('active', s === step);
        });
    }

    // -----------------------------------------------------------------------
    // Step 1: Upload & BLAST
    // -----------------------------------------------------------------------

    // The Functional Annotation actions require an uploaded peptidomic dataset.
    function enableFunctionalActions() {
        document.getElementById('start-blast-btn').disabled = false;
        document.getElementById('skip-blast-btn').disabled = false;
        var note = document.getElementById('functional-gate-note');
        if (note) note.classList.add('hidden');
    }

    document.getElementById('upload-form').addEventListener('submit', function(e) {
        e.preventDefault();
        var pepInput = document.getElementById('peptidomic_file');
        var pepFiles = pepInput ? pepInput.files : [];

        if (!pepFiles.length) {
            showError('Please select a peptidomic data file (add multiple files to merge them into one dataset).');
            return;
        }

        var btn = document.getElementById('upload-btn');
        btn.disabled = true;
        btn.innerHTML = '<span class="dt-spinner"></span> Uploading...';

        var formData = new FormData(this);
        // Two or more files → merge mode: send them all under 'merge_files'.
        // A single file is sent as 'peptidomic_file' (already in formData).
        if (pepFiles.length >= 2) {
            formData.delete('peptidomic_file');
            for (var i = 0; i < pepFiles.length; i++) {
                formData.append('merge_files', pepFiles[i]);
            }
        }
        ajax('POST', 'upload/', formData, function(resp) {
            btn.disabled = false;
            btn.innerHTML = '<i class="fas fa-upload"></i> Upload';

            var summary = '';
            if (resp.merged_from) {
                summary += '<div class="dt-alert dt-alert-info">' +
                    '<i class="fas fa-object-group"></i> Merged <strong>' + resp.merged_from +
                    '</strong> datasets into one.</div>';
            }
            summary += '<div class="dt-alert dt-alert-success">' +
                'Loaded <strong>' + resp.rows + '</strong> rows, <strong>' + resp.columns + '</strong> columns. ' +
                'Found <strong>' + resp.sequences + '</strong> unique sequences to search.</div>';
            if (resp.warning) {
                summary += '<div class="dt-alert dt-alert-warning">' + resp.warning + '</div>';
            }

            document.getElementById('upload-summary').innerHTML = summary;

            if (resp.has_mbpdb) {
                // Embedded function column detected — skip the Functional
                // Annotation step entirely.
                document.getElementById('functional-card').classList.add('hidden');
                ajax('POST', 'start-blast/', null, function(blastResp) {
                    document.getElementById('upload-summary').innerHTML = summary +
                        '<div class="dt-alert dt-alert-success">' +
                        '<i class="fas fa-check-circle"></i> Functional data loaded: <strong>' +
                        resp.mbpdb_rows + '</strong> records. Skipping database search.</div>' +
                        '<div class="dt-alert dt-alert-info"><i class="fas fa-database"></i> ' +
                        blastResp.message + ' (' + blastResp.count + ' records)</div>' +
                        '<div style="margin-top:12px;"><button class="dt-btn dt-btn-primary" ' +
                        'id="goto-step2-btn-embedded">Continue to Study Variables ' +
                        '<i class="fas fa-arrow-right"></i></button></div>';
                    document.getElementById('goto-step2-btn-embedded')
                        .addEventListener('click', function() { goToStep(2); loadStep2(); });
                });
            } else {
                document.getElementById('functional-card').classList.remove('hidden');
                enableFunctionalActions();
            }
        }, function(msg) {
            btn.disabled = false;
            btn.innerHTML = '<i class="fas fa-upload"></i> Upload';
            showError(msg);
        });
    });

    document.getElementById('start-blast-btn').addEventListener('click', function() {
        var btn = this;
        btn.disabled = true;

        var fnInput = document.getElementById('functional_file');
        var srcEl = document.getElementById('annotation_source');
        var threshEl = document.getElementById('similarity_threshold');
        var blastData = new FormData();
        if (fnInput && fnInput.files.length) { blastData.append('functional_file', fnInput.files[0]); }
        if (srcEl) { blastData.append('annotation_source', srcEl.value); }
        if (threshEl) { blastData.append('similarity_threshold', threshEl.value); }

        ajax('POST', 'start-blast/', blastData, function(resp) {
            document.getElementById('functional-actions').classList.add('hidden');
            if (resp.skipped) {
                document.getElementById('blast-summary').innerHTML =
                    '<div class="dt-alert dt-alert-info">' + resp.message +
                    ' (' + resp.count + ' rows)</div>';
                document.getElementById('blast-continue').classList.remove('hidden');
                return;
            }
            document.getElementById('blast-progress').classList.remove('hidden');
            pollCompletion(resp.task_id, function() {
                ajax('GET', 'blast-results/' + resp.task_id + '/', null, function(r) {
                    document.getElementById('blast-progress').classList.add('hidden');
                    document.getElementById('blast-summary').innerHTML =
                        '<div class="dt-alert dt-alert-success">Search complete. Found <strong>' +
                        r.count + '</strong> matches.</div>';
                    document.getElementById('blast-continue').classList.remove('hidden');
                });
            }, function(errMsg) {
                document.getElementById('blast-progress').classList.add('hidden');
                document.getElementById('functional-actions').classList.remove('hidden');
                btn.disabled = false;
                showError(errMsg || 'Database search failed. Check server logs.');
            });
        }, function(msg) {
            btn.disabled = false;
            showError(msg);
        });
    });

    document.getElementById('skip-blast-btn').addEventListener('click', function() {
        goToStep(2);
        loadStep2();
    });

    document.getElementById('reset-uploads-btn').addEventListener('click', function() {
        if (!confirm('This clears the peptidomic upload and the functional annotation ' +
                     'and starts over. Continue?')) return;
        ajax('POST', 'cleanup/', null, function() { location.reload(); },
             function() { location.reload(); });
    });

    document.getElementById('goto-step2-btn').addEventListener('click', function() {
        goToStep(2);
        loadStep2();
    });

    // -----------------------------------------------------------------------
    // Progress polling (reuses existing /check-progress/ endpoint)
    // -----------------------------------------------------------------------

    // Lightweight poll — no progress bar, just waits for completion/failure.
    function pollCompletion(taskId, onComplete, onFail) {
        var interval = setInterval(function() {
            ajax('GET', '/check-progress/' + taskId + '/', null, function(resp) {
                if (resp.status === 'complete') {
                    clearInterval(interval);
                    onComplete();
                } else if (resp.status === 'failed') {
                    clearInterval(interval);
                    if (onFail) onFail(resp.error || '');
                }
            });
        }, 2000);
    }

    function pollProgress(taskId, barId, textId, onComplete) {
        var bar = document.getElementById(barId);
        var text = document.getElementById(textId);
        var interval = setInterval(function() {
            ajax('GET', '/check-progress/' + taskId + '/', null, function(resp) {
                var pct = resp.percent_progress || 0;
                bar.style.width = pct + '%';

                if (resp.estimated_time_remaining != null) {
                    text.textContent = Math.round(pct) + '% - ~' +
                        Math.ceil(resp.estimated_time_remaining) + 's remaining';
                } else {
                    text.textContent = Math.round(pct) + '%';
                }

                if (resp.status === 'complete') {
                    clearInterval(interval);
                    bar.style.width = '100%';
                    text.textContent = 'Complete!';
                    onComplete();
                } else if (resp.status === 'failed') {
                    clearInterval(interval);
                    text.textContent = 'Failed.';
                    showError('Task failed');
                }
            });
        }, 2000);
    }

    // -----------------------------------------------------------------------
    // Step 2: Tech Replicate Assignment + Study Variable Grouping (combined)
    // -----------------------------------------------------------------------

    var rawAvailableColumns = [];  // pre-collapse columns (for tech rep panel), rename-applied
    var autoTechDupMapping = {};   // always {}: kept for computeBioRepColumns signature only
    var trAvailableColumns = [];   // raw columns available in tech rep panel
    var trSelectedColumns = [];    // columns selected for the current tech rep group being built
    var manualTechReps = [];       // [{name, columns}] — includes both auto-detected (pre-populated) and manually added groups
    var baseColumns = [];          // original (pre-rename) abundance columns
    var manualRenames = [];        // [{new, original}] — column-rename definitions

    // Compute bio rep columns from raw columns + manualTechReps mapping.
    // autoTechDupMapping is always empty; all groups live in manualTechReps.
    function computeBioRepColumns() {
        var colToBioRep = {};
        // (autoTechDupMapping is always {}, kept for historical compat)
        Object.keys(autoTechDupMapping).forEach(function(bioRep) {
            autoTechDupMapping[bioRep].forEach(function(col) { colToBioRep[col] = bioRep; });
        });
        // Manual
        manualTechReps.forEach(function(g) {
            g.columns.forEach(function(col) { colToBioRep[col] = g.name; });
        });

        var result = [];
        var inserted = {};
        rawAvailableColumns.forEach(function(col) {
            var bioRep = colToBioRep[col];
            if (bioRep) {
                if (!inserted[bioRep]) { result.push(bioRep); inserted[bioRep] = true; }
            } else {
                result.push(col);
            }
        });
        return result;
    }

    function loadStep2() {
        ajax('GET', 'step2/', null, function(resp) {
            // Base (original) columns + any saved renames drive the rename panel.
            baseColumns = resp.base_columns || resp.raw_columns || resp.columns || [];
            var savedRenames = resp.column_renames || {};   // {new: original}
            manualRenames = Object.keys(savedRenames).map(function(newName) {
                return {new: newName, original: savedRenames[newName]};
            });

            rawAvailableColumns = resp.raw_columns || resp.columns || [];
            var detectedMapping = resp.tech_dup_mapping || {};
            trAvailableColumns = rawAvailableColumns;
            trSelectedColumns = [];

            // Pre-populate manualTechReps with any auto-detected groups so the
            // user can review, remove, or add to them. autoTechDupMapping is
            // cleared to avoid double-counting in computeBioRepColumns().
            autoTechDupMapping = {};
            manualTechReps = Object.keys(detectedMapping).map(function(bioRep) {
                return {name: bioRep, columns: detectedMapping[bioRep]};
            });

            // Group panel uses bio rep columns (post-collapse)
            availableColumns = computeBioRepColumns();
            selectedColumns = [];

            // Restore previously saved groups (resume) or start fresh.
            definedGroups = (resp.saved_groups && resp.saved_groups.length)
                ? resp.saved_groups.slice()
                : [];

            // Reset inputs
            document.getElementById('rn-new-name-input').value = '';
            document.getElementById('tr-bio-rep-name-input').value = '';
            document.getElementById('tr-column-search').value = '';
            document.getElementById('group-name-input').value = '';
            document.getElementById('column-search').value = '';

            renderRenameOriginalSelect();
            renderRenameList();
            renderTrColumnPanels('');
            renderDefinedTechReps();
            renderGroups();
            renderColumnPanels('');
            document.getElementById('submit-groups-btn').disabled = (definedGroups.length === 0);
        });
    }

    // -----------------------------------------------------------------------
    // Column renames (applied before tech reps / grouping)
    // -----------------------------------------------------------------------

    // Apply the current manualRenames ({new, original}) to a list of columns.
    function applyRenamesToList(cols) {
        if (!manualRenames.length) return cols.slice();
        var originalToNew = {};
        manualRenames.forEach(function(r) { originalToNew[r.original] = r.new; });
        return cols.map(function(c) { return originalToNew.hasOwnProperty(c) ? originalToNew[c] : c; });
    }

    // Recompute rename-dependent column lists and refresh the tech-rep + group
    // panels so the simplified names appear everywhere downstream.
    function refreshAfterRenameChange() {
        rawAvailableColumns = applyRenamesToList(baseColumns);
        trAvailableColumns = rawAvailableColumns;
        availableColumns = computeBioRepColumns();
        renderTrColumnPanels(document.getElementById('tr-column-search').value);
        renderColumnPanels(document.getElementById('column-search').value);
    }

    // Populate the "Original Column" dropdown with base columns not yet renamed.
    function renderRenameOriginalSelect() {
        var sel = document.getElementById('rn-original-select');
        var used = {};
        manualRenames.forEach(function(r) { used[r.original] = true; });
        sel.innerHTML = '';
        var available = baseColumns.filter(function(c) { return !used[c]; });
        if (!available.length) {
            var opt = document.createElement('option');
            opt.value = '';
            opt.textContent = 'All columns already renamed';
            opt.disabled = true;
            sel.appendChild(opt);
            return;
        }
        var placeholder = document.createElement('option');
        placeholder.value = '';
        placeholder.textContent = 'Select a column…';
        sel.appendChild(placeholder);
        available.forEach(function(c) {
            var opt = document.createElement('option');
            opt.value = c;
            opt.textContent = c;
            sel.appendChild(opt);
        });
    }

    function renderRenameList() {
        var list = document.getElementById('rn-list');
        list.innerHTML = '';
        if (!manualRenames.length) {
            document.getElementById('rn-defined-list-wrap').classList.add('hidden');
            return;
        }
        document.getElementById('rn-defined-list-wrap').classList.remove('hidden');
        manualRenames.forEach(function(r, i) {
            var div = document.createElement('div');
            div.className = 'dt-group-item';
            div.innerHTML = '<span><strong>' + escHtml(r.new) + '</strong> &larr; <em>' +
                escHtml(r.original) + '</em></span>' +
                '<span class="remove-group" data-idx="' + i + '"><i class="fas fa-times"></i></span>';
            list.appendChild(div);
        });
        list.querySelectorAll('.remove-group').forEach(function(el) {
            el.addEventListener('click', function() {
                manualRenames.splice(parseInt(this.getAttribute('data-idx')), 1);
                persistRenames();
                renderRenameOriginalSelect();
                renderRenameList();
                refreshAfterRenameChange();
            });
        });
    }

    // Save the current rename mapping server-side so raw_columns stays in sync
    // for the tech-rep step and processing. Sends {new: original}.
    function persistRenames(onDone) {
        var mapping = {};
        manualRenames.forEach(function(r) { mapping[r.new] = r.original; });
        ajax('POST', 'submit-renames/', {renames: mapping}, function(resp) {
            if (onDone) onDone(resp);
        }, function(msg) { showStep2Error(msg, 'rename'); });
    }

    document.getElementById('rename-json-input').addEventListener('change', function() {
        if (!this.files.length) return;
        var formData = new FormData();
        formData.append('rename_file', this.files[0]);
        ajax('POST', 'upload-renames/', formData, function(resp) {
            var mapping = resp.renames || {};
            manualRenames = Object.keys(mapping).map(function(newName) {
                return {new: newName, original: mapping[newName]};
            });
            persistRenames();
            renderRenameOriginalSelect();
            renderRenameList();
            refreshAfterRenameChange();
        }, function(msg) { showStep2Error(msg, 'rename'); });
    });

    document.getElementById('rn-add-btn').addEventListener('click', function() {
        var original = document.getElementById('rn-original-select').value;
        var newName = document.getElementById('rn-new-name-input').value.trim();
        if (!original) { showStep2Error('Please select a column to rename', 'rename'); return; }
        if (!newName) { showStep2Error('Please enter a new name', 'rename'); return; }
        if (manualRenames.some(function(r) { return r.new === newName; })) {
            showStep2Error('The name "' + newName + '" is already used', 'rename'); return;
        }
        manualRenames.push({new: newName, original: original});
        document.getElementById('rn-new-name-input').value = '';
        persistRenames();
        renderRenameOriginalSelect();
        renderRenameList();
        refreshAfterRenameChange();
    });

    // Format a column name for display: _repN suffix → (repN)
    // e.g. "Sample_A_rep1" → "Sample_A (rep1)". Raw name is preserved as tooltip.
    function fmtRepCol(col) {
        return col.replace(/_rep(\d+)$/i, ' (rep$1)');
    }

    // Tech rep dual-panel: shows raw columns minus those already in a manual group
    function renderTrColumnPanels(filter) {
        var f = (filter || '').toLowerCase();
        var availList = document.getElementById('tr-column-available-list');
        var selList = document.getElementById('tr-column-selected-list');
        var availCount = document.getElementById('tr-avail-count');
        var selCount = document.getElementById('tr-selected-count');

        var usedInGroups = {};
        manualTechReps.forEach(function(g) { g.columns.forEach(function(c) { usedInGroups[c] = true; }); });

        availList.innerHTML = '';
        var shown = 0;
        trAvailableColumns.forEach(function(col) {
            if (trSelectedColumns.indexOf(col) !== -1) return;
            if (usedInGroups[col]) return;
            if (f && col.toLowerCase().indexOf(f) === -1) return;
            shown++;
            var div = document.createElement('div');
            div.className = 'dt-col-item';
            div.title = col;
            div.textContent = fmtRepCol(col);
            div.addEventListener('click', function() {
                if (trSelectedColumns.indexOf(col) === -1) {
                    trSelectedColumns.push(col);
                    renderTrColumnPanels(document.getElementById('tr-column-search').value);
                }
            });
            availList.appendChild(div);
        });
        if (!shown) availList.innerHTML = '<div class="dt-col-empty">No columns available</div>';
        availCount.textContent = '(' + shown + ')';

        selList.innerHTML = '';
        if (!trSelectedColumns.length) {
            selList.innerHTML = '<div class="dt-col-empty">Click columns to add</div>';
        } else {
            trSelectedColumns.forEach(function(col) {
                var div = document.createElement('div');
                div.className = 'dt-col-item selected-item';
                div.title = col;
                var txt = document.createElement('span');
                txt.textContent = fmtRepCol(col);
                txt.style.overflow = 'hidden';
                txt.style.textOverflow = 'ellipsis';
                var rm = document.createElement('span');
                rm.className = 'dt-col-remove';
                rm.textContent = '×';
                rm.title = 'Remove';
                rm.addEventListener('click', function(e) {
                    e.stopPropagation();
                    var idx = trSelectedColumns.indexOf(col);
                    if (idx !== -1) trSelectedColumns.splice(idx, 1);
                    renderTrColumnPanels(document.getElementById('tr-column-search').value);
                });
                div.appendChild(txt);
                div.appendChild(rm);
                selList.appendChild(div);
            });
        }
        selCount.textContent = '(' + trSelectedColumns.length + ')';
    }

    function renderDefinedTechReps() {
        var list = document.getElementById('tr-groups-list');
        list.innerHTML = '';
        if (!manualTechReps.length) {
            document.getElementById('tr-defined-groups').classList.add('hidden');
            return;
        }
        document.getElementById('tr-defined-groups').classList.remove('hidden');
        manualTechReps.forEach(function(g, i) {
            var div = document.createElement('div');
            div.className = 'dt-group-item';
            // Count occurrences so duplicate column names are numbered in the display
            var colCounts = {};
            g.columns.forEach(function(c) { colCounts[c] = (colCounts[c] || 0) + 1; });
            var colInstances = {};
            var colLabels = g.columns.map(function(c) {
                if (colCounts[c] > 1) {
                    colInstances[c] = (colInstances[c] || 0) + 1;
                    return escHtml(fmtRepCol(c)) + ' <span style="color:#888;font-size:0.8em;">[' + colInstances[c] + ']</span>';
                }
                return '<em>' + escHtml(fmtRepCol(c)) + '</em>';
            });
            div.innerHTML = '<span><strong>' + g.name + '</strong> &larr; ' +
                colLabels.join(', ') + '</span>' +
                '<span class="remove-group" data-idx="' + i + '"><i class="fas fa-times"></i></span>';
            list.appendChild(div);
        });
        list.querySelectorAll('.remove-group').forEach(function(el) {
            el.addEventListener('click', function() {
                manualTechReps.splice(parseInt(this.getAttribute('data-idx')), 1);
                renderDefinedTechReps();
                renderTrColumnPanels(document.getElementById('tr-column-search').value);
                // Refresh group panel to reflect removed tech rep group
                availableColumns = computeBioRepColumns();
                renderColumnPanels(document.getElementById('column-search').value);
            });
        });
    }

    document.getElementById('tech-rep-json-input').addEventListener('change', function() {
        if (!this.files.length) return;
        var formData = new FormData();
        formData.append('tech_rep_file', this.files[0]);
        ajax('POST', 'upload-tech-reps/', formData, function(resp) {
            manualTechReps = resp.groups || [];
            renderDefinedTechReps();
            renderTrColumnPanels(document.getElementById('tr-column-search').value);
            availableColumns = computeBioRepColumns();
            renderColumnPanels(document.getElementById('column-search').value);
        }, function(msg) { showStep2Error(msg, 'techrep'); });
    });

    document.getElementById('tr-column-search').addEventListener('input', function() {
        renderTrColumnPanels(this.value);
    });

    document.getElementById('tr-add-group-btn').addEventListener('click', function() {
        var name = document.getElementById('tr-bio-rep-name-input').value.trim();
        if (!name) { showStep2Error('Please enter a biological replicate name', 'techrep'); return; }
        if (trSelectedColumns.length < 2) { showStep2Error('Please select at least 2 columns to form a technical replicate group', 'techrep'); return; }
        manualTechReps.push({name: name, columns: trSelectedColumns.slice()});
        document.getElementById('tr-bio-rep-name-input').value = '';
        trSelectedColumns = [];
        renderDefinedTechReps();
        renderTrColumnPanels(document.getElementById('tr-column-search').value);
        // Refresh group panel so the new bio rep name appears in available columns
        availableColumns = computeBioRepColumns();
        renderColumnPanels(document.getElementById('column-search').value);
    });

    document.getElementById('back-to-step1-btn').addEventListener('click', function() {
        goToStep(1);
    });

    // Categorical group dual-panel
    function renderColumnPanels(filter) {
        var f = (filter || '').toLowerCase();
        var availList = document.getElementById('column-available-list');
        var selList = document.getElementById('column-selected-list');
        var availCount = document.getElementById('avail-count');
        var selCount = document.getElementById('selected-count');

        availList.innerHTML = '';
        var shown = 0;
        availableColumns.forEach(function(col) {
            if (selectedColumns.indexOf(col) !== -1) return;
            if (f && col.toLowerCase().indexOf(f) === -1) return;
            shown++;
            var div = document.createElement('div');
            div.className = 'dt-col-item';
            div.title = col;
            div.textContent = col;
            div.addEventListener('click', function() {
                if (selectedColumns.indexOf(col) === -1) {
                    selectedColumns.push(col);
                    renderColumnPanels(document.getElementById('column-search').value);
                }
            });
            availList.appendChild(div);
        });
        if (!shown) availList.innerHTML = '<div class="dt-col-empty">No columns available</div>';
        availCount.textContent = '(' + shown + ')';

        selList.innerHTML = '';
        if (!selectedColumns.length) {
            selList.innerHTML = '<div class="dt-col-empty">Click columns to add</div>';
        } else {
            selectedColumns.forEach(function(col) {
                var div = document.createElement('div');
                div.className = 'dt-col-item selected-item';
                div.title = col;
                var txt = document.createElement('span');
                txt.textContent = col;
                txt.style.overflow = 'hidden';
                txt.style.textOverflow = 'ellipsis';
                var rm = document.createElement('span');
                rm.className = 'dt-col-remove';
                rm.textContent = '×';
                rm.title = 'Remove';
                rm.addEventListener('click', function(e) {
                    e.stopPropagation();
                    var idx = selectedColumns.indexOf(col);
                    if (idx !== -1) selectedColumns.splice(idx, 1);
                    renderColumnPanels(document.getElementById('column-search').value);
                });
                div.appendChild(txt);
                div.appendChild(rm);
                selList.appendChild(div);
            });
        }
        selCount.textContent = '(' + selectedColumns.length + ')';
    }

    document.getElementById('column-search').addEventListener('input', function() {
        renderColumnPanels(this.value);
    });

    document.getElementById('group-json-input').addEventListener('change', function() {
        if (!this.files.length) return;
        var formData = new FormData();
        formData.append('group_file', this.files[0]);
        ajax('POST', 'upload-groups/', formData, function(resp) {
            definedGroups = [];
            Object.keys(resp.groups).forEach(function(name) {
                definedGroups.push({name: name, columns: resp.groups[name] || []});
            });
            renderGroups();
            document.getElementById('submit-groups-btn').disabled = false;
        }, function(msg) { showStep2Error(msg, 'groups'); });
    });

    document.getElementById('add-group-btn').addEventListener('click', function() {
        var name = document.getElementById('group-name-input').value.trim();
        if (!name) { showStep2Error('Please enter a group name', 'groups'); return; }
        if (!selectedColumns.length) { showStep2Error('Please select at least one column', 'groups'); return; }
        definedGroups.push({name: name, columns: selectedColumns.slice()});
        renderGroups();
        document.getElementById('group-name-input').value = '';
        selectedColumns = [];
        renderColumnPanels(document.getElementById('column-search').value);
        document.getElementById('submit-groups-btn').disabled = false;
    });

    function renderGroups() {
        var list = document.getElementById('groups-list');
        list.innerHTML = '';
        if (!definedGroups.length) {
            document.getElementById('defined-groups').classList.add('hidden');
            return;
        }
        document.getElementById('defined-groups').classList.remove('hidden');
        var dragSrcIdx = null;
        definedGroups.forEach(function(g, i) {
            var div = document.createElement('div');
            div.className = 'dt-group-item';
            div.draggable = true;
            div.innerHTML =
                '<span class="dt-group-drag-handle" title="Drag to reorder"><i class="fas fa-grip-vertical"></i></span>' +
                '<span style="flex:1;"><strong>' + escHtml(g.name) + '</strong> (' + g.columns.length + ' columns)</span>' +
                '<span class="remove-group" data-idx="' + i + '"><i class="fas fa-times"></i></span>';

            div.addEventListener('dragstart', function(e) {
                dragSrcIdx = i;
                e.dataTransfer.effectAllowed = 'move';
                div.classList.add('dt-group-dragging');
            });
            div.addEventListener('dragend', function() {
                div.classList.remove('dt-group-dragging');
                list.querySelectorAll('.dt-group-item').forEach(function(el) {
                    el.classList.remove('dt-group-drag-over');
                });
            });
            div.addEventListener('dragover', function(e) {
                e.preventDefault();
                e.dataTransfer.dropEffect = 'move';
                list.querySelectorAll('.dt-group-item').forEach(function(el) {
                    el.classList.remove('dt-group-drag-over');
                });
                if (dragSrcIdx !== i) div.classList.add('dt-group-drag-over');
            });
            div.addEventListener('drop', function(e) {
                e.preventDefault();
                if (dragSrcIdx !== null && dragSrcIdx !== i) {
                    var moved = definedGroups.splice(dragSrcIdx, 1)[0];
                    definedGroups.splice(i, 0, moved);
                    renderGroups();
                }
            });

            list.appendChild(div);
        });
        list.querySelectorAll('.remove-group').forEach(function(el) {
            el.addEventListener('click', function() {
                definedGroups.splice(parseInt(this.getAttribute('data-idx')), 1);
                renderGroups();
                if (!definedGroups.length) document.getElementById('submit-groups-btn').disabled = true;
            });
        });
    }

    // Reset a dropzone file input + its "drag & drop" visual back to empty state.
    function clearDropzone(inputId) {
        var input = document.getElementById(inputId);
        if (!input) return;
        input.value = '';
        var zone = input.closest('.dt-dropzone');
        if (!zone) return;
        var dt = zone.querySelector('.drop-text');
        if (dt) dt.classList.remove('hidden');
        var fn = zone.querySelector('.file-name');
        if (fn) { fn.textContent = ''; fn.classList.add('hidden'); }
    }

    document.getElementById('reset-groups-btn').addEventListener('click', function() {
        // Clear column renames (inputs, uploaded file, and the defined list).
        manualRenames = [];
        document.getElementById('rn-new-name-input').value = '';
        clearDropzone('rename-json-input');
        renderRenameOriginalSelect();
        renderRenameList();
        // Revert the tech-rep/group panels to the original (pre-rename) columns.
        refreshAfterRenameChange();
        // Persist the empty mapping so column_renames.json is cleared and the
        // raw column list reverts server-side too.
        persistRenames();

        // Clear grouping client-side state
        definedGroups = [];
        selectedColumns = [];
        document.getElementById('group-name-input').value = '';
        document.getElementById('column-search').value = '';
        renderGroups();
        renderColumnPanels('');
        document.getElementById('submit-groups-btn').disabled = true;

        // Clear the uploaded group-definition file input + its dropzone display
        clearDropzone('group-json-input');

        // Delete the server-side saved group definitions so a resume /
        // step-2 re-fetch does not restore the prior grouping variable.
        ajax('POST', 'reset-groups/', null, function() {}, function() {});
    });

    // Save Groups: submit tech reps first (to save mapping), then submit groups
    document.getElementById('submit-groups-btn').addEventListener('click', function() {
        var btn = this;
        btn.disabled = true;
        ajax('POST', 'submit-tech-reps/', {tech_reps: manualTechReps}, function() {
            ajax('POST', 'submit-groups/', {groups: definedGroups}, function() {
                goToStep(3);
                loadStep3();
            }, function(msg) { btn.disabled = false; showStep2Error(msg, 'actions'); });
        }, function(msg) { btn.disabled = false; showStep2Error(msg, 'actions'); });
    });

    // Skip Grouping: still save tech reps, then skip groups using bio rep columns
    document.getElementById('skip-groups-btn').addEventListener('click', function() {
        var btn = this;
        btn.disabled = true;
        ajax('POST', 'submit-tech-reps/', {tech_reps: manualTechReps}, function(resp) {
            var bioRepCols = resp.columns || availableColumns;
            ajax('POST', 'skip-groups/', {columns: bioRepCols}, function() {
                btn.disabled = false;
                goToStep(3);
                loadStep3();
            }, function(msg) { btn.disabled = false; showStep2Error(msg, 'actions'); });
        }, function(msg) { btn.disabled = false; showStep2Error(msg, 'actions'); });
    });

    // -----------------------------------------------------------------------
    // Step 3: Protein Mapping
    // -----------------------------------------------------------------------

    function loadStep3() {
        ajax('GET', 'step3/', null, function(resp) {
            renderStep3Info(resp);
        });
    }

    function renderStep3Info(resp) {
        var info = document.getElementById('protein-info');
        var html = '';

        // Known proteins summary
        if (resp.known_protein_count > 0) {
            html += '<div class="dt-alert dt-alert-success" style="margin-bottom: 8px;">' +
                '<strong>' + resp.known_protein_count + '</strong> protein(s) mapped from database.</div>';
            if (resp.known_proteins && resp.known_proteins.length > 0) {
                html += '<div style="max-height: 120px; overflow-y: auto; font-size: 0.82rem; background: #f8f9fa; border: 1px solid #dee2e6; border-radius: 4px; padding: 6px 10px; margin-bottom: 8px;">';
                resp.known_proteins.forEach(function(p) {
                    html += '<span style="display:inline-block; margin-right:16px; margin-bottom:2px;">' +
                        '<strong>' + p.id + '</strong>' +
                        (p.name ? ' — ' + p.name.split(' ').slice(0,4).join(' ') : '') +
                        (p.species ? ' <em>(' + p.species + ')</em>' : '') + '</span>';
                });
                if (resp.known_protein_count > 50) html += '...';
                html += '</div>';
            }
        }

        // Missing proteins
        if (resp.missing_protein_count > 0) {
            html += '<div class="dt-alert dt-alert-warning">' +
                '<strong>' + resp.missing_protein_count + '</strong> protein(s) not found in the dictionary. ' +
                'Fetch from UniProt to get names and species, or skip.</div>';
            if (resp.missing_ids && resp.missing_ids.length > 0) {
                html += '<div style="max-height: 100px; overflow-y: auto; font-size: 0.82rem; color: #666; margin-top: 4px;">' +
                    resp.missing_ids.join(', ') + (resp.missing_protein_count > 50 ? ', ...' : '') + '</div>';
            }
        } else if (resp.known_protein_count === 0) {
            html += '<div class="dt-alert dt-alert-info">No protein identifiers found in the data. ' +
                'Ensure your file has a protein/accession column.</div>';
        } else {
            html += '<div class="dt-alert dt-alert-success">All proteins are mapped.</div>';
        }

        info.innerHTML = html;

        if (resp.has_combinations && resp.combinations.length > 0) {
            document.getElementById('protein-combinations').classList.remove('hidden');
            renderCombinations(resp.combinations, resp.saved_decisions || {});
        } else {
            document.getElementById('protein-combinations').classList.add('hidden');
        }

        renderSourceRenames(resp.protein_sources || [], resp.saved_source_renames || []);
    }

    // -----------------------------------------------------------------------
    // Merge / Rename Protein Sources
    // -----------------------------------------------------------------------

    var _proteinSources = [];   // [{id, name, species, count}] from step3/
    var _mergeGroupSeq = 0;     // unique id per rendered card

    function renderSourceRenames(sources, savedGroups) {
        _proteinSources = sources || [];
        var panel = document.getElementById('source-renames');
        var container = document.getElementById('source-rename-groups');
        if (!panel || !container) return;

        container.innerHTML = '';
        if (_proteinSources.length === 0) {
            // Nothing to merge (e.g. no protein column) — hide the whole panel.
            panel.classList.add('hidden');
            return;
        }
        panel.classList.remove('hidden');

        (savedGroups || []).forEach(function(g) { addMergeGroupCard(g); });
    }

    function _sourceById(id) {
        for (var i = 0; i < _proteinSources.length; i++) {
            if (_proteinSources[i].id === id) return _proteinSources[i];
        }
        return null;
    }

    // IDs the combination step is actively changing (non "keep combined"),
    // read live from the on-screen controls, for the inline conflict hint.
    function getTouchedComboIds() {
        var touched = {};
        document.querySelectorAll('#combinations-list .dt-combo-row').forEach(function(row, idx) {
            var modeRadio = row.querySelector('input[name="combo_mode_' + idx + '"]:checked');
            var mode = modeRadio ? modeRadio.value : 'asis';
            if (mode === 'asis') return;
            var comboStrong = row.querySelector('strong');
            var combo = comboStrong ? comboStrong.textContent : '';
            combo.split(';').forEach(function(p) {
                var t = p.trim();
                if (t) touched[t] = true;
            });
            if (mode === 'custom') {
                var ci = row.querySelector('#custom_' + idx);
                if (ci && ci.value.trim()) touched[ci.value.trim()] = true;
            }
        });
        return touched;
    }

    function addMergeGroupCard(prefill) {
        prefill = prefill || {};
        var container = document.getElementById('source-rename-groups');
        if (!container) return;
        var idx = _mergeGroupSeq++;
        var savedSources = prefill.sources || [];

        var card = document.createElement('div');
        card.className = 'dt-merge-group';
        card.style.cssText = 'border:1px solid #dee2e6; border-radius:6px; padding:10px 12px; margin-bottom:10px; background:#fff;';

        var html = '';
        html += '<div style="display:flex; justify-content:space-between; align-items:center; margin-bottom:6px;">' +
            '<strong style="font-size:0.9rem; color:#555;">Merge group</strong>' +
            '<button type="button" class="dt-btn dt-btn-secondary merge-remove-btn" ' +
            'style="padding:2px 8px; font-size:0.8rem;"><i class="fas fa-times"></i> Remove</button></div>';

        html += '<div style="font-size:0.8rem; color:#888; margin-bottom:4px;">Protein sources to merge / rename:</div>';
        html += '<div style="position:relative; margin-bottom:6px;">' +
            '<i class="fas fa-search" style="position:absolute; left:8px; top:50%; transform:translateY(-50%); ' +
            'color:#adb5bd; font-size:0.8rem; pointer-events:none;"></i>' +
            '<input type="text" class="merge-source-search" placeholder="Filter by ID, name, or species…" ' +
            'style="width:100%; padding:4px 8px 4px 26px; box-sizing:border-box; font-size:0.85rem;"></div>';
        html += '<div class="merge-sources" style="max-height:150px; overflow-y:auto; background:#f8f9fa; ' +
            'border:1px solid #dee2e6; border-radius:4px; padding:6px 10px; margin-bottom:8px;">';
        _proteinSources.forEach(function(s) {
            var checked = savedSources.indexOf(s.id) !== -1 ? ' checked' : '';
            var haystack = [s.id, s.name, s.species].join(' ').toLowerCase();
            html += '<label class="merge-source-row" data-search="' + escHtml(haystack) + '" ' +
                'style="display:block; margin:3px 0; cursor:pointer;">' +
                '<input type="checkbox" class="merge-source-cb" value="' + escHtml(s.id) + '"' + checked + '> ' +
                '<strong>' + escHtml(s.id) + '</strong>' +
                (s.name ? ' — ' + escHtml(s.name) : '') +
                (s.species ? ' <em>(' + escHtml(s.species) + ')</em>' : '') +
                ' <span style="color:#999;">&middot; ' + s.count + ' row' + (s.count === 1 ? '' : 's') + '</span>' +
                '</label>';
        });
        html += '<div class="merge-sources-empty" style="display:none; color:#999; font-size:0.82rem; padding:4px 0;">' +
            'No proteins match your filter.</div>';
        html += '</div>';

        html += '<div style="display:flex; gap:10px; flex-wrap:wrap; align-items:flex-end;">';
        html += '<label style="font-size:0.85rem;">Merge into (target ID)<br>' +
            '<input type="text" class="merge-target-id" placeholder="e.g. P02666" ' +
            'value="' + escHtml(prefill.target_id || '') + '" style="width:160px;"></label>';
        html += '<label style="font-size:0.85rem;">Name <span style="color:#999;">(optional)</span><br>' +
            '<input type="text" class="merge-target-name" placeholder="e.g. Protein name" ' +
            'value="' + escHtml(prefill.target_name || '') + '" style="width:200px;"></label>';
        html += '</div>';

        html += '<div class="merge-conflict-hint" style="display:none; margin-top:8px;"></div>';

        card.innerHTML = html;
        container.appendChild(card);

        var targetId = card.querySelector('.merge-target-id');
        var targetName = card.querySelector('.merge-target-name');

        function autofillFromSource(id) {
            var s = _sourceById(id);
            if (!s) return;
            if (!targetName.value.trim() && s.name) targetName.value = s.name;
        }

        card.querySelectorAll('.merge-source-cb').forEach(function(cb) {
            cb.addEventListener('change', function() {
                if (this.checked && !targetId.value.trim()) {
                    // First pick becomes the default canonical target.
                    targetId.value = this.value;
                    autofillFromSource(this.value);
                }
                updateConflictHint(card);
            });
        });
        targetId.addEventListener('input', function() {
            autofillFromSource(this.value.trim());
            updateConflictHint(card);
        });

        card.querySelector('.merge-remove-btn').addEventListener('click', function() {
            card.parentNode.removeChild(card);
        });

        // Live search: hide rows that don't match the query. Checked rows always
        // stay visible so an active selection can't be filtered out of sight.
        var searchInput = card.querySelector('.merge-source-search');
        var emptyMsg = card.querySelector('.merge-sources-empty');
        searchInput.addEventListener('input', function() {
            var q = this.value.trim().toLowerCase();
            var anyVisible = false;
            card.querySelectorAll('.merge-source-row').forEach(function(row) {
                var cb = row.querySelector('.merge-source-cb');
                var match = !q || row.getAttribute('data-search').indexOf(q) !== -1 ||
                    (cb && cb.checked);
                row.style.display = match ? 'block' : 'none';
                if (match) anyVisible = true;
            });
            if (emptyMsg) emptyMsg.style.display = anyVisible ? 'none' : 'block';
        });

        updateConflictHint(card);
    }

    function _cardGroup(card) {
        var sources = Array.prototype.slice.call(
            card.querySelectorAll('.merge-source-cb:checked')).map(function(cb) { return cb.value; });
        return {
            sources: sources,
            target_id: card.querySelector('.merge-target-id').value.trim(),
            target_name: card.querySelector('.merge-target-name').value.trim()
        };
    }

    function updateConflictHint(card) {
        var hint = card.querySelector('.merge-conflict-hint');
        if (!hint) return;
        var g = _cardGroup(card);
        var ids = g.sources.slice();
        if (g.target_id) ids.push(g.target_id);
        var touched = getTouchedComboIds();
        var overlap = ids.filter(function(id) { return touched[id]; });
        // de-duplicate
        overlap = overlap.filter(function(v, i) { return overlap.indexOf(v) === i; });
        if (overlap.length) {
            hint.style.display = 'block';
            hint.innerHTML = '<div class="dt-alert dt-alert-warning" style="margin:0; font-size:0.82rem;">' +
                '<i class="fas fa-exclamation-triangle"></i> ' + escHtml(overlap.join(', ')) +
                ' also has a combined-protein decision above. This merge/rename is applied last and will override it.</div>';
        } else {
            hint.style.display = 'none';
            hint.innerHTML = '';
        }
    }

    function collectSourceRenames() {
        var groups = [];
        document.querySelectorAll('#source-rename-groups .dt-merge-group').forEach(function(card) {
            var g = _cardGroup(card);
            if (g.sources.length > 0 && g.target_id) groups.push(g);
        });
        return groups;
    }

    function renderCombinations(combos, savedDecisions) {
        savedDecisions = savedDecisions || {};
        savedDecisions = savedDecisions || {};
        var list = document.getElementById('combinations-list');
        list.innerHTML = '';
        combos.forEach(function(combo, idx) {
            // Restore a previously saved decision for this combination, if any.
            var saved = savedDecisions[combo.combo] || null;
            var mode = saved ? saved.mode : 'split';           // default: split
            var savedIds = (saved && saved.protein_ids) ? saved.protein_ids : null;
            var customVal = (saved && saved.mode === 'custom') ? (saved.protein_id || '') : '';
            var customNameVal = (saved && saved.mode === 'custom') ? (saved.protein_name || '') : '';
            var splitVisible = (mode === 'split');

            var div = document.createElement('div');
            div.className = 'dt-combo-row';
            var html = '<strong>' + escHtml(combo.combo) + '</strong> (' + combo.occurrences + ' rows)<br>';
            html += '<div class="combo-options" style="margin-top: 8px;">';

            // Mode 1: Keep combined
            html += '<label style="display: block; margin: 4px 0;">' +
                '<input type="radio" name="combo_mode_' + idx + '" value="asis"' +
                (mode === 'asis' ? ' checked' : '') + '> ' +
                'Keep combined ID</label>';

            // Mode 2: Split
            html += '<label style="display: block; margin: 4px 0;">' +
                '<input type="radio" name="combo_mode_' + idx + '" value="split"' +
                (mode === 'split' ? ' checked' : '') + '> ' +
                'Split into individual proteins:</label>';
            html += '<div class="combo-split-panel" id="split_panel_' + idx + '" ' +
                'style="display:' + (splitVisible ? 'block' : 'none') + '; margin-left:24px; margin-top:4px; padding:8px; ' +
                'background:#f8f9fa; border:1px solid #dee2e6; border-radius:4px;">';
            html += '<div style="font-size:0.8rem; color:#888; margin-bottom:6px;">' +
                'Checked proteins will each get their own row; unchecked proteins are removed.</div>';
            combo.proteins.forEach(function(p) {
                // When a split decision was saved, honor its selection; otherwise
                // fall back to the data-derived default.
                var isChecked = savedIds
                    ? (savedIds.indexOf(p.id) !== -1)
                    : (p.default_decision === 'new');
                var chk = isChecked ? ' checked' : '';
                html += '<label style="display:block; margin:3px 0; cursor:pointer;">' +
                    '<input type="checkbox" class="combo-protein-cb" ' +
                    'name="combo_proteins_' + idx + '" value="' + escHtml(p.id) + '"' + chk + '> ' +
                    '<strong>' + escHtml(p.id) + '</strong>' +
                    (p.name ? ' — ' + escHtml(p.name) : '') +
                    (p.species ? ' <em>(' + escHtml(p.species) + ')</em>' : '') + '</label>';
            });
            html += '</div>';

            // Mode 3: Custom ID — restore the saved custom ID + name when present.
            html += '<label style="display:block; margin:4px 0;">' +
                '<input type="radio" name="combo_mode_' + idx + '" id="mode_custom_' + idx + '" value="custom"' +
                (mode === 'custom' ? ' checked' : '') + '> ' +
                'Custom ID: <input type="text" id="custom_' + idx + '" ' +
                'placeholder="Enter protein ID" value="' + escHtml(customVal) + '" ' +
                'style="margin-left:6px; width:160px; display:inline;" ' +
                'onfocus="document.getElementById(\'mode_custom_' + idx + '\').checked=true;' +
                'document.getElementById(\'split_panel_' + idx + '\').style.display=\'none\';"> ' +
                'Name <span style="color:#999;">(optional)</span>: ' +
                '<input type="text" id="custom_name_' + idx + '" ' +
                'placeholder="e.g. Protein name" value="' + escHtml(customNameVal) + '" ' +
                'style="margin-left:4px; width:160px; display:inline;" ' +
                'onfocus="document.getElementById(\'mode_custom_' + idx + '\').checked=true;' +
                'document.getElementById(\'split_panel_' + idx + '\').style.display=\'none\';"></label>';

            html += '</div>';  // .combo-options
            div.innerHTML = html;

            // Show/hide split panel when mode changes
            div.querySelectorAll('input[name="combo_mode_' + idx + '"]').forEach(function(radio) {
                radio.addEventListener('change', function() {
                    div.querySelector('#split_panel_' + idx).style.display =
                        this.value === 'split' ? 'block' : 'none';
                });
            });

            // Checking any protein checkbox auto-activates "split" mode
            div.querySelectorAll('.combo-protein-cb').forEach(function(cb) {
                cb.addEventListener('change', function() {
                    var splitRadio = div.querySelector(
                        'input[name="combo_mode_' + idx + '"][value="split"]');
                    if (splitRadio) {
                        splitRadio.checked = true;
                        div.querySelector('#split_panel_' + idx).style.display = 'block';
                    }
                });
            });

            list.appendChild(div);
        });
    }

    document.getElementById('back-to-step2-btn').addEventListener('click', function() {
        goToStep(2);
        loadStep2();
    });

    document.getElementById('fetch-uniprot-btn').addEventListener('click', function() {
        var btn = this;
        btn.disabled = true;

        function afterFetch(savedCount) {
            // Persist any fetched results then fully refresh the protein panel
            ajax('POST', 'save-uniprot/', null, function(saveResp) {
                var n = (saveResp && saveResp.saved) ? saveResp.saved : (savedCount || 0);
                ajax('GET', 'step3/', null, function(r) {
                    renderStep3Info(r);
                    if (n > 0) {
                        document.getElementById('protein-info').insertAdjacentHTML('afterbegin',
                            '<div class="dt-alert dt-alert-success" style="margin-bottom:8px;">' +
                            '<i class="fas fa-check-circle"></i> UniProt identified <strong>' +
                            n + '</strong> protein(s). Protein map updated.</div>');
                    }
                    btn.disabled = false;
                });
            }, function() {
                ajax('GET', 'step3/', null, function(r) { renderStep3Info(r); });
                btn.disabled = false;
            });
        }

        ajax('POST', 'start-uniprot/', null, function(resp) {
            if (resp.skipped) {
                // No missing proteins — still refresh the panel in case a previous
                // fetch left results that haven't been rendered yet
                afterFetch(0);
                return;
            }
            document.getElementById('uniprot-progress').classList.remove('hidden');
            pollProgress(resp.task_id, 'uniprot-progress-bar', 'uniprot-progress-text', function() {
                document.getElementById('uniprot-progress').classList.add('hidden');
                afterFetch(resp.count);
            });
        }, function(msg) { btn.disabled = false; showError(msg); });
    });

    // Custom FASTA file — auto-validate and merge into the protein database
    // on file selection or drop, then reprocess the mapping (same as a
    // UniProt fetch) so any previously-ambiguous/unknown proteins pick up
    // the newly-supplied species and name information.
    document.getElementById('fasta_file').addEventListener('change', function() {
        if (!this.files.length) return;
        var statusEl = document.getElementById('fasta-upload-status');
        statusEl.innerHTML = '<div class="dt-alert dt-alert-info" style="margin-top:6px;">' +
            '<span class="dt-spinner"></span> Validating & applying FASTA file...</div>';

        var formData = new FormData();
        formData.append('fasta_file', this.files[0]);
        ajax('POST', 'upload-fasta/', formData, function(resp) {
            statusEl.innerHTML = '<div class="dt-alert dt-alert-success" style="margin-top:6px;">' +
                '<i class="fas fa-check-circle"></i> Loaded <strong>' + resp.added +
                '</strong> protein sequence(s) from FASTA file.</div>';
            ajax('GET', 'step3/', null, function(r) {
                renderStep3Info(r);
                document.getElementById('protein-info').insertAdjacentHTML('afterbegin',
                    '<div class="dt-alert dt-alert-success" style="margin-bottom:8px;">' +
                    '<i class="fas fa-check-circle"></i> Protein mapping re-processed using the ' +
                    'uploaded FASTA file.</div>');
            });
        }, function(msg) {
            statusEl.innerHTML = '<div class="dt-alert dt-alert-danger" style="margin-top:6px;">' +
                msg + '</div>';
        });
    });

    // Upload mapping key — auto-apply on file selection or drop
    document.getElementById('protein-map-input').addEventListener('change', function() {
        if (!this.files.length) return;
        var statusEl = document.getElementById('protein-map-upload-status');
        statusEl.innerHTML = '<div class="dt-alert dt-alert-info" style="margin-top:6px;">' +
            '<span class="dt-spinner"></span> Applying mapping key...</div>';

        var formData = new FormData();
        formData.append('map_file', this.files[0]);
        ajax('POST', 'upload-protein-map/', formData, function(resp) {
            statusEl.innerHTML = '<div class="dt-alert dt-alert-success" style="margin-top:6px;">' +
                '<i class="fas fa-check-circle"></i> Mapping key applied — ' +
                resp.applied + ' combination(s) processed. Continuing...</div>';
            setTimeout(function() { goToStep(4); }, 800);
        }, function(msg) {
            statusEl.innerHTML = '<div class="dt-alert dt-alert-danger" style="margin-top:6px;">' +
                msg + '</div>';
        });
    });

    document.getElementById('submit-proteins-btn').addEventListener('click', function() {
        var decisions = {};
        document.querySelectorAll('#combinations-list .dt-combo-row').forEach(function(row, idx) {
            var combo = row.querySelector('strong').textContent;
            var modeRadio = row.querySelector('input[name="combo_mode_' + idx + '"]:checked');
            var mode = modeRadio ? modeRadio.value : 'asis';

            if (mode === 'asis') {
                decisions[combo] = {action: 'ASIS'};
            } else if (mode === 'split') {
                var checkedBoxes = Array.prototype.slice.call(
                    row.querySelectorAll('input[name="combo_proteins_' + idx + '"]:checked')
                );
                var selectedIds = checkedBoxes.map(function(cb) { return cb.value; });
                if (selectedIds.length === 0) {
                    // Nothing selected in split mode → treat as keep combined
                    decisions[combo] = {action: 'ASIS'};
                } else {
                    decisions[combo] = {action: 'SPLIT', protein_ids: selectedIds};
                }
            } else if (mode === 'custom') {
                var customInput = row.querySelector('#custom_' + idx);
                var customNameInput = row.querySelector('#custom_name_' + idx);
                decisions[combo] = {action: 'CUSTOM', protein_id: customInput ? customInput.value.trim() : ''};
                var customNameVal2 = customNameInput ? customNameInput.value.trim() : '';
                if (customNameVal2) decisions[combo].protein_name = customNameVal2;
            }
        });
        var sourceRenames = collectSourceRenames();
        ajax('POST', 'submit-proteins/', {decisions: decisions, source_renames: sourceRenames}, function(resp) {
            if (resp && resp.warnings && resp.warnings.length) {
                var info = document.getElementById('protein-info');
                if (info) {
                    info.insertAdjacentHTML('afterbegin',
                        '<div class="dt-alert dt-alert-warning" style="margin-bottom:8px;">' +
                        '<i class="fas fa-exclamation-triangle"></i> ' +
                        resp.warnings.map(escHtml).join('<br>') + '</div>');
                }
            }
            goToStep(4);
        });
    });

    document.getElementById('reset-proteins-btn').addEventListener('click', function() {
        // Clear the uploaded mapping-key file input + its dropzone display and status.
        clearDropzone('protein-map-input');
        var statusEl = document.getElementById('protein-map-upload-status');
        if (statusEl) statusEl.innerHTML = '';

        // Clear the uploaded custom FASTA file input + its dropzone display and status.
        clearDropzone('fasta_file');
        var fastaStatusEl = document.getElementById('fasta-upload-status');
        if (fastaStatusEl) fastaStatusEl.innerHTML = '';

        // Delete the server-side saved decisions and merge/rename groups, then
        // re-fetch step 3 so combinations and merge groups redraw at their
        // untouched defaults.
        ajax('POST', 'reset-proteins/', null, function() { loadStep3(); }, function() {});
    });

    var addMergeBtn = document.getElementById('add-merge-group-btn');
    if (addMergeBtn) {
        addMergeBtn.addEventListener('click', function() { addMergeGroupCard(); });
    }

    document.getElementById('skip-proteins-btn').addEventListener('click', function() {
        ajax('POST', 'skip-proteins/', null, function() {
            goToStep(4);
        });
    });

    // -----------------------------------------------------------------------
    // Step 4: Process & Export
    // -----------------------------------------------------------------------

    document.getElementById('back-to-step3-btn').addEventListener('click', function() {
        goToStep(3);
        loadStep3();
    });

    document.getElementById('process-btn').addEventListener('click', function() {
        var btn = this;
        btn.disabled = true;
        document.getElementById('process-progress').classList.remove('hidden');

        ajax('POST', 'process/', null, function(resp) {
            btn.disabled = false;
            document.getElementById('process-progress').classList.add('hidden');
            document.getElementById('export-section').classList.remove('hidden');

            document.getElementById('process-summary').innerHTML =
                '<div class="dt-alert dt-alert-success">Processing complete. ' +
                'Result: <strong>' + resp.rows + '</strong> rows, <strong>' +
                resp.columns + '</strong> columns.</div>';

            renderExportButtons(resp.exports);
            document.getElementById('viz-links-section').classList.remove('hidden');
        }, function(msg, resp) {
            btn.disabled = false;
            document.getElementById('process-progress').classList.add('hidden');
            showError(msg);
            var detail = (resp && (resp.detail || resp.traceback)) ? (resp.detail || resp.traceback) : null;
            if (detail) {
                document.getElementById('export-section').classList.remove('hidden');
                document.getElementById('process-summary').innerHTML =
                    '<div class="dt-alert dt-alert-danger">' +
                    '<strong>Error:</strong> ' + escHtml(msg) +
                    '<pre style="margin-top:8px;font-size:0.75rem;white-space:pre-wrap;' +
                    'word-break:break-all;max-height:220px;overflow-y:auto;' +
                    'background:#fff3;border-radius:3px;padding:6px;">' +
                    escHtml(detail) + '</pre></div>';
            }
        });
    });

    var currentViewerSheets = [];

    function renderExportButtons(exports) {
        var container = document.getElementById('export-buttons');
        container.innerHTML = '';

        var items = [
            {key: 'merged_dataset', icon: 'fa-table', label: 'Merged Dataset',
             title: 'Full peptide-level dataset with abundances, functional annotations and study-variable groupings'},
            {key: 'mbpdb_results', icon: 'fa-database', label: 'Functional Database Results',
             title: 'Functional database search matches for each queried peptide, with the database sequence and reported bioactivity'},
            {key: 'summed_function', icon: 'fa-flask', label: 'Summed Functional Data',
             title: 'Peptide abundance and count summed by functional category for each group'},
            {key: 'column_rename_key', icon: 'fa-pen-to-square', label: 'Column Rename Key',
             title: 'Reusable key mapping instrument-generated column names to your simplified sample names'},
            {key: 'tech_rep_key', icon: 'fa-key', label: 'Technical Replicate Key',
             title: 'Reusable key mapping each biological replicate to its technical replicate columns'},
            {key: 'tech_rep_correlation', icon: 'fa-vials', label: 'Within-Sample Correlations (Technical Replicates)', corr: true,
             title: 'Within one sample: agreement between its technical replicate runs'},
            {key: 'group_definitions', icon: 'fa-layer-group', label: 'Group Definitions',
             title: 'Study-variable group assignments for each abundance column'},
            {key: 'replicate_correlation', icon: 'fa-chart-area', label: 'Within-Group Correlations (Biological Replicates)', corr: true,
             title: 'Within one group: agreement between its biological replicate samples'},
            {key: 'protein_map', icon: 'fa-file-export', label: 'Protein Mapping Key',
             download_url: 'download-protein-map/',
             title: 'Reusable key recording peptide-to-protein mapping and any merge/rename decisions'},
            {key: 'group_correlation', icon: 'fa-chart-line', label: 'Between-Group Correlations (Group Averages)', corr: true,
             title: 'Across groups: agreement between each pair of group averages'},
            {key: 'protein_analysis', icon: 'fa-dna', label: 'Protein Analysis',
             title: 'Protein-level abundance and peptide-count distributions across groups'},
            {key: 'summed_peptide', icon: 'fa-chart-bar', label: 'Summed Peptide Results',
             title: 'Total abundance and unique-peptide counts summed per group'},
            {key: 'sequence_list', icon: 'fa-list', label: 'Sequence List',
             title: 'Unique peptide sequences detected in each group'},
        ];

        items.forEach(function(item) {
            var enabled = exports[item.key];
            var row = document.createElement('div');
            row.className = 'dt-export-row' + (enabled ? '' : ' disabled');

            var lbl = document.createElement('span');
            lbl.className = 'dt-export-label';
            lbl.title = item.title || '';
            lbl.innerHTML = '<i class="fas ' + item.icon + '"></i> ' + item.label;
            row.appendChild(lbl);

            if (enabled) {
                var actions = document.createElement('div');
                actions.className = 'dt-export-actions';

                if (item.corr) {
                    actions.appendChild(buildCorrOpts(item.key));
                }

                if (!item.no_view) {
                    var viewBtn = document.createElement('button');
                    viewBtn.className = 'dt-export-view-btn';
                    viewBtn.innerHTML = '<i class="fas fa-eye"></i> View';
                    (function(key, label) {
                        viewBtn.addEventListener('click', function() { viewExport(key, label); });
                    })(item.key, item.label);
                    actions.appendChild(viewBtn);
                }

                var dlBtn = document.createElement('button');
                dlBtn.className = 'dt-export-dl-btn';
                dlBtn.innerHTML = '<i class="fas fa-download"></i> Download';
                if (item.download_url) {
                    (function(url) {
                        dlBtn.addEventListener('click', function() {
                            var iframe = document.createElement('iframe');
                            iframe.style.display = 'none';
                            iframe.src = url;
                            document.body.appendChild(iframe);
                            setTimeout(function() { document.body.removeChild(iframe); }, 10000);
                        });
                    })(item.download_url);
                } else {
                    (function(key) {
                        dlBtn.addEventListener('click', function() { downloadExport(key); });
                    })(item.key);
                }
                actions.appendChild(dlBtn);

                row.appendChild(actions);
            }
            container.appendChild(row);
        });
    }

    // The three correlation exports each carry their own Correlation Method +
    // Log10/Log2 controls (built inline per row in renderExportButtons), making
    // it explicit that the log transform is scoped to that correlation file and
    // is not applied to any other export.
    var CORR_EXPORTS = ['tech_rep_correlation', 'replicate_correlation', 'group_correlation'];

    // Build the compact inline controls for one correlation export row.
    function buildCorrOpts(key) {
        var wrap = document.createElement('div');
        wrap.className = 'dt-corr-opts';

        var sel = document.createElement('select');
        sel.id = 'ct--' + key;
        ['Pearson', 'Spearman'].forEach(function(m) {
            var o = document.createElement('option');
            o.value = m; o.textContent = m;
            sel.appendChild(o);
        });

        var logs = document.createElement('div');
        logs.className = 'dt-corr-logs';
        var l10 = _mkLogChk('l10--' + key, 'Log₁₀');
        var l2 = _mkLogChk('l2--' + key, 'Log₂');
        // Log10 / Log2 are mutually exclusive; "neither" is allowed.
        l10.input.addEventListener('change', function() { if (l10.input.checked) l2.input.checked = false; });
        l2.input.addEventListener('change', function() { if (l2.input.checked) l10.input.checked = false; });
        logs.appendChild(l10.label);
        logs.appendChild(l2.label);

        wrap.appendChild(sel);
        wrap.appendChild(logs);
        return wrap;
    }

    function _mkLogChk(id, text) {
        var label = document.createElement('label');
        var input = document.createElement('input');
        input.type = 'checkbox';
        input.id = id;
        label.appendChild(input);
        label.appendChild(document.createTextNode(' ' + text));
        return { label: label, input: input };
    }

    // Query fragment for a correlation export, read from its own inline controls.
    // Returns '' for non-correlation exports (they take no correlation params).
    function corrParams(key) {
        if (CORR_EXPORTS.indexOf(key) === -1) return '';
        var sel = document.getElementById('ct--' + key);
        var l2 = document.getElementById('l2--' + key);
        var l10 = document.getElementById('l10--' + key);
        if (!sel) return '';
        var base = l2 && l2.checked ? '2' : (l10 && l10.checked ? '10' : '');
        return '&correlation_type=' + encodeURIComponent(sel.value) +
               '&log_transform=' + (base !== '') + '&log_base=' + base;
    }

    function downloadExport(type) {
        var url = 'download/' + type + '/?_=' + Date.now() + corrParams(type);
        var iframe = document.createElement('iframe');
        iframe.style.display = 'none';
        iframe.src = url;
        document.body.appendChild(iframe);
        setTimeout(function() { document.body.removeChild(iframe); }, 30000);
    }

    document.getElementById('download-all-btn').addEventListener('click', function() {
        // Each correlation export's controls are passed namespaced by key so the
        // bundle honours the per-row selections.
        var url = 'download-all/?_=' + Date.now();
        CORR_EXPORTS.forEach(function(key) {
            var sel = document.getElementById('ct--' + key);
            if (!sel) return;
            var l2 = document.getElementById('l2--' + key);
            var l10 = document.getElementById('l10--' + key);
            var base = l2 && l2.checked ? '2' : (l10 && l10.checked ? '10' : '');
            url += '&' + key + '_correlation_type=' + encodeURIComponent(sel.value) +
                   '&' + key + '_log_transform=' + (base !== '') +
                   '&' + key + '_log_base=' + base;
        });
        var iframe = document.createElement('iframe');
        iframe.style.display = 'none';
        iframe.src = url;
        document.body.appendChild(iframe);
        setTimeout(function() { document.body.removeChild(iframe); }, 60000);
    });

    function viewExport(type, label) {
        var url = 'view/' + type + '/?_=' + Date.now() + corrParams(type);

        var panel = document.getElementById('dt-viewer-panel');
        var content = document.getElementById('dt-viewer-content');
        var wrap = document.getElementById('dt-viewer-wrap');
        document.getElementById('dt-viewer-title').textContent = label;
        document.getElementById('dt-viewer-tabs').innerHTML = '';
        document.getElementById('dt-viewer-info').textContent = '';
        content.innerHTML = '<div style="padding: 20px; text-align: center;"><span class="dt-spinner"></span> Loading...</div>';
        panel.classList.remove('hidden');
        panel.scrollIntoView({behavior: 'smooth', block: 'nearest'});

        ajax('GET', url, null, function(resp) {
            currentViewerSheets = resp.sheets || [];
            renderViewerSheet(0);
        }, function(msg) {
            content.innerHTML = '<div class="dt-alert dt-alert-danger" style="margin: 8px;">' + msg + '</div>';
        });
    }

    function renderViewerSheet(idx) {
        var sheets = currentViewerSheets;
        var tabsEl = document.getElementById('dt-viewer-tabs');
        var content = document.getElementById('dt-viewer-content');
        var info = document.getElementById('dt-viewer-info');

        // Render tabs if multiple sheets
        tabsEl.innerHTML = '';
        if (sheets.length > 1) {
            sheets.forEach(function(sheet, i) {
                var tab = document.createElement('div');
                tab.className = 'dt-viewer-tab' + (i === idx ? ' active' : '');
                tab.textContent = sheet.name;
                (function(i) {
                    tab.addEventListener('click', function() { renderViewerSheet(i); });
                })(i);
                tabsEl.appendChild(tab);
            });
        }

        var sheet = sheets[idx];
        if (!sheet) {
            content.innerHTML = '<div class="dt-alert dt-alert-info" style="margin: 8px;">No data</div>';
            return;
        }

        // Build table HTML
        var html = '<table class="dt-viewer-table"><thead><tr>';
        sheet.columns.forEach(function(col) {
            html += '<th>' + escHtml(String(col)) + '</th>';
        });
        html += '</tr></thead><tbody>';
        sheet.rows.forEach(function(row) {
            html += '<tr>';
            row.forEach(function(cell) {
                var val = cell === null || cell === undefined ? '' : String(cell);
                html += '<td title="' + escHtml(val) + '">' + escHtml(val.length > 80 ? val.substring(0, 80) + '…' : val) + '</td>';
            });
            html += '</tr>';
        });
        html += '</tbody></table>';
        content.innerHTML = html;

        var infoText = sheet.total_rows + ' row' + (sheet.total_rows !== 1 ? 's' : '');
        if (sheet.truncated) {
            infoText += ' — showing first ' + sheet.rows.length + ' rows';
        }
        infoText += ', ' + sheet.columns.length + ' column' + (sheet.columns.length !== 1 ? 's' : '');
        info.textContent = infoText;
    }

    function escHtml(str) {
        return String(str).replace(/&/g, '&amp;').replace(/</g, '&lt;').replace(/>/g, '&gt;').replace(/"/g, '&quot;');
    }

    document.getElementById('dt-viewer-close-btn').addEventListener('click', function() {
        document.getElementById('dt-viewer-panel').classList.add('hidden');
    });

    // -----------------------------------------------------------------------
    // Visualization transfer buttons
    // -----------------------------------------------------------------------

    function transferToViz(transferUrl, btnId, btnLabel, iconClass, redirectUrl, storageKey) {
        var btn = document.getElementById(btnId);
        btn.disabled = true;
        btn.innerHTML = '<span class="dt-spinner"></span> Transferring...';
        document.getElementById('viz-transfer-status').innerHTML = '';
        ajax('POST', transferUrl, null, function(resp) {
            try { sessionStorage.setItem(storageKey, JSON.stringify(resp)); } catch(e) {}
            window.location.href = redirectUrl;
        }, function(msg) {
            btn.disabled = false;
            btn.innerHTML = '<i class="fas ' + iconClass + '"></i> ' + btnLabel;
            document.getElementById('viz-transfer-status').innerHTML =
                '<div class="dt-alert dt-alert-danger">' + escHtml(msg) + '</div>';
        });
    }

    document.getElementById('open-in-da-btn').addEventListener('click', function() {
        transferToViz('/data_analysis/transfer-from-dt/', 'open-in-da-btn',
            'Open in Data Analysis', 'fa-chart-bar',
            '/data_analysis/?from_dt=1', 'da_from_dt');
    });
    document.getElementById('open-in-hm-btn').addEventListener('click', function() {
        transferToViz('/heatmap/transfer-from-dt/', 'open-in-hm-btn',
            'Open in Heatmap Visualization', 'fa-th',
            '/heatmap/?from_dt=1', 'hm_from_dt');
    });

    document.getElementById('start-over-btn').addEventListener('click', function() {
        if (!confirm('This will clear all data. Continue?')) return;
        ajax('POST', 'cleanup/', null, function() {
            location.reload();
        });
    });

    // Restore a dropzone's visual filename display (mirrors the internal showFile closure)
    function restoreDropzoneName(inputId, fileName) {
        if (!fileName) return;
        var input = document.getElementById(inputId);
        if (!input) return;
        var zone = input.closest('.dt-dropzone');
        if (!zone) return;
        zone.querySelector('.drop-text').classList.add('hidden');
        var el = zone.querySelector('.file-name');
        el.textContent = fileName;
        el.classList.remove('hidden');
    }

    // -----------------------------------------------------------------------
    // Resume banner — check for existing session on page load
    // -----------------------------------------------------------------------
    ajax('GET', 'session-state/', null, function(resp) {
        if (!resp.has_session || !resp.has_data) return;

        // Restore uploaded file names on the Step 1 dropzones
        var fn = resp.file_names || {};
        var savedFiles = resp.has_saved_files || {};

        if (fn.merge_files && fn.merge_files.length) {
            // Merge mode — show individual filenames in the peptidomic dropzone
            var mergeInput = document.getElementById('peptidomic_file');
            if (mergeInput) {
                var mergeZone = mergeInput.closest('.dt-dropzone');
                if (mergeZone) {
                    mergeZone.querySelector('.drop-text').classList.add('hidden');
                    var ml = mergeZone.querySelector('.merge-file-list');
                    if (ml) {
                        ml.innerHTML = '';
                        fn.merge_files.forEach(function(name) {
                            var item = document.createElement('div');
                            item.className = 'merge-file-item';
                            item.textContent = name;
                            ml.appendChild(item);
                        });
                        ml.classList.remove('hidden');
                    }
                }
            }
        } else {
            restoreDropzoneName('peptidomic_file', fn.peptidomic_file);
        }
        restoreDropzoneName('functional_file', fn.functional_file);
        restoreDropzoneName('fasta_file',      fn.fasta_file);

        // If saved file bytes are available, also reload them into the inputs
        // so the form can be resubmitted without the user re-selecting files.
        function reloadSavedFile(fileKey, inputId) {
            if (!savedFiles[fileKey]) return;
            fetch('reload-saved-file/' + fileKey + '/')
                .then(function(r) {
                    if (!r.ok) throw new Error('HTTP ' + r.status);
                    var origName = r.headers.get('X-Original-Filename') || fileKey;
                    return r.blob().then(function(blob) {
                        return {blob: blob, name: origName};
                    });
                })
                .then(function(obj) {
                    var file = new File([obj.blob], obj.name,
                        {type: obj.blob.type || 'application/octet-stream'});
                    var dt = new DataTransfer();
                    dt.items.add(file);
                    var input = document.getElementById(inputId);
                    if (input) {
                        input.files = dt.files;
                        // Don't dispatch 'change' here — we don't want to
                        // re-trigger the dropzone visual since names are already shown.
                    }
                })
                .catch(function() { /* silent — user can still pick manually */ });
        }

        reloadSavedFile('peptidomic_file', 'peptidomic_file');
        reloadSavedFile('functional_file', 'functional_file');
        reloadSavedFile('fasta_file',      'fasta_file');

        var banner = document.getElementById('dt-resume-banner');
        var btns = document.getElementById('dt-resume-btns');

        function makeResumeBtn(label, icon, onClick, extraClass) {
            var btn = document.createElement('button');
            btn.className = 'dt-resume-btn' + (extraClass ? ' ' + extraClass : '');
            btn.innerHTML = '<i class="fas ' + icon + '"></i> ' + label;
            btn.addEventListener('click', function() {
                banner.classList.add('hidden');
                onClick();
            });
            btns.appendChild(btn);
        }

        makeResumeBtn('Upload Data', 'fa-upload', function() {
            goToStep(1);
            if (resp.has_data) {
                document.getElementById('upload-summary').innerHTML =
                    '<div class="dt-alert dt-alert-info">' +
                    '<i class="fas fa-check-circle"></i> Peptidomic data already uploaded. ' +
                    'Re-upload above to replace it, or continue with the functional annotation below.</div>';
                document.getElementById('functional-card').classList.remove('hidden');
                enableFunctionalActions();
            }
        });
        makeResumeBtn('Study Variables', 'fa-layer-group', function() {
            goToStep(2); loadStep2();
        });
        makeResumeBtn('Protein Mapping', 'fa-dna', function() {
            goToStep(3); loadStep3();
        });
        if (resp.has_processed) {
            makeResumeBtn('Process &amp; Export', 'fa-cogs', function() {
                goToStep(4);
            });
        }
        makeResumeBtn('Reset', 'fa-trash-alt', function() {
            if (!confirm('This will clear all uploaded data and reset the session. Continue?')) {
                banner.classList.remove('hidden');
                return;
            }
            ajax('POST', 'cleanup/', null, function() {
                location.reload();
            });
        }, 'dt-resume-btn-reset');

        banner.classList.remove('hidden');
    });

})();
