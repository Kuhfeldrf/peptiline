# Splitting PeptiLine and MBPDB into Separate Apps/Containers

Status: draft plan, not yet executed
Owner: Russell Kuhfeld
Scope: split the single `mbpdbcontainer` Azure Container App (repo: `~/mbpdb`, one Django
project named `peptide`) into two independently deployed apps — MBPDB (existing) and
PeptiLine (this repo, `~/notebooks/peptiline`, branch `django_apps`) — each in its own
Docker container/Container App, with the MBPDB reference database (protein/peptide/
function/reference tables + BLAST index) **replicated** into PeptiLine rather than shared
live.

---

## 0. RESOLVED 2026-08-25 — PeptiLine is live and healthy on Azure

`peptilinecontainer` is now serving correctly at
`https://peptilinecontainer.lemonisland-71b15397.westus3.azurecontainerapps.io/` — revision
`peptilinecontainer--fix3`, `Healthy`/`Running`, 100% traffic. All 9 routes
(`/health/ / /data_transformation/ /data_analysis/ /heatmap/ /about_us/ /pepex/
/peptiline/supplementals/ /admin/`) verified 200/302.

**What actually happened, for the record** (the two code fixes below turned out to be right
the first time — the real blocker after that was infrastructure, not code):
1. The two code fixes described in the original "RESUME HERE" note below
   (`peptiline/middleware.py` health-check bypass, `peptiline/storage.py` lenient static
   storage) were correct and are what's running now.
2. **The actual reason redeploys kept failing wasn't the code — it was a masked push
   failure.** `docker push ... | tail -40` reports the exit code of `tail`, not of `docker
   push`; several pushes silently failed (`TLS handshake timeout`, `use of closed network
   connection` — this session's network was flaky) while still being reported as "completed
   (exit code 0)". Azure was correctly redeploying every time; it was just re-pulling the
   *same broken old image* each time because the registry never actually got the new one.
   Caught by comparing the pushed digest against what `docker manifest inspect` showed
   after a "successful" push — they matched the *original* image's digest, proving nothing
   had changed. Fixed by re-running the push with the exit code captured explicitly
   (`docker push ... | tee log; echo "PUSH_EXIT_CODE=${PIPESTATUS[0]}"`) until it actually
   returned 0, then deploying by immutable digest (`@sha256:...`) rather than `:latest`, to
   remove any ambiguity about which image content was actually being deployed.
3. Lesson for next time: **never trust a piped command's reported exit code** — pipe to a
   real destination and check `${PIPESTATUS[0]}`, or avoid the pipe for anything whose
   success/failure matters (pushes, deploys, migrations).

Not yet done: DNS/custom domain for this new hostname, Phase 5 (MBPDB redirect cutover),
Phase 2 (DB replica — MBPDB search still intentionally unavailable, correct for now).

---

## 0c. UPDATE 2026-08-25 — MBPDB-formatting removal, repo pushed to GitHub

Per user request: stripped the app down to just the landing page + the three PeptiLine
modules, and dropped MBPDB's shared nav/footer chrome in favor of PeptiLine's own design
system on every page (not just the landing page, which is all it covered before).

**Code changes:**
- `templates/peptide/base.html`: `<head>` now loads the `design-system.css`/Inter font/
  Font Awesome assets unconditionally (previously only `peptiline_landing.html` loaded
  them, so module pages never got the modern look). The old Bootstrap `menubar-nav` block
  is replaced with the same `pl-nav` markup the landing page used, trimmed to **About /
  Tools (Data Transformation, Data Analysis, Heatmap Visualization only) / Contact** —
  Home, Search, Help, and PepEx all dropped from the nav. Footer simplified to drop the
  unused `latest_peptides` loop (MBPDB-only context this app's views never populate).
- `templates/peptide/peptiline_landing.html`: its own duplicate head/nav/footer block
  overrides deleted entirely — now inherits the same chrome as every other page from
  `base.html`, so there's exactly one nav/footer implementation instead of two drifting
  copies.
- Removed the MBPDB-only `pepex`/`peptide_search` routes, view (`pepex_unavailable`), and
  template — this repo is now just the landing page + the three modules, nothing else.
- `static/peptide/hero.png` copied over from `~/mbpdb` (missing since the original vendoring
  pass, same class of bug as `MBPDB_Help.pdf` before it) — the landing page hero image now
  actually loads instead of silently degrading via the lenient static storage fallback.
- Contact email (`Contact-MBPDB@oregonstate.edu`) was left as-is — not asked to change, and
  a real replacement wasn't specified. Worth revisiting once PeptiLine has its own contact
  channel, since it's the one remaining visible MBPDB reference in the nav.

**Repo/deployment:**
- Committed and pushed to `origin/django_apps`.
- `main` was stale (last real content was the "Vendor v2" commit, `275fc15` merge); confirmed
  it was a strict ancestor of `django_apps` (clean fast-forward relationship), merged
  `django_apps` into local `main`, and pushed — `main` is now current at commit `2047fe5`.
  `origin/main` is `git@github.com-kuhfeldrf:Kuhfeldrf/peptiline.git` (the real GitHub repo
  — `.github/workflows/deploy.yml`'s `on: push: branches: [main]` trigger now points at
  actual content, though it still won't run successfully until the GitHub secrets it
  references are added, see Phase 4).
- Rebuilt the Docker image with these changes, pushed `mbpdb/peptiline:latest`
  (digest `sha256:8deb175e21068b5b638c3f144b02159511e8289609192afb1f3d3837cd9357c2`), and
  redeployed `peptilinecontainer` by digest (`--revision-suffix nav1`). Verified live:
  trimmed nav present (`About`/`Tools`/`Contact`, no `Home`/`Search`/`Help`), hero image
  200s, `/pepex/` 404s, all module routes still 200/302.

---

## 0b. Original "RESUME HERE" note (2026-08-24, superseded by the above — kept for context)

**Azure is currently in a broken state — read this before doing anything else.**

`peptilinecontainer` (Container App, resource group `COH_MBPDB_RG`) exists and has 100%
traffic on revision `peptilinecontainer--0000003`, which is **`ActivationFailed`** — it 500s/
was never healthy. There is an older, healthy-but-inactive revision `--dsgs7t7` at 0%
traffic. **The app is effectively down right now**, not serving real traffic to anyone
external yet (no DNS/custom domain points at it), but it's sitting broken. If anyone asks
"is peptilinecontainer up," the answer today is no.

**Why:** the image deployed to that revision (`mbpdb/peptiline:latest`, pushed earlier in
this session) has two now-fixed-locally-but-not-yet-repushed bugs:
1. Azure Container Apps' health probes hit the container over an internal cluster IP
   (e.g. `100.100.2.79`), not the public FQDN, so Django's `ALLOWED_HOSTS` check rejected
   every probe with `DisallowedHost` → container never passed its startup probe.
   **Fixed locally**: `peptiline/middleware.py` (`HealthCheckMiddleware`, first in
   `MIDDLEWARE`) answers `/health/` before the host-check middleware ever runs.
2. `django.test.Client`-independent bug: WhiteNoise's `CompressedManifestStaticFilesStorage`
   is strict by default — any `{% static %}` reference to a file that doesn't physically
   exist raises `ValueError` and 500s the *entire* page. This repo's "vendor v2" commit
   never copied over ~50 marketing/screenshot/demo/supplemental assets that
   `base.html`/`peptiline_landing.html`/`peptiline_supplementals.html` reference (hero
   image, module screenshots, example data files, demo plot HTML/PNGs, the whole
   `publications/zukaitis_2026/supplementals/` tree — grep the missing-file list via
   `grep -rhoE "\{% static '[^']+' %\}" templates/` and check each against `static/`).
   **Fixed locally**: `peptiline/storage.py` (`LenientManifestStaticFilesStorage`) falls back
   to the plain unhashed filename instead of raising — pages render with a broken
   link/image instead of a 500. (One file, `static/peptide/MBPDB_Help.pdf`, was copied over
   for real since it's small and clearly available; the rest are still genuinely missing —
   this only stops them from taking the whole page down.)

Both fixes were verified locally: rebuilt `peptiline:local`, ran it with
`DEBUG=False` (matching production) and a spoofed internal-IP `Host` header, and all of
`/ /health/ /data_transformation/ /data_analysis/ /heatmap/ /about_us/ /pepex/
/peptiline/supplementals/ /admin/` returned 200/302 as expected.

**Exact next steps to un-break this:**
1. `cd ~/notebooks/peptiline && docker build -t peptiline:local .` (picks up
   `peptiline/middleware.py` + `peptiline/storage.py`, already committed to disk, not yet
   pushed to Docker Hub).
2. `docker tag peptiline:local mbpdb/peptiline:latest && docker push mbpdb/peptiline:latest`
   (Docker Hub login already done this session, in the shell's own `~/.docker/config.json` —
   should still be valid; re-run `docker login -u mbpdb` yourself if it's expired).
3. Force a new revision so Azure actually re-pulls the `:latest` tag (same tag won't
   auto-redeploy):
   ```bash
   az containerapp update --name peptilinecontainer --resource-group COH_MBPDB_RG \
     --image docker.io/mbpdb/peptiline:latest --revision-suffix fix2
   ```
4. Watch it come up:
   ```bash
   az containerapp revision list --name peptilinecontainer --resource-group COH_MBPDB_RG -o json
   FQDN=$(az containerapp show --name peptilinecontainer --resource-group COH_MBPDB_RG \
     --query "properties.configuration.ingress.fqdn" -o tsv)
   curl -s -o /dev/null -w "%{http_code}\n" "https://$FQDN/health/"
   ```
5. Once healthy, smoke-test the same 9 routes listed above against `https://$FQDN/...`.
6. **Known gotcha already hit once**: `az containerapp update --set-env-vars` *replaces* the
   entire env var list, it does not merge — if you need to add/change one env var, re-supply
   all of them together (`DJANGO_SETTINGS_MODULE`, `WEBSITES_PORT`,
   `WEBSITES_CONTAINER_START_TIME_LIMIT`, `DJANGO_SECRET_KEY`, `DJANGO_ALLOWED_HOSTS`) in one
   call. This bit us once already mid-session (wiped `DJANGO_SECRET_KEY`, crashed the app)
   before being caught and fixed.
7. `git status` in this repo currently shows a large set of modified/new files, all
   uncommitted (nothing has been committed this session — only asked-for by the user when
   they explicitly request a commit). Review before committing.

**Not yet done at all** (unstarted, no partial state to worry about): Phase 2 (DB replica —
MBPDB search still raises a clear "not available" error, this is expected/by design for now),
DNS/custom domain for PeptiLine's new hostname, GitHub Actions secrets (workflow file exists
at `.github/workflows/deploy.yml` but references secrets that don't exist yet in the repo, so
don't expect CI pushes to work), Phase 5 (MBPDB redirect cutover), Phase 6 (README/docs
updates).

---

## 1. Current-state audit (as of 2026-08-24)

**Single Azure Container App**: `mbpdbcontainer`, resource group `COH_MBPDB_RG`, deployed
from `~/mbpdb` via `.github/workflows/mbpdbcontainer-AutoDeployTrigger-*.yml` on every push
to `main`. One Docker image (`mbpdb/mbpdb:latest`), one Django project (`peptide`), running
nginx → gunicorn (127.0.0.1:8001) → Django, plus a local Redis and one Celery worker, all in
the same container (`start.sh`).

**One Django project, one app namespace.** `INSTALLED_APPS` has a single app, `peptide`,
which is really the whole codebase: MBPDB's own views (`peptide_search`, `pepex_tool`,
BLAST-backed search, admin) live directly under `peptide/`, and PeptiLine's three modules
(`data_transformation`, `data_analysis`, `heatmap_viz`) are **sub-packages of `peptide`**,
mounted in the same `urls.py` under `/data_transformation/`, `/data_analysis/`, `/heatmap/`,
with a landing page at `/peptiline/`. There is no separation at the Python-package,
settings, or deployment level today — it's one monolith serving two products from one
origin (`mbpdb.nws.oregonstate.edu`).

**Database is SQLite, baked into the image, not a client-server DB.**
```python
DATABASES = {'default': {'ENGINE': 'django.db.backends.sqlite3', 'NAME': .../db.sqlite3}}
```
The Dockerfile does `touch db.sqlite3` at build time and chowns it — there is no evidence
of a mounted persistent volume in what's in the repo (no Azure Files/volume mount found in
the workflow or probes yaml). This needs to be confirmed directly against the live Container
App config (`az containerapp show`) before anything else — **if the running container's
`db.sqlite3` is not actually persisted between deployments/restarts today, "replication" is
moot until persistence itself is fixed.** This is the single most important fact to verify
first and is flagged as an open question in §7.

**Core MBPDB schema** (`peptide/models.py`): `ProteinInfo`, `ProteinVariant`, `PeptideInfo`,
`Function`, `Reference`, `Submission`, `Counter`, `GitHubActions`. `Submission`/`Counter`/
`GitHubActions` look like MBPDB-admin/site-metadata tables that PeptiLine has no business
touching.

**Exactly one live coupling point between PeptiLine's code and MBPDB's database:**
`data_transformation/services/blast_search.py`:
```python
from peptide.models import PeptideInfo, ProteinInfo, Function, Reference
```
This is where a user's uploaded peptides get matched against MBPDB (exact match or BLAST
homology search) to pull in bioactive-function annotations. Every other cross-import found
(`peptide.utils.uniprot_client`, `peptide.utils.lazy_import`, `peptide.toolbox`,
`peptide.data_analysis.services.stats` from `heatmap_viz`) is PeptiLine-internal code that
merely lives under the `peptide.` namespace today — not a dependency on MBPDB's data. **This
means the "shared SQL-like database" the two apps need is narrow**: four read-mostly tables
(`PeptideInfo`, `ProteinInfo`, `Function`, `Reference`) plus the BLAST protein FASTA/index
built from them (`peptide/blast_db/`, built via `makeblastdb` in `blast_search.py`).

**A standalone PeptiLine repo already exists and is partway extracted**
(`~/notebooks/peptiline`, branch `django_apps`, 2 commits: legacy notebook import, then
"Vendor v2 Django app (data_transformation, data_analysis, heatmap_viz, utils) into
standalone PeptiLine repo"). It has `requirements.txt`, `templates/`, `static/`,
`tests/`, `examples/`, `supplementals/`, and the three app directories — **but no
`manage.py`, no Django project package (`settings.py`, `urls.py`, `wsgi.py`), no
`Dockerfile`, no CI workflow.** The README already documents a `python manage.py migrate
&& python manage.py runserver` workflow that doesn't work yet, because that scaffolding
hasn't been created. This is the first gap to close, independent of the database question.

**The standalone repo's README already anticipates the coupling problem** and currently
documents a *deliberately decoupled* answer: "MBPDB search will not run locally... upload
your own functional annotation table, or upload a pre-downloaded MBPDB TSV." That's a valid
fallback UX, but it's not what was asked for — the user wants live-feeling MBPDB search to
keep working in the split-out PeptiLine deployment, backed by a **replica** of MBPDB's data
rather than the shared live DB. The README/docs will need updating once replication is
wired up, since right now they undersell what the split version can do.

**Publication URL stability matters.** The current README advertises
`https://mbpdb.nws.oregonstate.edu/peptiline/` as *the* live application URL, and it's cited
in a manuscript under review at JPR (MS pr-2025-01102w) and will be archived on Zenodo with
that URL baked into the citation trail. **Whatever new hostname PeptiLine gets, the old
`/peptiline/` path must keep working (redirect, not 404) indefinitely.**

**Other shared/duplicated concerns found:**
- Celery + Redis: one broker/worker per container today; PeptiLine's `tasks.py` uses Celery
  for BLAST search, UniProt fetch, exports — it will need its own Redis+Celery in its own
  container.
- `uploads/temp` directory, `WORK_DIRECTORY` for BLAST scratch space — per-container local
  disk, fine to duplicate, not shared state.
- `CORS_ALLOWED_ORIGINS` / `CSRF_TRUSTED_ORIGINS` are already environment-driven
  (`DJANGO_CORS_ORIGINS`, `DJANGO_CSRF_ORIGINS`) — good, this generalizes to a second app
  without code changes, just new env values.
- Static files: nginx serves `/static/` from `static_files/` collected at build time —
  duplicate this setup per container, no shared state needed.
- `ncbi-blast+`, `redis-server`, `gosu`, etc. are installed at the OS level in the Dockerfile
  — the new PeptiLine image will need `ncbi-blast+` too (for the BLAST-backed search) and
  `redis-server` (for its own Celery broker), i.e. it can't just drop those dependencies.

---

## 2. Target architecture

Two independently deployable units, each its own git repo, Docker image, and Azure
Container App:

1. **MBPDB app** (existing `~/mbpdb` repo, existing `mbpdbcontainer`) — becomes the
   canonical/source-of-truth database owner. Minimal changes: strip PeptiLine's three app
   modules and the `/peptiline/`, `/data_transformation/`, `/data_analysis/`, `/heatmap/`
   URL includes out once PeptiLine is confirmed working standalone; keep `peptide_search`,
   admin, BLAST-search-against-MBPDB, and add a replication export mechanism (§3).
2. **PeptiLine app** (this repo) — gets real Django project scaffolding, its own
   Dockerfile/nginx/start.sh modeled on MBPDB's, its own Azure Container App, its own Redis
   + Celery, and a **local read replica** of the four MBPDB tables + BLAST index, refreshed
   on a schedule (§3), used by `blast_search.py` unchanged (it already only depends on the
   Django ORM, not on which physical DB backs it).

Both apps sit behind the same Azure resource group (`COH_MBPDB_RG`, or split if there's a
billing/ownership reason not surfaced here — flagged as a question) but as two separate
Container Apps with two separate ingress hostnames, e.g. `mbpdb.nws.oregonstate.edu`
(unchanged) and a new `peptiline.<domain>` (or a path-based reverse proxy in front of both —
see §7 open question on whether `/peptiline/` must stay under the *same* hostname to satisfy
the "don't break the citation URL" requirement, which pushes toward a proxy/redirect
approach rather than a bare new hostname).

---

## 3. Database replication strategy

The user asked to keep this simple for now ("just replicate the SQL-like database"), so the
recommendation is the smallest thing that works, not a distributed-database rebuild:

**Recommended: one-way, scheduled export/import, MBPDB → PeptiLine, of the four read-only
tables, plus a rebuilt BLAST index.**

- MBPDB adds a management command (e.g. `dumpreplica`) that serializes `ProteinInfo`,
  `ProteinVariant`, `PeptideInfo`, `Function`, `Reference` (Django's `dumpdata` is enough
  for this volume) to a single artifact — either pushed to Azure Blob Storage on a timer, or
  exposed via an authenticated endpoint PeptiLine can pull from.
- PeptiLine adds a matching `loadreplica` management command, run on a schedule (Celery
  beat, or Azure Container Apps' own cron/scheduled job feature) and once at container
  startup, that loads the dump into its own local SQLite tables and re-runs `makeblastdb`
  to rebuild its BLAST index from the refreshed `ProteinInfo`/`ProteinVariant` FASTA.
- Because SQLite is single-writer and file-based, PeptiLine's replica tables should live in
  a **separate SQLite file** from PeptiLine's own operational tables (Celery result
  backend if any, sessions, uploads bookkeeping) so a replica reload can't lock out or
  corrupt PeptiLine's own writes — use Django's multi-database support
  (`DATABASES['mbpdb_replica']` + a router) rather than mixing replica and app tables in one
  file.
- This keeps `blast_search.py` completely unchanged in PeptiLine's copy — it queries
  `PeptideInfo`/`ProteinInfo`/`Function`/`Reference` from the ORM exactly as it does today;
  only the underlying `DATABASES` alias/router changes.

**Explicitly deferred (call out as future work, not now):** two-way sync, live
cross-container queries, migrating to a real client/server DB (Postgres/MySQL) shared over
the network, or event-driven replication (webhooks/CDC). All of these are more correct
long-term (SQLite-over-file-copy is a known operational headache — see risks below) but are
more work than "replicate the database" as scoped for this pass. Worth a one-line note back
to the user that Postgres + a real read replica is the natural next step once both apps are
live and the manual/scheduled dump proves too coarse (e.g. new MBPDB peptides not showing up
in PeptiLine fast enough).

**Open question to confirm before building this:** how often does MBPDB's underlying data
actually change (new submissions approved, corrections)? That determines whether "refresh on
container start + nightly" is enough, or whether something tighter is needed.

---

## 4. Coupling points to resolve in code

1. `data_transformation/services/blast_search.py` — only needs the `DATABASES` alias/router
   change described above; no logic change.
2. Everything else under `data_transformation/`, `data_analysis/`, `heatmap_viz/`, `utils/`
   currently does `from peptide.X import ...` — once PeptiLine is its own Django project,
   these need a namespace decision: keep the top-level package named `peptide` in the new
   repo too (least code churn, but confusing to have two unrelated repos both ship a
   `peptide` Python package) vs. rename to something like `peptiline` (clearer, but touches
   every import across `data_transformation/`, `data_analysis/`, `heatmap_viz/`,
   `heatmap_renderer.py`'s cross-import of `data_analysis.services.stats`, and templates
   that reference `{% static %}`/app namespaces). **Recommend renaming to avoid the
   two-different-`peptide`-packages footgun**, but this is a judgment call worth confirming
   since it's the single largest mechanical diff in this whole project.
3. `heatmap_viz/services/heatmap_renderer.py` imports from
   `peptide.data_analysis.services.stats` — an intra-PeptiLine cross-module import, fine to
   keep once the namespace is settled, just needs the same rename applied consistently.
4. Templates and any hardcoded `/data_transformation/`, `/heatmap/`, `/peptiline/` URLs —
   audit for `{% url %}` usage vs. hardcoded paths; hardcoded absolute paths will break if
   PeptiLine's URL prefixes change on the new host (they may not need to change at all if a
   redirect/proxy preserves `/peptiline/...` — see §7).
5. `settings.WORK_DIRECTORY`, `uploads/temp`, static file collection, `django_celery_progress`
   — all need to exist as first-class settings in PeptiLine's new `settings.py`, currently
   inherited implicitly from MBPDB's.

---

## 5. Infra / deployment changes

- New git-hosted CI: a `deploy.yml` in this repo modeled directly on
  `~/mbpdb/.github/workflows/mbpdbcontainer-AutoDeployTrigger-*.yml`, with its own Azure
  service principal credentials (new GitHub secret, e.g. `PEPTILINECONTAINER_AZURE_CREDENTIALS`),
  its own registry username/password, and its own `containerAppName`
  (e.g. `peptilinecontainer`).
- New Dockerfile for PeptiLine, based on MBPDB's: same OS deps
  (`ncbi-blast+`, `redis-server`, `gosu`, `nginx`, `dos2unix`, `sqlite3`, `build-essential`),
  same `celery_user` pattern, same `collectstatic` step, same `start.sh` shape (nginx +
  gunicorn + celery, plus the new replica-load step on boot).
- New `container-app-probes.yaml` — reuse the MBPDB pattern (`/health/` httpGet on all three
  probe types, not bare TCP) with a PeptiLine `health_check` view.
- New environment variables/secrets mirroring MBPDB's: `DJANGO_SECRET_KEY` (own, don't
  reuse MBPDB's), `DJANGO_ALLOWED_HOSTS`, `DJANGO_CORS_ORIGINS`, `DJANGO_CSRF_ORIGINS`,
  superuser vars if PeptiLine gets its own admin, plus whatever the replication pull step
  needs (a shared secret or SAS token to fetch MBPDB's dump — do **not** give PeptiLine
  direct network access to MBPDB's Django admin/DB).
- Decide Container App scaling/resources for PeptiLine independently — it runs the same
  BLAST/Celery workload pattern as MBPDB today, so start from the same CPU/memory allocation
  and tune from there.

---

## 6. Rollout / URL continuity plan

Given the citation-URL constraint (§1), the recommended sequence is:

1. Stand up PeptiLine as its own container at a **new** hostname/path first, fully working
   (including the DB replica), verified independently — no traffic cut over yet.
2. On MBPDB's side, replace the in-process `/peptiline/`, `/data_transformation/`,
   `/data_analysis/`, `/heatmap/` Django views with **HTTP redirects** (301) to the
   equivalent paths on the new PeptiLine host, rather than deleting them outright. This
   preserves every existing external link/citation without needing a shared-hostname proxy.
3. Only after redirects are confirmed working end-to-end should the old in-process
   PeptiLine app code be deleted from the MBPDB repo.
4. Update this repo's README (`docs/INSTALL.md`, `docs/REPRODUCIBILITY.md`) to reflect the
   new hosted URL and the fact that live MBPDB search now works locally-hosted too (once
   §3 is done) — the current "MBPDB search will not run locally" language becomes
   inaccurate for the new deployment and should be corrected.

---

## 7. Open questions / risks to confirm before or during implementation

- **Is `db.sqlite3` actually persisted today**, or does every MBPDB container restart start
  from an empty/build-time DB? This determines whether "replication" is even meaningful yet
  and should be checked first (`az containerapp show`, or check for a volume mount / init
  container not visible in this repo).
- SQLite has no built-in replication and is a poor fit for "keep two independent write
  paths in sync" long-term; this plan deliberately keeps PeptiLine's replica **read-only**
  and MBPDB as sole writer to sidestep that, but if PeptiLine ever needs to write back to
  MBPDB data (e.g. user-submitted corrections), this design breaks and Postgres becomes
  necessary sooner.
- Package-name collision (`peptide` in both repos) — needs an explicit decision (§4.2).
- Whether `COH_MBPDB_RG` / Azure subscription has room/budget for a second Container App,
  and who owns provisioning the new Azure service principal + registry credentials — this is
  an access/ops question outside the codebase.
- DNS/custom domain + TLS cert provisioning for a new PeptiLine hostname (if going the
  new-hostname route rather than a shared reverse proxy) is an Azure/DNS-admin task not
  visible from the repos and should be scheduled early since cert issuance can take time.
- Whether MBPDB's `Submission`/`Counter`/`GitHubActions` tables have any PeptiLine relevance
  — audit confirms no code reference from PeptiLine's three modules today, but worth a second
  look once the replica dump/load commands are actually written, so they're excluded
  deliberately rather than by omission.
- Test coverage: both repos have `tests/` with overlapping filenames (e.g.
  `test_heatmap_differential.py` exists in both `~/mbpdb/tests/` and this repo's `tests/`,
  and again nested at `~/mbpdb/include/peptide/tests/`) — reconcile which is authoritative
  before the split so PeptiLine's CI runs the real, current suite rather than a stale copy.

---

## 8. Detailed step-by-step implementation checklist

### Phase 0 — Verify assumptions
- [ ] Confirm whether MBPDB's `db.sqlite3` is persisted across container restarts/deploys
      today (check Azure Container App volume config directly).
- [x] Diff `~/mbpdb/include/peptide/peptide/{data_transformation,data_analysis,heatmap_viz,utils}`
      against `~/notebooks/peptiline`'s copies to find drift since the "Vendor v2" commit,
      and decide which is authoritative going forward. *(2026-08-24: the vendored copies in
      this repo are authoritative going forward — they already had `from peptide.X` imports
      rewritten to bare top-level imports (`utils.X`, `data_transformation.services.X`, etc.)
      except the one deliberate MBPDB coupling point, `blast_search.py`.)*
- [x] Decide the package-namespace question (`peptide` vs. rename) and the
      shared-hostname-vs-new-hostname question (§7) — both are prerequisites for concrete
      code changes below. *(2026-08-24: kept the Django **project** package named `peptiline`
      (new), but did NOT rename the `templates/peptide/` and `static/peptide/` namespace
      folders already vendored into this repo — app code already imports bare module names,
      not `peptide.X`, so there was no `peptide`-vs-`peptiline` collision to resolve; only
      `blast_search.py`'s single MBPDB-model import needed handling, see Phase 1 below.
      Shared-hostname-vs-new-hostname is still open, deferred to Phase 4/6.)*

### Phase 1 — Scaffold PeptiLine as a runnable standalone Django project — **done 2026-08-24**
- [x] Add `manage.py`, project package (`settings.py`, `urls.py`, `wsgi.py`, `celery.py`,
      `views.py`, `management/commands/bootstrap.py`) — see `peptiline/`.
- [x] Point `INSTALLED_APPS` at PeptiLine's three apps (`data_transformation`,
      `data_analysis`, `heatmap_viz`) plus `django_celery_progress` and `peptiline` itself
      (needed so the `bootstrap` management command is discovered — omitting it silently
      broke `manage.py bootstrap` inside the container on the first Docker test run).
- [x] `blast_search.py`'s module-level `from peptide.models import ...` (the one real MBPDB
      coupling point, already flagged in-repo as "KNOWN BLOCKER... L-7") now imports lazily
      inside a try/except, with a `_require_mbpdb_models()` guard raising a clear
      `RuntimeError` only when MBPDB search is actually invoked. Before this fix the whole
      Data Transformation dashboard 500'd at import time; now everything else works and only
      MBPDB search itself is unavailable until Phase 2 (replica) lands.
- [x] Added two stub views/templates (`about_us`, `pepex_unavailable`) and duplicate
      `peptiline_landing`/`peptide_search` URL names — `base.html`/`peptiline_landing.html`
      already reference these via `{% url %}` and would otherwise raise `NoReverseMatch`.
- [x] `python manage.py check`, `migrate --run-syncdb`, and `runserver` all confirmed working
      locally against a fresh SQLite DB with no MBPDB dependency; all app routes (`/`,
      `/health/`, `/data_transformation/`, `/data_analysis/`, `/heatmap/`, `/about_us/`,
      `/pepex/`, `/peptiline/supplementals/`, `/admin/`) return 200/302 as expected.
- [x] Ran the existing `tests/` suite against the new scaffolding: found and fixed three
      **pre-existing bugs in the vendored repo**, unrelated to the rename itself —
      `test_data_transformation.py`/`test_case_insensitive_mapping.py` pointed at a stale
      `peptide/static/peptide/examples/...` path from the old mbpdb layout (fixed to
      `examples/legacy_v1/...`); `test_data_transformation.py`'s `ARCHIVE` constant pointed
      at `peptide/notebooks/archive/`, whose 5 fixture CSVs were never vendored into this
      repo at all (copied over from `~/mbpdb/include/peptide/peptide/notebooks/archive/`
      into `examples/legacy_v1/archive/`); `matplotlib` was missing from `requirements.txt`
      even though `heatmap_renderer.py` genuinely imports it (added). Result: 239→314 of 318
      tests passing (the `not_in_use_examples` docstring in
      `test_data_transformation.py` also implied `PEAKS_example.csv`/`spectronaut.tsv`
      should live one level up from `needs_transformation/`; copied them up to match). The
      remaining 2 failures (`TestTransferFromDtFasta`) need a real `django.test.Client`,
      which needs `sessions`/`contenttypes` in a full `INSTALLED_APPS` + DB — the
      intentionally minimal `conftest.py` (`INSTALLED_APPS=[]`) doesn't provide that; a
      pre-existing test-harness gap, not a split-related regression, left as known future
      work rather than rebuilding the harness now.

### Phase 2 — Database replica plumbing — **done 2026-08-25** (simple-copy version)

Implemented per explicit user direction to keep this simple: a literal SQLite table copy
via `ATTACH DATABASE`, not the JSON dumpdata/loaddata round-trip originally drafted below
(tried first, worked, replaced anyway per "we can just make a simple copy of the DB SQLite
database"). Two independent DB files/deployments are accepted for now — no live sync between
MBPDB and PeptiLine, "figure out later" per user.

- `mbpdb_replica/` (new Django app): `ProteinInfo`, `ProteinVariant`, `PeptideInfo`,
  `Function`, `Reference` — field-for-field copies of MBPDB's `peptide.models`, managed
  models with their own migration.
- `peptiline/db_router.py` (`MBPDBReplicaRouter`): pins every `mbpdb_replica` model to the
  `mbpdb_replica` DB alias (`mbpdb_replica.sqlite3`, physically separate from PeptiLine's own
  `db.sqlite3`) for both reads and migrations; everything else stays on `default`.
- `mbpdb_replica/management/commands/loadreplica.py`: given a path to a copy of MBPDB's
  `db.sqlite3`, `ATTACH DATABASE`s it and does explicit-column `INSERT ... SELECT` for all
  five tables (children deleted/inserted in FK-safe order). No Django serialization step, so
  the replica is byte-identical to the source — verified: `run_blast_search` against the
  replica returns results identical to querying MBPDB's live `db.sqlite3` directly for the
  same peptides (both exact-match and BLAST-homology paths tested).
- `data_transformation/services/blast_search.py` now imports directly from
  `mbpdb_replica.models` — the Phase-1 lazy-import/`RuntimeError` guard is gone; an empty
  (not-yet-loaded) replica just returns no results rather than erroring.
- `mbpdb_seed.sqlite3` — a checked-in snapshot of MBPDB's `db.sqlite3` (2MB), loaded into
  `mbpdb_replica` at every container boot via `start.sh` (after `bootstrap`). Refreshing it
  is a manual "copy a newer MBPDB db.sqlite3 over this file, commit, rebuild" step for now.
  **Explicitly confirmed with the user before committing**: `kuhfeldrf/peptiline` is a public
  GitHub repo and the README states MBPDB's data is "privately maintained" and "not part of
  this codebase" — bundling a real data snapshot contradicts that line. User's call: commit
  it anyway (the bioactivity annotations aren't actually sensitive). README's "MBPDB
  integration" section should be updated to reflect that this is no longer accurate (Phase 6).
- Dockerfile: added `ncbi-blast+` back (needed now that search is real) and
  `touch`/chown/chmod for `mbpdb_replica.sqlite3` alongside `db.sqlite3`.
- `peptiline/management/commands/bootstrap.py`: also runs
  `migrate mbpdb_replica --database=mbpdb_replica` so the replica schema exists on first boot,
  independent of whether `loadreplica` succeeds.

**Verified locally end-to-end** (venv + full Docker container): migrations create both DBs
correctly, `loadreplica` against a copy of MBPDB's real `db.sqlite3` loads 2806 rows across
5 tables (47/26/727/912/1094 — matches MBPDB's live counts exactly), exact-match search
(`REKVLASS`, `YLGSRY`) and homology search (`REKVLASA` @ 80% → hits `REKVLASS` @ 87.5%) both
return correct results, cross-checked against direct SQL against MBPDB's own `db.sqlite3`.
All app routes still 200/302 in the container with the new boot step added.

Rebuilt, pushed, and redeployed to Azure (`peptilinecontainer--replica1`) — verified live:
all routes 200/302, example files and corrected supplementals load, search returns real
results. A real cross-container refresh mechanism (Blob storage push/pull, scheduled job)
remains future work — this is a one-time-per-image-build snapshot, not live sync, by design
for this pass.

---

## 0d. UPDATE 2026-08-25 — repo cleanup: legacy_v1/ folder, README refresh

Per user request: created top-level `legacy_v1/` (matching the naming convention already
used by `examples/legacy_v1/` and `utils/legacy_v1/`) and moved everything from the
pre-Django notebook era into it: the five `.ipynb` files, `_settings.py`, and
`notebook_requirements.txt`. Also moved the root `supplementals/` folder (60MB) in there —
confirmed via `0c` above that it's a stale, superseded set (S1-S13/S1-S6, `zukaitis_2026`-
style numbering) that the live app never actually served; the real, current supplementals
(S1-S15/S1-S5, `kuhfeld_2026`) already live under `static/peptide/publications/`. Added
`legacy_v1/README.md` explaining what's there and why. Updated `.dockerignore` (one
`legacy_v1` entry replaces the old separate `supplementals`/`notebook_requirements.txt`
lines) and the top-level `README.md`'s "Project structure", "Data and reproducibility"
(fixed the now-corrected supplementals claim), and "Legacy notebook implementation" sections
to match.

**Caught a real regression while verifying**: moving files was safe, but the pytest suite's
`conftest.py` broke — 13 new failures, all `AppRegistryNotReady`. Root cause was unrelated to
the move itself: earlier in this session (`0c`), `blast_search.py`'s MBPDB-model import was
changed from a try/except-guarded `peptide.models` import to a direct, unguarded
`mbpdb_replica.models` import. `mbpdb_replica` is a real local app, so defining its model
classes now requires the Django app registry to be populated — but `conftest.py` only ever
called `settings.configure(INSTALLED_APPS=[], ...)`, never `django.setup()`. Fixed by adding
`mbpdb_replica` to the test `INSTALLED_APPS`, giving it an in-memory `DATABASES` entry, and
calling `django.setup()`. Also dropped the `conftest.py`'s now-fully-dead `peptide.celery`
sys.modules stub (nothing has imported `peptide.*` since the Phase 1 rename). Back to the
same 314/318 passing (the 2 pre-existing, documented `TestTransferFromDtFasta` failures) after
the fix. Verified Docker build still produces the same 184 static files (nothing needed by
the app was in the moved directories) and all routes still 200/302.
- [ ] Add a second `DATABASES` alias (`mbpdb_replica`) + a Django database router to
      PeptiLine's settings so replica tables are physically separate from PeptiLine's own
      operational SQLite file.
- [ ] On the MBPDB side, add a `dumpreplica` management command exporting `ProteinInfo`,
      `ProteinVariant`, `PeptideInfo`, `Function`, `Reference`.
- [ ] Decide and implement the transport (push to Blob Storage on a schedule, vs. an
      authenticated pull endpoint) — favor push-to-Blob so PeptiLine never needs inbound
      network access to MBPDB.
- [ ] On the PeptiLine side, add a `loadreplica` command that reads the dump and loads it
      into the `mbpdb_replica` alias, then rebuilds the BLAST FASTA/index
      (`makeblastdb`) from the refreshed protein data.
- [ ] Wire `loadreplica` into PeptiLine's `start.sh` (run once at boot) and into a recurring
      schedule (Celery beat or Container Apps scheduled job).
- [ ] Point `blast_search.py`'s queries at the `mbpdb_replica` alias explicitly (via the
      router or `.using('mbpdb_replica')`).
- [ ] End-to-end test: run a real MBPDB search from PeptiLine's UI against replica data and
      confirm results match what MBPDB itself would return for the same peptide.

### Phase 3 — Containerize PeptiLine — **done 2026-08-24**
- [x] Wrote PeptiLine's `Dockerfile`, `nginx.conf`, `start.sh`, modeled on MBPDB's, but
      **without** `ncbi-blast+` or the PEPEX perl scripts (no BLAST-backed search exists in
      this deployment yet — MBPDB's OS deps that are only used for that are skipped; adding
      `ncbi-blast+` back is part of Phase 2 when the replica is built).
- [x] Added a `health_check` view + `/health/` route (`peptiline/views.py`,
      `peptiline/urls.py`) for probe compatibility; wrote the (correct, httpGet-based)
      `lib/container-app-probes.yaml` from the start rather than the plain-TCP probes MBPDB's
      live app currently runs (see Phase 4 notes below).
- [x] Added `.gitignore` and `.dockerignore` (repo had neither). The first Docker build
      without a `.dockerignore` came out to 2.09GB, ~246MB of which was `supplementals/`
      (60MB) and `examples/` (32MB) — manuscript figures/sample datasets the running app
      never reads — getting `COPY`'d in and then duplicated again by the `chown -R` layer.
      Adding `.dockerignore` (excluding those, `.git`, notebooks, `tests/`, `docs/`) dropped
      the build context from 96MB to 2.6MB and the final image to 1.85GB (the remaining size
      is ~1.1GB from the `python:3.10` base image itself — identical to what MBPDB's own
      image already carries — plus ~716MB of legitimately-needed numpy/pandas/scipy/plotly/
      kaleido/biopython/matplotlib).
- [x] Built and ran the image locally (`docker build`, `docker run -p 8124:8000`); confirmed
      nginx → gunicorn → Django → Celery all start, `manage.py bootstrap` runs migrations at
      boot, and all app routes return 200/302 as expected. (Replica loading is Phase 2, not
      yet applicable.)
- [x] Fixed one bug found only by the containerized run (not caught by local `runserver`
      testing): `peptiline`'s custom `bootstrap` management command wasn't discovered because
      the `peptiline` project package itself wasn't in `INSTALLED_APPS` — `start.sh` logged
      `Unknown command: 'bootstrap'` and silently skipped migrations, so `/` and
      `/data_transformation/` 500'd. Fixed by adding `"peptiline"` to `INSTALLED_APPS`;
      re-verified all routes 200/302 after the fix.
- [x] Tagged and pushed `mbpdb/peptiline:latest` to Docker Hub (same account/namespace as
      MBPDB's existing image, per user decision) — login done manually by the user outside
      this session (password/PAT never entered into the assistant conversation).

### Phase 4 — Azure + CI

**Confirmed MBPDB Container App settings (`az containerapp show`, 2026-08-24)** —
mirror these for `peptilinecontainer` unless there's a specific reason to diverge:
- Resource group: `COH_MBPDB_RG`, environment: `managedEnvironment-mbpdb-a6f4` (West US 3)
- Scale: `minReplicas: 0`, `maxReplicas: 3`, HTTP scale rule at 10 concurrent
  requests/replica, `cooldownPeriod: 300`, `pollingInterval: 30` — this is the scale-to-zero
  behavior assumed in the cost discussion above; PeptiLine should use the same.
- Resources per replica: `cpu: 1`, `memory: 2Gi`, `ephemeralStorage: 4Gi`
- Ingress: external, `targetPort: 8000`, `activeRevisionsMode: Single`
- Probes currently configured are plain `tcpSocket` on port 8000 (Liveness/Readiness/
  Startup), **not** the httpGet-`/health/` probes described in
  `lib/container-app-probes.yaml` in the MBPDB repo — that file exists but does not appear
  to have been applied to the live app. Decide whether PeptiLine should ship with the
  (better) httpGet probes from day one rather than copying the currently-live TCP probes.
- Registry: Docker Hub (`index.docker.io` / `docker.io`), same pattern as MBPDB
  (`registryUrl: docker.io` in the GitHub Actions workflow) — PeptiLine needs its own image
  repo (e.g. `mbpdb/peptiline:latest` or a new Docker Hub namespace) and its own registry
  credentials secret, not reuse of MBPDB's `dockerio-mbpdb` secret.
- Custom domain: MBPDB's `mbpdb.nws.oregonstate.edu` is bound directly to the Container App
  via an Azure-managed certificate (`SniEnabled`) — PeptiLine's new hostname will need the
  same managed-certificate binding set up once DNS is decided (§7).
- Note: `DJANGO_SECRET_KEY`, `DJANGO_SUPERUSER_PASSWORD`, and `GITHUB_PAT` are stored as
  **plain env vars** on MBPDB's container (visible via `az containerapp show`), not as
  Container Apps secrets referenced by name. Worth moving these to real `secretRef`s for
  both apps while touching this configuration, rather than reproducing the same pattern on
  PeptiLine.

- [x] `mbpdb/peptiline:latest` pushed to Docker Hub (see Phase 3) — ready to reference below.
- [x] **Provisioned 2026-08-25.** `peptilinecontainer` created in `COH_MBPDB_RG` /
      `managedEnvironment-mbpdb-a6f4`, mirroring MBPDB's settings (1 vCPU/2Gi,
      `minReplicas: 0`/`maxReplicas: 3`, http scale rule at 10 concurrent/replica). No
      registry credentials needed — `mbpdb/peptiline` is a public Docker Hub repo. Actual
      command used (differs slightly from the draft below — no registry flags needed, and
      `DJANGO_ALLOWED_HOSTS` had to be added after creation, see the incident log in §0):
      ```bash
      az containerapp create \
        --name peptilinecontainer --resource-group COH_MBPDB_RG \
        --environment managedEnvironment-mbpdb-a6f4 \
        --image docker.io/mbpdb/peptiline:latest \
        --target-port 8000 --ingress external \
        --cpu 1 --memory 2Gi --min-replicas 0 --max-replicas 3 \
        --scale-rule-name http-scaler --scale-rule-type http --scale-rule-http-concurrency 10 \
        --env-vars DJANGO_SETTINGS_MODULE=peptiline.settings WEBSITES_PORT=8000 \
          WEBSITES_CONTAINER_START_TIME_LIMIT=1800 DJANGO_SECRET_KEY="..."
      az containerapp update --name peptilinecontainer --resource-group COH_MBPDB_RG \
        --set-env-vars DJANGO_SETTINGS_MODULE=peptiline.settings WEBSITES_PORT=8000 \
          WEBSITES_CONTAINER_START_TIME_LIMIT=1800 DJANGO_SECRET_KEY="..." \
          DJANGO_ALLOWED_HOSTS="<fqdn>,localhost,127.0.0.1"   # --set-env-vars REPLACES the list, always resupply all of them
      az containerapp update --name peptilinecontainer --resource-group COH_MBPDB_RG \
        --yaml lib/container-app-probes.yaml
      ```
      Live at `https://peptilinecontainer.lemonisland-71b15397.westus3.azurecontainerapps.io/`.
      Original draft command, kept for reference (assumed a private registry, which turned
      out not to apply):
      ```bash
      az containerapp create \
        --name peptilinecontainer \
        --resource-group COH_MBPDB_RG \
        --environment managedEnvironment-mbpdb-a6f4 \
        --image docker.io/mbpdb/peptiline:latest \
        --registry-server index.docker.io \
        --registry-username mbpdb \
        --registry-password "$DOCKERHUB_PASSWORD_OR_TOKEN" \
        --target-port 8000 \
        --ingress external \
        --cpu 1 --memory 2Gi \
        --min-replicas 0 --max-replicas 3 \
        --scale-rule-name http-scaler \
        --scale-rule-type http \
        --scale-rule-http-concurrency 10 \
        --env-vars \
          DJANGO_SETTINGS_MODULE=peptiline.settings \
          WEBSITES_PORT=8000 \
          WEBSITES_CONTAINER_START_TIME_LIMIT=1800 \
          DJANGO_SECRET_KEY="$PEPTILINE_DJANGO_SECRET_KEY" \
          DJANGO_ALLOWED_HOSTS="$PEPTILINE_ALLOWED_HOSTS"
      # Then apply the httpGet probes (az containerapp create has no probe flags):
      az containerapp update \
        --name peptilinecontainer --resource-group COH_MBPDB_RG \
        --yaml lib/container-app-probes.yaml
      ```
      Per user decision (2026-08-24): leave `DJANGO_SECRET_KEY`/superuser password/etc. as
      plain env vars, matching MBPDB's current pattern, rather than Container Apps
      `secretRef`s — only the probe mechanism is being upgraded, not the secrets handling.
- [ ] Provision a new service principal / Azure credentials, registry credentials, and add
      them as GitHub secrets in this repo (`PEPTILINECONTAINER_AZURE_CREDENTIALS`,
      `PEPTILINECONTAINER_REGISTRY_USERNAME`/`_PASSWORD`, `PEPTILINE_DJANGO_SECRET_KEY`, etc.
      — see `.github/workflows/deploy.yml`, already written and mirrors MBPDB's workflow).
- [ ] Provision the Blob Storage container (or equivalent) used for the replica dump, with
      MBPDB granted write and PeptiLine granted read. (Phase 2 prerequisite, not yet needed.)
- [x] `.github/workflows/deploy.yml` added to this repo, modeled on MBPDB's deploy workflow
      (`containerAppName: peptilinecontainer`, `mbpdb/peptiline:latest`, its own secret names
      — not yet wired to real GitHub secrets, so this workflow will fail if triggered until
      the secrets above are added).
- [ ] Set up DNS/custom domain + TLS for PeptiLine's new hostname (or configure the
      reverse-proxy/redirect approach chosen in Phase 0 instead). Currently only reachable
      at the default `*.azurecontainerapps.io` hostname.
- [x] Deploy PeptiLine standalone; smoke-test all three modules — **done 2026-08-25**, all 9
      routes (`/ /health/ /data_transformation/ /data_analysis/ /heatmap/ /about_us/ /pepex/
      /peptiline/supplementals/ /admin/`) return 200/302 against the live deployment. MBPDB
      search itself also verified live (see §0c) — the replica loads at boot and returns
      real, correct results.

### Phase 5 — Cut over MBPDB's in-process PeptiLine routes
- [ ] Replace MBPDB's `/peptiline/`, `/data_transformation/`, `/data_analysis/`,
      `/heatmap/` views with 301 redirects to the new PeptiLine host's equivalent paths.
- [ ] Verify old links (including any already published/cached) resolve correctly through
      the redirect.
- [ ] After a confirmation period, delete PeptiLine's in-process app code from the MBPDB
      repo (`data_transformation/`, `data_analysis/`, `heatmap_viz/` under
      `include/peptide/peptide/`) and the corresponding `INSTALLED_APPS`/`urls.py` entries.
- [ ] Update MBPDB's `Dockerfile`/`requirements.txt` to drop now-unused dependencies
      (BLAST is still needed for MBPDB's own search, so likely nothing to drop there, but
      recheck `requirements.txt` for PeptiLine-only packages).

---

## 0e. UPDATE 2026-08-25 — synced with MBPDB feature commits since vendoring

The 2026-08-14 vendoring pass was a point-in-time snapshot; MBPDB's `main` kept moving.
Diffed `~/mbpdb/include/peptide/peptide/{data_transformation,data_analysis,heatmap_viz,utils}`
against this repo's copies (the Phase-0 diff task, done properly this time — the earlier
"done" checkmark only verified the namespace-rewrite, not actual content drift) and found six
commits' worth of real changes to port: Data Analysis's Plot Filter dropdown removal
(derive the mode from selections instead), Data Transformation's FASTA-upload move to
Protein Mapping + iterrows() vectorization, and Heatmap's axis-lock/zoom fix. Ported each,
reapplying PeptiLine's standalone import adjustments on top (see commit `fa313b2`). Also
brought over `utils/lazy_import.py` (LazyModule), previously not ported — it's directly
relevant to PeptiLine too, not just MBPDB, since both run at `minReplicas: 0`.

**Also found and fixed**: the landing-page gallery plots (`static/peptide/plots/`, 27 files)
were never vendored at all — same class of gap as `hero.png`/screenshots/examples before.

**Found, did not fix**: porting `heatmap_renderer.py` surfaced that MBPDB's own `main` has
had 6 failing tests in `test_heatmap_legend_position.py` since commit `6f1b606d` (changed
`_below_legend_geometry`'s return signature without updating the test) — confirmed by running
MBPDB's own test suite directly, not just inferring from the port. This is upstream, not
something to silently paper over from this repo; ported the same (currently broken) test file
for accurate parity and flagged it to the user. **Someone should fix this in MBPDB directly**
and PeptiLine will need the corresponding test-file re-sync when that happens.

Verified: `manage.py check` clean, full suite gives the same 308 passed/8 failed as MBPDB's
own working tree, Docker build confirms the actual UI changes (plot-filter dropdown gone,
FASTA upload under Protein Mapping) and all routes/gallery assets 200.

---

### Phase 6 — Documentation and cleanup
- [ ] Update this repo's `README.md` (hosted URL, MBPDB integration section — live search
      now works in the standalone deployment via replica, not just TSV upload/local-only).
- [ ] Write `docs/INSTALL.md` and `docs/REPRODUCIBILITY.md` referenced by the README but not
      yet confirmed to exist — verify/create.
- [ ] Document the replication mechanism itself (this file, or a new
      `docs/DATABASE_REPLICATION.md`) so the next person understands the one-way,
      scheduled-dump design and its known limitations (§7).
- [ ] Note Postgres/real-replica migration as explicit future work, not silently forgotten.
