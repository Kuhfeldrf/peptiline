#!/bin/bash
# Apply deploy/probes.json to the live Container App's probe config.
#
# Probes are NOT covered by the deploy workflow's `az containerapp update
# --set-env-vars`, and there is no other IaC for them -- so without this script
# they only ever exist as hand-edited Azure state. That has already burned us
# once: MBPDB shipped a /health/ view and a commit claiming to have repointed
# the probes at it, but the live probes stayed on bare tcpSocket:8000 (which
# passes as soon as nginx binds :8000, ~7s before gunicorn is listening behind
# it on :8001, so a replica reports ready while requests still 502).
#
# az containerapp update takes only --yaml, but YAML is a superset of JSON, so
# we hand it JSON and skip needing PyYAML on the runner.
#
# Idempotent: if the live probes already match, ACA sees an unchanged template
# and does not cut a new revision.
set -euo pipefail

APP="${APP:-peptilinecontainer}"
RG="${RG:-COH_MBPDB_RG}"
PROBES="$(dirname "$0")/probes.json"

workdir=$(mktemp -d)
trap 'rm -rf "$workdir"' EXIT

echo "Fetching current config for $APP..."
az containerapp show -n "$APP" -g "$RG" -o json --only-show-errors > "$workdir/current.json"

# Inject probes, then drop server-managed fields that `update` rejects.
jq --slurpfile probes "$PROBES" '
    .properties.template.containers[0].probes = $probes[0]
  | del(.id, .name, .type, .systemData, .location, .identity)
  | del(.properties.provisioningState,
        .properties.latestRevisionName,
        .properties.latestRevisionFqdn,
        .properties.latestReadyRevisionName,
        .properties.customDomainVerificationId,
        .properties.outboundIpAddresses,
        .properties.eventStreamEndpoint,
        .properties.runningStatus,
        .properties.workloadProfileName,
        .properties.environmentId,
        .properties.managedEnvironmentId,
        .properties.delegatedIdentities)
  | del(.properties.configuration.ingress.fqdn)
' "$workdir/current.json" > "$workdir/desired.yaml"

echo "Applying probes:"
jq -r '.[] | "  \(.type): \(.httpGet.path // "tcp:\(.tcpSocket.port)") every \(.periodSeconds)s x\(.failureThreshold)"' "$PROBES"

az containerapp update -n "$APP" -g "$RG" --yaml "$workdir/desired.yaml" --only-show-errors > /dev/null

echo "Verifying live probe config..."
az containerapp show -n "$APP" -g "$RG" \
  --query "properties.template.containers[0].probes" -o json --only-show-errors \
  | jq -S 'sort_by(.type)' > "$workdir/live.json"
jq -S 'sort_by(.type)' "$PROBES" > "$workdir/want.json"

if diff -q "$workdir/want.json" "$workdir/live.json" > /dev/null; then
  echo "OK: live probes match deploy/probes.json"
else
  echo "ERROR: live probes do not match deploy/probes.json" >&2
  diff "$workdir/want.json" "$workdir/live.json" >&2 || true
  exit 1
fi
