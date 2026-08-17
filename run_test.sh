#!/usr/bin/env bash
# One-command test setup: install R deps into renv, then start Shiny.
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$ROOT"

echo "==> Installing dependencies (skip if already present)..."
Rscript scripts/install_deps.R

echo "==> Starting Shiny (Ctrl+C to stop)..."
Rscript scripts/run_app.R
