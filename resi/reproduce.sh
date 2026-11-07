#!/bin/sh
# Rerun every analysis and figure from the cached public data.
# Usage: sh reproduce.sh CACHE_DIR   (fetch first: python fetch_cd20.py CACHE_DIR; python fetch_cd20.py CACHE_DIR npc)
set -e
CACHE="$1"
cd "$(dirname "$0")"
PY="${PYTHON:-python}"
$PY verify_exact.py
$PY verify_model.py
$PY verify_triplets.py
$PY verify_angles.py
$PY pairs_vs_triplets.py
$PY run_cd20.py "$CACHE"
$PY registration_cd20.py "$CACHE"
$PY summarize_cd20.py "$CACHE"
$PY run_cd20_molecules.py "$CACHE"
$PY triplets_cd20.py "$CACHE"
$PY permanence_cd20.py "$CACHE"
$PY design_cd20.py "$CACHE"
$PY budget_cd20.py "$CACHE"
$PY shape_cd20.py "$CACHE"
$PY npc_validation.py "$CACHE"
$PY zaza_npc.py "$CACHE"
$PY robustness.py "$CACHE"
$PY figures.py
