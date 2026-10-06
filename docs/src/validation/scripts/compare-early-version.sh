#!/usr/bin/env bash
# Run compare-early-version.jl for an early release (default v0.1.6) and for
# the current checkout, at t = N and t = 10N. Run from the repository root:
#
#   docs/src/validation/scripts/compare-early-version.sh [old_tag=v0.1.6] [reps=30] [seed=2026] [threads=8] \
#       > docs/src/validation/early-version-comparison.txt
#
# The early release is checked out into a temporary Git worktree, which is
# removed on exit.
set -euo pipefail

OLD=${1:-v0.1.6}
REPS=${2:-30}
SEED=${3:-2026}
THREADS=${4:-8}

ROOT=$(git rev-parse --show-toplevel)
SCRIPT="$ROOT/docs/src/validation/scripts/compare-early-version.jl"
WT=$(mktemp -d)
trap 'git -C "$ROOT" worktree remove --force "$WT" >/dev/null 2>&1 || true; rm -rf "$WT"' EXIT

git -C "$ROOT" worktree add --detach -f "$WT" "$OLD" >/dev/null 2>&1
julia --startup-file=no --project="$WT" -e 'using Pkg; Pkg.instantiate()' >/dev/null 2>&1
julia --startup-file=no --project="$ROOT" -e 'using Pkg; Pkg.instantiate()' >/dev/null 2>&1

TREE=$(git -C "$ROOT" diff --quiet HEAD && echo clean || echo dirty)
echo "# Date: $(date -u +%Y-%m-%dT%H:%MZ)  Commit: $(git -C "$ROOT" rev-parse --short HEAD) ($(git -C "$ROOT" describe --tags --always)); tree: $TREE"
echo "# Early release: $OLD ($(git -C "$ROOT" rev-parse --short "$OLD^{commit}"))"
echo "# CPU: $(lscpu 2>/dev/null | sed -n 's/^Model name:[ \t]*//p' | head -n1); $(julia --version)"
echo "# \$ docs/src/validation/scripts/compare-early-version.sh $OLD $REPS $SEED $THREADS"
echo

for T in 1 10; do
    for PROJ in "$WT" "$ROOT"; do
        julia --startup-file=no -t "$THREADS" --project="$PROJ" "$SCRIPT" "$T" "$REPS" "$SEED" 2>/dev/null
    done
done
