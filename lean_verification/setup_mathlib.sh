#!/bin/bash
# Robust Mathlib setup for lean_verification/
# Handles intermittent wifi dropouts with retries.
# Run from: c:\irreducible-twistor-geometry\lean_verification\
#
# Usage: bash setup_mathlib.sh

set -e
export PATH="$HOME/.elan/bin:$PATH"

MAX_RETRIES=10
RETRY_DELAY=15

retry() {
  local n=0
  local cmd="$@"
  until [ $n -ge $MAX_RETRIES ]; do
    echo "[Attempt $((n+1))/$MAX_RETRIES] $cmd"
    if eval "$cmd"; then
      return 0
    fi
    n=$((n+1))
    if [ $n -lt $MAX_RETRIES ]; then
      echo "  Failed. Retrying in ${RETRY_DELAY}s..."
      sleep $RETRY_DELAY
    fi
  done
  echo "  FAILED after $MAX_RETRIES attempts."
  return 1
}

echo "============================================"
echo "Mathlib Setup for TwistorVerified"
echo "============================================"
echo ""

# Step 1: Clean any previous failed state
echo "[1/5] Cleaning previous state..."
rm -rf .lake lake-manifest.json lake-packages 2>/dev/null || true

# Step 2: Write the Mathlib-compatible lakefile
echo "[2/5] Writing lakefile with Mathlib dependency..."
cat > lakefile.lean << 'LAKEFILE'
import Lake
open Lake DSL

package «TwistorVerified» where
  leanOptions := #[
    ⟨`autoImplicit, false⟩
  ]

require mathlib from git
  "https://github.com/leanprover-community/mathlib4"

@[default_target]
lean_lib «TwistorVerified» where
  roots := #[`TwistorVerified]
LAKEFILE

# Step 3: Resolve dependencies (clones Mathlib — this is the big download)
echo "[3/5] Resolving dependencies (cloning Mathlib, ~2GB)..."
echo "       This may take 10-30 minutes on first run."
retry "lake update 2>&1 | tail -5"

# Step 4: Download prebuilt Mathlib oleans (avoids building from source)
echo ""
echo "[4/5] Downloading prebuilt Mathlib cache (~4GB)..."
echo "       This avoids a multi-hour build from source."
retry "lake exe cache get 2>&1 | tail -5"

# Step 5: Test build
echo ""
echo "[5/5] Building TwistorVerified..."
if lake build 2>&1; then
  echo ""
  echo "============================================"
  echo "SUCCESS! Mathlib is ready."
  echo "You can now use norm_num, interval_cases, ℚ, etc."
  echo "============================================"
else
  echo ""
  echo "============================================"
  echo "Build failed. Check errors above."
  echo "The Mathlib download may have succeeded —"
  echo "errors might be in the .lean files themselves."
  echo "Try: cd lean_verification && lake build"
  echo "============================================"
  exit 1
fi
