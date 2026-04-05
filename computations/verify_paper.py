"""Verify structural integrity of the paper LaTeX file — updated for current state."""
import re

PAPER = r"paper/Independently Meandering on Through Entangling Informational And Computational Geometry for The Explicit Purpose of Having Quite a Lot of Fun by J. R. Manuel.tex"

with open(PAPER, "r", encoding="utf-8") as f:
    content = f.read()
    lines = content.split("\n")

print(f"Total lines: {len(lines)}")
print()

# =========================================================================
# 1. LaTeX environment matching
# =========================================================================
envs = [
    "theorem", "proposition", "remark", "definition", "conjecture",
    "construction", "lemma", "corollary", "proof", "document",
    "abstract", "enumerate", "itemize", "center", "tabular",
    "equation", "align",
]
env_ok = True
for env in envs:
    opens = content.count(r"\begin{" + env + "}")
    closes = content.count(r"\end{" + env + "}")
    if opens != closes:
        print(f"MISMATCH: {env} has {opens} opens and {closes} closes")
        env_ok = False
if env_ok:
    print("OK: All LaTeX environments matched")

# Count specific types
for env in ["theorem", "proposition", "remark", "definition", "lemma", "conjecture", "proof"]:
    n = content.count(r"\begin{" + env + "}")
    if n > 0:
        print(f"  {env}: {n}")

print()

# =========================================================================
# 2. Label/ref consistency
# =========================================================================
labels = set(re.findall(r"\\label\{([^}]+)\}", content))
refs = re.findall(r"\\ref\{([^}]+)\}", content)
missing_refs = sorted(set(r for r in refs if r not in labels))
if missing_refs:
    print(f"WARNING: {len(missing_refs)} refs to undefined labels:")
    for m in missing_refs:
        print(f"  \\ref{{{m}}}")
else:
    print(f"OK: All {len(refs)} refs resolve ({len(labels)} labels defined)")

# Unused labels
used_labels = set(refs)
unused = sorted(labels - used_labels)
if unused:
    print(f"INFO: {len(unused)} labels never referenced:")
    for u in unused:
        print(f"  {u}")

print()

# =========================================================================
# 3. Bibliography: cited vs defined
# =========================================================================
bibitems = set(re.findall(r"\\bibitem\{([^}]+)\}", content))
cites_raw = re.findall(r"\\cite\{([^}]+)\}", content)
cited = set()
for c in cites_raw:
    for key in c.split(","):
        cited.add(key.strip())

uncited = sorted(bibitems - cited)
undefined_cites = sorted(cited - bibitems)

if uncited:
    print(f"UNCITED bibliography entries ({len(uncited)}):")
    for u in uncited:
        print(f"  {u}")
else:
    print(f"OK: All {len(bibitems)} bibliography entries are cited")

if undefined_cites:
    print(f"WARNING: {len(undefined_cites)} citations to undefined bibitems:")
    for u in undefined_cites:
        print(f"  {u}")

print()

# =========================================================================
# 4. Document structure
# =========================================================================
structure = [
    (r"\begin{document}", "\\begin{document}"),
    (r"\end{document}", "\\end{document}"),
    (r"\maketitle", "\\maketitle"),
    (r"\tableofcontents", "\\tableofcontents"),
    (r"\begin{thebibliography}", "bibliography"),
    (r"\end{thebibliography}", "end bibliography"),
]
for pattern, name in structure:
    if pattern in content:
        print(f"OK: {name}")
    else:
        print(f"MISSING: {name}")

# Check parts
parts = re.findall(r"\\part\{(.+?)\}", content)
print(f"\nParts found ({len(parts)}):")
for p in parts:
    print(f"  {p}")

print()

# =========================================================================
# 5. Key numerical claims
# =========================================================================
print("=== Key Numerical Claims ===")
numerics = [
    ("K^*=7/30", "K* = 7/30"),
    ("1/27", "Born floor 1/27"),
    ("0.2829", "lambda_0"),
    ("0.2813", "lambda_1"),
    ("1.263", "spectral gap Delta"),
    ("1.2626", "Delta precise value"),
    ("0.00084", "min |det(M)|"),
    ("810/7", "instanton action S"),
    ("0.011672", "proton/W ratio"),
    ("80{,}385", "M_W prediction"),
    ("0.2320", "K* empirical"),
    ("0.0033", "K* CV"),
]
for val, name in numerics:
    count = content.count(val)
    if count > 0:
        print(f"  OK: {name} ({val}) appears {count} time(s)")
    else:
        print(f"  MISSING: {name} ({val})")

print()

# =========================================================================
# 6. Check for known bad patterns
# =========================================================================
print("=== Integrity Checks ===")

# No 57.5sigma
if "57.5" in content:
    for i, line in enumerate(lines):
        if "57.5" in line:
            print(f"WARNING: 57.5 still appears at line {i+1}: {line.strip()[:80]}")
else:
    print("OK: No 57.5sigma claims remain")

# No false 50-digit Delta
if "1.26596590" in content:
    print("WARNING: False 50-digit Delta value still present")
else:
    print("OK: No false 50-digit Delta")

# No degree-32 claim
if "degree 32" in content.lower() or "degree-32" in content.lower():
    print("WARNING: Degree-32 claim still present")
else:
    print("OK: No degree-32 claim")

# No |K| < 1 for complex masses
if "|K|<1" in content or "|K| < 1" in content:
    for i, line in enumerate(lines):
        if "|K|<1" in line or "|K| < 1" in line:
            print(f"WARNING: |K|<1 at line {i+1}: {line.strip()[:80]}")
else:
    print("OK: No |K|<1 claims for complex masses")

# Check the phase example sums to 1
if "0.5e^{i\\pi/4}" in content:
    print("WARNING: Old broken phase example still present")
else:
    print("OK: Old broken phase example removed")

# H^4/K* vs H^3/K*
h4_count = content.count("H^4/K")
h3_count = content.count("H^3/K")
if h4_count > 0:
    print(f"CHECK: H^4/K appears {h4_count} time(s), H^3/K appears {h3_count} time(s)")
else:
    print(f"OK: H^3/K* used consistently ({h3_count} occurrences)")

# No conjectures
conj_count = content.count(r"\begin{conjecture}")
if conj_count > 0:
    print(f"WARNING: {conj_count} conjecture(s) remain")
else:
    print("OK: 0 conjectures")

# Ward→Mason correction checks
if "rem:floor_inactive" in content:
    print("WARNING: rem:floor_inactive still exists (should be rem:floor_active)")
else:
    print("OK: rem:floor_inactive removed")

if "rem:floor_active" in content:
    print("OK: rem:floor_active present")
else:
    print("WARNING: rem:floor_active missing")

# Check no false claims about dbar=0 at fixed point
import re as _re
dbar_zero_fp = _re.findall(r'dbar.*Phi.*=.*0.*fixed.point|floor.*inactive.*fixed', content, _re.IGNORECASE)
if dbar_zero_fp:
    print(f"WARNING: Possible false dbar=0 at fixed point claim found")
else:
    print("OK: No false dbar=0-at-fixed-point claims")

# Check Mason spectral identification exists
if "Mason Spectral Identification" in content:
    print("OK: Mason Spectral Identification theorem present")
elif "Ward Spectral Equivalence" in content:
    print("WARNING: Old Ward Spectral Equivalence still present")
else:
    print("CHECK: Neither Mason nor Ward spectral theorem found")

# Check ||dbar Phi|| = 1.12 is stated
if "1.12" in content:
    print("OK: ||dbar Phi|| = 1.12 value present")
else:
    print("WARNING: ||dbar Phi|| = 1.12 not found in paper")

# Check chain rule error is documented
if "chain rule" in content.lower():
    print("OK: Chain rule correction documented")
else:
    print("CHECK: Chain rule correction not explicitly mentioned")

# Check Mason uniform applicability
if "rem:mason_uniform" in content:
    print("OK: Mason uniform applicability remark present")
else:
    print("CHECK: Mason uniform remark not found")

# Non-triviality: proper numbers not fabricated
if "57.5" in content:
    print("WARNING: 57.5sigma still present")
elif "9.0" in content and "14.5" in content:
    print("OK: Non-triviality has correct numbers (9.0sigma, 14.5sigma)")
else:
    print("CHECK: Non-triviality numbers may need verification")

print()

# =========================================================================
# 7. Theorem/proof pairing
# =========================================================================
print("=== Theorem/Proof Pairing ===")
theorem_labels = re.findall(r"\\begin\{theorem\}\[.*?\]\\label\{([^}]+)\}", content)
theorem_labels += re.findall(r"\\begin\{theorem\}\\label\{([^}]+)\}", content)
# Also multi-line: label on next line
for i, line in enumerate(lines):
    if r"\begin{theorem}" in line and r"\label{" not in line:
        if i + 1 < len(lines) and r"\label{" in lines[i+1]:
            m = re.search(r"\\label\{([^}]+)\}", lines[i+1])
            if m and m.group(1) not in theorem_labels:
                theorem_labels.append(m.group(1))

# Count proofs
proof_count = content.count(r"\begin{proof}")
print(f"Theorems with labels: {len(theorem_labels)}")
print(f"Proofs: {proof_count}")

print()
print("=== SUMMARY ===")
print(f"Lines: {len(lines)}")
print(f"Theorems: {content.count(chr(10) + chr(92) + 'begin{theorem}')}")
print(f"Propositions: {content.count(chr(10) + chr(92) + 'begin{proposition}')}")
print(f"Lemmas: {content.count(chr(10) + chr(92) + 'begin{lemma}')}")
print(f"Remarks: {content.count(chr(10) + chr(92) + 'begin{remark}')}")
print(f"Definitions: {content.count(chr(10) + chr(92) + 'begin{definition}')}")
print(f"Conjectures: {conj_count}")
print(f"Bibliography: {len(bibitems)} entries, {len(cited)} cited, {len(uncited)} uncited")
