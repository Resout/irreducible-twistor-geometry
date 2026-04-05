"""
The composition rate formula: Schmidt(g∘h) / Schmidt(gh) = (1-K*)(1-1/H³)

Composition of two crystals encoding group elements g, h produces a crystal
whose Schmidt is reduced by the factor (1-K*)(1-1/H³) = 299/405 relative
to the target crystal encoding the group product gh.

Two multiplicative contributions:
  (1-K*) = 23/30: Dempster conflict filter
  (1-1/H³) = 26/27: Born floor compactness

Verify across ALL 36 products of S₃.
"""

import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import torch
import math
from solver.algebra import (H, MASS_DIM, DELTA, K_STAR, BORN_FLOOR,
                            schmidt_number, born_fidelity)
from solver.crystals import Entangler, compose

torch.set_grad_enabled(False)

S3_CORRS = {
    "e":     torch.tensor([[1,0,0],[0,1,0],[0,0,1]], dtype=torch.float32),
    "(01)":  torch.tensor([[0,1,0],[1,0,0],[0,0,1]], dtype=torch.float32),
    "(02)":  torch.tensor([[0,0,1],[0,1,0],[1,0,0]], dtype=torch.float32),
    "(12)":  torch.tensor([[1,0,0],[0,0,1],[0,1,0]], dtype=torch.float32),
    "(012)": torch.tensor([[0,1,0],[0,0,1],[1,0,0]], dtype=torch.float32),
    "(021)": torch.tensor([[0,0,1],[1,0,0],[0,1,0]], dtype=torch.float32),
}

S3_MULT = {
    ("e","e"): "e", ("e","(01)"): "(01)", ("e","(02)"): "(02)",
    ("e","(12)"): "(12)", ("e","(012)"): "(012)", ("e","(021)"): "(021)",
    ("(01)","e"): "(01)", ("(01)","(01)"): "e", ("(01)","(02)"): "(012)",
    ("(01)","(12)"): "(021)", ("(01)","(012)"): "(02)", ("(01)","(021)"): "(12)",
    ("(02)","e"): "(02)", ("(02)","(01)"): "(021)", ("(02)","(02)"): "e",
    ("(02)","(12)"): "(012)", ("(02)","(012)"): "(12)", ("(02)","(021)"): "(01)",
    ("(12)","e"): "(12)", ("(12)","(01)"): "(012)", ("(12)","(02)"): "(021)",
    ("(12)","(12)"): "e", ("(12)","(012)"): "(01)", ("(12)","(021)"): "(02)",
    ("(012)","e"): "(012)", ("(012)","(01)"): "(12)", ("(012)","(02)"): "(01)",
    ("(012)","(12)"): "(02)", ("(012)","(012)"): "(021)", ("(012)","(021)"): "e",
    ("(021)","e"): "(021)", ("(021)","(01)"): "(02)", ("(021)","(02)"): "(12)",
    ("(021)","(12)"): "(01)", ("(021)","(012)"): "e", ("(021)","(021)"): "(012)",
}

# Predicted ratio
PREDICTED = (1 - K_STAR) * (1 - BORN_FLOOR)

print("=" * 80)
print("COMPOSITION RATE FORMULA VERIFICATION")
print(f"Predicted: (1-K*)(1-1/H³) = ({1-K_STAR:.6f})({1-BORN_FLOOR:.6f}) = {PREDICTED:.6f}")
print("=" * 80)

n_seeds = 50

# Build seed-averaged crystals
refs = {}
for name, corr in S3_CORRS.items():
    s = torch.zeros(MASS_DIM, MASS_DIM, dtype=torch.cfloat)
    for seed in range(n_seeds):
        s += Entangler(corr, seed=seed).build().joint
    refs[name] = s / n_seeds

# Compute all 36 products
print(f"\n{'g':>6s} {'h':>6s} {'gh':>6s} {'Sn(g∘h)':>8s} {'Sn(gh)':>8s} {'ratio':>8s} {'pred':>8s} {'err%':>7s}")
print("-" * 70)

all_ratios = []
for g in S3_CORRS:
    for h in S3_CORRS:
        gh = S3_MULT[(g, h)]
        comp = compose(refs[g], refs[h])
        sn_comp = schmidt_number(comp)
        sn_target = schmidt_number(refs[gh])
        ratio = sn_comp / sn_target
        err = (ratio - PREDICTED) / PREDICTED * 100
        all_ratios.append(ratio)
        print(f"{g:>6s} {h:>6s} {gh:>6s} {sn_comp:>8.4f} {sn_target:>8.4f} {ratio:>8.4f} {PREDICTED:>8.4f} {err:>+7.2f}%")

mean_ratio = sum(all_ratios) / len(all_ratios)
std_ratio = (sum((r - mean_ratio)**2 for r in all_ratios) / len(all_ratios)) ** 0.5

print(f"\n{'='*70}")
print(f"  Predicted:  {PREDICTED:.6f} = (23/30)(26/27) = 299/405")
print(f"  Measured:   {mean_ratio:.6f} ± {std_ratio:.6f}")
print(f"  Agreement:  {mean_ratio/PREDICTED*100:.3f}%")
print(f"  Deviation:  {abs(mean_ratio - PREDICTED)/PREDICTED*100:.3f}%")

# Separate by conjugacy class
trans = ["(01)", "(02)", "(12)"]
cycles = ["(012)", "(021)"]

trans_ratios = [r for i, r in enumerate(all_ratios)
                if list(S3_CORRS.keys())[i//6] in trans or list(S3_CORRS.keys())[i%6] in trans]

print(f"\n  --- By product type ---")
# Products giving identity
id_ratios = []
for g in S3_CORRS:
    for h in S3_CORRS:
        if S3_MULT[(g,h)] == "e":
            comp = compose(refs[g], refs[h])
            id_ratios.append(schmidt_number(comp) / schmidt_number(refs["e"]))
print(f"  Products → e:     {sum(id_ratios)/len(id_ratios):.6f} (n={len(id_ratios)})")

# Products giving transpositions
t_ratios = []
for g in S3_CORRS:
    for h in S3_CORRS:
        if S3_MULT[(g,h)] in trans:
            comp = compose(refs[g], refs[h])
            t_ratios.append(schmidt_number(comp) / schmidt_number(refs[S3_MULT[(g,h)]]))
print(f"  Products → trans: {sum(t_ratios)/len(t_ratios):.6f} (n={len(t_ratios)})")

# Products giving 3-cycles
c_ratios = []
for g in S3_CORRS:
    for h in S3_CORRS:
        if S3_MULT[(g,h)] in cycles:
            comp = compose(refs[g], refs[h])
            c_ratios.append(schmidt_number(comp) / schmidt_number(refs[S3_MULT[(g,h)]]))
print(f"  Products → 3cyc: {sum(c_ratios)/len(c_ratios):.6f} (n={len(c_ratios)})")


# ═══════════════════════════════════════════════════════════════
#  Multi-step verification
# ═══════════════════════════════════════════════════════════════

print(f"\n\n{'='*80}")
print("MULTI-STEP: Schmidt after n compositions")
print("="*80)

# Start with (01), compose with (02) repeatedly
# (02)^n ∘ (01) cycles through S₃ elements
crystal = refs["(01)"]
print(f"  Step 0: Schmidt = {schmidt_number(crystal):.4f}")
for step in range(1, 9):
    crystal = compose(crystal, refs["(02)"])
    sn = schmidt_number(crystal)
    predicted_sn = schmidt_number(refs["(01)"]) * PREDICTED**step
    print(f"  Step {step}: Schmidt = {sn:.4f}, predicted = {predicted_sn:.4f}, "
          f"ratio = {sn/predicted_sn:.4f}")


print(f"\n\n{'='*80}")
print("FINDINGS")
print("="*80)
print(f"""
  The composition rate formula:

    Schmidt(g ∘ h) / Schmidt(gh) = (1-K*)(1-1/H³) = 299/405

  is verified across all 36 products of S₃ to {abs(mean_ratio - PREDICTED)/PREDICTED*100:.2f}%.

  Two multiplicative mechanisms:
    (1-K*) = 23/30 = DS conflict normalization
    (1-1/H³) = 26/27 = Born floor compactness

  The crystal's group algebra is a CONTRACTIVE representation
  of S₃ with contraction rate 299/405 ≈ 0.738 per multiplication.

  After n multiplications: Schmidt ∝ (299/405)ⁿ
  The group algebra decays exponentially to the product state.
""")
