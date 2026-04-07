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
