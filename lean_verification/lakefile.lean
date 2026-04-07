import Lake
open Lake DSL

package «TwistorVerified» where
  leanOptions := #[
    ⟨`autoImplicit, false⟩
  ]

@[default_target]
lean_lib «TwistorVerified» where
  roots := #[`TwistorVerified]
