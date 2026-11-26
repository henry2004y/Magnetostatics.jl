# MagnetostaticsCore.jl

Lightweight, pure-Julia core package for analytical magnetostatic and inductive
electromagnetic fields with zero heavy solver dependencies (`LinearAlgebra`,
`StaticArrays`, and `PrecompileTools` only).

Provides:
- Abstract base types (`AbstractMagneticField`, `AbstractCurrentSource`, `AbstractSolver`)
- Field algebra (`NullField`, `CompositeField`, `ScaledField`, `TransformedField`,
  `translate_field`, `rotate_field`)
- Analytical primitives supporting both SI (`μ₀ / 4π`) and normalized code units (`k`):
  `UniformField`, `Dipole`, `LineDipole2D`, `InfiniteWire`, `WireSegment`,
  `CurrentLoop`, `CurrentLoopAnalytic`, `ToroidalField`, and
  `SphericalHarmonicField` with thread-safe `SphericalHarmonicWorkspace`.
