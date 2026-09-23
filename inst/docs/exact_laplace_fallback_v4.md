# Exact discrete-Laplace fallback (1.4.0)

The v4 production fallback gives each of two independent noise authorities one
complete global-epsilon, global-L1-sensitivity discrete-Laplace draw per
coordinate. Each draw is an arbitrary-precision, unbounded integer. It is
reduced modulo `2^128` when added to its Ring128 share. The final sum uses the
existing fixed signed decoding and coordinate clamp.

Artifacts state `guarantee: "pure-dp-under-ideal-bits"` and
`randomness: "keyed-stream-computational"`. Under ideal independent bits,
exact geometric sampling with `q = exp(-epsilon/S)` gives a peer-noise PMF
proportional to `q^abs(z)`. Sensitivity at most S gives the epsilon likelihood
ratio bound. Conditional on either designated peer's simulatable source-share
view, the other independent complete draw supplies that guarantee. The other
draw, modular reduction, fixed signed decoding and clamp are post-processing.
Keyed HKDF/ChaCha20 replay additionally assumes pseudorandomness; the
variable-time sampler does not claim timing or transcript DP.

The mechanism and implementation deltas are zero. The Gaussian route retains
its positive declared delta. Staged grouped/LMM/GLMM/GEE/Cox source paths still
require the existing finite exact-GC validity gate and refuse zero-delta plans
before private material is accessed; their positive-delta route is unchanged.
The exact v4 fallback does not replace that separate source-validity gate.
Plans enforce `0 < epsilon <= 10000`, positive
integer S, `epsilon/S >= 2^-100`, `1 <= d <= 1000000` and W=128.

`representability_bound` and `wrap_bound` carry a positive outward scientific
decimal upper bound for `min(1, 2*d*q^B)`, B=`2^127`.
`wrap_bound_certified: true` certifies the event named by `wrap_bound_event`:
at least one peer draw is outside signed Ring128. It is a utility event and
is never charged to privacy delta. The actual statistic-plus-two-draw sum has
its separate `sum_wrap_bound` using T=`floor((B-U)/2)`, where U is the signed
public maximum scaled source coordinate. Neither certificate asserts that
wrapping is impossible. Accuracy readers reserve this actual sum-wrap bound
in addition to the exact convolution tail, preserving a positive outward
numeric bound when its decimal value is below double precision.

The backend is `independent_full_global_draw_convolution_ring128_v4`; the
sampler is `hkdf-sha256-chacha20-independent-full-draw-exact-geometric-v4`;
the stream domain is `dsVert/joint-dp/vector-convolution-private-stream/v4/`.
Plan, input, share and finalizer identifiers carry v4. Readers retain the
finite v3 identifiers and their positive implementation-delta semantics;
existing artifacts are never relabelled.

The seven-family `validate_dp_statistical_methods.R --quick` battery uses the
production v4 oracle with delta zero. `--ideal-sampler` remains an explicit
base-R comparison. The production record includes complete public synthetic
seeds, inputs, plans, exact draws and released coordinates for replay. This
short battery is not a DSLite or multi-peer deployment evaluation.
