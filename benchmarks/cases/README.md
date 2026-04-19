# Executable Benchmark Cases

Each case file defines one benchmark scenario derived from `benchmarks/use_case_catalog.yaml`.

## Required Top-Level Fields

- `id`
- `title`
- `source_refs`
- `family`
- `input_fixtures`
- `prompt`
- `expected_route`
- `forbidden_routes`
- `forbidden_fallbacks`
- `required_steps`
- `required_outputs`
- `safety_requirements`
- `provenance_requirements`

## Scorer Contract

`benchmarks/scoring/use_case_gate_scorer.py` expects run observations with:

- `observed_route`: selected route family for this prompt.
- `fallback_targets`: list of fallback route targets used (if any).
- `observed_steps`: normalized list of steps actually executed.
- `observed_output_types`: normalized output classes produced.
- `observed_safety_flags`: safety flags surfaced in response.
- `observed_provenance_signals`: provenance/reproducibility signals emitted.
- `replay_status`: one of `replay_match`, `replay_within_tolerance`, or failure value.

The gate fails if a forbidden route/fallback occurs, provenance is incomplete,
replay fails, or the weighted total score drops below threshold.

