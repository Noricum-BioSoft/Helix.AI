# Sandbox Baseline (P0)

This baseline defines minimum sandbox expectations for untrusted execution paths.

## Required Controls

- Isolated execution context for untrusted code/tool runs.
- Resource limits (CPU, memory, wall-clock timeout).
- Restricted filesystem scope to session-local working directories.
- Restricted network egress by default; explicit allowlist required.
- Structured execution audit logs with session and tool attribution.

## Validation Checklist

- [ ] Untrusted execution paths run via sandbox executor in staging.
- [ ] Timeout and memory limits are enforced and tested.
- [ ] Writes outside allowed workspace are blocked.
- [ ] Network-restricted mode is validated for benchmark scenarios.
- [ ] Audit logs include tool name, session id, and policy decision.

## P0 Decision

Release readiness is blocked if sandbox baseline controls are not met for
execution paths handling untrusted user-provided data or generated code.
