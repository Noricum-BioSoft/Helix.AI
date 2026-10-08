"""
Secure Science pre-routing assessment (Phase 2).

Answers "can this *class* of action proceed?" before any provider is chosen.
Provider-specific checks ("may it proceed via THIS provider under THESE
conditions?") live in the Phase 4 ProviderAuthorizer, not here.

Outcomes combine as DENY > REQUIRE_REVIEW > ALLOW_WITH_APPROVAL > ALLOW; the
PolicyEnvelope is the intersection of every check's envelope. Fail closed:
a missing screening adapter is REQUIRE_REVIEW, never ALLOW.
"""
