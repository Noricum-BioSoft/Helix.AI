from __future__ import annotations

import uuid

from backend.contracts.ids import canonical_json, is_trace_id, new_id, new_trace_id, stable_hash, uuid7


def test_uuid7_is_version_7_variant_rfc_and_time_ordered():
    a, b = uuid7(), uuid7()
    assert a.version == 7 and b.version == 7
    assert a.variant == uuid.RFC_4122
    assert a.int >> 80 <= b.int >> 80  # timestamp prefix non-decreasing


def test_new_id_and_trace_id_formats():
    uuid.UUID(new_id())
    t = new_trace_id()
    assert is_trace_id(t) and t.startswith("orch_")
    assert not is_trace_id("orch_short") and not is_trace_id(new_id())


def test_stable_hash_is_order_independent_and_deterministic():
    assert stable_hash({"a": 1, "b": [1, 2]}) == stable_hash({"b": [1, 2], "a": 1})
    assert stable_hash({"a": 1}) != stable_hash({"a": 2})
    assert canonical_json({"b": 1, "a": "é"}) == '{"a":"é","b":1}'
