from types import SimpleNamespace

import pytest
from fastapi import HTTPException

import backend.main as main


def _request_with_ip(ip: str):
    return SimpleNamespace(client=SimpleNamespace(host=ip))


def test_session_create_rate_limit_blocks_after_threshold(monkeypatch):
    monkeypatch.setattr(main, "MAX_SESSIONS_PER_HOUR", 2)
    main._hourly_session_create_counters.clear()
    req = _request_with_ip("203.0.113.10")

    main._check_and_increment_session_create_counter(req)
    main._check_and_increment_session_create_counter(req)

    with pytest.raises(HTTPException) as exc:
        main._check_and_increment_session_create_counter(req)
    assert exc.value.status_code == 429


def test_session_create_rate_limit_skips_localhost(monkeypatch):
    monkeypatch.setattr(main, "MAX_SESSIONS_PER_HOUR", 1)
    main._hourly_session_create_counters.clear()
    req = _request_with_ip("127.0.0.1")

    # Should never raise for localhost.
    main._check_and_increment_session_create_counter(req)
    main._check_and_increment_session_create_counter(req)
