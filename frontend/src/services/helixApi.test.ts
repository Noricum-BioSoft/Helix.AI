import { describe, expect, it, vi } from "vitest";

const axiosMocks = vi.hoisted(() => ({
  get: vi.fn(),
  post: vi.fn(),
}));

vi.mock("axios", () => ({
  default: axiosMocks,
}));

import { API_BASE_URL, helixApi, normalizeBaseUrl } from "./helixApi";

describe("helixApi regression", () => {
  describe("normalizeBaseUrl", () => {
    it("strips trailing slashes", () => {
      expect(normalizeBaseUrl("http://localhost:8001/")).toBe("http://localhost:8001");
      expect(normalizeBaseUrl("http://localhost:8001///")).toBe("http://localhost:8001");
    });

    it("returns undefined for empty / whitespace-only input", () => {
      expect(normalizeBaseUrl("")).toBeUndefined();
      expect(normalizeBaseUrl("   ")).toBeUndefined();
      expect(normalizeBaseUrl(undefined)).toBeUndefined();
    });

    it("preserves a clean URL unchanged", () => {
      expect(normalizeBaseUrl("http://localhost:8001")).toBe("http://localhost:8001");
      expect(normalizeBaseUrl("https://api.example.com")).toBe("https://api.example.com");
    });
  });

  it("API_BASE_URL is a non-empty string", () => {
    // API_BASE_URL is environment-dependent (localhost in dev, AWS in deployed
    // environments).  We only assert it is a non-empty string so the module
    // loaded correctly; normalizeBaseUrl unit tests above cover the logic.
    expect(typeof API_BASE_URL).toBe("string");
    expect(API_BASE_URL.length).toBeGreaterThan(0);
  });

  it("executePipelinePlan sends execute_plan payload", async () => {
    axiosMocks.post.mockResolvedValueOnce({ data: { ok: true } });

    const result = await helixApi.executePipelinePlan("run this workflow", "session-123");

    expect(axiosMocks.post).toHaveBeenCalledWith(
      `${API_BASE_URL}/execute`,
      {
        command: "run this workflow",
        session_id: "session-123",
        execute_plan: true,
      },
    );
    expect(result).toEqual({ ok: true });
  });
});
