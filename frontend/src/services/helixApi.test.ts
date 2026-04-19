import { describe, expect, it, vi } from "vitest";

const axiosMocks = vi.hoisted(() => ({
  get: vi.fn(),
  post: vi.fn(),
}));

vi.mock("axios", () => ({
  default: axiosMocks,
}));

import { API_BASE_URL, helixApi } from "./helixApi";

describe("helixApi regression", () => {
  it("uses local backend base URL by default in tests", () => {
    expect(API_BASE_URL).toBe("http://localhost:8001");
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
