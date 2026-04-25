import { defineConfig } from "vitest/config";
import react from "@vitejs/plugin-react";

export default defineConfig({
  plugins: [react()],
  test: {
    // Node environment for pure TS/API tests; jsdom for React component tests.
    // Individual test files can override with @vitest-environment docblock.
    environment: "jsdom",
    include: ["src/**/*.test.{ts,tsx}"],
    globals: true,
    clearMocks: true,
    setupFiles: ["./src/test-setup.ts"],
  },
});
