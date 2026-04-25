import { defineConfig, devices } from "@playwright/test";

export default defineConfig({
  testDir: "./e2e",
  fullyParallel: false,
  retries: 1,
  timeout: 30_000,
  use: {
    baseURL: "http://localhost:5173",
    // Backend URL injected as env var when starting the dev server
    headless: true,
    screenshot: "only-on-failure",
    video: "off",
  },
  projects: [
    {
      name: "chromium",
      use: { ...devices["Desktop Chrome"] },
    },
  ],
  // Start the Vite dev server automatically for E2E runs
  webServer: {
    command: "VITE_API_BASE_URL=http://localhost:8001 npm run dev",
    url: "http://localhost:5173",
    reuseExistingServer: true,
    timeout: 20_000,
  },
  // Visual snapshot directory
  snapshotDir: "./e2e/snapshots",
});
