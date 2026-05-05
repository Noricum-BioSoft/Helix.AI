/**
 * Visual regression tests — Helix.AI frontend
 *
 * Uses Playwright's built-in screenshot comparison.  On first run these tests
 * create baseline snapshots in e2e/snapshots/.  On subsequent runs they diff
 * against the baseline and fail if any pixel threshold is exceeded.
 *
 * Update baselines:  npx playwright test --update-snapshots
 * Run:               npx playwright test e2e/visual.spec.ts
 */
import { test, expect } from "@playwright/test";

test.describe("Visual regression", () => {
  test("main page matches baseline", async ({ page }) => {
    await page.goto("/");
    // Wait for the UI to settle — server badge must be visible
    await expect(page.getByText(/server.*healthy/i)).toBeVisible({
      timeout: 10_000,
    });
    // Small extra wait for fonts / animations to settle
    await page.waitForTimeout(500);

    await expect(page).toHaveScreenshot("main-page.png", {
      fullPage: true,
      // Allow up to 1% pixel difference (anti-aliasing, font rendering)
      maxDiffPixelRatio: 0.01,
    });
  });

  test("empty conversation pane matches baseline", async ({ page }) => {
    await page.goto("/");
    await expect(page.getByText(/server.*healthy/i)).toBeVisible({
      timeout: 10_000,
    });
    await page.waitForTimeout(300);

    // The empty-state conversation pane is anchored by the CapabilityGrid
    // component (see frontend/src/components/CapabilityGrid.tsx) which
    // replaced the legacy "Run a command to see your prompts here" placeholder.
    const pane = page.getByText(
      /here's what helix\.ai can do/i,
    ).locator("..");
    await expect(pane).toHaveScreenshot("empty-conversation-pane.png", {
      maxDiffPixelRatio: 0.01,
    });
  });

  test("input area matches baseline", async ({ page }) => {
    await page.goto("/");
    await expect(page.getByText(/server.*healthy/i)).toBeVisible({
      timeout: 10_000,
    });
    await page.waitForTimeout(300);

    const inputArea = page
      .getByPlaceholder(/ask anything/i)
      .locator("../..");
    await expect(inputArea).toHaveScreenshot("input-area.png", {
      maxDiffPixelRatio: 0.01,
    });
  });
});
