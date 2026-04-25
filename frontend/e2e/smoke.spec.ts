/**
 * E2E smoke tests — Helix.AI frontend
 *
 * These tests run against the live Vite dev server (started automatically by
 * playwright.config.ts webServer).  They require the backend to be running on
 * http://localhost:8001.
 *
 * Run:  npx playwright test e2e/smoke.spec.ts
 */
import { test, expect } from "@playwright/test";

test.describe("Page load", () => {
  test("renders the main prompt heading", async ({ page }) => {
    await page.goto("/");
    await expect(
      page.getByText("What data do you want to process today?"),
    ).toBeVisible({ timeout: 10_000 });
  });

  test("shows Server: healthy badge when backend is up", async ({ page }) => {
    await page.goto("/");
    await expect(page.getByText(/server.*healthy/i)).toBeVisible({
      timeout: 10_000,
    });
  });

  test("action buttons are present", async ({ page }) => {
    await page.goto("/");
    await expect(page.getByRole("button", { name: /examples/i })).toBeVisible();
    await expect(page.getByRole("button", { name: /demo/i })).toBeVisible();
    await expect(page.getByRole("button", { name: /jobs/i })).toBeVisible();
  });
});

test.describe("Command submission", () => {
  test("submits a command and receives a non-empty response", async ({
    page,
  }) => {
    await page.goto("/");

    // Wait for the server healthy badge before submitting
    await expect(page.getByText(/server.*healthy/i)).toBeVisible({
      timeout: 10_000,
    });

    const input = page.getByPlaceholder(/ask anything/i);
    await input.fill("What tools do you have available?");

    // Run button should become enabled once there is text
    const runBtn = page.getByRole("button", { name: /^run$/i });
    await expect(runBtn).toBeEnabled({ timeout: 3_000 });
    await runBtn.click();

    // A response card should appear within 30 s
    await expect(
      page.locator(".card, [data-testid='response']").first(),
    ).toBeVisible({ timeout: 30_000 });
  });

  test("conversation section shows user prompt after submission", async ({
    page,
  }) => {
    await page.goto("/");
    await expect(page.getByText(/server.*healthy/i)).toBeVisible({
      timeout: 10_000,
    });

    const input = page.getByPlaceholder(/ask anything/i);
    await input.fill("List available tools");
    await page.getByRole("button", { name: /^run$/i }).click();

    // The user's message should appear in the conversation pane
    await expect(page.getByText("List available tools")).toBeVisible({
      timeout: 5_000,
    });
  });
});

test.describe("File upload UI", () => {
  test("file input accepts .csv files", async ({ page }) => {
    await page.goto("/");
    await expect(page.getByText(/server.*healthy/i)).toBeVisible({
      timeout: 10_000,
    });

    // Trigger file chooser via the attachment button
    const [fileChooser] = await Promise.all([
      page.waitForEvent("filechooser"),
      page.locator("input[type='file']").first().evaluate((el) => {
        (el as HTMLInputElement).dispatchEvent(new MouseEvent("click"));
      }),
    ]);

    // Set a valid CSV file
    await fileChooser.setFiles({
      name: "test_data.csv",
      mimeType: "text/csv",
      buffer: Buffer.from("gene,logFC,padj\nGENE1,2.1,0.001\nGENE2,-1.3,0.04"),
    });

    // File should appear in the upload list
    await expect(page.getByText("test_data.csv")).toBeVisible({
      timeout: 5_000,
    });
  });
});
