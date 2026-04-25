/**
 * Unit tests for file upload validation logic.
 *
 * Tests the pure validateSelectedFiles function and the associated constants
 * without rendering any React components.
 */
import { describe, expect, it } from "vitest";
import {
  validateSelectedFiles,
  MAX_UPLOAD_BYTES,
  MAX_UPLOAD_MB,
  ALLOWED_UPLOAD_EXTENSIONS,
} from "./uploadValidation";

// ---------------------------------------------------------------------------
// Helpers
// ---------------------------------------------------------------------------

function makeFile(name: string, sizeBytes = 100): File {
  const buf = new Uint8Array(sizeBytes);
  return new File([buf], name, { type: "application/octet-stream" });
}

// ---------------------------------------------------------------------------
// Constants
// ---------------------------------------------------------------------------

describe("upload constants", () => {
  it("MAX_UPLOAD_BYTES is 20 MiB", () => {
    expect(MAX_UPLOAD_BYTES).toBe(20 * 1024 * 1024);
  });

  it("MAX_UPLOAD_MB is 20", () => {
    expect(MAX_UPLOAD_MB).toBe(20);
  });

  it("ALLOWED_UPLOAD_EXTENSIONS includes all expected types", () => {
    const required = [
      ".fasta", ".fa", ".fas",
      ".fastq", ".fq", ".gz",
      ".csv", ".tsv",
      ".xlsx", ".xls",
      ".txt",
    ];
    for (const ext of required) {
      expect(ALLOWED_UPLOAD_EXTENSIONS).toContain(ext);
    }
  });
});

// ---------------------------------------------------------------------------
// validateSelectedFiles — accepted files
// ---------------------------------------------------------------------------

describe("validateSelectedFiles — accepted files", () => {
  it.each([
    "sequences.fasta",
    "reads.fastq",
    "reads.fq",
    "archive.gz",
    "data.csv",
    "data.tsv",
    "workbook.xlsx",
    "workbook.xls",
    "notes.txt",
    "sequences.fa",
    "sequences.fas",
    "reads.fastq.gz",
    "reads.fq.gz",
  ])("accepts %s", (filename) => {
    const { accepted, rejected } = validateSelectedFiles([makeFile(filename)]);
    expect(accepted).toHaveLength(1);
    expect(rejected).toHaveLength(0);
  });

  it("accepts multiple valid files at once", () => {
    const files = [
      makeFile("counts.csv"),
      makeFile("metadata.txt"),
      makeFile("sequences.fasta"),
    ];
    const { accepted, rejected } = validateSelectedFiles(files);
    expect(accepted).toHaveLength(3);
    expect(rejected).toHaveLength(0);
  });

  it("accepts file at exactly the size limit", () => {
    const { accepted } = validateSelectedFiles([
      makeFile("edge.csv", MAX_UPLOAD_BYTES),
    ]);
    expect(accepted).toHaveLength(1);
  });
});

// ---------------------------------------------------------------------------
// validateSelectedFiles — rejected: unsupported extension
// ---------------------------------------------------------------------------

describe("validateSelectedFiles — rejected: unsupported extension", () => {
  it.each([
    "image.png",
    "document.pdf",
    "script.py",
    "archive.zip",
    "data.json",
    "config.yaml",
    "video.mp4",
  ])("rejects %s with unsupported-type error", (filename) => {
    const { accepted, rejected } = validateSelectedFiles([makeFile(filename)]);
    expect(accepted).toHaveLength(0);
    expect(rejected).toHaveLength(1);
    expect(rejected[0].status).toBe("failed");
    expect(rejected[0].error).toMatch(/unsupported type/i);
  });

  it("returns the original filename in the rejected entry", () => {
    const { rejected } = validateSelectedFiles([makeFile("bad.pdf")]);
    expect(rejected[0].name).toBe("bad.pdf");
  });
});

// ---------------------------------------------------------------------------
// validateSelectedFiles — rejected: oversized
// ---------------------------------------------------------------------------

describe("validateSelectedFiles — rejected: oversized", () => {
  it("rejects a file that exceeds MAX_UPLOAD_BYTES", () => {
    const oversized = makeFile("big.csv", MAX_UPLOAD_BYTES + 1);
    const { accepted, rejected } = validateSelectedFiles([oversized]);
    expect(accepted).toHaveLength(0);
    expect(rejected).toHaveLength(1);
    expect(rejected[0].status).toBe("failed");
    expect(rejected[0].error).toMatch(/exceeds/i);
    expect(rejected[0].error).toContain(`${MAX_UPLOAD_MB}MB`);
  });
});

// ---------------------------------------------------------------------------
// validateSelectedFiles — mixed batch
// ---------------------------------------------------------------------------

describe("validateSelectedFiles — mixed batch", () => {
  it("correctly partitions a mixed list of files", () => {
    const files = [
      makeFile("good.csv", 1000),
      makeFile("bad.pdf", 1000),
      makeFile("good.fasta", 1000),
      makeFile("toobig.csv", MAX_UPLOAD_BYTES + 1),
    ];
    const { accepted, rejected } = validateSelectedFiles(files);
    expect(accepted).toHaveLength(2);
    expect(accepted.map((f) => f.name)).toEqual(["good.csv", "good.fasta"]);
    expect(rejected).toHaveLength(2);
    expect(rejected.map((f) => f.name)).toEqual(["bad.pdf", "toobig.csv"]);
  });

  it("returns empty accepted and rejected for an empty input array", () => {
    const { accepted, rejected } = validateSelectedFiles([]);
    expect(accepted).toHaveLength(0);
    expect(rejected).toHaveLength(0);
  });
});
