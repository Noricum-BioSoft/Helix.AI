/**
 * File upload validation utilities.
 *
 * Extracted from App.tsx so they can be unit-tested independently of the
 * React component tree.
 */

import type { UploadedSessionFile } from "../services/helixApi";

export const MAX_UPLOAD_BYTES = 20 * 1024 * 1024;
export const MAX_UPLOAD_MB = 20;

export const ALLOWED_UPLOAD_EXTENSIONS = [
  // Sequence / FASTA / FASTQ
  ".fasta",
  ".fa",
  ".fas",
  ".fastq",
  ".fq",
  ".gz",
  // Tabular / structured data
  ".csv",
  ".tsv",
  ".xlsx",
  ".xls",
  // Generic text
  ".txt",
];

export interface ValidationResult {
  accepted: File[];
  rejected: UploadedSessionFile[];
}

/**
 * Partition *files* into accepted (pass validation) and rejected (fail) lists.
 * Rejection reasons: unsupported extension or file size exceeds MAX_UPLOAD_BYTES.
 */
export function validateSelectedFiles(files: File[]): ValidationResult {
  const accepted: File[] = [];
  const rejected: UploadedSessionFile[] = [];

  files.forEach((file) => {
    const lowerName = file.name.toLowerCase();
    const hasValidExt =
      ALLOWED_UPLOAD_EXTENSIONS.some((ext) => lowerName.endsWith(ext)) ||
      lowerName.endsWith(".fastq.gz") ||
      lowerName.endsWith(".fq.gz");

    if (!hasValidExt) {
      rejected.push({
        name: file.name,
        size: file.size,
        status: "failed",
        error:
          "Unsupported type — accepted: FASTA/FASTQ (.fa .fas .fasta .fastq .fq .gz), tabular (.csv .tsv .xlsx .xls), or .txt",
      });
      return;
    }

    if (file.size > MAX_UPLOAD_BYTES) {
      rejected.push({
        name: file.name,
        size: file.size,
        status: "failed",
        error: `Exceeds ${MAX_UPLOAD_MB}MB`,
      });
      return;
    }

    accepted.push(file);
  });

  return { accepted, rejected };
}
