#!/usr/bin/env nextflow
/*
 * Helix.AI — ChIP-seq / CUT&RUN peak calling pipeline (PoC)
 *
 * Runs MACS3 peak calling inside a BioContainers Docker image.
 * This is a minimal workflow designed to prove the Nextflow execution tier.
 * For production use, migrate to nf-core/chipseq.
 *
 * Parameters
 * ----------
 * --treatment_bam   Path to IP/treatment BAM file
 * --control_bam     Path to IgG/input control BAM file (optional)
 * --genome_size     MACS3 genome size string: hs, mm, ce, dm, or integer (default: hs)
 * --peak_type       narrow | broad (default: narrow)
 * --outdir          Output directory (default: ./results)
 * --helix_job_id    Internal Helix job ID injected for weblog correlation
 */

nextflow.enable.dsl = 2

// ---------------------------------------------------------------------------
// Parameters with defaults
// ---------------------------------------------------------------------------
params.treatment_bam = ""
params.control_bam   = ""
params.genome_size   = "hs"
params.peak_type     = "narrow"
params.outdir        = "results"
params.helix_job_id  = ""

// ---------------------------------------------------------------------------
// Processes
// ---------------------------------------------------------------------------

process QC_FLAGSTAT {
    tag "$treatment_bam.simpleName"
    container "quay.io/biocontainers/samtools:1.18--h50ea8bc_1"
    cpus 2
    memory "2 GB"

    input:
    path treatment_bam

    output:
    path "flagstat.txt"

    script:
    """
    samtools flagstat ${treatment_bam} > flagstat.txt
    """
}

process PEAK_CALLING {
    tag "macs3"
    container "quay.io/biocontainers/macs3:3.0.2--py310h4ccef5b_0"
    cpus 2
    memory "4 GB"
    publishDir "${params.outdir}/peaks", mode: "copy"

    input:
    path treatment_bam
    path control_bam
    val  genome_size
    val  peak_type

    output:
    path "peaks/*"
    path "peaks/*_peaks.xls",  emit: xls
    path "peaks/*_peaks.bed",  emit: bed,  optional: true
    path "peaks/*_summits.bed", emit: summits, optional: true

    script:
    def control_arg = control_bam.name != "NO_FILE" ? "-c ${control_bam}" : ""
    def broad_arg   = peak_type == "broad" ? "--broad" : ""
    """
    mkdir -p peaks
    macs3 callpeak \\
        -t ${treatment_bam} \\
        ${control_arg} \\
        -g ${genome_size} \\
        -n peaks \\
        --outdir peaks \\
        ${broad_arg} \\
        --nomodel --extsize 200 \\
        -q 0.05
    """
}

process SUMMARISE_RESULTS {
    tag "summary"
    container "quay.io/biocontainers/python:3.11"
    publishDir "${params.outdir}", mode: "copy"

    input:
    path flagstat
    path peaks_xls

    output:
    path "summary.json"

    script:
    """
    python3 - <<'EOF'
import json, re, pathlib

summary = {}

# Parse flagstat
text = pathlib.Path("${flagstat}").read_text()
m = re.search(r"(\\d+) \\+ \\d+ mapped", text)
summary["mapped_reads"] = int(m.group(1)) if m else None

# Parse MACS3 XLS header for peak count
xls_text = pathlib.Path("${peaks_xls}").read_text()
data_lines = [l for l in xls_text.splitlines() if l and not l.startswith("#")]
# First line is header
summary["peak_count"] = max(0, len(data_lines) - 1)

summary["genome_size"] = "${params.genome_size}"
summary["peak_type"]   = "${params.peak_type}"
summary["helix_job_id"] = "${params.helix_job_id}"

pathlib.Path("summary.json").write_text(json.dumps(summary, indent=2))
print(json.dumps(summary, indent=2))
EOF
    """
}

// ---------------------------------------------------------------------------
// Workflow
// ---------------------------------------------------------------------------

workflow {
    if (!params.treatment_bam) {
        error "Please provide --treatment_bam <path>"
    }

    treatment_ch = Channel.fromPath(params.treatment_bam, checkIfExists: true)

    // Control BAM is optional — use a dummy file path when absent
    if (params.control_bam) {
        control_ch = Channel.fromPath(params.control_bam, checkIfExists: true)
    } else {
        control_ch = Channel.fromPath("NO_FILE", checkIfExists: false)
            .ifEmpty { file("NO_FILE") }
    }

    // QC
    flagstat_ch = QC_FLAGSTAT(treatment_ch)

    // Peak calling
    peaks_out = PEAK_CALLING(
        treatment_ch,
        control_ch.first(),
        params.genome_size,
        params.peak_type,
    )

    // Summary JSON
    SUMMARISE_RESULTS(flagstat_ch, peaks_out.xls)
}
