from backend.session_param_extractor import get_fastqc_inputs_from_session


def test_fastqc_explicit_r1_r2_fastq_not_misclassified_as_fasta():
    command = (
        "Run FastQC. Forward reads: s3://b/r1.fastq "
        "Reverse reads: s3://b/r2.fastq Output: s3://b/out/"
    )
    r1, r2, out, err = get_fastqc_inputs_from_session({}, command)
    assert err is None
    assert r1 is None and r2 is None and out is None


def test_fastqc_merged_fasta_is_rejected_with_clear_error():
    command = "Run FastQC on merged.fasta"
    r1, r2, out, err = get_fastqc_inputs_from_session({}, command)
    assert r1 is None and r2 is None and out is None
    assert isinstance(err, str) and "cannot run on merged reads" in err.lower()
