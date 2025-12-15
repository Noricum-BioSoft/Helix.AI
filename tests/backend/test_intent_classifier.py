from backend.intent_classifier import classify_intent


def test_intent_classifier_detects_qa():
    d = classify_intent("What is a FASTQ file?")
    assert d.intent == "qa"


def test_intent_classifier_detects_execute_from_verbs():
    d = classify_intent("Run FastQC on my reads")
    assert d.intent == "execute"


def test_intent_classifier_detects_execute_from_s3_uri():
    d = classify_intent("analyze s3://my-bucket/sample_R1.fastq and s3://my-bucket/sample_R2.fastq")
    assert d.intent == "execute"


