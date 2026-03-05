## Demo data generators

This folder contains small, reproducible dataset generators used by the Helix.AI demo scenarios.

### Bulk RNA-seq: Infection × Time (2×2 factorial)

Generates:
- `tgondii_counts.csv`: raw-ish integer counts (genes × samples)
- `tgondii_metadata.csv`: sample sheet with `infection_status` and `time_point`

Run:

```bash
python3 scripts/demo_data/generate_tgondii_factorial_demo.py --out-dir /tmp/tgondii_demo
```

Optional upload to the public demo prefix used by the frontend:

```bash
aws s3 cp /tmp/tgondii_demo/tgondii_counts.csv s3://noricum-ngs-data/demo/rnaseq/tgondii_counts.csv
aws s3 cp /tmp/tgondii_demo/tgondii_metadata.csv s3://noricum-ngs-data/demo/rnaseq/tgondii_metadata.csv
```

