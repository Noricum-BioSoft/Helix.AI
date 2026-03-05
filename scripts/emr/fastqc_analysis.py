#!/usr/bin/env python3
"""
PySpark FASTQ Quality Analysis
Generates FastQC-like statistics for paired-end FASTQ files.

Usage:
    spark-submit fastqc_analysis.py \
        --input-r1 s3://bucket/path/mate_R1.fq \
        --input-r2 s3://bucket/path/mate_R2.fq \
        --output s3://bucket/output/
"""

import argparse
import json
from collections import Counter
from pyspark.sql import SparkSession
from pyspark.sql.functions import col, length, regexp_extract, when, lit
from pyspark.sql.types import StructType, StructField, StringType, IntegerType, DoubleType
import sys


def parse_fastq_record(lines):
    """Parse a 4-line FASTQ record into components."""
    if len(lines) != 4:
        return None
    header = lines[0]
    sequence = lines[1]
    separator = lines[2]
    quality = lines[3]
    
    if not header.startswith('@') or not separator.startswith('+'):
        return None
    
    return {
        'header': header[1:],  # Remove @
        'sequence': sequence,
        'quality': quality,
        'length': len(sequence)
    }


def calculate_phred_score(quality_char, offset=33):
    """Convert quality character to Phred score."""
    return ord(quality_char) - offset


def calculate_gc_content(sequence):
    """Calculate GC content percentage."""
    if not sequence:
        return 0.0
    gc_count = sequence.upper().count('G') + sequence.upper().count('C')
    return (gc_count / len(sequence)) * 100.0


def main():
    parser = argparse.ArgumentParser(description='FASTQ Quality Analysis with PySpark')
    parser.add_argument('--input-r1', required=True, help='S3 path to R1 FASTQ file')
    parser.add_argument('--input-r2', required=True, help='S3 path to R2 FASTQ file')
    parser.add_argument('--output', required=True, help='S3 output path for results')
    parser.add_argument('--sample-size', type=int, default=1000000, 
                       help='Sample size for some statistics (default: 1M reads)')
    
    args = parser.parse_args()
    
    # Initialize Spark
    spark = SparkSession.builder \
        .appName("FASTQ Quality Analysis") \
        .config("spark.sql.adaptive.enabled", "true") \
        .config("spark.sql.adaptive.coalescePartitions.enabled", "true") \
        .getOrCreate()
    
    sc = spark.sparkContext
    sc.setLogLevel("WARN")
    
    print(f"Processing R1: {args.input_r1}")
    print(f"Processing R2: {args.input_r2}")
    
    # Read FASTQ files - use wholeTextFiles for better handling of large files
    print("Reading FASTQ files from S3...")
    r1_text = sc.wholeTextFiles(args.input_r1).values().first()
    r2_text = sc.wholeTextFiles(args.input_r2).values().first()
    
    # Split into lines and process
    def process_fastq_text(text):
        """Process FASTQ text into records."""
        lines = text.split('\n')
        records = []
        for i in range(0, len(lines) - 3, 4):
            record = parse_fastq_record(lines[i:i+4])
            if record:
                records.append(record)
        return records
    
    # Process R1 and R2
    print("Parsing R1 reads...")
    r1_records_list = process_fastq_text(r1_text)
    print(f"Parsed {len(r1_records_list)} R1 records")
    
    print("Parsing R2 reads...")
    r2_records_list = process_fastq_text(r2_text)
    print(f"Parsed {len(r2_records_list)} R2 records")
    
    # Convert to RDDs for distributed processing
    r1_records = sc.parallelize(r1_records_list)
    r2_records = sc.parallelize(r2_records_list)
    
    # Convert to DataFrame for easier analysis
    schema = StructType([
        StructField("header", StringType(), True),
        StructField("sequence", StringType(), True),
        StructField("quality", StringType(), True),
        StructField("length", IntegerType(), True)
    ])
    
    # Create DataFrames
    r1_df = spark.createDataFrame(r1_records, schema)
    r2_df = spark.createDataFrame(r2_records, schema)
    
    print(f"R1 reads parsed: {r1_df.count()}")
    print(f"R2 reads parsed: {r2_df.count()}")
    
    # Calculate statistics
    results = {}
    
    # Basic Statistics
    print("Calculating basic statistics...")
    r1_total = r1_df.count()
    r2_total = r2_df.count()
    r1_total_bases = r1_df.agg({"length": "sum"}).collect()[0][0]
    r2_total_bases = r2_df.agg({"length": "sum"}).collect()[0][0]
    
    results['basic_statistics'] = {
        'r1': {
            'total_sequences': r1_total,
            'total_bases': r1_total_bases,
            'average_length': r1_total_bases / r1_total if r1_total > 0 else 0
        },
        'r2': {
            'total_sequences': r2_total,
            'total_bases': r2_total_bases,
            'average_length': r2_total_bases / r2_total if r2_total > 0 else 0
        }
    }
    
    # Sequence Length Distribution
    print("Calculating length distribution...")
    r1_length_dist = r1_df.groupBy("length").count().orderBy("length").collect()
    r2_length_dist = r2_df.groupBy("length").count().orderBy("length").collect()
    
    results['length_distribution'] = {
        'r1': {row.length: row['count'] for row in r1_length_dist},
        'r2': {row.length: row['count'] for row in r2_length_dist}
    }
    
    # GC Content (sample-based for performance)
    print("Calculating GC content...")
    r1_sample = r1_df.sample(False, min(1.0, args.sample_size / r1_total)) if r1_total > 0 else r1_df
    r2_sample = r2_df.sample(False, min(1.0, args.sample_size / r2_total)) if r2_total > 0 else r2_df
    
    def calc_gc_udf(seq):
        return calculate_gc_content(seq) if seq else 0.0
    
    from pyspark.sql.functions import udf
    from pyspark.sql.types import DoubleType
    
    gc_udf = udf(calc_gc_udf, DoubleType())
    r1_gc = r1_sample.withColumn("gc_content", gc_udf(col("sequence"))).select("gc_content")
    r2_gc = r2_sample.withColumn("gc_content", gc_udf(col("sequence"))).select("gc_content")
    
    r1_gc_stats = r1_gc.describe("gc_content").collect()
    r2_gc_stats = r2_gc.describe("gc_content").collect()
    
    results['gc_content'] = {
        'r1': {
            'mean': float(r1_gc_stats[1][1]) if len(r1_gc_stats) > 1 else 0.0,
            'std': float(r1_gc_stats[2][1]) if len(r1_gc_stats) > 2 else 0.0
        },
        'r2': {
            'mean': float(r2_gc_stats[1][1]) if len(r2_gc_stats) > 1 else 0.0,
            'std': float(r2_gc_stats[2][1]) if len(r2_gc_stats) > 2 else 0.0
        }
    }
    
    # Per-base quality scores (sample-based)
    print("Calculating per-base quality scores...")
    # This is simplified - full implementation would calculate quality per position
    results['per_base_quality'] = {
        'note': 'Per-base quality requires position-by-position analysis',
        'summary': 'See full FastQC report for detailed per-base metrics'
    }
    
    # Save results
    print(f"Saving results to {args.output}...")
    results_json = json.dumps(results, indent=2)
    
    # Write to S3
    sc.parallelize([results_json]).saveAsTextFile(args.output)
    
    print("\n" + "="*60)
    print("FASTQ Quality Analysis Complete!")
    print("="*60)
    print(f"\nResults saved to: {args.output}")
    print(f"\nSummary:")
    print(f"  R1: {r1_total:,} reads, {r1_total_bases:,} bases")
    print(f"  R2: {r2_total:,} reads, {r2_total_bases:,} bases")
    print(f"  R1 Average Length: {results['basic_statistics']['r1']['average_length']:.1f} bp")
    print(f"  R2 Average Length: {results['basic_statistics']['r2']['average_length']:.1f} bp")
    print(f"  R1 GC Content: {results['gc_content']['r1']['mean']:.2f}%")
    print(f"  R2 GC Content: {results['gc_content']['r2']['mean']:.2f}%")
    
    spark.stop()


if __name__ == "__main__":
    main()

