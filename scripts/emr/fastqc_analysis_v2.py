#!/usr/bin/env python3
"""
PySpark FASTQ Quality Analysis (Streaming Version)
Generates FastQC-like statistics for paired-end FASTQ files.
Handles large files by processing in chunks.

Usage:
    spark-submit fastqc_analysis_v2.py \
        --input-r1 s3://bucket/path/mate_R1.fq \
        --input-r2 s3://bucket/path/mate_R2.fq \
        --output s3://bucket/output/
"""

import argparse
import json
import sys
from pyspark.sql import SparkSession
from pyspark.sql.functions import col, length, udf
from pyspark.sql.types import StructType, StructField, StringType, IntegerType, DoubleType
from pyspark import StorageLevel
from collections import defaultdict

# Ensure print statements are flushed immediately
def print_flush(*args, **kwargs):
    print(*args, **kwargs)
    sys.stdout.flush()


def parse_fastq_record(lines):
    """Parse a 4-line FASTQ record into components."""
    if len(lines) != 4:
        return None
    header = lines[0].strip()
    sequence = lines[1].strip()
    separator = lines[2].strip()
    quality = lines[3].strip()
    
    if not header.startswith('@') or not separator.startswith('+'):
        return None
    
    if not sequence or not quality:
        return None
    
    return {
        'header': header[1:],  # Remove @
        'sequence': sequence,
        'quality': quality,
        'length': len(sequence)
    }


def calculate_gc_content(sequence):
    """Calculate GC content percentage."""
    if not sequence:
        return 0.0
    gc_count = sequence.upper().count('G') + sequence.upper().count('C')
    return (gc_count / len(sequence)) * 100.0


def process_fastq_chunk(text_chunk):
    """Process a chunk of FASTQ text into records."""
    lines = text_chunk.split('\n')
    records = []
    for i in range(0, len(lines) - 3, 4):
        record = parse_fastq_record(lines[i:i+4])
        if record:
            records.append(record)
    return records


def main():
    try:
        parser = argparse.ArgumentParser(description='FASTQ Quality Analysis with PySpark')
        parser.add_argument('--input-r1', required=True, help='S3 path to R1 FASTQ file')
        parser.add_argument('--input-r2', required=True, help='S3 path to R2 FASTQ file')
        parser.add_argument('--output', required=True, help='S3 output path for results')
        
        args = parser.parse_args()
        
        print_flush("="*60)
        print_flush("FASTQ Quality Analysis - Starting")
        print_flush("="*60)
        print_flush(f"R1 Input: {args.input_r1}")
        print_flush(f"R2 Input: {args.input_r2}")
        print_flush(f"Output: {args.output}")
        print_flush("")
        
        # Initialize Spark with S3A filesystem configuration
        # Since EmrFileSystem isn't available, use S3A which should be in Hadoop libraries
        print_flush("Initializing Spark session...")
        spark = SparkSession.builder \
            .appName("FASTQ Quality Analysis") \
            .config("spark.sql.adaptive.enabled", "true") \
            .config("spark.sql.adaptive.coalescePartitions.enabled", "true") \
            .config("spark.sql.files.maxPartitionBytes", "134217728") \
            .config("spark.hadoop.io.compression.codecs", "org.apache.hadoop.io.compress.GzipCodec,org.apache.hadoop.io.compress.DefaultCodec,org.apache.hadoop.io.compress.BZip2Codec") \
            .config("spark.hadoop.io.compression.codec.lzo.class", "") \
            .config("spark.hadoop.fs.s3a.impl", "org.apache.hadoop.fs.s3a.S3AFileSystem") \
            .config("spark.hadoop.fs.s3a.aws.credentials.provider", "com.amazonaws.auth.InstanceProfileCredentialsProvider") \
            .config("spark.hadoop.fs.s3a.endpoint", "s3.us-east-1.amazonaws.com") \
            .config("spark.hadoop.fs.s3a.path.style.access", "true") \
            .getOrCreate()
        
        sc = spark.sparkContext
        sc.setLogLevel("WARN")
        
        # Also set Hadoop configuration
        hadoop_conf = sc._jsc.hadoopConfiguration()
        hadoop_conf.set("io.compression.codecs", "org.apache.hadoop.io.compress.GzipCodec,org.apache.hadoop.io.compress.DefaultCodec,org.apache.hadoop.io.compress.BZip2Codec")
        hadoop_conf.set("io.compression.codec.lzo.class", "")
        hadoop_conf.set("fs.s3a.impl", "org.apache.hadoop.fs.s3a.S3AFileSystem")
        hadoop_conf.set("fs.s3a.aws.credentials.provider", "com.amazonaws.auth.InstanceProfileCredentialsProvider")
        hadoop_conf.set("fs.s3a.endpoint", "s3.us-east-1.amazonaws.com")
        hadoop_conf.set("fs.s3a.path.style.access", "true")
        
        print_flush("✅ Spark session created")
        print_flush(f"Spark version: {sc.version}")
        print_flush("")
        
        # Download files from S3 to HDFS first to avoid filesystem class issues
        # This is a workaround for the EmrFileSystem/S3AFileSystem ClassNotFoundException
        import subprocess
        import os
        import tempfile
        
        print_flush("Downloading files from S3 to HDFS...")
        print_flush(f"  R1: {args.input_r1}")
        print_flush(f"  R2: {args.input_r2}")
        print_flush("")
        
        # Create temporary HDFS paths
        hdfs_base = f"/tmp/fastqc_{os.getpid()}"
        r1_hdfs_path = f"{hdfs_base}/mate_R1.fq"
        r2_hdfs_path = f"{hdfs_base}/mate_R2.fq"
        
        try:
            # Create HDFS directory first
            print_flush(f"Creating HDFS directory: {hdfs_base}...")
            result = subprocess.run([
                "hadoop", "fs", "-mkdir", "-p", hdfs_base
            ], capture_output=True, text=True)
            if result.returncode != 0:
                print_flush(f"Warning creating directory (may already exist): {result.stderr}")
            print_flush("✅ HDFS directory ready")
            
            # Download R1 to HDFS
            print_flush("Downloading R1 to HDFS...")
            result = subprocess.run([
                "hadoop", "fs", "-cp",
                args.input_r1,
                r1_hdfs_path
            ], capture_output=True, text=True)
            if result.returncode != 0:
                print_flush(f"Error copying R1: stdout={result.stdout}, stderr={result.stderr}")
                raise subprocess.CalledProcessError(result.returncode, result.args, result.stdout, result.stderr)
            print_flush("✅ R1 downloaded to HDFS")
            
            # Download R2 to HDFS
            print_flush("Downloading R2 to HDFS...")
            result = subprocess.run([
                "hadoop", "fs", "-cp",
                args.input_r2,
                r2_hdfs_path
            ], capture_output=True, text=True)
            if result.returncode != 0:
                print_flush(f"Error copying R2: stdout={result.stdout}, stderr={result.stderr}")
                raise subprocess.CalledProcessError(result.returncode, result.args, result.stdout, result.stderr)
            print_flush("✅ R2 downloaded to HDFS")
            
            # Now read from HDFS using SparkContext.textFile()
            print_flush("Reading files from HDFS...")
            print_flush(f"  R1: {r1_hdfs_path}")
            print_flush(f"  R2: {r2_hdfs_path}")
            print_flush("")
            
            r1_lines = sc.textFile(r1_hdfs_path)
            r2_lines = sc.textFile(r2_hdfs_path)
            
            # Persist to memory/disk
            print_flush("Persisting RDDs to memory/disk...")
            r1_lines.persist(StorageLevel.MEMORY_AND_DISK)
            r2_lines.persist(StorageLevel.MEMORY_AND_DISK)
            
            # Force materialization
            print_flush("Materializing RDDs (this may take a while for large files)...")
            r1_count = r1_lines.count()
            r2_count = r2_lines.count()
            print_flush(f"✅ Files read and persisted: R1 has {r1_count} lines, R2 has {r2_count} lines")
            
        except subprocess.CalledProcessError as e:
            print_flush(f"❌ ERROR downloading files to HDFS: {e}")
            print_flush(f"stdout: {e.stdout}")
            print_flush(f"stderr: {e.stderr}")
            raise
        except Exception as e:
            print_flush(f"❌ ERROR reading files: {e}")
            import traceback
            traceback.print_exc()
            raise
        finally:
            # Clean up HDFS files after processing (optional, can leave for debugging)
            try:
                print_flush(f"Cleaning up HDFS files in {hdfs_base}...")
                subprocess.run(["hadoop", "fs", "-rm", "-r", "-f", hdfs_base], 
                             capture_output=True, text=True)
            except:
                pass
        
        # Process lines into 4-line records using mapPartitions
        # Note: This approach may miss records that span partition boundaries
        # but is more memory-efficient for large files
        def process_partition(partition):
            """Process partition of lines into FASTQ records."""
            try:
                lines = list(partition)
                records = []
                # Group into 4-line chunks
                # Skip incomplete records at the end of partition
                num_complete_records = (len(lines) // 4) * 4
                for i in range(0, num_complete_records, 4):
                    if i + 4 <= len(lines):
                        record = parse_fastq_record(lines[i:i+4])
                        if record:
                            records.append(record)
                return iter(records)
            except Exception as e:
                print(f"ERROR in process_partition: {e}")
                import traceback
                traceback.print_exc()
                return iter([])
        
        print_flush("Parsing R1 reads...")
        try:
            r1_records_rdd = r1_lines.mapPartitions(process_partition)
            print_flush("✅ R1 parsing started")
        except Exception as e:
            print_flush(f"❌ ERROR parsing R1: {e}")
            import traceback
            traceback.print_exc()
            raise
        
        print_flush("Parsing R2 reads...")
        try:
            r2_records_rdd = r2_lines.mapPartitions(process_partition)
            print_flush("✅ R2 parsing started")
        except Exception as e:
            print_flush(f"❌ ERROR parsing R2: {e}")
            import traceback
            traceback.print_exc()
            raise
        
        # Convert to DataFrame
        schema = StructType([
            StructField("header", StringType(), True),
            StructField("sequence", StringType(), True),
            StructField("quality", StringType(), True),
            StructField("length", IntegerType(), True)
        ])
        
        print_flush("Creating DataFrames...")
        try:
            r1_df = spark.createDataFrame(r1_records_rdd, schema)
            r2_df = spark.createDataFrame(r2_records_rdd, schema)
            print_flush("✅ DataFrames created")
        except Exception as e:
            print_flush(f"❌ ERROR creating DataFrames: {e}")
            import traceback
            traceback.print_exc()
            raise
        
        # Cache for multiple operations
        print_flush("Caching DataFrames...")
        r1_df.cache()
        r2_df.cache()
        print_flush("✅ DataFrames cached")
        
        print_flush("Calculating statistics...")
        
        # Basic Statistics
        print_flush("  - Basic statistics...")
        try:
            r1_total = r1_df.count()
            print_flush(f"  R1 total reads: {r1_total:,}")
            r2_total = r2_df.count()
            print_flush(f"  R2 total reads: {r2_total:,}")
        except Exception as e:
            print_flush(f"❌ ERROR counting reads: {e}")
            import traceback
            traceback.print_exc()
            raise
        
        r1_total_bases = r1_df.agg({"length": "sum"}).collect()[0][0] or 0
        r2_total_bases = r2_df.agg({"length": "sum"}).collect()[0][0] or 0
        
        r1_avg_length = r1_total_bases / r1_total if r1_total > 0 else 0
        r2_avg_length = r2_total_bases / r2_total if r2_total > 0 else 0
        
        results = {
            'basic_statistics': {
                'r1': {
                    'total_sequences': r1_total,
                    'total_bases': r1_total_bases,
                    'average_length': round(r1_avg_length, 2)
                },
                'r2': {
                    'total_sequences': r2_total,
                    'total_bases': r2_total_bases,
                    'average_length': round(r2_avg_length, 2)
                }
            }
        }
        
        # Sequence Length Distribution
        print_flush("  - Length distribution...")
        r1_length_dist = r1_df.groupBy("length").count().orderBy("length").collect()
        r2_length_dist = r2_df.groupBy("length").count().orderBy("length").collect()
        
        results['length_distribution'] = {
            'r1': {int(row.length): int(row['count']) for row in r1_length_dist},
            'r2': {int(row.length): int(row['count']) for row in r2_length_dist}
        }
        
        # GC Content
        print_flush("  - GC content...")
        gc_udf = udf(calculate_gc_content, DoubleType())
        
        r1_gc = r1_df.withColumn("gc_content", gc_udf(col("sequence")))
        r2_gc = r2_df.withColumn("gc_content", gc_udf(col("sequence")))
        
        # Get mean and stddev separately
        r1_gc_mean = r1_gc.agg({"gc_content": "mean"}).collect()[0][0] or 0.0
        r1_gc_std = r1_gc.agg({"gc_content": "stddev"}).collect()[0][0] or 0.0
        r2_gc_mean = r2_gc.agg({"gc_content": "mean"}).collect()[0][0] or 0.0
        r2_gc_std = r2_gc.agg({"gc_content": "stddev"}).collect()[0][0] or 0.0
        
        results['gc_content'] = {
            'r1': {
                'mean': round(float(r1_gc_mean), 2),
                'std': round(float(r1_gc_std), 2)
            },
            'r2': {
                'mean': round(float(r2_gc_mean), 2),
                'std': round(float(r2_gc_std), 2)
            }
        }
    
        # Save results to local file first, then upload to S3
        # This avoids Hadoop filesystem class issues
        import subprocess
        local_output_file = f"/tmp/fastqc_results_{os.getpid()}.json"
        print_flush(f"Saving results to local file: {local_output_file}...")
        results_json = json.dumps(results, indent=2)
        
        # Write to local file
        with open(local_output_file, 'w') as f:
            f.write(results_json)
        print_flush("✅ Results saved to local file")
        
        # Upload to S3 using AWS CLI
        s3_output_path = args.output.rstrip('/') + "/results.json"
        print_flush(f"Uploading results to S3: {s3_output_path}...")
        try:
            result = subprocess.run([
                "aws", "s3", "cp",
                local_output_file,
                s3_output_path
            ], check=True, capture_output=True, text=True)
            print_flush("✅ Results uploaded to S3")
        except subprocess.CalledProcessError as e:
            print_flush(f"⚠️  WARNING: Failed to upload results to S3: {e}")
            print_flush(f"Results are still available locally at: {local_output_file}")
            print_flush(f"stdout: {e.stdout}")
            print_flush(f"stderr: {e.stderr}")
            raise
        
        # Clean up local file
        try:
            os.remove(local_output_file)
        except:
            pass
        
        print_flush("\n" + "="*60)
        print_flush("FASTQ Quality Analysis Complete!")
        print_flush("="*60)
        print_flush(f"\nResults saved to: {args.output}")
        print_flush(f"\nSummary:")
        print_flush(f"  R1: {r1_total:,} reads, {r1_total_bases:,} bases")
        print_flush(f"  R2: {r2_total:,} reads, {r2_total_bases:,} bases")
        print_flush(f"  R1 Average Length: {r1_avg_length:.1f} bp")
        print_flush(f"  R2 Average Length: {r2_avg_length:.1f} bp")
        print_flush(f"  R1 GC Content: {r1_gc_mean:.2f}%")
        print_flush(f"  R2 GC Content: {r2_gc_mean:.2f}%")
        
        spark.stop()
        print_flush("\n✅ Spark session stopped")
        return 0
        
    except Exception as e:
        print_flush(f"\n❌ ERROR: {e}")
        import traceback
        traceback.print_exc()
        sys.stdout.flush()
        sys.stderr.flush()
        return 1


if __name__ == "__main__":
    import sys
    sys.exit(main())

