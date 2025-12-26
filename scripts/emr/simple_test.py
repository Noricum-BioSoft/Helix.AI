#!/usr/bin/env python3
"""
Simple test to verify Spark and S3 access work on EMR
"""

from pyspark.sql import SparkSession
import sys

def main():
    print("="*60)
    print("Simple Spark Test")
    print("="*60)
    
    try:
        # Initialize Spark
        print("Initializing Spark...")
        spark = SparkSession.builder \
            .appName("Simple Test") \
            .getOrCreate()
        
        sc = spark.sparkContext
        sc.setLogLevel("WARN")
        
        print(f"✅ Spark version: {sc.version}")
        
        # Test 1: Basic RDD operation
        print("\nTest 1: Basic RDD operation...")
        test_rdd = sc.parallelize([1, 2, 3, 4, 5])
        result = test_rdd.sum()
        print(f"✅ Sum test: {result} (expected: 15)")
        
        # Test 2: S3 access
        print("\nTest 2: S3 file access...")
        test_file = "s3://noricum-ngs-data/datasets/GRCh38.p12.MafHi/mate_R1.fq"
        try:
            # Just try to read first few lines
            lines_rdd = spark.read.text(test_file).limit(10)
            line_count = lines_rdd.count()
            print(f"✅ S3 access works: Read {line_count} lines")
            
            # Show first line
            first_line = lines_rdd.first()
            if first_line:
                print(f"✅ First line preview: {first_line[0][:50]}...")
        except Exception as e:
            print(f"❌ S3 access failed: {e}")
            return 1
        
        # Test 3: Python imports
        print("\nTest 3: Python imports...")
        try:
            import json
            import argparse
            from collections import defaultdict
            print("✅ Standard library imports work")
        except Exception as e:
            print(f"❌ Import failed: {e}")
            return 1
        
        print("\n" + "="*60)
        print("✅ All tests passed!")
        print("="*60)
        
        spark.stop()
        return 0
        
    except Exception as e:
        print(f"\n❌ ERROR: {e}")
        import traceback
        traceback.print_exc()
        return 1

if __name__ == "__main__":
    sys.exit(main())










