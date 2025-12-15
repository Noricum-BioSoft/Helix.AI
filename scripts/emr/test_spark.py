#!/usr/bin/env python3
"""
Simple test script to verify Spark is working on EMR
"""

from pyspark.sql import SparkSession
import sys

def main():
    print("Starting Spark test...")
    
    try:
        spark = SparkSession.builder \
            .appName("Spark Test") \
            .getOrCreate()
        
        sc = spark.sparkContext
        sc.setLogLevel("WARN")
        
        print("✅ Spark session created successfully")
        
        # Simple test
        test_data = sc.parallelize([1, 2, 3, 4, 5])
        result = test_data.sum()
        
        print(f"✅ Test calculation: sum([1,2,3,4,5]) = {result}")
        print("✅ Spark is working correctly!")
        
        spark.stop()
        return 0
        
    except Exception as e:
        print(f"❌ Error: {e}")
        import traceback
        traceback.print_exc()
        return 1

if __name__ == "__main__":
    sys.exit(main())







