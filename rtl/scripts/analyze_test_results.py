#!/usr/bin/env python3
"""
Script Name: analyze_test_results.py
Purpose: Analyze unit test results from output/test_results/ directory
Usage: python scripts/analyze_test_results.py
Author: AV2 RTL Team
Date: 2026-02-12
"""

import os
import sys
import re
from pathlib import Path

def analyze_test_log(log_file):
    """Analyze a single test log file"""
    with open(log_file, 'r') as f:
        content = f.read()
    
    results = {
        'passed': 0,
        'failed': 0,
        'errors': [],
        'warnings': []
    }
    
    # Check for test completion
    if 'Test Summary' in content:
        # Extract total tests
        match = re.search(r'Total tests:\s*(\d+)', content)
        if match:
            results['total_tests'] = int(match.group(1))
        
        # Extract errors
        match = re.search(r'Errors:\s*(\d+)', content)
        if match:
            results['failed'] = int(match.group(1))
            results['passed'] = results['total_tests'] - results['failed']
        
        # Check for pass/fail message
        if 'ALL TESTS PASSED' in content:
            results['status'] = 'PASSED'
        elif 'SOME TESTS FAILED' in content:
            results['status'] = 'FAILED'
        else:
            results['status'] = 'UNKNOWN'
    
    # Extract errors
    error_pattern = re.compile(r'✗\s*Error:.+', re.MULTILINE)
    errors = error_pattern.findall(content)
    results['errors'] = errors
    
    # Extract warnings
    warning_pattern = re.compile(r'⚠\s*.+', re.MULTILINE)
    warnings = warning_pattern.findall(content)
    results['warnings'] = warnings
    
    return results

def main():
    results_dir = Path('output/test_results')
    
    if not results_dir.exists():
        print(f"❌ Test results directory not found: {results_dir}")
        return 1
    
    print("=" * 60)
    print("Unit Test Results Analysis")
    print("=" * 60)
    print(f"Results directory: {results_dir.absolute()}")
    print()
    
    # Find all result files (both .log and .txt)
    log_files = list(results_dir.glob('*.log')) + list(results_dir.glob('*.txt'))
    
    if not log_files:
        print("❌ No test log files found")
        return 1
    
    # Analyze each test
    all_passed = True
    total_tests = 0
    total_passed = 0
    total_failed = 0
    
    test_results = {}
    
    for log_file in sorted(log_files):
        test_name = log_file.stem
        print(f"📊 Test: {test_name}")
        print("-" * 60)
        
        try:
            results = analyze_test_log(log_file)
            test_results[test_name] = results
            
            if 'status' in results:
                status_icon = '✅' if results['status'] == 'PASSED' else '❌'
                print(f"   Status: {status_icon} {results['status']}")
                
                if results['status'] == 'PASSED':
                    print(f"   Tests: {results.get('total_tests', 'N/A')} passed")
                else:
                    all_passed = False
                    print(f"   Tests: {results.get('total_tests', 'N/A')} total, "
                          f"{results.get('passed', 0)} passed, "
                          f"{results.get('failed', 0)} failed")
                
                total_tests += results.get('total_tests', 0)
                total_passed += results.get('passed', 0)
                total_failed += results.get('failed', 0)
            else:
                all_passed = False
                print(f"   Status: ⚠️  INCOMPLETE")
            
            if results['errors']:
                print(f"   Errors ({len(results['errors'])}):")
                for error in results['errors'][:5]:  # Show first 5 errors
                    print(f"      - {error}")
                if len(results['errors']) > 5:
                    print(f"      ... and {len(results['errors']) - 5} more errors")
            
            if results['warnings']:
                print(f"   Warnings ({len(results['warnings'])}):")
                for warning in results['warnings'][:3]:  # Show first 3 warnings
                    print(f"      - {warning}")
                if len(results['warnings']) > 3:
                    print(f"      ... and {len(results['warnings']) - 3} more warnings")
            
        except Exception as e:
            all_passed = False
            print(f"   ❌ Error analyzing log: {e}")
        
        print()
    
    # Summary
    print("=" * 60)
    print("Summary")
    print("=" * 60)
    print(f"Total test suites: {len(log_files)}")
    print(f"Total test cases: {total_tests}")
    print(f"Passed: {total_passed} ({total_tests/total_passed*100 if total_passed > 0 else 0:.1f}%)")
    print(f"Failed: {total_failed}")
    
    if all_passed:
        print()
        print("✅ ALL TESTS PASSED!")
        return 0
    else:
        print()
        print("❌ SOME TESTS FAILED")
        return 1

if __name__ == '__main__':
    sys.exit(main())