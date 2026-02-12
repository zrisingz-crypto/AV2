#!/usr/bin/env python3
"""
Run unit tests with timeout protection
"""
import subprocess
import sys
import time

def run_with_timeout(cmd, timeout=90):
    """Run command with timeout"""
    print(f"Running: {' '.join(cmd)}")
    print(f"Timeout: {timeout} seconds")
    print("-" * 50)
    
    start = time.time()
    p = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, 
                        text=True, bufsize=1, universal_newlines=True)
    
    while True:
        # Check if process completed
        if p.poll() is not None:
            # Print remaining output
            for line in p.stdout:
                print(line, end='')
            break
        
        # Check timeout
        if time.time() - start > timeout:
            p.terminate()
            p.wait()
            print(f"\n{'=' * 50}")
            print(f"ERROR: Test timed out after {timeout} seconds")
            print(f"{'=' * 50}")
            return False
        
        # Read and print output
        try:
            line = p.stdout.readline()
            if line:
                print(line, end='')
        except:
            pass
    
    # Return exit code
    if p.returncode == 0:
        print(f"\n{'=' * 50}")
        print("SUCCESS: All tests passed")
        print(f"{'=' * 50}")
        return True
    else:
        print(f"\n{'=' * 50}")
        print(f"ERROR: Tests failed with exit code {p.returncode}")
        print(f"{'=' * 50}")
        return False

if __name__ == "__main__":
    cmd = ['bash', 'scripts/run_unit_tests.sh']
    success = run_with_timeout(cmd, timeout=90)
    sys.exit(0 if success else 1)