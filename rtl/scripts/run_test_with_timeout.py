#!/usr/bin/env python3
import subprocess
import sys
import time

def main():
    # Run simulation with timeout
    try:
        process = subprocess.Popen(
            ["vvp sim_real_tile_decoder"],
            shell=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True
        )
        
        start_time = time.time()
        timeout = 120  # seconds - increased for larger frames
        
        while process.poll() is None:
            if time.time() - start_time > timeout:
                print("Simulation timed out after", timeout, "seconds")
                process.kill()
                return 1
            time.sleep(0.1)
        
        # Get output
        stdout, stderr = process.communicate()
        print(stdout[:5000])  # Print first 5000 characters
        if stderr:
            print("\nError output:")
            print(stderr)
            
        return process.returncode
        
    except Exception as e:
        print(f"Error: {e}")
        return 1

if __name__ == "__main__":
    sys.exit(main())
