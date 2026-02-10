#!/usr/bin/env python3
"""
Convert IVF file to Verilog testbench format
Extracts bitstream data for RTL simulation
"""

import struct
import sys
from pathlib import Path

# IVF file header format (32 bytes)
# DKIF[0:3], version[4], header_size[5], fourcc[6:9], width[10:11], height[12:13]
# timebase_denom[14:17], timebase_numer[18:21], num_frames[22:25], unused[26:29]

def parse_ivf_header(data):
    """Parse IVF file header"""
    if len(data) < 32:
        raise ValueError("IVF file too small")
    
    magic = data[0:4]
    if magic != b'DKIF':
        raise ValueError(f"Invalid IVF magic: {magic}")
    
    version = struct.unpack('<H', data[4:6])[0]
    header_size = struct.unpack('<H', data[6:8])[0]
    fourcc = data[8:12]
    width = struct.unpack('<H', data[12:14])[0]
    height = struct.unpack('<H', data[14:16])[0]
    timebase_denom = struct.unpack('<I', data[16:20])[0]
    timebase_numer = struct.unpack('<I', data[20:24])[0]
    num_frames = struct.unpack('<I', data[24:28])[0]
    
    return {
        'version': version,
        'header_size': header_size,
        'fourcc': fourcc,
        'width': width,
        'height': height,
        'timebase_denom': timebase_denom,
        'timebase_numer': timebase_numer,
        'num_frames': num_frames
    }

def parse_ivf_frame(data, offset):
    """Parse a single IVF frame header"""
    if offset + 12 > len(data):
        return None
    
    frame_header = struct.unpack('<III', data[offset:offset+12])
    
    return {
        'frame_size': frame_header[0],
        'timestamp': frame_header[1],
        'frame_number': frame_header[2],
        'data_offset': offset + 12,
        'data': data[offset+12:offset+12+frame_header[0]]
    }

def generate_verilog_mem(ivf_file, output_file, max_frames=2):
    """Generate Verilog memory file from IVF"""
    
    print(f"Reading IVF file: {ivf_file}")
    with open(ivf_file, 'rb') as f:
        ivf_data = f.read()
    
    # Parse header
    header = parse_ivf_header(ivf_data)
    print(f"IVF Header:")
    print(f"  FourCC: {header['fourcc']}")
    print(f"  Resolution: {header['width']}x{header['height']}")
    print(f"  Frames: {header['num_frames']}")
    
    # Parse frames
    offset = 32  # Skip IVF header
    frames = []
    frame_count = 0
    
    while offset < len(ivf_data) and frame_count < max_frames:
        frame = parse_ivf_frame(ivf_data, offset)
        if frame is None:
            break
        
        frames.append(frame)
        print(f"\nFrame {frame_count}:")
        print(f"  Size: {frame['frame_size']} bytes")
        print(f"  Timestamp: {frame['timestamp']}")
        print(f"  First bytes: {frame['data'][:16].hex()}")
        
        offset += 12 + frame['frame_size']
        frame_count += 1
    
    if not frames:
        raise ValueError("No frames found in IVF file")
    
    # Generate Verilog code
    print(f"\nGenerating Verilog test data: {output_file}")
    
    with open(output_file, 'w') as f:
        f.write("//==============================================================================\n")
        f.write("// IVF Bitstream Data\n")
        f.write(f"// Source: {ivf_file}\n")
        f.write(f"// Resolution: {header['width']}x{header['height']}\n")
        f.write(f"// Frames extracted: {len(frames)}\n")
        f.write("//==============================================================================\n\n")
        
        f.write(f"`define IVF_WIDTH {header['width']}\n")
        f.write(f"`define IVF_HEIGHT {header['height']}\n")
        f.write(f"`define IVF_NUM_FRAMES {len(frames)}\n")
        f.write(f"`define IVF_HEADER_SIZE 12\n\n")
        
        # Generate frame data
        for i, frame in enumerate(frames):
            f.write(f"// Frame {i}: {frame['frame_size']} bytes\n")
            f.write(f"localparam FRAME_{i}_SIZE = {frame['frame_size']};\n")
            f.write(f"reg [7:0] frame_{i}_data [0:{frame['frame_size']}-1];\n")
            f.write(f"initial $readmemh(\"frame_{i}.hex\", frame_{i}_data);\n\n")
            
            # Also write hex file for this frame
            hex_file = output_file.rsplit('.', 1)[0] + f"_frame_{i}.hex"
            with open(hex_file, 'w') as hf:
                for byte in frame['data']:
                    hf.write(f"{byte:02x}\n")
            print(f"  Wrote {hex_file}")
        
        # Generate size array
        f.write(f"\n// Frame sizes array\n")
        f.write(f"reg [31:0] frame_sizes [0:{len(frames)-1}];\n")
        f.write(f"initial begin\n")
        for i, frame in enumerate(frames):
            f.write(f"    frame_sizes[{i}] = FRAME_{i}_SIZE;\n")
        f.write(f"end\n\n")
        
        # Generate total data
        total_size = sum(f['frame_size'] for f in frames)
        f.write(f"// Total bitstream size\n")
        f.write(f"localparam BITSTREAM_SIZE = {total_size};\n")
        
    print(f"\n✓ Generated {output_file}")
    print(f"✓ Total bitstream size: {total_size} bytes")
    return header, frames

if __name__ == '__main__':
    ivf_file = 'crowd_64x64p30f32_av2.ivf'
    output_file = 'ivf_bitstream_data.v'
    
    try:
        header, frames = generate_verilog_mem(ivf_file, output_file, max_frames=2)
        print("\n✓ Success!")
    except Exception as e:
        print(f"\n✗ Error: {e}", file=sys.stderr)
        sys.exit(1)