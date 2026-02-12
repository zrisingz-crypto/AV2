#!/usr/bin/env python3
with open('output/sw_output.yuv', 'rb') as f:
    content = f.read()

# Find the start of Y plane (after "FRAME" header)
frame_header_pos = content.find(b'FRAME')
if frame_header_pos == -1:
    raise Exception("FRAME header not found")

# Skip FRAME header and newline
y_start = frame_header_pos + len(b'FRAME') + 1

# For 64x64 frame (only Y plane for testing)
frame_size = 64 * 64
y_plane = content[y_start:y_start + frame_size]

if len(y_plane) < frame_size:
    raise Exception(f"Y plane is only {len(y_plane)} bytes, expected {frame_size}")

# Generate hex format for $readmemh
hex_lines = []

for i in range(len(y_plane)):
    pixel_value = y_plane[i]
    hex_str = f"{pixel_value:02x}"
    hex_lines.append(hex_str)

# Write to ROM file
with open('output/sw_pixel_rom.txt', 'w') as f:
    for line in hex_lines:
        f.write(line + '\n')

print(f"Generated ROM for {len(y_plane)} pixels")
print("Output written to output/sw_pixel_rom.txt")
print(f"First 10 pixels (hex): {hex_lines[:10]}")
