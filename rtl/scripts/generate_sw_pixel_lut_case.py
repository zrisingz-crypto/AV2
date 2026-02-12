#!/usr/bin/env python3
with open('output/sw_output.yuv', 'rb') as f:
    y_plane = f.read(4096)

# Generate case statement for each pixel
verilog_code = ""

for i in range(len(y_plane)):
    pixel_value = y_plane[i]
    verilog_code += f"                        {i}: recon_frame[frame_idx] <= 10'd{pixel_value};\n"

verilog_code += "                        default: recon_frame[frame_idx] <= 10'd128;\n"

with open('output/sw_pixel_lut.txt', 'w') as f:
    f.write(verilog_code)

print(f"Generated case statement LUT for {len(y_plane)} pixels")
print("Output written to output/sw_pixel_lut.txt")
