#!/usr/bin/env python3
import sys
import struct

def hex_to_yuv(input_hex, output_yuv, width=64, height=64, bit_depth=10, to_8bit=True):
    print(f"[*] Converting {input_hex} to {output_yuv}...")
    print(f"[*] Resolution: {width}x{height}, Bit Depth: {bit_depth}-bit")

    try:
        with open(input_hex, 'r') as f_in, open(output_yuv, 'wb') as f_out:
            pixels_all = []
            
            for line in f_in:
                line = line.strip()
                if not line or line.startswith('/'): continue
                
                # 将 128-bit 十六进制字符串转换为整数
                val = int(line, 16)
                
                # 根据 Verilog 里的打包逻辑提取像素 (每个 128-bit 存了 12 个 10-bit 像素)
                for i in range(12):
                    pixel = (val >> (i * 10)) & 0x3FF
                    
                    if to_8bit:
                        # 转换到 8-bit: 除以 4 (右移 2 位)
                        f_out.write(struct.pack('B', pixel >> 2))
                    else:
                        # 保持 10-bit: 以 Little Endian 写入 2 字节
                        f_out.write(struct.pack('<H', pixel))
                    
                    pixels_all.append(pixel)
                    if len(pixels_all) >= width * height:
                        break
                
                if len(pixels_all) >= width * height:
                    break

            print(f"[+] Success! Total pixels processed: {len(pixels_all)}")
            if to_8bit:
                print(f"[!] View using: ffplay -f rawvideo -pixel_format gray -video_size {width}x{height} {output_yuv}")
            else:
                print(f"[!] View using: ffplay -f rawvideo -pixel_format gray10le -video_size {width}x{height} {output_yuv}")

    except Exception as e:
        print(f"[X] Error: {e}")

if __name__ == "__main__":
    if len(sys.argv) < 3:
        print("Usage: python3 hex2yuv.py <input.hex> <output.yuv> [width] [height]")
        sys.exit(1)
    
    in_file = sys.argv[1]
    out_file = sys.argv[2]
    w = int(sys.argv[3]) if len(sys.argv) > 3 else 64
    h = int(sys.argv[4]) if len(sys.argv) > 4 else 64
    
    hex_to_yuv(in_file, out_file, w, h)
