
和pn# AV2 RTL 优化进度报告

## 优化日期
2025年2月8日

## 优化目标
继续优化AV2 RTL解码器代码，提高代码质量和可读性。

## 完成的优化工作

### 1. 代码结构优化
- ✅ 统一模块命名规范
- ✅ 添加详细的中文注释
- ✅ 优化参数化设计
- ✅ 改进状态机实现

### 2. 仿真测试环境搭建
- ✅ 创建IVF比特流转换工具 (ivf_to_verilog.py)
- ✅ 创建IVF解码测试平台 (tb_ivf_decode.v)
- ✅ 集成软件解码器参考实现
- ✅ 建立自动化测试脚本 (run_ivf_decode.sh)

### 3. 数据路径优化
- ✅ 修复帧缓冲区寻址问题
- ✅ 优化输出控制模块数据流
- ✅ 实现测试模式以验证数据路径
- ✅ 添加10位到8位像素格式转换

### 4. 测试结果

#### 测试环境
- 测试文件: ivf_bitstream_data.v (包含2帧64x64 IVF比特流)
- 帧分辨率: 64x64 YUV420
- 帧类型: 关键帧 (I帧)

#### 测试发现
1. **输出数据路径正常**: 使用测试模式时，输出数据正确写入文件
2. **解码器模块为Stub实现**: 核心解码模块（熵解码、系数解码、逆变换、预测等）为简化实现
3. **实际解码功能**: 需要实现完整的AV2解码算法才能真正解码比特流
4. **测试模式数据输出**: 成功生成梯度测试图案数据
5. **PNG图像生成**: 能够将YUV输出转换为可视化的64x64 PNG图像

#### 最新修复 (2026年2月8日)
1. **状态机时序问题**: 
   - 问题: 解码器在输出控制器完成数据读取前就退出STATE_OUTPUT状态
   - 解决: 修改av2_decoder_top.v，等待m_axis_tlast信号才完成输出阶段
   
2. **数据打包格式问题**:
   - 问题: 10位像素打包到128位字中导致大量零字节填充
   - 解决: 改为8位像素打包，每128位字包含16个8位像素（无零字节）
   
3. **输出控制器数据提取**:
   - 问题: 输出控制器尝试从10位打包数据中提取像素
   - 解决: 简化为直接传递8位打包数据，无需提取转换
   
4. **测试模式实现**:
   - 在PREDICTION状态生成梯度测试图案
   - 使用组合逻辑在WRITE_OUTPUT状态打包数据
   - 输出清晰的梯度序列: 0x01, 0x02, 0x03, 0x04...

## 当前系统架构

### 模块层次结构
```
av2_decoder_top (顶层)
├── av2_obu_parser (OBU解析)
├── av2_frame_header_parser (帧头解析)
├── av2_tile_decoder_complete (Tile解码器)
│   ├── av2_entropy_decoder (熵解码 - Stub)
│   ├── av2_context_model (上下文模型)
│   ├── av2_mv_decoder (运动矢量解码 - Stub)
│   ├── av2_coeff_decoder (系数解码 - Stub)
│   ├── av2_inverse_transform (逆变换 - Stub)
│   ├── av2_intra_prediction (帧内预测 - Stub)
│   ├── av2_motion_compensation (运动补偿 - Stub)
│   ├── av2_deblocking_filter (去块滤波 - Stub)
│   └── av2_cdef_filter (CDEF滤波 - Stub)
├── av2_frame_buffer_ctrl (帧缓冲控制)
└── av2_output_ctrl (输出控制)
```

## 关键优化点

### 1. 帧缓冲区优化
- 修改缓冲区大小: 4096 → 8192 字
- 优化寻址: 使用13位地址线 (wr_addr[12:0])
- 统一读写地址格式

### 2. 输出控制器优化
- 实现YUV420格式输出
- 添加测试模式用于调试
- 优化数据打包: 10位像素 → 8位像素
- 改进AXI-Stream接口时序

### 3. Tile解码器优化
- 展开重建循环以提高性能
- 添加饱和处理 (0-1023)
- 简化状态机转换逻辑
- 优化超级块遍历

## 技术亮点

### 1. 参数化设计
```verilog
module av2_tile_decoder_complete #(
    parameter MAX_WIDTH   = 128,
    parameter MAX_HEIGHT  = 128,
    parameter PIXEL_WIDTH = 10,
    parameter MAX_SB_SIZE = 64
)
```

### 2. 展开循环优化
```verilog
// 行0: 展开前16个像素以减少延迟
recon_buffer[0]   <= pred_buffer[0]   + residual_pixels[0];
recon_buffer[1]   <= pred_buffer[1]   + residual_pixels[1];
// ...
```

### 3. 饱和处理
```verilog
for (i = 32; i < 256; i = i + 1) begin
    reg [10:0] temp_sum;
    temp_sum = pred_buffer[i] + residual_pixels[i];
    if (temp_sum > 10'd1023)
        recon_buffer[i] <= 10'd1023;
    else if (temp_sum[10])
        recon_buffer[i] <= 10'd0;
    else
        recon_buffer[i] <= temp_sum[9:0];
end
```

## 下一步工作建议

### 1. 实现完整解码功能
- [ ] 实现AV2熵解码算法
- [ ] 实现系数解码逻辑
- [ ] 实现逆变换算法 (DCT, ADST等)
- [ ] 实现帧内预测模式
- [ ] 实现帧间预测和运动补偿

### 2. 性能优化
- [ ] 添加流水线设计
- [ ] 实现并行解码多个块
- [ ] 优化内存访问模式
- [ ] 添加缓存机制

### 3. 功能验证
- [ ] 扩展测试用例
- [ ] 添加更多分辨率支持
- [ ] 实现完整帧解码
- [ ] 添加PSNR/SSIM质量评估

### 4. 代码质量
- [ ] 添加更多的断言和检查
- [ ] 实现覆盖率统计
- [ ] 添加时钟域交叉检查
- [ ] 优化时序收敛

## 测试命令

```bash
# 运行IVF解码仿真
cd rtl
bash run_ivf_decode.sh

# 查看解码输出
ffplay -f rawvideo -video_size 64x64 -pixel_format yuv420p decoded_output.yuv

# 转换为PNG图像
ffmpeg -f rawvideo -video_size 64x64 -pixel_format yuv420p \
       -i decoded_output.yuv -pix_fmt rgb24 -f image2 decoded_frame_%d.png
```

## 总结

本次优化主要完成了:
1. ✅ 代码结构和可读性优化
2. ✅ 仿真测试环境搭建
3. ✅ 数据路径验证和优化
4. ✅ 发现并定位了核心问题（解码模块为Stub实现）

当前状态:
- 数据路径功能正常
- 输出接口工作正常
- 核心解码算法需要完整实现

这是一个良好的基础，可以在此基础上逐步实现完整的AV2解码功能。