# Scripts 技能文档

本文档说明rtl/scripts/目录下各个脚本所需的技能和使用方法。

## 脚本分类

### 测试运行脚本

#### run_unit_tests.sh
**用途**: 运行所有单元测试  
**所需技能**:
- Shell脚本编写
- Icarus Verilog仿真
- 日志分析
- 进程管理

**使用方法**:
```bash
bash scripts/run_unit_tests.sh
```

**输出**:
- 测试结果摘要
- 详细日志文件: output/test_results/

**依赖**:
- iverilog
- vvp
- 各个模块的testbench文件

---

#### run_real_decoder.sh
**用途**: 运行真实IVF解码器测试  
**所需技能**:
- Shell脚本编写
- IVF bitstream处理
- YUV文件操作
- 日志管理

**使用方法**:
```bash
bash scripts/run_real_decoder.sh [ivf_file]
```

**输出**:
- 解码输出YUV文件
- 仿真日志

**依赖**:
- iverilog
- vvp
- ivf_bitstream_data.v
- tb_ivf_decode_real.v

---

#### compile_ivf_real_decoder.sh
**用途**: 编译IVF真实解码器  
**所需技能**:
- Icarus Verilog编译
- 模块链接
- 文件路径管理

**使用方法**:
```bash
bash scripts/compile_ivf_real_decoder.sh
```

**输出**:
- 可执行文件: sim_real_decoder

**依赖**:
- iverilog
- 所有相关RTL模块

---

#### compile_real_decoder.sh
**用途**: 编译真实解码器  
**所需技能**:
- Icarus Verilog编译
- 模块链接
- 文件路径管理

**使用方法**:
```bash
bash scripts/compile_real_decoder.sh
```

**输出**:
- 可执行文件

**依赖**:
- iverilog
- RTL模块文件

---

#### run_ivf_decode.sh
**用途**: 运行IVF解码测试  
**所需技能**:
- Shell脚本编写
- IVF文件解析
- 解码流程控制

**使用方法**:
```bash
bash scripts/run_ivf_decode.sh [ivf_file]
```

**输出**:
- 解码结果
- 性能统计

---

#### run_full_decode.sh
**用途**: 运行完整解码流程  
**所需技能**:
- 流程编排
- 多步骤协调
- 错误处理

**使用方法**:
```bash
bash scripts/run_full_decode.sh [ivf_file]
```

---

#### run_real_decoder_test.sh
**用途**: 运行真实解码器测试  
**所需技能**:
- 测试用例管理
- 结果验证
- 性能分析

**使用方法**:
```bash
bash scripts/run_real_decoder_test.sh [ivf_file]
```

---

#### run_test_with_timeout.py
**用途**: 带超时的测试运行器  
**所需技能**:
- Python多进程
- 超时控制
- 进程管理

**使用方法**:
```bash
python scripts/run_test_with_timeout.py [command] [timeout]
```

---

### 数据分析脚本

#### analyze_test_results.py
**用途**: 分析测试结果  
**所需技能**:
- Python数据分析
- 日志解析
- 统计计算
- 报告生成

**使用方法**:
```bash
python scripts/analyze_test_results.py [log_file]
```

**输出**:
- 测试通过率
- 性能统计
- 错误汇总

---

#### analyze_sw_decoder.py
**用途**: 分析软件解码器输出  
**所需技能**:
- Python数据分析
- 解码器原理
- 数据格式理解

**使用方法**:
```bash
python scripts/analyze_sw_decoder.py [output_file]
```

**输出**:
- 解码统计信息
- 帧信息
- 质量指标

---

#### analyze_sw_output.py
**用途**: 分析软件输出  
**所需技能**:
- Python文本处理
- 日志分析
- 模式匹配

**使用方法**:
```bash
python scripts/analyze_sw_output.py [log_file]
```

---

### 对比脚本

#### compare_yuv.py
**用途**: 对比两个YUV文件  
**所需技能**:
- Python图像处理
- YUV格式理解
- 像素级对比
- PSNR/SSIM计算

**使用方法**:
```bash
python scripts/compare_yuv.py [file1.yuv] [file2.yuv] [width] [height]
```

**输出**:
- PSNR值
- SSIM值
- 差异统计
- 差异可视化（可选）

---

#### compare_yuv_files.py
**用途**: 对比YUV文件（批量）  
**所需技能**:
- Python批量处理
- 文件管理
- 结果汇总

**使用方法**:
```bash
python scripts/compare_yuv_files.py [dir1] [dir2]
```

---

#### compare_decoders.py
**用途**: 对比硬件和软件解码器  
**所需技能**:
- Python数据分析
- 解码器对比
- 性能分析
- 差异定位

**使用方法**:
```bash
python scripts/compare_decoders.py [hw_output] [sw_output]
```

---

### 数据提取脚本

#### extract_rtl_output.py
**用途**: 从RTL仿真输出提取数据  
**所需技能**:
- Python正则表达式
- 日志解析
- 数据格式转换
- YUV文件生成

**使用方法**:
```bash
python scripts/extract_rtl_output.py [sim_log] [output_yuv]
```

**输出**:
- YUV格式文件
- 元数据文件

---

### 生成脚本

#### generate_test_data.sh
**用途**: 生成测试数据  
**所需技能**:
- Shell脚本编写
- 测试数据生成
- 随机数据生成
- 数据格式转换

**使用方法**:
```bash
bash scripts/generate_test_data.sh [count]
```

**输出**:
- 测试数据文件
- 测试向量

---

#### generate_sw_pixel_lut.py
**用途**: 生成软件像素LUT  
**所需技能**:
- Python数组操作
- 查找表生成
- 像素处理算法

**使用方法**:
```bash
python scripts/generate_sw_pixel_lut.py
```

**输出**:
- LUT数据文件
- Verilog ROM文件

---

#### generate_sw_pixel_lut_case.py
**用途**: 生成软件像素LUT的case语句  
**所需技能**:
- Python代码生成
- Verilog case语句
- 查找表映射

**使用方法**:
```bash
python scripts/generate_sw_pixel_lut_case.py
```

**输出**:
- Verilog case语句代码

---

#### generate_sw_pixel_lut_rom.py
**用途**: 生成软件像素LUT的ROM文件  
**所需技能**:
- Python文件操作
- Verilog ROM格式
- 内存初始化

**使用方法**:
```bash
python scripts/generate_sw_pixel_lut_rom.py
```

**输出**:
- Verilog ROM文件

---

### 转换脚本

#### hex2yuv.py
**用途**: 十六进制转YUV格式  
**所需技能**:
- Python数据转换
- 十六进制解析
- YUV格式理解

**使用方法**:
```bash
python scripts/hex2yuv.py [hex_file] [yuv_file]
```

**输出**:
- YUV格式文件

---

#### ivf_to_verilog.py
**用途**: IVF转Verilog数据文件  
**所需技能**:
- Python文件解析
- IVF格式理解
- Verilog数组生成
- 位操作

**使用方法**:
```bash
python scripts/ivf_to_verilog.py [ivf_file] [output_file]
```

**输出**:
- Verilog数据文件
- 可用于仿真

---

#### yuv_to_png.py
**用途**: YUV转PNG图像  
**所需技能**:
- Python图像处理
- YUV到RGB转换
- PNG编码
- 图像库使用

**使用方法**:
```bash
python scripts/yuv_to_png.py [yuv_file] [width] [height] [output.png]
```

**输出**:
- PNG图像文件

---

## 开发技能要求

### Shell脚本开发
- [ ] Bash编程基础
- [ ] 条件判断和循环
- [ ] 函数定义
- [ ] 参数处理
- [ ] 文件操作
- [ ] 进程管理
- [ ] 错误处理

### Python脚本开发
- [ ] Python基础语法
- [ ] 文件I/O
- [ ] 正则表达式
- [ ] 数据处理
- [ ] 图像处理（PIL/Pillow）
- [ ] 进程和线程
- [ ] 日志记录

### 视频编解码知识
- [ ] YUV格式理解
- [ ] IVF格式解析
- [ ] 像素处理
- [ ] PSNR/SSIM计算
- [ ] 解码流程
- [ ] 比特流处理

### RTL仿真知识
- [ ] Icarus Verilog使用
- [ ] iverilog编译
- [ ] vvp仿真
- [ ] 日志解析
- [ ] 波形分析

## 最佳实践

### 脚本编写
1. **添加文件头注释**: 说明用途、作者、日期
2. **参数验证**: 检查输入参数合法性
3. **错误处理**: 捕获和处理异常
4. **日志记录**: 记录重要操作和错误
5. **文档化**: 添加使用说明和示例

### 测试脚本
1. **自动化**: 尽量自动化测试流程
2. **超时保护**: 添加超时机制防止死锁
3. **结果收集**: 收集和汇总测试结果
4. **错误报告**: 清晰的错误信息
5. **可重复**: 确保测试可重复执行

### 分析脚本
1. **可配置**: 支持配置参数
2. **可扩展**: 易于添加新分析功能
3. **可视化**: 提供直观的输出
4. **性能优化**: 处理大文件时优化性能
5. **文档完整**: 清晰的使用说明

## 常见问题

### Q: 如何添加新的测试脚本？
A: 遵循命名规范（run_*.sh或test_*.sh），添加文件头注释，使用统一的输出格式。

### Q: 如何调试脚本？
A: 使用bash -x运行Shell脚本，使用pdb调试Python脚本，添加详细的日志输出。

### Q: 如何处理大文件？
A: 使用流式处理，避免一次性加载整个文件，考虑使用多进程加速。

### Q: 如何集成到CI/CD？
A: 创建CI配置文件，定义测试流程，设置超时和失败处理。

---

**版本**: V1.0  
**创建日期**: 2026-02-12  
**最后更新**: 2026-02-12