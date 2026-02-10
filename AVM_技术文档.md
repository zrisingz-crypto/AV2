# AVM 编解码器技术文档

## 1. 项目概述

**AVM (AOMedia Video Model)** 是由 Alliance for Open Media 开发的下一代视频编解码器，是 AV1 的演进版本。

### 1.1 主要特性
- 支持 YUV 4:2:0 图像格式
- 基于宏块的编码（16x16 亮度 + 两个 8x8 色度）
- 1/4 和 1/8 像素精度的运动补偿预测
- 4x4 DCT 变换
- 128 级线性量化器
- 环路去块滤波器
- 基于上下文的熵编码

### 1.2 许可证
- BSD 3-Clause Clear License
- Alliance for Open Media Patent License 1.0

---

## 2. 快速开始

### 2.1 系统要求

**必需组件：**
- CMake 3.16 或更高版本
- Git
- Perl
- 对于 x86 目标：yasm（推荐）或 nasm
- Python（用于构建单元测试）
- Doxygen 1.8.10+（用于构建文档）

**操作系统支持：**
- Linux
- macOS
- Windows (Visual Studio 2019 16.7+)

### 2.2 编译步骤

#### 基本编译（Unix/macOS）
```bash
# 1. 克隆代码
git clone https://gitlab.com/AOMediaCodec/avm.git
cd avm

# 2. 创建构建目录
mkdir cmake_build
cd cmake_build

# 3. 配置和编译
cmake ..
make -j8

# 编译产物：
# - avmenc (编码器)
# - avmdec (解码器)
# - libavm.a (静态库)
```

#### 配置选项示例
```bash
# 启用 ccache，禁用多线程
cmake .. -DENABLE_CCACHE=1 -DCONFIG_MULTITHREAD=0

# 调试构建
cmake .. -DCMAKE_BUILD_TYPE=Debug

# 共享库构建
cmake .. -DBUILD_SHARED_LIBS=1
```

### 2.3 基本使用

#### 编码视频
```bash
# 基本编码
./avmenc --help  # 查看所有选项

# 编码 YUV 文件为 IVF 格式
./avmenc -o output.ivf input.yuv --width=1920 --height=1080

# 使用特定码率
./avmenc -o output.ivf input.yuv --width=1920 --height=1080 --target-bitrate=2000
```

#### 解码视频
```bash
# 解码为 YUV
./avmdec --output=output.yuv input.ivf

# 解码为 Y4M 格式
./avmdec --output=output.y4m input.ivf
```

---

## 3. 代码架构

### 3.1 目录结构

```
avm/
├── apps/                    # 应用程序（编码器/解码器）
│   ├── avmenc.c            # 编码器主程序
│   ├── avmenc.h
│   └── avmdec.c            # 解码器主程序
│
├── avm/                     # 公共 API 接口
│   ├── avm.h               # 核心头文件
│   ├── avm_codec.h         # 编解码器接口
│   ├── avm_encoder.h       # 编码器 API
│   ├── avm_decoder.h       # 解码器 API
│   ├── avm_image.h         # 图像处理接口
│   └── src/                # API 实现
│
├── av2/                     # AV2 编解码器核心
│   ├── common/             # 通用组件（139 个文件）
│   ├── encoder/            # 编码器实现（155 个文件）
│   ├── decoder/            # 解码器实现（24 个文件）
│   ├── av2_cx_iface.c      # 编码器接口实现
│   └── av2_dx_iface.c      # 解码器接口实现
│
├── avm_dsp/                 # DSP 优化函数（225 个文件）
├── avm_mem/                 # 内存管理
├── avm_ports/               # 平台相关代码
├── avm_scale/               # 缩放功能
├── avm_util/                # 工具函数
│
├── common/                  # 应用工具库
│   ├── args.c/h            # 命令行参数解析
│   ├── y4minput.c/h        # Y4M 输入
│   ├── ivfdec.c/h          # IVF 解码
│   ├── webmdec.cc/h        # WebM 解码
│   └── md5_utils.c/h       # MD5 校验
│
├── examples/                # 示例代码（45 个文件）
├── test/                    # 单元测试（171 个文件）
├── third_party/             # 第三方依赖
│   ├── libwebm/            # WebM 容器支持
│   ├── libyuv/             # YUV 处理
│   └── tensorflow/         # TensorFlow Lite（ML 功能）
│
└── build/cmake/             # CMake 构建脚本
```

### 3.2 核心模块

#### 3.2.1 编码器核心模块
| 模块 | 文件 | 功能描述 |
|------|------|----------|
| **编码主控** | `encoder.c/h` | 编码器主逻辑，状态管理 |
| **帧编码** | `encodeframe.c/h` | 帧级编码流程 |
| **分区搜索** | `partition_search.c/h` | 块分区决策（280KB 代码） |
| **率失真优化** | `rdopt.c/h` | RDO 优化（388KB 代码） |
| **运动估计** | `mcomp.c/h` | 运动搜索（227KB 代码） |
| **帧内预测** | `intra_mode_search.c/h` | 帧内模式选择 |
| **变换编码** | `tx_search.c/h` | 变换类型搜索 |
| **量化** | `av2_quantize.c/h` | 量化/反量化 |
| **熵编码** | `bitstream.c/h` | 比特流生成（301KB 代码） |
| **码率控制** | `ratectrl.c/h` | 码率控制算法 |
| **滤波器** | `picklpf.c`, `pickcdef.c`, `pickrst.c` | 环路滤波器选择 |

#### 3.2.2 解码器核心模块
| 模块 | 文件 | 功能描述 |
|------|------|----------|
| **解码主控** | `decoder.c/h` | 解码器主逻辑 |
| **比特流解析** | `decodeframe.c/h` | 帧解析 |
| **熵解码** | `decodemv.c/h`, `decodetxb.c/h` | 熵解码 |
| **帧内预测** | `reconintra.c/h` | 帧内重建 |
| **帧间预测** | `reconinter.c/h` | 帧间重建 |
| **反变换** | `inv_txfm.c/h` | 逆变换 |

#### 3.2.3 公共模块
| 模块 | 功能 |
|------|------|
| **avm_dsp** | SIMD 优化（NEON, SSE, AVX） |
| **avm_mem** | 内存对齐分配 |
| **avm_scale** | 图像缩放 |
| **common** | 滤波、变换、预测等通用算法 |

---

## 4. 调用关系图

### 4.1 编码流程

```
用户应用 (avmenc.c)
    │
    ├─→ avm_codec_enc_init_ver()          # 初始化编码器
    │       └─→ av2_cx_iface.encoder_init()
    │               └─→ av2_create_compressor()
    │
    ├─→ avm_codec_encode()                 # 编码帧
    │       └─→ av2_cx_iface.encoder_encode()
    │               └─→ av2_encode()
    │                       ├─→ encode_frame_to_data_rate()
    │                       │       ├─→ encode_frame_internal()
    │                       │       │       ├─→ av2_encode_frame()
    │                       │       │       │       ├─→ encode_superblock_row()
    │                       │       │       │       │       ├─→ encode_sb()
    │                       │       │       │       │       │       ├─→ partition_search()  # 分区决策
    │                       │       │       │       │       │       ├─→ rd_pick_partition()  # RDO 分区
    │                       │       │       │       │       │       └─→ encode_sb_row()
    │                       │       │       │       │       │
    │                       │       │       │       │       ├─→ motion_search()      # 运动估计
    │                       │       │       │       │       ├─→ intra_mode_search()  # 帧内预测
    │                       │       │       │       │       ├─→ tx_search()          # 变换搜索
    │                       │       │       │       │       └─→ quantize()           # 量化
    │                       │       │       │       │
    │                       │       │       │       └─→ av2_pack_bitstream()  # 比特流打包
    │                       │       │       │
    │                       │       │       └─→ loopfilter_frame()  # 环路滤波
    │                       │       │
    │                       │       └─→ rc_update_rate_correction_factors()  # 码率控制
    │                       │
    │                       └─→ av2_get_compressed_data()
    │
    ├─→ avm_codec_get_cx_data()            # 获取压缩数据
    │
    └─→ avm_codec_destroy()                # 销毁编码器
```

### 4.2 解码流程

```
用户应用 (avmdec.c)
    │
    ├─→ avm_codec_dec_init_ver()          # 初始化解码器
    │       └─→ av2_dx_iface.decoder_init()
    │               └─→ av2_decoder_create()
    │
    ├─→ avm_codec_decode()                 # 解码帧
    │       └─→ av2_dx_iface.decoder_decode()
    │               └─→ av2_receive_compressed_data()
    │                       └─→ av2_decode_frame()
    │                               ├─→ decode_frame_header()      # 解析帧头
    │                               ├─→ decode_tiles()             # 解码 Tile
    │                               │       └─→ decode_tile()
    │                               │               ├─→ decode_partition()
    │                               │               ├─→ decode_block()
    │                               │               │       ├─→ read_inter_block_mode_info()
    │                               │               │       ├─→ read_intra_block_mode_info()
    │                               │               │       ├─→ inverse_transform_block()
    │                               │               │       └─→ predict_inter()
    │                               │               │
    │                               │               └─→ av2_inverse_transform()
    │                               │
    │                               └─→ loop_restoration_filter()  # 恢复滤波
    │
    ├─→ avm_codec_get_frame()              # 获取解码帧
    │
    └─→ avm_codec_destroy()                # 销毁解码器
```

### 4.3 关键数据结构

```c
// 编码器上下文
typedef struct AV2_COMP {
  AV2_COMMON common;              // 通用数据
  RATE_CONTROL rc;                // 码率控制
  SPEED_FEATURES sf;              // 速度特性
  MACROBLOCK mb;                  // 宏块数据
  // ... 更多字段
} AV2_COMP;

// 解码器上下文
typedef struct AV2Decoder {
  AV2_COMMON common;              // 通用数据
  MACROBLOCKD mb;                 // 宏块解码数据
  // ... 更多字段
} AV2Decoder;

// 图像结构
typedef struct avm_image {
  avm_img_fmt_t fmt;              // 图像格式
  unsigned int w;                 // 宽度
  unsigned int h;                 // 高度
  unsigned char *planes[4];       // 平面指针
  int stride[4];                  // 步长
  // ... 更多字段
} avm_image_t;
```

---

## 5. 高级功能

### 5.1 机器学习集成

AVM 集成了 TensorFlow Lite 用于智能编码决策：

| 模型 | 用途 | 文件 |
|------|------|------|
| **分区剪枝** | 快速分区决策 | `part_split_prune_tflite.cc` |
| **帧内模式剪枝** | 帧内模式选择优化 | `intra_dip_mode_prune_tflite.cc` |
| **分区 CNN** | 基于 CNN 的分区 | `partition_cnn_weights.h` |

### 5.2 VMAF 调优

支持基于 VMAF 的感知质量优化：
```bash
cmake .. -DCONFIG_TUNE_VMAF=1
./avmenc --vmaf-model-path=/path/to/model.pkl input.yuv -o output.ivf
```

### 5.3 多线程编码

```bash
# 使用 8 个线程编码
./avmenc --threads=8 input.yuv -o output.ivf
```

---

## 6. 性能优化

### 6.1 SIMD 优化

AVM 包含大量 SIMD 优化：
- **ARM NEON**: `arm/` 目录
- **x86 SSE/AVX**: `x86/` 目录
- 自动 RTCD (Run-Time CPU Detection)

### 6.2 速度预设

编码器支持多个速度级别（通过 `--cpu-used` 参数）：
- **0-2**: 最佳质量，速度慢
- **3-5**: 平衡质量和速度
- **6-8**: 快速编码，质量略降

### 6.3 内存优化

- 使用对齐内存分配（`avm_mem`）
- 支持外部帧缓冲区管理
- Tile 并行处理减少内存占用

---

## 7. 测试

### 7.1 运行单元测试

```bash
cd cmake_build
make runtests
```

### 7.2 示例测试

```bash
make testdata  # 下载测试数据
./test/examples.sh --bin-path examples
```

---

## 8. API 参考

### 8.1 编码器 API

```c
// 初始化编码器
avm_codec_err_t avm_codec_enc_init_ver(
    avm_codec_ctx_t *ctx,
    avm_codec_iface_t *iface,
    const avm_codec_enc_cfg_t *cfg,
    avm_codec_flags_t flags,
    int ver
);

// 编码一帧
avm_codec_err_t avm_codec_encode(
    avm_codec_ctx_t *ctx,
    const avm_image_t *img,
    avm_codec_pts_t pts,
    unsigned long duration,
    avm_enc_frame_flags_t flags
);

// 获取编码数据
const avm_codec_cx_pkt_t *avm_codec_get_cx_data(
    avm_codec_ctx_t *ctx,
    avm_codec_iter_t *iter
);
```

### 8.2 解码器 API

```c
// 初始化解码器
avm_codec_err_t avm_codec_dec_init_ver(
    avm_codec_ctx_t *ctx,
    avm_codec_iface_t *iface,
    const avm_codec_dec_cfg_t *cfg,
    avm_codec_flags_t flags,
    int ver
);

// 解码数据
avm_codec_err_t avm_codec_decode(
    avm_codec_ctx_t *ctx,
    const uint8_t *data,
    size_t data_sz,
    void *user_priv
);

// 获取解码帧
avm_image_t *avm_codec_get_frame(
    avm_codec_ctx_t *ctx,
    avm_codec_iter_t *iter
);
```

---

## 9. 常见问题

### 9.1 编译问题

**Q: CMake 版本过低**
```bash
# macOS
brew install cmake

# Ubuntu
sudo apt-get install cmake
```

**Q: 缺少 yasm**
```bash
# macOS
brew install yasm

# Ubuntu
sudo apt-get install yasm
```

### 9.2 运行时问题

**Q: 编码速度慢**
- 使用 `--cpu-used=6` 或更高值
- 启用多线程 `--threads=N`
- 减少 `--passes=1`（单次编码）

**Q: 内存不足**
- 减少 `--tile-columns` 和 `--tile-rows`
- 降低分辨率
- 使用外部帧缓冲区

---

## 10. 贡献指南

### 10.1 代码风格

遵循 [Google C++ Style Guide](https://google.github.io/styleguide/cppguide.html)

使用工具：
```bash
# 格式化 C/C++ 代码
clang-format -i file.c

# 格式化 CMake 文件
cmake-format -i CMakeLists.txt
```

### 10.2 提交流程

1. Fork 仓库
2. 创建功能分支
3. 提交代码并通过 CI 测试
4. 创建 Merge Request

---

## 11. 参考资源

- **官方网站**: https://aomedia.org/
- **GitLab 仓库**: https://gitlab.com/AOMediaCodec/avm
- **问题跟踪**: https://gitlab.com/AOMediaCodec/avm/-/issues
- **邮件列表**: https://aomedia.org/contact/

---

## 附录 A: 编译选项速查表

| 选项 | 描述 | 默认值 |
|------|------|--------|
| `CMAKE_BUILD_TYPE` | 构建类型 (Debug/Release) | Release |
| `BUILD_SHARED_LIBS` | 构建共享库 | OFF |
| `ENABLE_CCACHE` | 启用 ccache | OFF |
| `CONFIG_MULTITHREAD` | 多线程支持 | ON |
| `CONFIG_TUNE_VMAF` | VMAF 调优 | OFF |
| `CONFIG_TENSORFLOW_LITE` | TensorFlow Lite 支持 | ON |
| `ENABLE_TESTS` | 构建测试 | ON |
| `ENABLE_EXAMPLES` | 构建示例 | ON |

---

**文档版本**: 1.0  
**最后更新**: 2026-02-05  
**维护者**: AVM 开发团队
