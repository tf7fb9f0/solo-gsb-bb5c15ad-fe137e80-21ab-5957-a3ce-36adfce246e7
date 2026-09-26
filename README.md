# 同位素峰包络计算 API

使用 Python + FastAPI 实现的分子同位素峰包络（isotopic envelope）计算服务。
输入分子中各元素的原子数以及每种同位素的**整数质量数、十进制精确质量和丰度**，
服务以各元素最小质量数为基准，对各原子的同位素分布进行**独立抽样的完整卷积**，
输出按质量偏移升序排列的非零概率峰及概率加权精确质量质心。

## 计算规则

1. 以每种元素**最小质量数**同位素作为该元素的整数质量偏移基准（偏移 = 0）。
2. 同一元素内，相同质量数的多个同位素先合并（概率相加，质量按概率加权）。
3. 每种元素的 `atom_count` 个原子视为独立抽样，用该元素的单原子分布做完整
   n 重卷积（平方倍增实现）；不同元素之间再做卷积。
4. 卷积过程中**相同偏移合并概率**，并累加概率加权质量；**不会提前丢弃任何
   低概率中间项**（内部全程使用 `Decimal` 高精度运算，精度 80 位有效数字）。
5. 每个偏移的质心 = 该偏移下所有原子配置的概率加权平均精确质量；仅输出概率
   严格非零的峰，按 `offset` 升序排列。
6. 丰度允许存在 JSON 浮点舍入误差：合计与 1 的绝对偏差不超过 `1e-9` 即视为
   合法，进入卷积前自动归一化，因此响应中 `total_probability` 恒为 `1.0`。

## 一条命令启动（Docker Compose）

需要已安装 Docker 与 Docker Compose 插件：

```bash
docker compose up --build
```

该命令会构建镜像、在构建阶段安装 `requirements.txt` 中的依赖并启动服务。
服务默认监听宿主机 **8000** 端口。

后台运行：

```bash
docker compose up --build -d
```

停止：

```bash
docker compose down
```

## 配置

唯一需要的配置是**宿主机端口**。容器内固定监听 8000，宿主机端口默认也为
8000，可用环境变量 `API_PORT` 覆盖：

```bash
API_PORT=9000 docker compose up --build
```

其余无需任何环境变量或外部服务。

## 接口地址与端口

服务启动后（以默认端口为例）：

| 用途 | 地址 |
| --- | --- |
| 健康检查 | <http://localhost:8000/health> |
| 峰包络计算 | `POST http://localhost:8000/envelope` |
| Swagger 交互文档 | <http://localhost:8000/docs> |
| OpenAPI Schema | <http://localhost:8000/openapi.json> |

## 请求

`POST /envelope`，`Content-Type: application/json`：

```json
{
  "elements": [
    {
      "symbol": "C",
      "atom_count": 2,
      "isotopes": [
        {"mass_number": 12, "exact_mass": "12.000000", "abundance": 0.9893},
        {"mass_number": 13, "exact_mass": "13.003355", "abundance": 0.0107}
      ]
    }
  ]
}
```

### 字段说明

| 字段 | 类型 | 约束 |
| --- | --- | --- |
| `elements` | 数组 | 至少 1 个元素；不允许未知字段 |
| `elements[].symbol` | 字符串 | 可选，仅用于阅读，最长 16 |
| `elements[].atom_count` | 非负整数 | 所有元素原子数合计必须为 **1～60** |
| `elements[].isotopes` | 数组 | 至少 1 个同位素 |
| `isotopes[].mass_number` | 正整数 | 同元素内最大质量数 − 最小质量数 ≤ **4** |
| `isotopes[].exact_mass` | 十进制数 | 有限正数，建议以字符串传入保证十进制精度 |
| `isotopes[].abundance` | 十进制数 | `0 ≤ abundance ≤ 1`，同元素合计为 1（容差 `1e-9`） |

## 响应

`200 OK`：

```json
{
  "baseline_mass": 24.0,
  "total_probability": 1.0,
  "peaks": [
    {"offset": 0, "probability": 0.97871449, "centroid_mass": 24.0},
    {"offset": 1, "probability": 0.02117102, "centroid_mass": 25.003355},
    {"offset": 2, "probability": 0.00011449, "centroid_mass": 26.00671}
  ]
}
```

| 字段 | 说明 |
| --- | --- |
| `baseline_mass` | 所有原子均取各元素最小质量数同位素时的精确质量（u），即偏移 0 基准 |
| `total_probability` | 全部峰概率之和（归一化后为 1.0） |
| `peaks[].offset` | 相对基准的整数质量偏移，从 0 开始升序 |
| `peaks[].probability` | 该偏移合并后的概率（JSON number，IEEE-754 双精度；极小概率如低于约 `5e-324` 会在**输出序列化**时下溢为 0，但计算过程中不会被丢弃，该峰仍保留） |
| `peaks[].centroid_mass` | 该偏移的概率加权绝对精确质量质心（u） |

### 调用示例

```bash
curl -X POST http://localhost:8000/envelope \
  -H 'Content-Type: application/json' \
  -d '{
    "elements": [
      {"symbol":"C","atom_count":1,"isotopes":[
        {"mass_number":12,"exact_mass":"12.000000","abundance":0.9893},
        {"mass_number":13,"exact_mass":"13.003355","abundance":0.0107}]},
      {"symbol":"O","atom_count":2,"isotopes":[
        {"mass_number":16,"exact_mass":"15.994915","abundance":0.99757},
        {"mass_number":17,"exact_mass":"16.999132","abundance":0.00038},
        {"mass_number":18,"exact_mass":"17.999161","abundance":0.00205}]}
    ]
  }'
```

## 错误响应（422）

任何非法输入均返回 HTTP `422 Unprocessable Entity`，响应体 `detail` 数组中
每项的 `loc` 字段给出精确的错误字段位置（从 `body` 起的路径），例如：

* 原子数合计为 0 或超过 60：`loc: ["body", ""]`（模型级错误）
* 某元素丰度合计不为 1：`loc: ["body", "elements", 0, "isotopes"]`
* 质量数跨度超过 4：`loc: ["body", "elements", 0, "isotopes"]`
* 某同位素丰度为负：`loc: ["body", "elements", 0, "isotopes", 0, "abundance"]`
* 原子数为负：`loc: ["body", "elements", 0, "atom_count"]`
* 缺少必填字段 / 类型错误：定位到对应字段

## 本地开发（不使用 Docker）

需要 Python 3.11+：

```bash
pip install -r requirements.txt
uvicorn app.main:app --reload --host 0.0.0.0 --port 8000
```

运行测试：

```bash
pip install -r requirements-dev.txt
pytest -q
```

## 项目结构

```
.
├── app/
│   ├── __init__.py
│   ├── main.py        # FastAPI 路由与序列化
│   ├── schemas.py     # 请求/响应模型与 422 校验
│   └── envelope.py    # Decimal 完整卷积引擎
├── tests/
│   └── test_envelope.py
├── requirements.txt
├── requirements-dev.txt
├── Dockerfile
├── docker-compose.yml
└── README.md
```
