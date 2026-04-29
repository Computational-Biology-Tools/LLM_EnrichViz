# LLM_EnrichViz

> **LLM-Powered Multi-Method Gene Set Enrichment Analysis and Visualization**  
> **基于大语言模型的多方法基因集富集分析与可视化**

---

## Table of Contents · 目录

- [Overview · 概述](#overview--概述)
- [Features · 功能特性](#features--功能特性)
- [Installation · 安装](#installation--安装)
- [Quick Start · 快速开始](#quick-start--快速开始)
- [Methods · 分析方法](#methods--分析方法)
- [LLM Providers · 大语言模型提供商](#llm-providers--大语言模型提供商)
- [Output · 输出文件](#output--输出文件)
- [Contributing · 贡献指南](#contributing--贡献指南)
- [License · 许可证](#license--许可证)

---

## Overview · 概述

**English**

`LLM_EnrichViz` is an R package that combines four gold-standard gene set
enrichment methods — **ORA**, **GSEA**, **GSVA**, and **ROAST** — with the
interpretive intelligence of state-of-the-art Large Language Models.
Each method × database pair is handled by an independent LLM *agent*
running in parallel.  A final synthesis call integrates all agent outputs into
a unified mechanistic hypothesis, a confidence-ranked pathway list, and an
interactive network visualisation.

**中文**

`LLM_EnrichViz` 是一个 R 软件包，将四种主流基因集富集分析方法——**ORA**、
**GSEA**、**GSVA** 和 **ROAST**——与最先进大语言模型的解释能力相结合。  
每种方法与数据库的组合由独立的 LLM *智能体*并行处理。最终，一次综合调用将
所有智能体的输出整合为统一的机制假说、置信度排序的通路列表以及交互式网络可视化。

---

## Features · 功能特性

**English** | **中文**

| Feature | 功能 |
|---------|------|
| ✅ Four enrichment methods (ORA, GSEA, GSVA, ROAST) | ✅ 四种富集分析方法 |
| ✅ Four pathway databases (KEGG, GO, Reactome, WikiPathways) | ✅ 四个通路数据库 |
| ✅ Four LLM providers (Gemini, OpenAI, Claude, DeepSeek) | ✅ 四大语言模型提供商 |
| ✅ Parallel multi-agent LLM execution | ✅ 并行多智能体执行 |
| ✅ Cross-method confidence tiers | ✅ 跨方法置信度分级 |
| ✅ Interactive `visNetwork` HTML network | ✅ 交互式 HTML 网络图 |
| ✅ Automated HTML report via RMarkdown | ✅ 自动化 HTML 报告 |
| ✅ Multi-species support (14 organisms) | ✅ 多物种支持（14种生物） |
| ✅ SYMBOL / ENSEMBL / ENTREZID input | ✅ 多种基因 ID 格式输入 |

---

## Installation · 安装

### Dependencies · 依赖项

**English:** Install Bioconductor and CRAN dependencies first.  
**中文：** 请先安装 Bioconductor 和 CRAN 依赖项。

```r
# Bioconductor packages · Bioconductor 包
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install(c(
  "clusterProfiler",
  "ReactomePA",
  "GSVA",
  "limma"
))

# Install organism-specific package based on your species
# 根据您的物种安装相应的 org.db 包
# Examples: org.Hs.eg.db (human), org.Mm.eg.db (mouse), org.Rn.eg.db (rat), etc.
# 示例：org.Hs.eg.db（人类）、org.Mm.eg.db（小鼠）、org.Rn.eg.db（大鼠）等

# CRAN packages · CRAN 包
install.packages(c(
  "msigdbr", "httr", "R6", "future", "furrr",
  "progressr", "dplyr", "tidyr", "ggplot2",
  "visNetwork", "stringr", "rmarkdown", "tibble"
))
```

### Install LLM_EnrichViz · 安装软件包

```r
# From GitHub · 从 GitHub 安装
devtools::install_github("Computational-Biology-Tools/LLM_EnrichViz")
```

---

## Quick Start · 快速开始

### Set up your API key · 配置 API 密钥

**English:** Store your API key in a plain-text file.  
**中文：** 将 API 密钥存储在纯文本文件中。

```bash
# Gemini
echo "AIza..." > ~/.gemini_key

# OpenAI / ChatGPT
echo "sk-..."  > ~/.openai_key

# Claude (Anthropic)
echo "sk-ant-..." > ~/.claude_key

# DeepSeek
echo "sk-..."  > ~/.deepseek_key
```

### ORA only · 仅运行 ORA

```r
library(LLM_EnrichViz)

LLM_EnrichViz(
  deg_list            = c("TP53", "MYC", "BRCA1", "EGFR", "VEGFA",
                           "CDK2", "MDM2", "BCL2", "PTEN", "RB1"),
  gene_type           = "SYMBOL",
  organism            = "human",  # Supported: human, mouse, rat, chicken, zebrafish, fly, worm, yeast, arabidopsis, pig, dog, cow, ecoli
  analyses            = "ORA",
  databases           = c("KEGG", "GO"),
  llm_provider        = "openai",
  api_key_file        = "~/.openai_key",
  experimental_design = "Breast cancer vs normal, bulk RNA-Seq",
  output_dir          = "ora_results"
)
```

### Full pipeline · 完整流程

```r
library(LLM_EnrichViz)

LLM_EnrichViz(
  expression_data     = expr_matrix,        # genes × samples
  ranked_list         = sort(wald_stats, decreasing = TRUE),
  deg_list            = sig_genes,
  gene_type           = "SYMBOL",
  organism            = "human",  # Supported: human, mouse, rat, chicken, zebrafish, fly, worm, yeast, arabidopsis, pig, dog, cow, ecoli
  design_matrix       = model.matrix(~ condition, data = sinfo),
  contrast_vector     = c(0, 1),
  analyses            = c("ORA", "GSEA", "GSVA", "ROAST"),
  databases           = c("KEGG", "GO", "Reactome"),
  llm_provider        = "claude",
  llm_model           = "claude-sonnet-4-6",
  api_key_file        = "~/.claude_key",
  experimental_design = "Lupus PBMC vs healthy controls, RNA-Seq, n=80",
  output_dir          = "lupus_enrichment"
)
```

---

## Methods · 分析方法

### ORA — Over-Representation Analysis · 过度代表性分析

**English:** Tests whether DEGs are statistically over-represented in gene sets
using a hypergeometric test.  **Input required:** `deg_list`.

**中文：** 使用超几何检验测试 DEG 是否在基因集中统计显著过度代表。  
**所需输入：** `deg_list`

| Pros · 优点 | Cons · 缺点 |
|-------------|-------------|
| Fast and simple · 快速简单 | Ignores expression magnitude · 忽略表达量 |
| Widely used standard · 广泛使用的标准 | Threshold-dependent · 依赖阈值 |

---

### GSEA — Gene Set Enrichment Analysis · 基因集富集分析

**English:** Uses a full ranked gene list to compute enrichment scores.
NES > 0 = activated; NES < 0 = repressed.  **Input required:** `ranked_list`
(named, sorted numeric vector).

**中文：** 使用完整排序基因列表计算富集评分。NES > 0 = 激活；NES < 0 = 抑制。  
**所需输入：** `ranked_list`（命名的排序数值向量）

| Pros · 优点 | Cons · 缺点 |
|-------------|-------------|
| Uses all genes · 使用所有基因 | Sensitive to ranking metric · 对排序指标敏感 |
| Captures subtle signals · 捕获微弱信号 | Slower than ORA · 比 ORA 慢 |

---

### GSVA — Gene Set Variation Analysis · 基因集变异分析

**English:** Computes per-sample pathway activity scores, enabling study of
inter-sample heterogeneity.  **Input required:** `expression_data` (matrix).

**中文：** 计算每个样本的通路活性评分，用于研究样本间异质性。  
**所需输入：** `expression_data`（表达矩阵）

| Pros · 优点 | Cons · 缺点 |
|-------------|-------------|
| Sample-level resolution · 样本级分辨率 | Relative scores · 相对评分 |
| Detects sub-populations · 检测亚群 | Computationally heavy · 计算密集 |

---

### ROAST — Rotation-based Gene Set Testing · 基于旋转的基因集检验

**English:** Tests for coordinated directional expression using random rotations
of residuals.  Outputs Up / Down / Mixed direction.
**Inputs required:** `expression_data`, `design_matrix`, `contrast_vector`.

**中文：** 通过随机旋转残差检验协调的方向性表达。  
输出方向：Up（上调）/ Down（下调）/ Mixed（混合）。  
**所需输入：** `expression_data`、`design_matrix`、`contrast_vector`

| Pros · 优点 | Cons · 缺点 |
|-------------|-------------|
| Accounts for gene correlations · 考虑基因相关性 | Needs design matrix · 需要设计矩阵 |
| Directional output · 方向性输出 | Computationally intensive · 计算密集 |

---

## LLM Providers · 大语言模型提供商

| Provider · 提供商 | Default model · 默认模型 | API endpoint · 接口 |
|-------------------|--------------------------|---------------------|
| **Gemini** (Google) | `gemini-1.5-flash-latest` | generativelanguage.googleapis.com |
| **OpenAI** (ChatGPT) | `gpt-4o` | api.openai.com |
| **Claude** (Anthropic) | `claude-sonnet-4-6` | api.anthropic.com |
| **DeepSeek** | `deepseek-chat` | api.deepseek.com |

**English:** All providers use a unified interface — switch providers by
changing `llm_provider` and `api_key_file` with no other code changes required.

**中文：** 所有提供商使用统一接口——仅需更改 `llm_provider` 和 `api_key_file`，
无需修改其他代码即可切换模型。

---

## Output · 输出文件

**English** | **中文**

| File · 文件 | Description · 描述 |
|-------------|-------------------|
| `ORA_{DB}_results.txt` | ORA enrichment table · ORA 富集结果表 |
| `GSEA_{DB}_results.txt` | GSEA table with NES · 含 NES 的 GSEA 结果 |
| `GSVA_{DB}_summary.txt` | GSVA score summary · GSVA 评分摘要 |
| `GSVA_{DB}_scores_matrix.txt` | Full per-sample GSVA matrix · 样本级 GSVA 矩阵 |
| `ROAST_{DB}_results.txt` | ROAST directional table · ROAST 方向性结果 |
| `ORA_{DB}_barplot.pdf` | ORA bar chart · ORA 柱状图 |
| `GSEA_{DB}_bubble.pdf` | GSEA bubble plot · GSEA 气泡图 |
| `GSVA_{DB}_heatmap.pdf` | GSVA heatmap · GSVA 热图 |
| `ROAST_{DB}_barplot.pdf` | ROAST bar chart · ROAST 柱状图 |
| `agent_responses.txt` | All LLM agent outputs · 所有智能体输出 |
| `final_response.txt` | Integrated LLM synthesis · 综合分析结果 |
| `integrated_network.html` | Interactive network · 交互式网络图 |
| `analysis_parameters.txt` | Run configuration · 运行配置 |
| `report_template.html` | Full HTML report · 完整 HTML 报告 |

---

## Pipeline Architecture · 流程架构

```
Input data · 输入数据
     │
     ▼
Gene ID mapping (bitr) · 基因 ID 转换
     │
     ├──── ORA  (enrichKEGG / enrichGO / enrichPathway / enrichWP)
     ├──── GSEA (gseKEGG / gseGO / gsePathway)
     ├──── GSVA (gsva)
     └──── ROAST (mroast)
               │
               ▼
     LLM Prompts (method × database) · 提示词构建
               │
               ▼
     Parallel LLM Agents (future/furrr) · 并行智能体
               │
               ▼
     Integration Prompt · 综合分析提示词
               │
               ▼
     Final LLM Response · 最终综合响应
               │
        ┌──────┴──────┐
        ▼             ▼
  Visualizations   HTML Report
  可视化图表      HTML 报告
```

---

## Contributing · 贡献指南

**English:** Contributions are welcome!  Please open an issue to discuss
proposed changes before submitting a pull request.  Run `devtools::check()`
and ensure all tests pass before submitting.

**中文：** 欢迎贡献代码！在提交拉取请求之前，请先开 Issue 讨论您的修改方案。  
提交前请运行 `devtools::check()` 并确保所有测试通过。

```r
# Run tests locally · 本地运行测试
devtools::test()

# Check package · 检查包
devtools::check()
```

---

## Citation · 引用

**English:** If you use LLM_EnrichViz in your research, please cite:

**中文：** 如果您在研究中使用了 LLM_EnrichViz，请引用：

```
Your Name (2024). LLM_EnrichViz: LLM-Powered Multi-Method Gene Set
Enrichment Analysis and Visualization. R package version 0.1.0.
https://github.com/yourname/LLM_EnrichViz
```

---

## License · 许可证

MIT © 2024 Your Name

**English:** This software is provided "as is", without warranty of any kind.
The authors are not responsible for any conclusions drawn from LLM outputs.
Always verify LLM-generated biological interpretations against primary
literature.

**中文：** 本软件按"现状"提供，不提供任何形式的保证。  
作者对 LLM 输出所得出的任何结论不承担责任。请始终通过原始文献验证 LLM
生成的生物学解释。

---

<div align="center">

**LLM_EnrichViz** — *From statistics to mechanistic insight, automatically.*  
**LLM_EnrichViz** — *从统计数据到机制洞察，全程自动化。*

</div>
