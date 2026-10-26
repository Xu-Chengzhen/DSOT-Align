# DSOT-Align

**基于 D/S 有序树比对的 RNA 二级结构聚类与分类方法**  
*D/S Ordered-Tree Alignment for RNA Secondary Structure Clustering and Classification*

DSOT-Align 将 RNA 二级结构表示为由 **D 节点（多分枝结构）** 和 **S 节点（茎区结构）** 组成的有序树，并在此基础上计算结构相似度。项目提供四项主要功能：RNA 点括号结构转换、有序子树动态规划（DP）比对、家族内部完全连接层次聚类及 Medoid/共有拓扑提取、基于代表模板和浅层 Wrapper 修正的 RNA 分类。
DSOT Align represents the RNA secondary structure as an ordered tree consisting of * * D nodes (multi branched structure) * * and * * S nodes (stem structure) * *, and calculates structural similarity based on this. The project provides four main functions: RNA bracket structure conversion, ordered subtree dynamic programming (DP) alignment, family fully connected hierarchical clustering and Medoid/shared topology extraction, and RNA classification based on representative templates and shallow wrapper correction.

本仓库是面向**核心算法使用**包含的评分与聚类，分类方法。
This code is designed for use with * * core algorithms, including scoring, clustering, and classification methods

## 1. 项目目录 / 1. Project Catalog

```text
DSOT-Align/
├── README.md
├── README_CN.md                  # 中文说明（本文件）
├── pyproject.toml
├── .gitignore
├── src/treernaalign_ds/
│   ├── __init__.py
│   ├── model.py                 # 有序树与节点数据模型 /Ordered Tree and Node Data Model
│   ├── io.py                    # RNA CSV 数据读取 /RNA CSV data reading
│   ├── conversion.py            # 点括号结构转换 /Point bracket structure conversion
│   ├── alignment.py             # 有序子树 DP 比对 /Ordered subtree DP comparison
│   ├── wrapper.py               # 浅层 Wrapper 修正 /Shallow Wrapper Correction
│   ├── clustering.py            # 聚类、Medoid 与共有拓扑提取 /Clustering, Medoid, and Common Topology Extraction
│   ├── classification.py        # 基于模板的分类 /Template based classification
│   └── cli.py                   # 统一命令行入口 /Unified Command Line Entrance
└── data/
    ├── rnastralign/
    │   ├── 16SrRNA.csv
    │   ├── group_I_intron.csv
    │   ├── RNaseP.csv
    │   ├── SRP.csv
    │   └── tmRNA.csv
    └── medoids_251.csv          # 251 条 Medoid 分类模板 /251 Medoid classification templates
```


## 2. 环境要求与安装 /2. Environmental requirements and installation

- Python **3.10 或更高版本**/3.10 or higher version
- NumPy
- SciPy

在项目根目录运行： /Run in the root directory of the project:

```bash
python -m pip install -e .
dsot-align --help
```

程序提供四个子命令(The program provides four sub commands)：`convert`、`align`、`cluster` 和 `classify`。

在 Windows PowerShell 中，如果通过项目根目录的 `.venv` 安装，可使用以下可执行文件路径： 
/In Windows PowerShell, if installed through `. venv ` in the project root directory, the following executable file path can be used:

```powershell
.\.venv\Scripts\dsot-align.exe --help
```

## 3. 基本使用方法  /3. Basic usage methods

### 3.1 `convert`：将点括号结构转换为 D/S 有序树  /Convert: Convert the bracket structure into a D/S ordered tree

```bash
dsot-align convert --structure '((..((...))..))'
```

程序以 JSON 格式输出结构转换结果，包括 D/S 节点、节点区间、有序子节点、是否存在虚拟根节点，以及是否采用特殊 Single-S 表示。
/The program outputs the structural transformation results in JSON format, including D/S nodes, node intervals, ordered child nodes, whether there is a virtual root node, and whether a special Single-S representation is used.


**注意：** 节点区间采用**从 0 开始的闭区间索引**（zero-based，包含起止位置）。
/Note: The node interval uses a closed interval index starting from 0 (zero based, including start and end positions).

### 3.2 `align`：计算两条 RNA 二级结构的相似度 /Align: Calculate the similarity between two RNA secondary structures

```bash
dsot-align align --query '((..((...))..))' --template '(((..((...))..)))' --trace
```

默认执行**严格、根节点锚定的有序子树动态规划比对**（Default execution of * * strict, root node anchored ordered subtree dynamic programming comparison * *.）。 
可添加 `--wrapper`，允许对一侧满足条件的浅层 D-under-D 层级进行一次候选收缩；Wrapper 默认罚分为 **0.10**。
/Add '-- wrapper' to allow for a candidate contraction of shallow D-under-D layers that meet the conditions on one side; The default penalty for the Wrapper is * * 0.10 * *.


输出包括归一化结构相似度（范围为 `[0,1]`）。使用 `--trace` 时，还会输出匹配的节点索引对。
/The output includes normalized structural similarity (within the range of '[0,1]'). When using ` -- trace `, matching node index pairs will also be output.


### 3.3 `cluster`：家族内部聚类及 Medoid 选择 /`Cluster: Internal family clustering and Medoid selection

对 SRP 家族执行聚类： /Perform clustering on the SRP family:

```bash
dsot-align cluster --data-dir data/rnastralign --family SRP --cutoff 0.75 --out-dir out/srp
```

主要输出文件： /Main output files:

- `clusters.csv`：聚类簇信息； /Cluster information;
- `assignments.csv`：各 RNA 的簇归属； /Cluster assignment of each RNA;
- `medoids.csv`：满足保留条件的簇代表 RNA。 /Clusters that meet the retention criteria represent RNA.

添加 `--consensus` 可进一步根据 **DP 回溯结果**提取共有拓扑的节点与边的支持度，生成 `consensus_nodes.csv`。默认对**至少包含 10 条 RNA** 的簇进行共有拓扑提取，支持度阈值为 **0.80**。
/Adding '-- consensus' can further extract the support of nodes and edges with a shared topology based on the DP backtracking results, and generate' consensus_nodes. csv '. By default, common topology extraction is performed on clusters containing at least 10 RNAs, with a support threshold of 0.80.


该功能**仅提取结构数据，不生成图片**。/This feature only extracts structural data and does not generate images.

不指定 `--family SRP` 时，程序将分别对五个已知 RNA 家族执行内部聚类。对于样本数量较多的家族，两两动态规划比对可能消耗较多计算资源。如果只需要分类，可以直接使用仓库中已提供的 251 条 Medoid，无需重新聚类。
/When the '-- family SRP' is not specified, the program will perform internal clustering on each of the five known RNA families. For families with a large sample size, pairwise dynamic programming comparisons may consume more computational resources. If only classification is needed, the 251 Medoids already provided in the warehouse can be used directly without the need for re clustering.


### 3.4 `classify`：RNA 家族分类  /`Classify `: RNA family classification

分类单条 RNA： /Classification of single RNA:

```bash
dsot-align classify --structure '((..((...))..))'
```

批量分类 CSV 文件： /Batch classification CSV file:

```bash
dsot-align classify --csv my_queries.csv --out out/predictions.csv
```

分类器默认使用 `data/medoids_251.csv` 中的 **251 条 Medoid 模板**，并沿用现有 Wrapper 分类逻辑。
/The classifier defaults to using * * 251 Medoid templates * * from 'data/medoids_251. csv' and follows the existing Wrapper classification logic.

- 使用 `--templates PATH` 可以指定其他模板库；
- 使用 `--no-wrapper` 可以关闭普通 D/S 树分类中的 Wrapper；
- 对触发特殊 **Single-S 表示**的查询结构，现有实现使用专门的 S-tree 严格 DP 路径，**不执行 Wrapper 候选评估**。
-Using ` -- templates PATH ` allows you to specify other template libraries;
-Using '-- no wrapper' can disable the wrapper in regular D/S tree classification;
-For query structures that trigger special Single-S representations, existing implementations use specialized S-tree strict DP paths and do not perform wrapper candidate evaluation.

预测结果包含 RNA 家族标签、最佳 Medoid 模板 ID、结构相似度和 Wrapper 相关信息。
The predicted results include RNA family tags, optimal Medoid template ID, structural similarity, and information related to the wrapper.

## 4. CSV 数据格式 /4. CSV data format

输入 RNA 数据与 Medoid 模板使用以下基本列名： /Enter RNA data and Medoid template using the following basic column names:

```csv
id,RNA,Structure,file_name
```

| 列名 | 含义 |
|---|---|
| `id` | RNA 记录标识符 |
| `RNA` | 核苷酸序列；未知家族 RNA 分类时可留空 |
| `Structure` | RNA 二级结构的点括号表示 |
| `file_name` | RNA 家族标签；待分类 RNA 可不填写 |

| Column Name | Meaning
| `id` | RNA record identifier|
| `RNA` | nucleotide sequence; When classifying unknown family RNA, it can be left blank|
| `Structure` | Point bracket representation of RNA secondary structure|
| `file_name` | RNA family tags; RNA to be classified can be left blank|

`Structure` 使用 `(`、`)` 和 `.` 表示 RNA 二级结构，**不支持假结**。
/`Structure` uses (,) and  `.`Indicates RNA secondary structure, * * does not support false knots * *.

Medoid 模板文件还包含以下字段： /The Medoid template file also includes the following fields:

| 字段 | 含义 |
|---|---|
| `cluster` | 所属聚类簇 |
| `cluster_size` | 簇内 RNA 数量 |
| `source_row` | 对应原始输入 CSV 的行号 |
| `medoid_avg_similarity` | Medoid 与簇内其他成员的平均结构相似度 |

|Field | Meaning|
|---|---|
|`cluster`  | Belonging cluster|
|`cluster_size` | Number of RNAs within the cluster|
|` source_row ` | corresponds to the line number of the original input CSV|
|`medoid_avg_similarity`  | The average structural similarity between Medoid and other members within the cluster|

`source_row` 采用**从 1 开始、且包含表头的 CSV 实际行号**：第 2 行表示第一条数据记录。保留此字段是为了区分部分 RNA 家族中可能重复出现的 RNA ID。
/`source_row 'uses * * CSV actual line numbers starting from 1 and including the header * *: the second line represents the first data record. Retaining this field is to distinguish RNA IDs that may be duplicated in some RNA families.


## 5. 核心算法与默认参数 /5. Core algorithms and default parameters

**D/S 有序树构建。** 保留节点类型、节点跨度、直接父子关系以及同级子节点的 5′→3′ 顺序。若存在多个互不包含的顶层结构块，则构建虚拟 D 根节点；对于 Single-S 表示，使用特殊浅层 S-tree 处理。
**D/S ordered tree construction. **Retain node types, node spans, direct parent-child relationships, and the 5 '→ 3' order of same level child nodes. If there are multiple top-level structural blocks that do not contain each other, construct a virtual D root node; For Single-S representation, use special shallow S-tree processing.


**节点相似度。** 仅允许相同 D/S 类型的节点直接匹配。局部节点得分为：
/**Node similarity. **Only nodes of the same D/S type are allowed to directly match. The local node score is:

\[
s(u,v)=2\exp\left(-2\frac{|L_u-L_v|}{\max(L_u,L_v,1)}\right)
\]

其中，\(L_u\) 和 \(L_v\) 分别表示两个待比较节点的跨度。/Among them\ (L_u \) and \ (L_v \) respectively represent the spans of two nodes to be compared.

**有序子树 DP。** 比对过程以根节点为锚点，保持直接父子关系及兄弟节点顺序。跳过某个子节点时，将跳过该节点所对应的整棵子树；被跳过节点的后代不会提升到其父节点之上继续匹配。
/**Ordered subtree DP. * * The comparison process takes the root node as the anchor point, maintaining a direct parent-child relationship and sibling node order. When skipping a sub node, the entire subtree corresponding to that node will be skipped; Descendants of skipped nodes will not be promoted above their parent nodes to continue matching.

**全局结构相似度。** 将根节点累计得分除以两棵**原始树**的节点数量之和，并将结果限制在 `[0,1]`。
/**Global structural similarity. **Divide the cumulative score of the root node by the sum of the number of nodes in two original trees, and limit the result to '[0,1]'.

**聚类。** 使用严格 DSOT-Align 结构相似度构建距离矩阵：
/**Clustering. **Constructing a distance matrix using strict DSOT Align structural similarity:

\[
d(Q,T)=1-S(Q,T)
\]

随后执行**完全连接层次聚类（complete linkage）**，默认结构相似度阈值为 **0.75**。程序根据已知家族标签分别执行**家族内部聚类**，并非将五个家族混合后进行无监督家族发现。
/Subsequently, complete linkage clustering is performed, with a default structural similarity threshold of 0.75. The program performs intra family clustering based on known family labels, rather than mixing five families for unsupervised family discovery.

**Medoid 选择。** 对样本量**不少于 4 条（\(n\geq4\)）**的簇，选取簇内与其他成员的平均结构相似度最高的真实 RNA 作为代表模板。
/**Medoid selection. **For clusters with a sample size of no less than 4 (n \ geq4 \) * *, select the real RNA with the highest average structural similarity to other members within the cluster as the representative template.

**共有拓扑提取。** 可选使用 DP 回溯结果，提取节点及直接连接边的支持度；默认对 \(n\geq10\) 的簇使用 **0.80** 支持度阈值。
/**Common topology extraction. **Optional use of DP backtracking results to extract support for nodes and directly connected edges; Default support threshold of * * 0.80 * * is used for clusters of \ (n \ geq10 \).

**分类与 Wrapper。** 将待分类 RNA 与保留的 Medoid 模板逐一比对。Wrapper 最多对其中一侧的浅层 D-under-D 结构执行候选收缩，候选得分扣除 **0.10** 罚分，并始终使用收缩前两棵树的节点总数归一化。选择最终相似度最高的模板，将该模板对应的家族作为预测结果。特殊 Single-S 查询沿用无需 Wrapper 的 S-tree 严格 DP 路径。
/**Classification and Wrapper. **Compare the RNA to be classified with the retained Medoid template one by one. The Wrapper can perform candidate shrinkage on one side of the shallow D-under-D structure at most, deduct a penalty of * * 0.10 * * from the candidate score, and always normalize the total number of nodes in the first two trees before shrinkage. Select the template with the highest similarity and use the corresponding family as the prediction result. The special Single-S query follows the strict DP path of S-tree without the need for a wrapper.


## 6. 附带数据与 Medoid 模板 /6. Attached data and Medoid templates

仓库中的五个 CSV 是从原项目整理得到的 **RNAStrAlign 五家族预处理数据**，对应当前 **4,280 条 RNA** 的聚类版本。整理过程中仅统一了列名 `id,RNA,Structure,file_name`，没有改变 RNA 序列、点括号结构或记录顺序，也没有额外添加预测生成的 RNA。
/The five CSV files in the warehouse are pre processed data from the RNAStrAlign five family, which were organized from the original project and correspond to the clustering versions of the current 4280 RNA clusters. During the sorting process, only the column names' id ',' RNA ',' Structure ', and' file_name 'were unified, without changing the RNA sequence, bracket structure, or recording order, nor adding any additional predicted generated RNA.

| RNA 家族文件 | 预处理后 RNA 数量 | 保留 Medoid 数（\(n\geq4\)） |
|---|---:|---:|
| `16SrRNA.csv` | 1,243 | 45 |
| `group_I_intron.csv` | 1,681 | 135 |
| `RNaseP.csv` | 381 | 24 |
| `SRP.csv` | 466 | 26 |
| `tmRNA.csv` | 509 | 21 |
| **合计** | **4,280** | **251** |

|RNA Family Files | Pre processed RNA Quantity | Retained Medoid Number (\ (n \ geq4 \))|
|---|---:|---:|
| `16SrRNA.csv` | 1,243 | 45 |
| `group_I_intron.csv` | 1,681 | 135 |
| `RNaseP.csv` | 381 | 24 |
| `SRP.csv` | 466 | 26 |
| `tmRNA.csv` | 509 | 21 |
| **Total ** | **4,280** | **251** |


`data/medoids_251.csv` 根据原项目最新的聚类结果整理：在相似度阈值 **0.75** 下，五个家族共形成 **587 个簇**，其中 **251 个簇**满足 \(n\geq4\) 的模板保留条件。每条 Medoid 均根据所属家族的输入 CSV 行号及聚类成员关系进行了核对。
/`According to the latest clustering results of the original project, 'data/medoids-251. csv' is organized as follows: under the similarity threshold of * * 0.75 * *, a total of * * 587 clusters * * are formed among the five families, of which * * 251 clusters * * meet the template retention condition of \ (n \ geq4 \). Each Medoid has been verified based on the input CSV line number and cluster membership relationship of its family.





## 7. 仓库范围

本仓库包含：

- 核心算法；
- 五个 RNA 家族的预处理数据；
- 251 条 Medoid 代表模板；
- 一个统一的命令行入口。
## 7. Warehouse scope
This warehouse contains:
-Core algorithm;
-Preprocessed data of five RNA families;
-251 Medoid representative templates;
-A unified command-line entry point.
