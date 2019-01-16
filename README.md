## 流程介绍

### KEGG/GO Enrichment
KEGG/GO Analysis for hsa with R Code

#### 第一步:准备文件

###### (1)单样本输入文件

样本名字 |基因名字
:-------|:------
SampleName | GeneName
PP_PN | PLEKHB1
PP_PN | NELL2
PP_PN | DPP4
PP_PN | LEF1
PP_PN | LRRK1
PP_PN | ZNF683

######  (2)多样本输入文件

样本名字 |基因名字
:-------|:------
SampleName | GeneName
F17092477270-KY409 | CREBBP
F17092477270-KY409 | CYP2D6
F17092477270-KY409 | DUSP2
F17092477270-KY409 | NF2
F17092477270-KY409 | PIK3CA
F17092477270-KY409 | PMS1
....               |...
F17092477271-KY409 | CREBBP
F17092477271-KY409 | CYP2D6
F17092477271-KY409 | DUSP2
F17092477271-KY409 | NF2

##### 第二步: 代码参数说明

######  2.1`plot_go_kegg_plot_2.R` 参数说明
```text
1. inputName: 输入文件（'\t'分隔的文件）
2. outdir: 输出目录
3. pmodule: 1 为p.ajust作为筛选条件 0为pvalue为筛选条件
4. pvalue：p值小于0.05为筛选条件，筛选显著性通路
5. width：图片的宽度
6. height：突变的高度
```

###### 2.2 `plot_go_kegg_plot_2.R` 参数说明

```text
如果需要对*.no.filter文件进行重新筛选显著性通路与作图
1. prefix: *go_*.no.filter.xlsx与*pathway.no.filter.xlsx文件的前缀
2. pmodule: 1 为p.ajust作为筛选条件 0为pvalue为筛选条件
3. pvalue：p值小于0.05为筛选条件，筛选显著性通路（默认是0.05）
4. width：图片的宽度
5. height：突变的高度
6. kegg_go_plot:1代表画所有显著性的GO term和pathway；0为画显著性的pathway
7. topnum:画top number的GO term和或pathway（默认是20）
```


`注意：如果需要重新筛选显著性通路并绘图，才需用脚本plot_go_kegg_plot_2.R`

##### 第三步：运行脚本

```shell
plot_go_kegg_enrich_plot.R inputName outdir 1 0.05 10 30
```

```shell
plot_go_kegg_plot_2.R prefix 1 0.05 10 30 1 20 

```
`注：plot_go_kegg_enrich_plot.R代码运行完成之后会输出workLOG.txt文件,该文件
打印出通过p.ajust或pvalue筛选后显著性通路/GO term数量`
