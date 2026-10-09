```bash
project_name="Acer"

ln -s \
    "/home/liuzhiyu/Projects/neo_caste/downstream_analysis_scVI/species/${project_name}/data/6_cluster-output/concat.h5ad" \
    "/home/liuzhiyu/Projects/neo_caste/anno_cells/${project_name}/matrix/cluster_unannotated.h5ad"
```
```bash
ln -s "/data/share/data/Zhou_lab_seq_data/20251022_scRNA_lzy/1_LZU/ori_h5ad/acer_qw_concat_clustered_annotation.h5ad" \
    "/home/liuzhiyu/Projects/neo_caste/anno_cells/${project_name}/matrix/Acer_annobyHY.h5ad"
```

```bash
ln -sf "/home/liuzhiyu/Projects/neo_caste/ApisSCENIC/01_Apis_mellifera/01_create_cistarget_db/gtf/Apis_mellifera.tsv" "./metadata"
```