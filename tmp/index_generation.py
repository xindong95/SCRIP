QC = pd.read_csv("/mnt/Storage/home/dongxin/temp/cistrome/human_factor_full_QC.txt", sep = "\t")
factor_list = QC.Factor.unique().tolist()
for f in factor_list:
    if len(QC[(QC["Factor"] == f)]) > 2:
        need_remove_table_factor = QC[(QC["Factor"] == f)]
        cell_line_list = need_remove_table_factor.Cell_line.unique().tolist()
        for cl in cell_line_list:
            if len(need_remove_table_factor[need_remove_table_factor["Cell_line"] == cl]) > 2:
                need_remove_table_cl = need_remove_table_factor[need_remove_table_factor["Cell_line"] == cl]
                remove_id = list(set(need_remove_table_cl.index.tolist()) - set(need_remove_table_cl.sort_values(by="PeaksFoldChangeAbove10", ascending=False).head(2).index.tolist()))
                QC = QC.drop(index=remove_id)
            else:
                continue
    else:
        continue


# rename index bed file
for i in os.listdir("/mnt/Storage/home/dongxin/temp/cistrome/filter_bed/fiveFoldPeak/"):
    dcid = int(i.split(".")[0])
    factor = QC[QC["DCid"] == dcid]["Factor"].to_string().split()[1]
    if "/" in factor or "_" in factor:
        factor = factor.replace("/", "or").replace("_", "")
#     Cell_line = QC[QC["DCid"] == dcid]["Cell_line"].to_string().split()[1]
#     Cell_type = 
#     Tissue_type
    os.system("mv /mnt/Storage/home/dongxin/temp/cistrome/filter_bed/fiveFoldPeak/%s.bed.gz /mnt/Storage/home/dongxin/temp/cistrome/filter_bed/fiveFoldPeak/%s_%s.bed.gz" % (dcid,dcid,factor))



annotation_path = '/mnt/Storage2/home/dongxin/Projects/scATAC/reference/giggle.all/CistromeDB.sample.annotation.txt'
anno_table = pd.read_csv(annotation_path,sep="\t")
factorid = list(map(str,anno_table[(anno_table["FactorType"]=='tf')]["ID"].tolist()))
factor_dict = anno_table[(anno_table["FactorType"]=='tf')][["ID","Factor"]].set_index("ID").to_dict()["Factor"]
# map factor by id "_"
factor_index_list = []
for i in result_p.index:
    factor_name = i.split("_")
    if len(factor_name) != 2:
        try:
            factor_index_list.append(factor_dict[int(factor_name[0])])
        except:
            print(factor_name)
            factor_index_list.append("None")
    else:
        factor_index_list.append(factor_name[1])