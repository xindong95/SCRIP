def search_giggle(bed_path, result_path, index_path):
    bed = bed_path.split("/")[-1]
    if bed_path.endswith('gz'):
        cmd = 'giggle search -i {index_path} -s -q {bed_path} > {result_path}\n'.format(index_path=index_path, result_path=result_path, bed_path=bed_path)
    else:
        cmd = 'sort --buffer-size 2G -k1,1 -k2,2n -k3,3n {bed_path} | bgzip -c > {bed_path}.gz\n'.format(bed_path=bed_path)
        cmd += 'giggle search -i {index_path} -s -q {bed_path}.gz > {result_path}\n'.format(index_path=index_path, result_path=result_path, bed_path=bed_path)
        cmd += 'rm {bed_path}'.format(bed_path=bed_path)
    subprocess.run(cmd, shell=True, check=True)

def search_giggle_batch(bed_folder, result_folder, index_path, thread=8):
    if not os.path.exists(result_folder):
        os.makedirs(result_folder)
    beds = os.listdir(bed_folder)
    args = []
    for bed in beds:
        barcodes = bed.split('.')[0]
        args.append((bed_folder + '/' + bed, result_folder + '/' + barcodes + '.txt', index_path))
    with Pool(thread) as p:
        p.starmap(search_giggle, args)

def read_giggle_result(path, filename_split="_", index_prefix="top10k/"):
    """For giggle stored path, return a table, col is cell cluster / cell and row is factor"""
    for i in range(len(os.listdir(path))):
        giggle_result = os.listdir(path)[i]
        cell_bc = giggle_result.split(filename_split)[0]
#         print(cell_bc)
        dtframe = pd.read_csv(os.path.join(path, giggle_result), sep="\t", index_col=None, comment="#", header=None)
        dtframe.columns = ["file", "file_size","overlaps","odds_ratio","fishers_two_tail","fishers_left_tail","fishers_right_tail","combo_score","NA"]
        if i == 0:
            dtframe = dtframe.loc[:,["file", "combo_score"]]
            total = dtframe.rename(columns={'combo_score':cell_bc}).copy()
        else:
            newcol = dtframe[["combo_score"]]
#             print(newcol)
            total[cell_bc] = newcol
    idList=[]
    for i in total['file']:
        idList.append(i.replace(index_prefix,"").replace(".bed.gz",""))
    total = total.rename(columns={"file":"id"})
    total["id"] = idList
    total = total.set_index("id")
    return total
