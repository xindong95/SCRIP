
python setup.py install

# SCRIPT -h

# SCRIPT enrich -h

# SCRIPT setting --show

SCRIPT setting --human_chip_index '/fs/home/dongxin/Projects/SCRIPT/indices/index/human_ChIP_index' --human_motif_index '/fs/home/dongxin/Projects/SCRIPT/indices/index/human_motif_index'
SCRIPT setting --mouse_chip_index '/fs/home/dongxin/Projects/SCRIPT/indices/index/mouse_ChIP_index' --mouse_motif_index '/fs/home/dongxin/Projects/SCRIPT/indices/index/mouse_motif_index'
SCRIPT setting --show

SCRIPT enrich -e ../test/sc_experiment.h5ad -i ../test/peak_count.h5 -s hs -p ../test/project --cell_number auto --peak_confidence auto --bg_iter 1000 --reference integration -t 8 -y

ls ../test/project/*

