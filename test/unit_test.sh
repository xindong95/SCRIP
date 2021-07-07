
python setup.py install

# SCRIPT -h

# SCRIPT enrich -h

# SCRIPT setting --show

SCRIPT setting --human_chip_index '/fs/home/dongxin/Projects/SCRIPT/indices/human/ChIP_top10k' --human_motif_index '/fs/home/dongxin/Projects/SCRIPT/indices/human/motif_top10k_resize'
SCRIPT setting --mouse_chip_index '/fs/home/dongxin/Projects/SCRIPT/indices/mouse/ChIP_top10k' --mouse_motif_index '/fs/home/dongxin/Projects/SCRIPT/indices/mouse/motif_top10k_resize'
SCRIPT setting --show

SCRIPT enrich -e ../test/sc_experiment.h5ad -i ../test/peak_count.h5 -s hs -p ../test/project --cell_number auto --peak_confidence auto --reference integration -t 8 -y

ls ../test/project/*

