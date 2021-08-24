
python setup.py install

# SCRIPT -h

# SCRIPT enrich -h

# SCRIPT setting --show

SCRIPT setting --human_tf_index '/fs/home/dongxin/Projects/SCRIPT/indices/human/tf_chip_qc_1_10k'
SCRIPT setting --mouse_tf_index '/fs/home/dongxin/Projects/SCRIPT/indices/mouse/tf_chip_qc_1_10k'
SCRIPT setting --show


SCRIPT enrich -i ../test/peak_count.h5 -s hs -p ../test/project -t 24 -y

ls ../test/project/*

