import os
import sys
import ruamel.yaml
from SCRIP.utilities.utils import read_config


def update_setting( args ):
    yaml = ruamel.yaml.YAML()
    yaml.default_flow_style = False
    CONFIG, CONFIG_PATH = read_config()
    if args.show:
        sys.stdout.write('The reference indices you set:\n')
        yaml.dump(CONFIG['index'], sys.stdout)
        sys.stdout.write('\n')
    else:
        if args.human_tf_index:
            CONFIG['index']['human_tf_index'] = os.path.abspath(args.human_tf_index)
        if args.human_hm_index:
            CONFIG['index']['human_hm_index'] = os.path.abspath(args.human_hm_index)
        if args.mouse_tf_index:
            CONFIG['index']['mouse_tf_index'] = os.path.abspath(args.mouse_tf_index)
        if args.mouse_hm_index:
            CONFIG['index']['mouse_hm_index'] = os.path.abspath(args.mouse_hm_index)
        with open(os.path.join(CONFIG_PATH, 'config.yml'), 'w+') as config_file:
            yaml.dump(CONFIG, config_file)
        print('Update sucessfully!')
