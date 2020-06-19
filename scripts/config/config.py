import os
import yaml


def get_config():
    this_dirname = os.path.dirname(os.path.realpath(__file__))
    config_fn = os.path.join(this_dirname, 'config.yaml')
    with open(config_fn) as f:
        config = yaml.safe_load(f)
    return config
