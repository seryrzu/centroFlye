import os
import errno


def list_only_dirs(path):
    for file in os.listdir(path):
        if os.path.isdir(os.path.join(path, file)):
            yield file


def list_only_files(path):
    for file in os.listdir(path):
        if os.path.isfile(os.path.join(path, file)):
            yield file


def smart_mkdir(dirname):
    try:
        os.mkdir(dirname)
    except OSError as exc:
        if exc.errno != errno.EEXIST:
            raise exc


def smart_makedirs(dirname):
    try:
        os.makedirs(dirname)
    except OSError as exc:
        if exc.errno != errno.EEXIST:
            raise exc
