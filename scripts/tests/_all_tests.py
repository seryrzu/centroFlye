# (c) 2020 by Authors
# This file is a part of centroFlye program.
# Released under the BSD license (see LICENSE file)

import logging

import os
import sys
import unittest

this_dirname = os.path.dirname(os.path.realpath(__file__))
sys.path.append(os.path.join(this_dirname, os.path.pardir))


from monomers_tests import MonomersTests
from sd_parser_tests import SDParserTests, SDParserWOHPCTests
from graph_tests import DBGraphTests
from standard_logger import get_logger
from utils.git import get_git_revision_short_hash


def main():
    logger = get_logger('test.log',
                        logger_name='centroFlye',
                        level=logging.DEBUG,
                        filemode='w',
                        stdout=False)
    logger.info(f'cmd: {sys.argv}')
    logger.info(f'git hash: {get_git_revision_short_hash()}')

    unittest.main()


if __name__ == "__main__":
   main()
