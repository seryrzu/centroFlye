import logging

import os
import sys
import unittest

this_dirname = os.path.dirname(os.path.realpath(__file__))
sys.path.append(os.path.join(this_dirname, os.path.pardir))


from monomers_tests import MonomersTests
from standard_logger import get_logger


def main():
    logger = get_logger('test.log',
                        logger_name='centroFlye',
                        level=logging.DEBUG,
                        filemode='w',
                        stdout=False)
    unittest.main()


if __name__ == "__main__":
    main()
