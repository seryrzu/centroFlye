import standard_logger

from test_module import add

def main():
    logger = standard_logger.get_logger('test')
    logger.info('info message')
    logger.critical('critical message')
    logger.debug('debug message')
    add(1, 2)

main()
