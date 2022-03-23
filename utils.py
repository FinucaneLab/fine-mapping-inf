import numpy as np
import pandas as pd
import os
import sys
import logging
from tqdm import tqdm

class TqdmHandler(logging.StreamHandler):
    def __init__(self):
        logging.StreamHandler.__init__(self)

    def emit(self, record):
        msg = self.format(record)
        tqdm.write(msg)

def configure_logger(out_log):

    logFormatter = logging.Formatter("[%(levelname)s]  %(message)s")
    logger = logging.getLogger()
    logger.setLevel(logging.NOTSET)

    consoleHandler = TqdmHandler()
    consoleHandler.setFormatter(logFormatter)
    logger.addHandler(consoleHandler)

    fileHandler = logging.FileHandler(out_log+'.log')
    fileHandler.setFormatter(logFormatter)
    logger.addHandler(fileHandler)

