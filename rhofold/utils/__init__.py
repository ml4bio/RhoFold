import torch
from inspect import isfunction
import contextlib
import tempfile
import time
from typing import Optional
import shutil
import logging

from .ss_utils import *

def exists(val):
    return val is not None

def default(val, d):
    if exists(val):
        return val
    return d() if isfunction(d) else d

def get_device(device) -> str:
    """
    """
    if device is None:
        if torch.cuda.is_available():
            return "cuda"
        else:
            return 'cpu'
    elif device == 'cpu':
        return device
    elif device.startswith('cuda'):
        if torch.cuda.is_available():
            return device
        else:
            raise ValueError(f"Cuda is not available")
    else:
        raise ValueError(f"Device{device} is not available")

@contextlib.contextmanager
def tmpdir(base_dir: Optional[str] = None):
  """Context manager that deletes a temporary directory on exit."""
  tmpdir = tempfile.mkdtemp(dir=base_dir)
  try:
    yield tmpdir
  finally:
    shutil.rmtree(tmpdir, ignore_errors=True)

@contextlib.contextmanager
def timing(msg: str, logger: logging.Logger):
  logger.info('Started %s', msg)
  tic = time.time()
  yield
  toc = time.time()
  logger.info('Finished %s in %.3f seconds', msg, toc - tic)

