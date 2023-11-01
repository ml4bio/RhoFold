from setuptools import setup, find_packages

setup(
  name = 'ribofold',
  packages = find_packages(),
  version = '0.0.1',
  license='',
  description = 'ribofold',
  author = '',
  author_email = '',
  url = '',
  keywords = [
    'artificial intelligence',
    'attention mechanism',
    'rna structure prediction'
  ],
  install_requires=[
  ],
  setup_requires=[
    'pytest-runner',
  ],
  tests_require=[
    'pytest'
  ],
  classifiers=[
    'Development Status :: 4 - Beta',
    'Intended Audience :: Developers',
    'Topic :: Scientific/Engineering :: Artificial Intelligence',
    'License :: OSI Approved :: MIT License',
    'Programming Language :: Python :: 3.7',
  ],
)
