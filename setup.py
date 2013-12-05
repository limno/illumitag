from distutils.core import setup

setup(
      name             = 'illumitag',
      version          = '0.2.0',
      description      = 'Pipeline for 16S Illumina paired-end sequencing',
      long_description = open('README.txt').read(),
      license          = 'MIT',
      url              = 'http://github.com/limno/illumitag/',
      author           = 'Lucas Sinclair',
      author_email     = 'lucas.sinclair@me.com',
      classifiers      = ['Topic :: Scientific/Engineering :: Bio-Informatics'],
      packages         = ['illumitag'],
      requires         = ['biopython', 'matplotlib', 'threadpool', 'sh', 'patsy', 'pandas', 'statsmodels', 'fastqident', 'rpy2'],
)