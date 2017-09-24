from distutils.core import setup
setup(
  name = 'trapes',
  packages = ['trapes'], # this must be the same as the name above
  version = '0.1',
  description = 'TCR Reconstruction Algorithm for Paired-End Single cell',
  author = 'Shaked Afik',
  author_email = 'safik@berkeley.edu',
  url = 'https://github.com/gunjanbaid/TRAPeS', # use the URL to the github repo
  download_url = 'https://github.com/gunjanbaid/TRAPeS/archive/0.1.tar.gz', # TODO Need to set up tags for original repo
  keywords = ['RNA-seq', 'T-cell'], # arbitrary keywords
  classifiers = [],
  entry_points={ 
    'console_scripts': [
      'trapes=trapes.trapes:main',
    ],
  },
)
