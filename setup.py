#!/usr/bin/env python

from distutils.core import setup

setup(name='REDItools',
	version='1.3',
	description='Python Scripts for RNA editing detection by RNA-Seq data',
	author='Ernesto Picardi',
	author_email='ernesto.picardi@gmail.com',
	url='https://github.com/BioinfoUNIBA/REDItools',
	scripts=['main/REDItoolDenovo.py',
	'main/REDItoolDnaRna.py',
	'main/REDItoolKnown.py',
	'accessory/AnnotateTable.py',
	'accessory/FilterTable.py',
	'accessory/SearchInTable.py',
	'accessory/selectPositions.py',
	'accessory/GFFtoTabix.py',
	'accessory/SortGFF.py',
	'accessory/SortTable.py',
	'accessory/TableToGFF.py',
	'accessory/tableToTabix.py',
	'accessory/readPsl.py',
	'accessory/subCount.py',
	'accessory/subCount2.py',
	'accessory/rediportal2recoding.py'
	],
	license='LICENSE.txt',
	classifiers=[
          'Intended Audience :: Computational biologists',
          'License :: OSI Approved :: MIT',
          'Operating System :: MacOS :: MacOS X',
          'Operating System :: POSIX',
          'Programming Language :: Python',
          ],
	long_description=open('README').read(),
	platforms=['Linux','Unix','MacOS']
)

