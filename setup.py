#!/usr/bin/python3
from setuptools import setup

setup(
   name='ExonSurfer',
   version='1.3',
   author='Pablo Monfort-Lanzas, Cristina Rusu',
   author_email='pablo.monfort@i-med.ac.at',
   packages=['ExonSurfer','ExonSurfer.specificity','ExonSurfer.ensembl',\
      'ExonSurfer.primerDesign','ExonSurfer.resources','ExonSurfer.visualization',\
       'ExonSurfer.readFiles','ExonSurfer.dimers'],
   scripts=['bin/exon_surfer.py',],
   url='https://github.com/CrisRu95/ExonPrimerSurfer/',
   license='LICENSE.txt',
   description='',
   long_description=open('README.md').read(),
   long_description_content_type='text/markdown',
   install_requires=[
   ],
   include_package_data=True,
   package_data={},
)
