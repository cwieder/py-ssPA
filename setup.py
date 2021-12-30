from setuptools import setup

setup(
    name='sspa',
    version='0.1.0',
    packages=['sspa'],
    package_dir={'':'src'},
    url='https://github.com/cwieder/sspa',
    license='GNU 3.0',
    author='Cecilia Wieder',
    author_email='cw2019@ic.ac.uk',
    description='Single sample pathway analysis tools for omics data',
    package_data={'sspa': ['example_data/*', 'pathway_databases/*']},
    classifiers=[
        'Programming Language :: Python',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.10',
    ]
)
