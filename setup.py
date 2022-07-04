from setuptools import setup
from pathlib import Path
this_directory = Path(__file__).parent
long_description = (this_directory / "README.md").read_text()


setup(
    name='sspa',
    version='0.1.2',
    packages=['sspa'],
    package_dir={'':'src'},
    url='https://github.com/cwieder/sspa',
    license='GNU 3.0',
    author='Cecilia Wieder',
    author_email='cw2019@ic.ac.uk',
    description='Single sample pathway analysis tools for omics data',
    long_description=long_description,
    long_description_content_type='text/markdown',
    package_data={'sspa': ['example_data/*', 'pathway_databases/*']},
    classifiers=[
        'Programming Language :: Python',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.10',
    ]
)
