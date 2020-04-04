import setuptools

with open('README.md', 'r') as fh:
    long_description = fh.read()

setuptools.setup(
    name='mutyper',
    version='0.1',
    author='William DeWitt',
    author_email='wsdewitt@gmail.com',
    description='Ancestral k-mer mutation types for SNP data',
    long_description=long_description,
    long_description_content_type='text/markdown',
    url='https://github.com/harrispopgen/mutyper',
    packages=setuptools.find_packages(),
    classifiers=[
        'Programming Language :: Python :: 3',
        'License :: OSI Approved :: MIT License',
        'Operating System :: OS Independent',
    ],
    python_requires='>=3.6',
    scripts=['bin/mutyper'],
    install_requires=[
        'cyvcf2',
        'pyfaidx',
        'biopython',
        'pandas'
    ]
)
