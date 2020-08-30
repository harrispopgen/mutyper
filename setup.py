import setuptools
import versioneer


with open('README.md', 'r') as fh:
    long_description = fh.read()

setuptools.setup(
    name='mutyper',
    version=versioneer.get_version(),
    cmdclass=versioneer.get_cmdclass(),
    author='William DeWitt',
    author_email='wsdewitt@gmail.com',
    description='ancestral k-mer mutation types for SNP data',
    long_description=long_description,
    long_description_content_type='text/markdown',
    url='https://github.com/harrispopgen/mutyper',
    entry_points={'console_scripts': ['mutyper=mutyper.cli:main']},
    # packages=setuptools.find_packages(exclude=['tests', 'docs', 'docsrc']),
    packages=['mutyper'],
    classifiers=[
        'Programming Language :: Python :: 3',
        'License :: OSI Approved :: MIT License',
        'Operating System :: OS Independent',
    ],
    python_requires='>=3.7',
    install_requires=[
        'cyvcf2',
        'pyfaidx',
        'biopython',
        'pandas',
        'pyliftover'
    ]
)
