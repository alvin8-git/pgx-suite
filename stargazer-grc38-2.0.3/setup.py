from setuptools import setup, find_packages

exec(open("stargazer/version.py").read())

setup(
    name="stargazer",
    version=__version__,
    author='Seung-been "Steven" Lee',
    author_email="sbstevenlee@gmail.com",
    description="A tool for calling star alleles from pharmacogenomic data",
    long_description=open("README.md").read(),
    packages=find_packages(),
    package_data={
        "stargazer": [
            "sv_table.tsv",
            "star_table.tsv",
            "snp_table.tsv",
            "gene_table.tsv",
            "beagle.05May22.33a.jar",
            "1kgp_vcf/*.vcf.gz",
        ],
    },
    entry_points={"console_scripts": ["stargazer=stargazer.__main__:main"]}
)
