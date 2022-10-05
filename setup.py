from setuptools import setup

with open("pypi.md", "r") as fh:
    long_description = fh.read()

setup(name="uorf4u",
      version="0.5.0",
      description="A tool for short uORF annotation.",
      url="https://art-egorov.github.io/uorf4u/",
      author="Artyom Egorov",
      author_email="artem.egorov@med.lu.se",
      license="WTFPL",
      packages=["uorf4u"],
      install_requires=["biopython", "configs", "argparse", "statistics", "logomaker", "matplotlib", "pandas",
                        "reportlab"],
      long_description=long_description,
      long_description_content_type="text/markdown",
      include_package_data=True,
      scripts=["bin/uorf4u"],
      zip_safe=True,
      package_data={"uorf4u": ["uorf4u_data/*", "uorf4u_data/fonts/*", "uorf4u_data/bin/*"]})