from setuptools import setup

with open("uorf4u/docs/pypi.md", "r") as fh:
    long_description = fh.read()

setup(name="uorf4u",
      version="0.5.1",
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
      package_data={"uorf4u": ["docs/pypi.md", "uorf4u_data/*", "uorf4u_data/fonts/*", "uorf4u_data/bin/*"]},
      include_package_data=True,
      scripts=["bin/uorf4u"],
      zip_safe=True)
