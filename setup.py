import setuptools
import os

with open("uorf4u/docs/pypi.md", "r") as fh:
    long_description = fh.read()


def package_files(directory):
    paths = []
    for (path, directories, filenames) in os.walk(directory):
        for filename in filenames:
            paths.append(os.path.join("..", path, filename))
    return paths


extra_files = package_files("uorf4u/uorf4u_data")
extra_files.append("docs/pypi.md")

setuptools.setup(name="uorf4u",
      version="0.9.6.1",
      description="A tool for short uORF annotation.",
      python_requires='>=3.7',
      url="https://art-egorov.github.io/uorf4u/",
      author="Artyom Egorov",
      author_email="artem.egorov@med.lu.se",
      license="WTFPL",
      packages=["uorf4u"],
      package_data={"uorf4u": extra_files},
      install_requires=["biopython", "configs", "statistics", "logomaker==0.8", "matplotlib", "pandas==1.4.0",
                        "numpy==1.26.4", "reportlab", "msa4u"],
      long_description=long_description,
      long_description_content_type="text/markdown",
      scripts=["bin/uorf4u"],
      zip_safe=False)
