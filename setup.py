from setuptools import setup

with open("README.md", "r") as fh:
    long_description = fh.read()

setup(name='uorf4u',
      version='0.3.0',
      description='A tool for short uORF annotation.',
      url='https://art-egorov.github.io/uorf4u/',
      author='Artyom Egorov',
      author_email='artyom.egorov@med.lu.se',
      license='WTFPL',
      packages=['uorf4u'],
      install_requires=['biopython', 'configs', 'argparse', 'statistics', 'logomaker', 'matplotlib', 'pandas'],
      long_description=long_description,
      long_description_content_type="text/markdown",
      scripts=['bin/uorf4u'],
      zip_safe=False,
      package_data={'uorf4u': ['uorf4u_data/*']},
      include_package_data=True)
