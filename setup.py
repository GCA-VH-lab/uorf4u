from setuptools import setup

with open("README.md", "r") as fh:
    long_description = fh.read()

setup(name='uorf4u',
      version='0.1',
      description='A tool for short uORF annotation.',
      url='-',
      author='Artyom Egorov',
      author_email='artyom.egorov@hotmail.com',
      license='WTFPL',
      packages=['uorf4u'],
      install_requires=['biopython', 'configs', 'argparse', 'statistics'],
      long_description=long_description,
      long_description_content_type="text/markdown",
      scripts=['bin/uorf4u'],
      zip_safe=False,
      package_data={'uorf4u': ['uorf4u_data/*']},
      include_package_data=True)
