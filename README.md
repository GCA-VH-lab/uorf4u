## Description

uorf4u is a simple bioinformatics tool for retrieving protein's homologous, corresponding genes upstream sequences and annotation conserved short upstream ORFs.  

**Programming languages:** python3 (main), R;  
**OS:** MacOS, Linux;  
**Python dependencies:** biopython, configs, argparse, statistics;  
**R dependencies:** ggmsa, ggplot2, optparse;  
**OS-level dependencies:** muscle;  
**License:** [WTFPL](http://www.wtfpl.net);  
**Version:** 0.1 (July 2022)

[**Detailed documentation**](https://art-egorov.github.io/uorf4u)

### Data analysis pipeline:

<img  src="docs/img/pipeline.svg" width="550"/>


## Installation


- The most stable release of uorf4u can be installed directly from pypi:

```
python3 -m pip install uorf4u
```

- The development version is available at github :

```
git clone https://github.com/art-egorov/uorf4u.git
cd uorf4u
python3 -m pip install --upgrade pip
python3 -m pip install wheel
python3 setup.py sdist bdist_wheel
python3 -m pip install -e .
```


## Reference

If you find uorf4u useful, please cite:


Artyom A. Egorov, Gemma C. Atkinson **uorf4u: ...,** *---, [doi]()*


## Contact

Please contact us by e-mail _artem**dot**egorov**AT**med**dot**lu**dot**se_ or use Issues to report any technical problems.  


## Authors

uorf4u is developed by Artyom Egorov at the [Atkinson Lab](https://atkinson-lab.com), Department of Experimental Medical Science, Lund University, Sweden. We are open for suggestions to extend and improve svist4get functionality. Please don't hesitate to share your ideas or feature requests.

