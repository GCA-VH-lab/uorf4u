# Сommand-line parameters

	
**POST-INSTALL EXAMPLE DATA**

- `--sampledata`  
Creates the *uorf4u_data* folder in the current working directory.
The folder will contain an adjustable configuration file template, SD-aSD energies table and helper scripts.


**MANDATORY ARGUMENTS**

- `-an` *accession_number*  
Protein's RefSeq accession number.

	OR

- `-hl` *accession_number1 [accession_number2, ...]*  
Space separated list of proteins accession numbers which will be used as list of homologous.

	OR

- `-hlf` *file.txt*
Path to a file with list of accession numbers. File format: one accession number per line, no header.


**OPTIONAL ARGUMENTS**

- `-sao`  
Save information about annotated ORFs as a set of tsv tables (one per each upstream sequence) [default: false].


- `-o` *dirname*  
Output dirname. It will be created if it's not exist. All output dirs will be then created in this folder [default: uorf4u_{current_date}; e.g. uorf4u_2022_07_25-20_41].


- `-c` *file.cfg*  
Path to a configuration file [default: internal].


**MISCELLANEOUS ARGUMENTS**

- `-h`, `--help`  
Show this help message and exit.

- `-v`, `--version`  
Show program version.

- `--debug`  #to be added  
Provide detailed stack trace for debugging purposes.

- `--verbose`  # to be added  
Show all progress messages [default: False]
