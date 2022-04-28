# VCFannotator
A bioinformatics tool to annotate VCF files

## Requirements
* Linux or MacOS
* python >= 3.5

## Installation

* [download](https://github.com/open-projects/VCFannotator/zipball/master) latest stable VCFannotator build from this page
* unzip the archive
* add resulting folder to your ``PATH`` variable
  * or add symbolic link for ``VCFannotator`` script to your ``bin`` folder
  * or use VCFannotator directly by specifying full path to the executable script

## Examples of usage

VCFannotator.py -i test_data/input.vcf.gz -a test_data/dbSNP.vcf.gz -o annotation.csv

## License
Copyright (c) 2022, D. Malko
All Rights Reserved

VCFannotator is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see [http://www.gnu.org/licenses/](http://www.gnu.org/licenses/).
