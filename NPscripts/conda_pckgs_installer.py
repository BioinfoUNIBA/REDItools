#!/usr/bin/python 
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
# THE SOFTWARE.

import os, subprocess

def install_conda_packages(conda_bin):
	"""Installs conda packages required by the protocol"""
	install_cmd = os.system(cmd + ' install -n nature_protocol bcftools==1.9 bedtools==2.28.0 \
        bzip2==1.0.6 bwa==0.7.17 bx-python==0.8.2 fastp==0.20.0 fastqc==0.11.8 \
        fisher==0.1.4 git==2.21.0 gmap==2018.07.04 htslib==1.9 libdeflate==1.0 \
        numpy==1.16.2 pysam==0.15.2 rseqc==2.6.4 samtools==1.9 scipy==1.2.1 \
        star==2.7.0f wget==1.20.1')
	return install_cmd

if subprocess.getstatusoutput('conda')[0] != 0:
	cwd = os.getcwd()
	installation_path = cwd + '/opt'
	if not os.path.exists(installation_path):
		os.mkdir(installation_path)
	os.chdir('./opt')
	conda_url = 'wget https://repo.anaconda.com/miniconda/Miniconda2-latest-Linux-x86_64.sh'
	i = 0
	while os.system(conda_url) != 0 and i <= 5:
		os.system(conda_url)
		i+=1
	os.system('chmod +x Miniconda2-latest-Linux-x86_64.sh')
	os.system('bash Miniconda2-latest-Linux-x86_64.sh')
	home_folder = os.path.expanduser('~')
	cmd = home_folder + '/miniconda2/bin/conda'
	os.system(cmd + ' config --add channels defaults')
	os.system(cmd + ' config --add channels bioconda')
	os.system(cmd + ' config --add channels conda-forge')
	os.system(cmd + ' create -n nature_protocol python=2.7 anaconda')
	install_conda_packages(cmd)
	print("Your conda environment has been succesfully created, now close your terminal and open a new one." + "\n" + \
	      "Type in order:" + "\n" + \
	      "source " + home_folder + "/.bashrc" + "\n" + \
	      "conda activate nature_protocol")
else:
        home_folder = os.path.expanduser('~')
        cmd = home_folder + '/miniconda2/bin/conda'
        os.system(cmd + ' config --add channels defaults')
        os.system(cmd + ' config --add channels bioconda')
        os.system(cmd + ' config --add channels conda-forge')
        os.system(cmd + ' create -n nature_protocol python=2.7 anaconda')
        install_conda_packages(cmd)
        print("Your conda environment has been succesfully created, now close your terminal and open a new one." + "\n" + 
		"Type in order:" + "\n" + \
		"source " + home_folder + "/.bashrc" + "\n" + \
		"conda activate nature_protocol")

