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


#GFF structure
#chr/tvalue/tfeature/tstart/tend/t./tstrand/./t gene_id ed_numb; transcript_id ed_numb;

import sys

try:
	in_table = sys.argv[1]
except:
	sys.exit('<RDPRTL_table>') 

i=0
with open('table1_hg19.tsv','r') as e:
	e.readline()
	for line in e:
		line = line.split('\t')
		if line[7] == 'NONREP' and line[10] == 'exonic':
			if ('\t'.join(line).count('nonsynonymous')) == 3:
				i+=1
				valore = line[13].split(':')[0] + '_' +  line[13].split('.')[-1]
				gff_row = line[1] + '\t'+ valore + '\t' + 'ed' + '\t' + line[2] + \
				'\t' + line[2] + '\t' + '.' + '\t' + line[5] + '\t' + '.' + '\t' + \
				'gene_id' + ' '  + '"ed_%s";' %(i) + ' ' + 'transcript_id' + ' ' + '"ed_%s";' %(i)
				print gff_row
