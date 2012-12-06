#! /usr/bin/env python
#coding=utf-8
import subprocess
from zipfile import ZipFile, ZIP_DEFLATED


class FIMO(object):
    """docstring for FIMO"""

    # params = {
    #     '--output-pthresh': None,
    #     '--output-qthresh': None,
    #     '--text': False,
    #     '--max-stored-scores': '100000',
    #     '--verbosity': '1',
    # }
    meme_motif = ''
    result = ''

    def __init__(self, motif_file, fasta_file, ouput_path):
        super(FIMO, self).__init__()
        self.motif_file = motif_file
        self.fasta_file = fasta_file
        # self.params['-oc'] = ouput_path
        self.params = {
            '--output-pthresh': None,
            '--output-qthresh': None,
            '--text': False,
            '--max-stored-scores': '100000',
            '--verbosity': '1',
            '-oc': ouput_path}

    def run(self):
        args = ['/home/gahoo/meme/bin/fimo']
        for arg in self.params:
            value = self.params[arg]
            if value == True:
                args.append(arg)
            elif value:
                args.extend([arg, self.params[arg]])
        args += [self.motif_file, self.fasta_file]
        #print args

        if self.params['--text'] and self.motif_file == '-':
            #args.insert(1, '--no-qvalue')
            fimo = subprocess.Popen(args, stdin=subprocess.PIPE, stdout=subprocess.PIPE)
            (self.result, err) = fimo.communicate(self.meme_motif)
            if err:
                print err
        elif not self.params['--text'] and self.motif_file == '-':
            fimo = subprocess.Popen(args, stdin=subprocess.PIPE)
            fimo.communicate(self.meme_motif)
        elif self.params['--text'] and not self.motif_file == '-':
            #args.insert(1, '--no-qvalue')
            fimo = subprocess.Popen(args, stdout=subprocess.PIPE)
            (self.result, err) = fimo.communicate(self.meme_motif)
            if err:
                print err
        else:
            fimo = subprocess.Popen(args)
        return self.parseResult()

    def zip(self):
        self.params['--text'] = True
        self.run()
        for line in self.meme_motif.split('\n'):
            line = line.split(' ')
            if line[0] == "MOTIF":
                Motif = line[1]
                break
        dumpzfile = "%s/Scan/%s.zip" % (self.params['-oc'], Motif)
        zf = ZipFile(dumpzfile, "w", ZIP_DEFLATED)
        zf.writestr(Motif + ".scan", self.result)
        zf.close()

    def parseResult(self):
        if not self.result:
            self.result = open(self.params['-oc'] + '/fimo.txt').read()
        return [rs.split('\t') for rs in self.result.split('\n')[1:-1]]

if __name__ == '__main__':
    #fimo = FIMO('kk.meme', 'TAp73alpha.fa', 'out_dir')
    fimo = FIMO('-', 'TAp73alpha.fa', 'out_dir')
    #fimo.params['--text'] = True
    fimo.meme_motif = """MEME version 4.4

ALPHABET= ACGT

strands: + -

Background letter frequencies (from uniform background):
A 0.25000 C 0.25000 G 0.25000 T 0.25000 

MOTIF TTGCTAG 

letter-probability matrix: alength= 4 w= 7 nsites= 20 E= 0
  0.000000    0.000000    0.000000    1.000000  
  0.000000    0.000000    0.000000    1.000000  
  0.000000    0.000000    1.000000    0.000000  
  0.000000    1.000000    0.000000    0.000000  
  0.000000    0.000000    0.000000    1.000000  
  1.000000    0.000000    0.000000    0.000000  
  0.000000    0.000000    1.000000    0.000000
    """
    fimo.run()
    print fimo.parseResult()
