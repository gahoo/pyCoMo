#! /usr/bin/env python
#coding=utf-8
#TODO parallel execute fimo
#TODO motif by motif, store current status
from Queue import Queue
import shutil
import pp
from misc.meme_motif_parsing import join_meme
from misc.jaspar import JasparDB
from misc.plantcare import PlantCARE
from misc.place import Place
from misc.FIMO import FIMO
from hashlib import md5
from misc.united_db import UnitiedDB
from optparse import OptionParser
import os
import shelve
import zshelve
import sys
import time


def loadMEME_ID(subset):
    care = PlantCARE('./databases/plantcare.db')
    place = Place('./databases/place.db')
    jaspar = JasparDB('./databases/jaspar.db')
    meme_ids = []
    print subset
    if 'plantcare' in subset.keys() and subset['plantcare'] != None:
        meme_ids.extend(care.getMotifSeq(subset['plantcare']))
    if 'place' in subset.keys() and subset['place'] != None:
        tmp_seq = place.getMotifSeq(subset['place'])
        #for those contain U would change by iupac2meme to T
        #causing KeyError
        tmp_seq = [seq.replace("U", "T") for seq in tmp_seq]
        meme_ids.extend(tmp_seq)
    if 'jaspar' in subset.keys() and subset['jaspar'] != None:
        meme_ids.extend(jaspar.getBASE_ID(subset['jaspar']))
    if 'denovo' in subset.keys() and subset['denovo'] != None:
        meme_ids.extend(united.getMEME_ID(subset['denovo']))
    return list(set(meme_ids))


def newFIMO(options, meme_slv, meme_ids):
    folder = md5(str(meme_ids)).hexdigest()[:6]
    fimo = FIMO('-', options.fasta_file, options.output_path + "/%s/" % folder)
    fimo.meme_motif = join_meme(meme_slv, meme_ids)
    #A nature article "Unsupervised pattern discovery in human chromatin structure through genomic segmentation" use q<0.1
    #An other "The accessible chromatin landscape of the human genome" use p<10e-5
    #fimo.params['--output-qthresh'] = '0.6'
    fimo.params['--output-pthresh'] = '10e-5'
    #the more stored the more time comsumed
    #fimo.params['--max-stored-scores'] = '1000000'
    if 'bg_file' in dir(options):
        fimo.params['--bgfile'] = options.bg_file
    if 'pvalue' in dir(options):
        fimo.params['--output-pthresh'] = options.pvalue
    if 'qvalue' in dir(options):
        fimo.params['--output-qthresh'] = options.qvalue
    if 'params' in dir(options) and options.params:
        print options.params
        fimo.params = options.params
    return fimo


def scanMotifs(options, united, meme_slv, meme_ids):
    log = zshelve.open(options.dbname + ".log")
    total_num = len(meme_ids) / 100.0
    for i, meme_id in enumerate(meme_ids):
        if meme_id in log.keys():
            sys.stderr.write("Skipping:%d%%\t%d\t%s\r" % (i / total_num, len(log.keys()), meme_id))
            continue
        sys.stderr.write("Scanning:%d%%\t%d\t%s\r" % (i / total_num, len(log.keys()), meme_id))
        fimo = newFIMO(options, meme_slv, [meme_id])
        result = fimo.run()
        united.importScan(result)
        log[str(meme_id)] = result
    return


def scanMotifsParallel(options, united, meme_slv, meme_ids):
    log = zshelve.open(options.dbname + ".log")
    total_num = len(meme_ids) / 100.0
    qResults = Queue()
    jobs = []
    for i, meme_id in enumerate(meme_ids):
        fimo = newFIMO(options, meme_slv, [meme_id])
        jobs.append(options.parallel.submit(fimo.run, (), (), ('subprocess',)))

    for i, (meme_id, job) in enumerate(zip(meme_ids, jobs)):
        if meme_id in log.keys():
            sys.stderr.write("Skipping:%d%%\t%s\r" % (i / total_num, meme_id))
            continue
        sys.stderr.write("Scanning:%d%%\t%s\r" % (i / total_num, meme_id))
        result = job()
        qResults.put((meme_id, result))
        while qResults.qsize() > 0:
            (meme_id, result) = qResults.get()
            united.importScan(result)
            log[str(meme_id)] = result
            time.sleep(0.5)
            folder = md5(str([meme_id])).hexdigest()[:6]
            shutil.rmtree(options.output_path + "/%s/" % folder)
    shutil.rmtree(options.output_path)
    return


def getOptions():
    parser = OptionParser(usage="Usage: scanMotifs [options] -d dbname -f <fastafile>")
    parser.add_option("-d", "--dbname", dest="dbname",
                help="Database Name", metavar="FILE")
    parser.add_option("-f", "--fasta", dest="fasta_file",
                help="Fasta File", metavar="FILE")
    parser.add_option("-b", "--background", dest="bg_file",
                help="Background File can be created from any FASTA sequence file using the fasta-get-marko",
                metavar="FILE")
    parser.add_option("-p", "--pvalue", dest="pvalue",
                help="Pvalue Threshold to output",
                metavar="FLOAT")
    parser.add_option("-q", "--qvalue", dest="qvalue",
                help="Qvalue Threshold to output",
                metavar="FLOAT")
    parser.add_option("-s", "--subset", dest="subset",
                help="Subset Database to scan \t\t\t"
                """{'plantcare': {'Organism': [u'Antirrhinum majus']},
                    'place': {'Keyword': ['phyA', 'chloroplast']},
                    'jaspar': {"COLLECTION": ["CORE"],
                            "tax_group": ["vertebrates", "insect"],
                            "family": ["GATA"]},\t
                    'denovo': {'Database': ['Denovo'],
                                'Motif': ['AATTC']}}""",
                default="""{'plantcare': {},
                        'place': {},
                        'jaspar': {},
                        'denovo': {'Database': ['Denovo']}}""",
                metavar="DICT")
    parser.add_option("-l", "--params", dest="params",
                help="Set parameters of FIMO \t\t\t\t\t"
                """{'--output-pthresh': None,\t\t\t
                '--output-qthresh': None,\t\t\t
                '--text': False,\t\t\t
                '--max-stored-scores': '100000',\t\t\t
                '--verbosity': '4'}""",
                metavar="DICT")
    parser.add_option("-o", "--output", dest="output_path",
                help="FIMO result output path", default="", metavar="DIR")
    parser.add_option("-n", "--parallel", dest="parallel",
                default=False,
                action="store_true",
                help="Run in parallel mode")
    return parser

if __name__ == '__main__':
    parser = getOptions()
    (options, args) = parser.parse_args()
    if not options.dbname or not options.fasta_file:
        parser.print_help()
        sys.exit(1)
    if not options.output_path:
        options.output_path = options.dbname
    if not os.path.exists(options.output_path):
        os.mkdir(options.output_path)

    if options.parallel:
        options.parallel = pp.Server()

    united = UnitiedDB(options.dbname + '.db')
    meme_slv = shelve.open(options.dbname + ".slv")
    meme_ids = loadMEME_ID(eval(options.subset, {}, {}))
    meme_ids = list(set(united.getMEME_ID({})) & set(meme_ids))
    print len(meme_ids), "to be scanned"
    timestamp = time.time()
    if options.parallel:
        scanMotifsParallel(options, united, meme_slv, meme_ids)
    else:
        scanMotifs(options, united, meme_slv, meme_ids)
    print time.time() - timestamp
