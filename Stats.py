#! /usr/bin/env python
#coding=utf-8
from misc.united_db import UnitiedDB
from misc.equations import hypergeo_cdf_enrich, enrichmen_ratio, holm_adjustment, hypergeo_cdf, ppCacl, ssCacl
from optparse import OptionParser
from collections import defaultdict
from hashlib import md5
from meliae import scanner
import gzip
import zshelve
import sys
import os
import random
import time
import pp


def getDist(posA, posB, low, high):
    '''
    Get the distances of a pair of motifs on a single seqname
    posA=posB=[(start,stop),...]
    comb=[((startA,stopA),(startB,stopB)),...]
    CoDist={((startA,stopA),(startB,stopB)): d,...}
    '''
    CoDist = {}
    for a in posA:
        a_ = sorted(a)
        for b in posB:
            b_ = sorted(b)
            if not is_overlap(a_, b_):
                d = caclDist(a_, b_)
                abs_d = abs(d)
                if abs_d >= low and abs_d <= high:
                    CoDist[(a, b)] = d
    return CoDist


def is_overlap(posA, posB):
    '''
    '''
    #if A.Stop < B.Start or A.Start >B.Stop then it's not overlap
    (min_a, max_a) = posA
    (min_b, max_b) = posB
    if max_a < min_b or min_a > max_b:
        return False
    else:
        return True


def caclDist(posA, posB):
    (min_a, max_a) = posA
    (min_b, max_b) = posB
    #Start - Stop
    d1 = min_b - max_a
    d2 = min_a - max_b
    return {True: d1, False: -d2}[d1 > d2]


def loadList(list_file):
    return [l.strip() for l in open(list_file).readlines()]


def getBasename(filename):
    return ".".join(filename.split('.')[:-1])


def countSeqNames(Motifs):
    num_Motifs = len(Motifs) / 100 + 1
    Cnts = {}
    for i, uid in enumerate(Motifs.iterkeys()):
        sys.stderr.write("Progress:%d%%\t%s\r" % (i / num_Motifs, uid))
        Cnts[uid] = len(Motifs[uid])
    print "Progress:100%"
    return Cnts
    #return dict([(uid, len(Motifs[uid])) for uid in Motifs.keys()])


def countPairSeqNames(Motifs, MotifPairs):
    num_Pairs = len(MotifPairs) / 100 + 1
    Cnts = {}
    for i, pair in enumerate(MotifPairs):
        pair1 = "%s,%s" % pair
        pair2 = "%s,%s" % pair[::-1]
        sys.stderr.write("Progress:%d%%\t%s\r" % (i / num_Pairs, str(pair)))
        try:
            Cnts[pair1] = len(Motifs[pair1])
        except KeyError:
            try:
                Cnts[pair2] = len(Motifs[pair2])
            except KeyError:
                pass
    print "Progress:100%"
    return Cnts
    #return dict([(uid, len(Motifs[uid])) for uid in Motifs.keys()])


def getMotif_SeqName(Motifs):
    return dict([(uid, set(Motifs[uid].keys())) for uid in Motifs.iterkeys()])


def stats(options):
    united = UnitiedDB(options.dbname)
    REF_SeqNames = [united.SeqName2REF_SeqName[seqname] for seqname in loadList(options.list_file)]
    #For debug random choose 200
    #REF_SeqNames = random.sample(united.SeqName2REF_SeqName.values(), 1000)
    #REF_SeqNames = united.SeqName2REF_SeqName.values()[:2000]
    list_size = len(REF_SeqNames)
    db_size = len(united.SeqName2REF_SeqName)
    if options.motif_list_file:
        REF_uids = [united.Motif2uid[motif] for motif in loadList(options.motif_list_file)]
    else:
        REF_uids = []
    print "Load Database..."
    Motifs_List = loadMotifs(united, options, getBasename(options.list_file), REF_SeqNames, REF_uids)
    Motifs_DB = loadMotifs(united, options, "DB", [], REF_uids)
    file_suffix = getBasename(os.path.basename(options.dbname)) + getFileSuffix(options, getBasename(options.list_file), REF_SeqNames + REF_uids)
    Nodes = enrichment(united, options, Motifs_List, Motifs_DB, list_size, db_size)
    Edges = cooccur(united, options, Motifs_List, Motifs_DB, list_size, db_size)
    output = newZslv(options, "Result_%s" % getBasename(options.list_file), REF_SeqNames + REF_uids)
    output["Nodes"] = Nodes
    saveData(Nodes, "Nodes", "%s/%s.Nodes" % (options.output_path, file_suffix))
    output["Edges"] = Edges
    saveData(Edges, "Edges", "%s/%s.Edges" % (options.output_path, file_suffix))


def loadMotifs(united, options, zslv_type, REF_SeqNames, REF_uids):
    print "Loading Motifs in %s" % zslv_type
    (dist_min, dist_max) = map(int, options.dist_range)
    Motifs = newZslv(options, "MotifScan_%s" % zslv_type, REF_SeqNames + REF_uids)
    if not Motifs:
        Motifs.update(united.dumpScan2Dict(options.upstream, options.strand, options.pValue, options.qValue,
                                        dist_min, dist_max, REF_SeqNames, REF_uids))
    return Motifs


def enrichment(united, options, Motifs_List, Motifs_DB, list_size, db_size):
    print "Counting Motifs"
    SeqNames_Cnts_List = countSeqNames(Motifs_List)
    SeqNames_Cnts_DB = countSeqNames(Motifs_DB)

    print "%s Motifs to Cacl Enrichment" % len(SeqNames_Cnts_List.keys())
    inputs = [(SeqNames_Cnts_List[uid], list_size, SeqNames_Cnts_DB[uid], db_size) \
                for uid in SeqNames_Cnts_List.iterkeys()]
    if options.parallel:
        ERs = ppCacl(options.parallel, inputs, enrichmen_ratio)
        Ps = ppCacl(options.parallel, inputs, hypergeo_cdf_enrich)
    else:
        ERs = ssCacl(inputs, enrichmen_ratio)
        Ps = ssCacl(inputs, hypergeo_cdf_enrich)
    adjPs = holm_adjustment(Ps)
    Nodes = formatNodes(united, SeqNames_Cnts_List, list_size,
                        SeqNames_Cnts_DB, db_size, ERs, Ps, adjPs)
    return Nodes


def formatNodes(united, SeqNames_Cnts_List, list_size, SeqNames_Cnts_DB, db_size, ERs, Ps, adjPs):
    Nodes = []
    for i, uid in enumerate(SeqNames_Cnts_List.iterkeys()):
        #when we need all nodes information, the thresh should be drop
        #if (adjPs[i] < 0.05 and SeqNames_Cnts_List[uid] > 5):
            #add MEME_ID or the Motifseq
            Nodes.append((united.uid2Motif[int(uid)], united.uid2Desc[int(uid)], united.uid2MEME_ID[int(uid)],
                         SeqNames_Cnts_List[uid], list_size,
                         SeqNames_Cnts_DB[uid], db_size,
                         ERs[i], Ps[i], adjPs[i]))
    return Nodes


def cooccur(united, options, Motifs_List, Motifs_DB, list_size, db_size):
    PairedMotifs_Cnts_List = countPairs(united, options, Motifs_List, getBasename(options.list_file))
    PairedMotifs_Cnts_DB = countPairs(united, options, Motifs_DB, "DB")

    (inputs, Ps, Edges) = ([], [], [])
    num_MotifPairs = len(PairedMotifs_Cnts_List)
    print "%d Pairs to Cacl Co-occuring" % num_MotifPairs
    num_MotifPairs = num_MotifPairs / 100 + 1
    for pair in PairedMotifs_Cnts_List.iterkeys():
        (uidA, uidB) = map(int, pair.split(','))
        try:
            inputs.append((PairedMotifs_Cnts_List[pair], list_size, PairedMotifs_Cnts_DB[pair], db_size))
        except KeyError:
            pair2 = "%s,%s" % (uidB, uidA)
            inputs.append((PairedMotifs_Cnts_List[pair], list_size, PairedMotifs_Cnts_DB[pair2], db_size))
        #sys.stderr.write("Progress:%d%%\t%s\t%s\r" % (i / num_MotifPairs, uidA, uidB))
    if options.parallel:
        print "Cacling Enrichment Ratio"
        ERs = ppCacl(options.parallel, inputs, enrichmen_ratio)
        print "Cacling pValues"
        Ps = ppCacl(options.parallel, inputs, hypergeo_cdf)
    else:
        print "Cacling Enrichment Ratio"
        ERs = ssCacl(inputs, enrichmen_ratio)
        print "Cacling pValues"
        Ps = ssCacl(inputs, hypergeo_cdf)
    print "Adjusting pValues"
    adjPs = holm_adjustment(Ps)
    Edges = formatEdges(united, PairedMotifs_Cnts_List, inputs, ERs, Ps, adjPs)
    return Edges


def formatEdges(united, PairedMotifs_Cnts_List, inputs, ERs, Ps, adjPs):
    Edges = []
    for i, pair in enumerate(PairedMotifs_Cnts_List.iterkeys()):
        (uidA, uidB) = map(int, pair.split(','))
        if adjPs[i] < 0.05 and inputs[i][0] > 5:
            Edges.append([united.uid2Motif[uidA], united.uid2Motif[uidB]] + list(inputs[i]) + [ERs[i], Ps[i], adjPs[i]])
    return Edges


def saveData(data, datatype, filename=""):
    print "%s has %d records." % (datatype, len(data))
    if not data:
        return []
    if datatype == "Nodes":
        output = "\n".join(["%s\t%s\t%s\t%d\t%d\t%d\t%d\t%.2f\t%.2E\t%.2E" % tuple(node) for node in data])
    elif datatype == "Edges":
        output = "\n".join(["%s\t%s\t%d\t%d\t%d\t%d\t%.2f\t%.2E\t%.2E" % tuple(edge) for edge in data])
    if filename:
        outfile = open(filename, 'w')
        outfile.write(output)
        outfile.close()
    else:
        print output
    return output


def countPairs(united, options, Motifs_zslv, zslv_type):
    print "Loading Motif Pair Counts in %s" % zslv_type
    PairedMotifs_Cnts = newZslv(options, "MotifCnts_%s" % zslv_type, Motifs_zslv.keys())
    if not PairedMotifs_Cnts:
        print "Getting Motif Pairs in %s\t" % zslv_type
        PairedMotifs = getPairedSeqNames(united, options, Motifs_zslv, "MotifPaired_%s" % zslv_type)
        print "Counting Motif Pairs in %s\t" % zslv_type
        MotifPairs = getPairedList(Motifs_zslv.keys())
        PairedMotifs_Cnts.update(countPairSeqNames(PairedMotifs, MotifPairs))
    #load to RAM will not improve many here
    Cnts = {}
    Cnts.update(PairedMotifs_Cnts)
    return Cnts


def newZslv(options, zslv_type, tohash=None):
    zslv_file = getBasename(os.path.basename(options.dbname)) + getFileSuffix(options, zslv_type, tohash) + ".zslv"
    #print "Saving to file:%s" % zslv_file
    #return zshelve.open("%s/%s" % (options.output_path, zslv_file))
    if "DB" in zslv_type:
        zslv_path = os.path.abspath(os.path.dirname(__file__))
        return zshelve.open("%s/cache/%s" % (zslv_path, zslv_file))
    else:
        return zshelve.open("%s/%s" % (options.output_path, zslv_file))


def getFileSuffix(options, ftype, tohash=None):
    if not tohash or 'DB' in ftype:
        hashREF = ""
    else:
        hashREF = md5(str(sorted(tohash))).hexdigest()[:6]
    (dist_min, dist_max) = map(int, options.dist_range)
    file_suffix = "_%s" % "_".join(map(str,
        [ftype, hashREF, options.upstream, options.strand, options.pValue, options.qValue, dist_min, dist_max]))
    return file_suffix


def countPairs(united, options, Motifs_zslv, zslv_type):
    if options.low_RAM:
        Motifs = Motifs_zslv
    else:
        #load Motifs into RAM to boost speed.
        Motifs = {}
        Motifs.update(Motifs_zslv)
    if options.big_RAM:
        #load everything in RAM
        PairedMotifs = {}
        PairedMotifs.update(newZslv(options, "PairedMotifs_%s" % zslv_type, Motifs.keys()))
        PairedMotifs_Cnts = {}
        PairedMotifs_Cnts.update(newZslv(options, "PairedMotifsCnts_%s" % zslv_type, Motifs.keys()))
    else:
        PairedMotifs = newZslv(options, "PairedMotifs_%s" % zslv_type, Motifs.keys())
        PairedMotifs_Cnts = newZslv(options, "PairedMotifsCnts_%s" % zslv_type, Motifs.keys())

    (dist_min, dist_max) = map(int, options.dist_range)
    SeqNames = getMotif_SeqName(Motifs)
    # done = set([eval("('%s','%s')" % tuple(k.split(','))) for k in PairedMotifs.iterkeys()])
    #(done, log) = newLog(options, zslv_type, Motifs.keys())
    #MotifPairs = set(getPairedList(Motifs.keys())) - done
    MotifPairs = getPairedList(Motifs.keys())
    num_MotifPairs = len(MotifPairs)
    print "%s Pairs of Motif to Count in %s\t" % (num_MotifPairs, zslv_type)
    num_MotifPairs = num_MotifPairs / 100 + 1

    for i, (uidA, uidB) in enumerate(MotifPairs):
        key = "%s,%s" % (uidA, uidB)
        # if PairedMotifs_Cnts.has_key(key):
        if key in PairedMotifs_Cnts:
            sys.stderr.write("Skipping:%d%%\t%s\t%s\r" % (i / num_MotifPairs, uidA, uidB))
            continue
        intersect = list(SeqNames[uidA] & SeqNames[uidB])
        if intersect:
            PosA = Motifs[uidA]
            PosB = Motifs[uidB]
            sys.stderr.write("Progress:%d%%\t%s\t%s\r" % (i / num_MotifPairs, uidA, uidB))

            distances = {}
            for seqname in intersect:
                distance = getDist(PosA[seqname], PosB[seqname], dist_min, dist_max)
                if distance:
                    distances[seqname] = distance
            if distances:
                PairedMotifs[key] = distances
                PairedMotifs_Cnts[key] = len(distances)

    print "Progress:100%"
    return PairedMotifs_Cnts


def newLog(options, zslv_type, tohash=None):
    log_file = getBasename(os.path.basename(options.dbname)) + getFileSuffix(options, zslv_type, tohash) + ".log"
    if not os.path.exists(log_file):
        log = gzip.open(log_file, 'w')
        done = []
    else:
        log = gzip.open(log_file, 'rw')
        #print log.readlines()
        done = set([eval("('%s','%s')" % tuple(k.strip().split(','))) for k in log.readlines()])
    return (set(done), log)


def getPairedList(uids):
    #uids = map(int, uids)
    return [(id1, id2) for i, id1 in enumerate(uids) for id2 in uids[i:] if id1 != id2]


def getOptions():
    parser = OptionParser(usage="Usage: Stats [options] -d dbname -l <list_file>")
    parser.add_option("-d", "--dbname", dest="dbname",
                help="Database Name", metavar="FILE")
    parser.add_option("-l", "--list_file", dest="list_file",
                help="File contains SeqNames", metavar="FILE")
    parser.add_option("-m", "--motif_list_file", dest="motif_list_file",
                help="File contains motif name", metavar="FILE")
    parser.add_option("-p", "--pvalue", dest="pValue",
                default='0.05',
                help="Pvalue Threshold to output",
                metavar="FLOAT")
    parser.add_option("-q", "--qvalue", dest="qValue",
                default='1.0',
                help="Qvalue Threshold to output",
                metavar="FLOAT")
    parser.add_option("-u", "--upstream", dest="upstream",
                default=None,
                help="Only count upstream N bp",
                metavar="FLOAT")
    parser.add_option("-s", "--strand", dest="strand",
                default='+',
                help="Which strand to count",
                metavar="FLOAT")
    parser.add_option("-r", "--dist_range", dest="dist_range",
                nargs=2,
                default=(10, 250),
                help="The range of distance of motif pairs",
                metavar="INT")
    parser.add_option("-o", "--output", dest="output_path",
                help="Result output path",
                default="",
                metavar="DIR")
    parser.add_option("-a", "--low_RAM", dest="low_RAM",
                default=False,
                action="store_true",
                help="Run in low RAM mode to save Memory")
    parser.add_option("-b", "--big_RAM", dest="big_RAM",
                default=False,
                action="store_true",
                help="Keep All data in RAM to boost speed")
    parser.add_option("-n", "--parallel", dest="parallel",
                default=False,
                action="store_true",
                help="Run in parallel mode")
    return parser


if __name__ == '__main__':
    parser = getOptions()
    (options, args) = parser.parse_args()

    if not options.dbname or not options.list_file:
        parser.print_help()
        sys.exit(1)

    if options.parallel:
        options.parallel = pp.Server()

    if not options.output_path:
        abs_path = os.path.abspath(options.list_file)
        folder = getBasename(os.path.basename(abs_path))
        abs_path = os.path.dirname(abs_path)
        options.output_path = "%s/%s" % (abs_path, folder)

    if not os.path.exists(options.output_path):
        os.makedirs(options.output_path)

#    options.dbname = "RAP3k.db"
#    options.list_file = "RAP3k"

    timestamp = time.time()
    stats(options)
    print time.time() - timestamp

    #print united.dumpScan2Dict(2000, 1, 0.01, 0.1)
    #united.dumpScan2Zslv(200, "both", 0.05, 0.5)
    #print getDist([(11,13),(22,20)],[(18,16),(16,18)],0,200)
