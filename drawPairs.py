#! /usr/bin/env python
#coding=utf-8
import gzip
import os
import sys
import zshelve
import rpy2.robjects as robj
from optparse import OptionParser
from Stats import loadList, newZslv, getBasename, getFileSuffix
from rpy2.robjects.packages import importr
from rpy2.rinterface import RRuntimeError
from misc.united_db import UnitiedDB

def drawRGraph(motif_dist, enriched, outpath):
    '''
    motif_dist[motif][SeqName]=[(start, end),...]
    use in family mode only
    '''
    grdevices = importr('grDevices')
    graphics = importr('graphics')
    geneplotter = importr('geneplotter')
    RColorBrewer = importr('RColorBrewer')
    xlim = robj.IntVector([-3000, 0])
    ylim = robj.FloatVector([0, 0.002])
    for motif in enriched.keys():
        print motif
        if not motif_dist[motif]:
            print "empty"
            continue
        motif_seq = {}
        #画multidensity
        for MotifSeq in motif_dist[motif]:
            #merge+,-
            pos = mergeList(motif_dist[motif][MotifSeq])
            if len(pos) > 1:
                pos = map(lambda x: x - 3000, pos)
                motif_seq[MotifSeq] = robj.IntVector(pos)
        #merge MotifSeq
        dist = mergeList(motif_seq)
        #print motif_seq
        if not motif_seq:
            print "empty"
            continue
        motif_seq = robj.ListVector(motif_seq)
        names = robj.r('names(sort(sapply(%s, length), decreasing=T))' % motif_seq.r_repr())
        #print motif_seq.r_repr(),names.r_repr()
        motif = motif.replace('/', '_')
        grdevices.png(file="%s/%s.png" % (outpath, motif), width=512, height=512)
        dist = robj.IntVector(dist)
        density = robj.r.density(dist)
        if len(motif_seq) > 1:
            geneplotter.multidensity(motif_seq.rx(names), lwd=3, xlab="Distribution", main=motif, xlim=xlim, ylim=ylim)
        else:
            graphics.plot(density, lty='dashed', lwd=3, main=motif, xlab="Distribution", xlim=xlim, ylim=ylim)
        #graphics.hist(dist,add=True,breaks=50)
        graphics.rug(dist)
        graphics.lines(density, lty='dashed', lwd='4')
        grdevices.dev_off()
    return


def drawRboxplot(united, Motifs, enriched, REF_SeqNames, output):
    '''
    motif_dist[motif][SeqName]=[(start, end),...]
    '''
    grdevices = importr('grDevices')
    graphics = importr('graphics')
    #Merged_dist=dict([(k,list(v)) for k,v in Merged_dist.items()])
    Positions = {}
    for motif in enriched:
        Pos = []
        for seqname in REF_SeqNames:
            Pos.extend([start for start, stop in \
                Motifs[str(united.Motif2uid[motif])][seqname]])
        if Pos:
            Positions[motif] = robj.IntVector(map(lambda x: x - 3000, Pos))

    Positions = robj.ListVector(Positions)
    enriched_names = robj.StrVector(enriched)
    grdevices.png(file="%s/%s_box.png" % \
        (output, output), width=512, height=512)
    margin = robj.IntVector([3, 9, 4, 2])
    graphics.par(mar=margin)
    graphics.boxplot(Positions.rx(enriched_names),
        main="Boxplot of Motif Positions",
        horizontal=True, las=1, col='lightblue')
    grdevices.dev_off()


def drawRbarplot(united, Motifs, Nodes, output, top=25):
    if not enriched:
        return
    grdevices = importr('grDevices')
    graphics = importr('graphics')

    Counts = [len(seqnames) for seqnames in Motifs.itervalues()]
    Counts_names = [united.uid2Motif[int(uid)] for uid in Motifs.iterkeys()]
    Counts = robj.IntVector(Counts)
    Counts.names = robj.StrVector(Counts_names)

    Ps = [node[-1] for node in Result['Nodes'] if node[-1] < 0.05]
    Ps_names = [node[0] for node in Result['Nodes'] if node[-1] < 0.05]
    Ps = robj.FloatVector(Ps)
    Ps.names = robj.StrVector(Ps_names)
    Ps = robj.r.sort(Ps, decreasing=True)

    enriched_counts = Counts.rx(Ps.names)
    grdevices.png(file="%s/%s_bar.png" % (output, output),
        width=512, height=512)
    margin = robj.IntVector([3, 9, 4, 2])
    graphics.par(mar=margin)
    bar = graphics.barplot(enriched_counts,
        main="Enriched Motifs Counts", horiz=True, las=1, col='lightblue')
    Ps_lab = robj.r('format(signif(%s ,digits=2), scientific=T)' % Ps.r_repr())
    graphics.text(x=enriched_counts, y=bar,
        label=Ps_lab, po=2)
    #graphics.text(bar,labels=top_counts,pos=4,offset=10)
    grdevices.dev_off()


def drawRCoocur(co_Occur, Edges, outpath):
    '''
    co_Occur[(motifa,motifb)][REF_SeqName]={(motifa_pos,motifb_pos):distance}
    '''
    for edge in Edges:
        (uidA, uidB) = (united.Motif2uid[edge[0]], united.Motif2uid[edge[1]])
        motif_pair = "%s,%s" % (uidA, uidB)
        print edge[:2], len(co_Occur[motif_pair])
        drawRdistance(co_Occur[motif_pair], outpath, edge[:2])
        drawRCoDistibution(co_Occur[motif_pair], outpath, edge[:2])


def drawRdistance(motif_pair, outpath, pair_name):
    filename = "&".join(pair_name).replace('/', '_')
    distances = []
    for REF_SeqName in motif_pair.keys():
        distances.extend(motif_pair[REF_SeqName].values())
    dist = robj.IntVector(distances)
    grdevices = importr('grDevices')
    graphics = importr('graphics')
    xlim = robj.IntVector([-250, 250])
    grdevices.png(file="%s/%s_distance.png" % (outpath, filename), width=512, height=512)
    graphics.par(cex=2.3, mar=robj.FloatVector([2.0, 2.0, 0.2, 0.2]))
    graphics.hist(dist, lwd=5, breaks=50, border=0, col=1, main="", xlab="", xlim=xlim, font=2)
    graphics.axis(1, at=robj.IntVector([-200, -100, 0, 100, 200]))
    #graphics.hist(dist,breaks=50,border=4,main="Distance Histogram of \n%s" % "&".join(pair_name),xlab="Distances")
    grdevices.dev_off()


def drawRCoDistibution(motif_pair, outpath, pair_name):
    motifA_pos = []
    motifB_pos = []
    filename = "&".join(pair_name).replace('/', '_')
    (motifa, motifb) = pair_name
    #print motif_pair
    for REF_SeqName in motif_pair.keys():
        #print motif_pair[REF_SeqName].keys()
        for motifa_pos, motifb_pos in motif_pair[REF_SeqName].keys():
            motifA_pos.append(motifa_pos[0])
            motifB_pos.append(motifb_pos[0])
            #print motifa_pos,motifb_pos
    grdevices = importr('grDevices')
    graphics = importr('graphics')
    #!!!R must have geneplotter package installed
    geneplotter = importr('geneplotter')
    RColorBrewer = importr('RColorBrewer')
    motifA_pos = robj.IntVector(map(lambda x: x - 3000, motifA_pos))
    motifB_pos = robj.IntVector(map(lambda x: x - 3000, motifB_pos))
    Pos = {motifa: motifA_pos, motifb: motifB_pos}
    Pos = robj.ListVector(Pos)
    grdevices.png(file="%s/%s_distribution.png" % (outpath, filename), width=512, height=512)
    color = RColorBrewer.brewer_pal(n=9, name="Set1")
    xlim = robj.IntVector([-3000, 0])
    leg = robj.ListVector({'x': 'topleft', 'legend': robj.StrVector([motifa, motifb]), 'fill': color})
    graphics.par(cex=2.5, mar=robj.FloatVector([2.0, 2.0, 0.2, 0.2]))
    #geneplotter.multidensity(Pos.rx(),lwd=3,xlab="Distribution",main="Distribution of \n%s" % filename,legend=leg,font=2)
    geneplotter.multidensity(Pos.rx(), lwd=10, xlab="", main="", ylab="", legend=leg, font=2, xlim=xlim)
    graphics.axis(2, lwd=5)
    graphics.axis(1, lwd=5)
    #geneplotter.multidensity(Pos.rx(),lwd=3,xlab="Distribution",main="Distribution of \n%s" % filename)
    graphics.rug(motifA_pos, col=4, lwd=3)
    graphics.rug(motifB_pos, col=2, lwd=3)
    grdevices.dev_off()
    #Scatter plot
    grdevices.png(file="%s/%s_Scatter.png" % (outpath, filename), width=512, height=512)
    limit = robj.IntVector([-3000, 0])
    '''
    graphics.plot(motifA_pos,motifB_pos,main="Position Scatter Plot of\n%s&%s" % (motifa,motifb), \
                    xlab="Positions of %s" % motifa, \
                    ylab="Positions of %s" % motifb, \
                    xlim=limit,ylim=limit)
    '''
    graphics.par(cex=2, mar=robj.FloatVector([2.0, 2.0, 0.2, 0.2]))
    graphics.plot(motifA_pos, motifB_pos, main="", xlab="", ylab="", xlim=limit, ylim=limit, lwd=2, font=2)
    graphics.abline(1, 1, lwd=5)
    graphics.axis(2, lwd=5)
    graphics.axis(1, lwd=5)
    graphics.text(-1000, -2500, motifa, cex=1.2, font=2)
    graphics.text(-2000, -500, motifb, cex=1.2, font=2)
    grdevices.dev_off()


def loadZslv(options, zslv_type, tohash=None):
    zslv_file = getBasename(os.path.basename(options.dbname)) + getFileSuffix(options, zslv_type, tohash) + ".zslv"
    abs_zslv_file = "%s/%s" % (options.output_path, zslv_file)
    if os.path.exists(abs_zslv_file):
        return zshelve.open("%s/%s" % (options.output_path, zslv_file))
    else:
        print "%s not exists, run Stats.py first" % abs_zslv_file
        sys.exit(1)


def getOptions():
    parser = OptionParser(usage="Usage: drawRgraph [options] -d dbname -l <list_file>")
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

    return parser

if __name__ == '__main__':
    parser = getOptions()
    (options, args) = parser.parse_args()

    if not options.output_path:
        abs_path = os.path.abspath(options.list_file)
        options.output_path = getBasename(os.path.basename(abs_path))
        # abs_path = os.path.dirname(abs_path)
        # options.output_path = "%s/%s" % (abs_path, basename)

    if not os.path.exists(options.output_path):
        os.makedirs(options.output_path)

    united = UnitiedDB(options.dbname)
    REF_SeqNames = [united.SeqName2REF_SeqName[seqname] for seqname in loadList(options.list_file)]
    if options.motif_list_file:
        REF_uids = [united.Motif2uid[motif] for motif in loadList(options.motif_list_file)]
    else:
        REF_uids = []
    # Motifs = zshelve.open(options.motif_zslv)
    # PairedMotifs = zshelve.open(options.paired_zslv)
    # Result = zshelve.open(options.result_file)

    Motifs = loadZslv(options, "MotifScan_%s" % options.output_path, REF_SeqNames + REF_uids)
    PairedMotifs = loadZslv(options, "PairedMotifs_%s" % options.output_path, Motifs.keys())
    PairedMotifs_Cnts = loadZslv(options, "PairedMotifsCnts_%s" % options.output_path, Motifs.keys())
    Result = loadZslv(options, "Result_%s" % options.output_path, REF_SeqNames + REF_uids)

    enriched = [node[0] for node in Result['Nodes'] if node[-1] < 0.05]
    #drawRGraph(Merged_MotifSeq_dist, enriched, "%s/Enriched" % output)
#    drawRboxplot(united, Motifs, enriched, REF_SeqNames, options.output_path)
    drawRbarplot(united, Motifs, enriched, options.output_path, top=25)
    drawRCoocur(PairedMotifs, Result['Edges'], options.output_path)
