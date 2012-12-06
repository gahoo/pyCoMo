#! /usr/bin/env python
#coding=utf-8
from subprocess import Popen, PIPE
from collections import defaultdict
import shelve
from jaspar import JasparDB
from plantcare import PlantCARE
from place import Place
#import zshelve


def iupac2meme(sequences):
    iupac2meme = Popen(
        args=['/home/gahoo/meme/bin/iupac2meme'] + sequences,
        stdout=PIPE)
    (meme_motifs, err) = iupac2meme.communicate()
    if err:
        print err
    return meme_motifs


def jaspar2meme(pfmpath):
    jaspar2meme = Popen(
        args=['/home/gahoo/meme/bin/jaspar2meme', '-pfm', pfmpath],
        stdout=PIPE)
    (meme_motifs, err) = jaspar2meme.communicate()
    if err:
        print err
    return meme_motifs


def split_meme(meme_motifs):
    motifs = defaultdict(str)
    motif = meme_motifs.split('\n')[:9]
    motif_id = 'header'
    for line in meme_motifs.split('\n')[9:]:
        if line.count("MOTIF"):
            motifs[motif_id] = "\n".join(motif)
            motif_id = line.split(' ')[1]
            motif = [line]
        else:
            motif.append(line)
    motifs[motif_id] = "\n".join(motif)
    return motifs


def join_meme(meme_motifs, motif_ids):
    '''
    meme_motifs[motif_id]=""
    motif_ids=[]
    '''
    return "\n".join([meme_motifs['header']] +
                 [meme_motifs[str(motif_id)] for motif_id in motif_ids])


if __name__ == '__main__':
    p = PlantCARE('plantcare.db')
    pl = Place('place.db')
    meme_slv = shelve.open("meme_motifs.slv")
    # meme_zslv = zshelve.open("meme_motifs.zslv")
    meme_motifs = jaspar2meme(
        '/home/gahoo/Project/jaspar/all_data/FlatFileDir/')
    # meme_motifs = split_meme(meme_motifs)
    #meme_slv['jaspar'] = split_meme(meme_motifs)
    meme_slv.update(split_meme(meme_motifs))
    # meme_zslv['jaspar'] = split_meme(meme_motifs)
    # print join_meme(meme_motifs, meme_motifs.keys()[5:9])
    meme_motifs = iupac2meme(p.getMotifSeq())
    # meme_motifs = split_meme(meme_motifs)
    #print join_meme(meme_motifs, meme_motifs.keys())
    #meme_slv['plantcare'] = split_meme(meme_motifs)
    meme_slv.update(split_meme(meme_motifs))
    # meme_zslv['plantcare'] = split_meme(meme_motifs)
    #print meme_zslv['plantcare']['header']
    #meme_slv.close()
    meme_motifs = iupac2meme(pl.getMotifSeq())
    #meme_slv['place'] = split_meme(meme_motifs)
    meme_slv.update(split_meme(meme_motifs))
    j = JasparDB('Jaspar.db')
    print meme_slv.keys()
    # print join_meme(meme_slv['jaspar'], j.getBASE_ID({"COLLECTION": ["CORE"], "tax_group": ["vertebrates", "insect"], "family": ["GATA"]}))
    print join_meme(meme_slv, j.getBASE_ID({"COLLECTION": ["CORE"], "tax_group": ["vertebrates", "insect"], "family": ["GATA"]}))
    # print join_meme(meme_slv['plantcare'], p.getMotifSeq({'Organism': [u'Antirrhinum majus', u'Flaveria trinervia', u'Nicotiana plumbaginifolia']}))
    # print join_meme(meme_slv['place'], pl.getMotifSeq({'Keyword': ['phyA', 'chloroplast']}))
    # # meme_zslv.close()
