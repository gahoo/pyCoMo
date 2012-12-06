#! /usr/bin/env python
#coding=utf-8
from misc.jaspar import JasparDB
from misc.plantcare import PlantCARE
from misc.place import Place
from misc.united_db import UnitiedDB
from misc.meme_motif_parsing import jaspar2meme, iupac2meme, split_meme
from optparse import OptionParser
import shelve
import os
import sys

parser = OptionParser()
parser.add_option("-n", "--dbname", dest="dbname", help="Database Name", metavar="FILE")
parser.add_option("-j", "--jaspar", dest="jaspar", help="Import Jaspar", default=False, action="store_true")
parser.add_option("-c", "--plantcare", dest="plantcare", help="Import PlantCARE", default=False, action="store_true")
parser.add_option("-p", "--place", dest="place", help="Import PLACE", default=False, action="store_true")
parser.add_option("-d", "--denovo", dest="denovo_path", help="Import Denovo Motifs(MEME format only)", metavar="DIR")

(options, args) = parser.parse_args()

if not options.dbname:
    parser.print_help()
    sys.exit(1)

jaspar_path = './databases/jaspar.db'
plantcare_path = './databases/plantcare.db'
place_path = './databases/place.db'

meme_slv = shelve.open(options.dbname + ".slv")
united = UnitiedDB(options.dbname + ".db")
#initialize databases
if options.jaspar:
    if not os.path.exists(jaspar_path):
        jaspar = JasparDB(jaspar_path)
        jaspar.importTables('./databases/jaspar/')
    else:
        jaspar = JasparDB(jaspar_path)
    meme_motifs = jaspar2meme('./databases/jaspar_pfm/')
    meme_slv.update(split_meme(meme_motifs))
    united.importDB('Jaspar', jaspar_path)

if options.plantcare:
    if not os.path.exists(plantcare_path):
        plantcare = PlantCARE(plantcare_path)
        plantcare.importMotifs('./databases/plantcare.txt')
    else:
        plantcare = PlantCARE(plantcare_path)
    meme_motifs = iupac2meme(plantcare.getMotifSeq())
    meme_slv.update(split_meme(meme_motifs))
    united.importDB('PlantCARE', plantcare_path)

if options.place:
    if not os.path.exists(place_path):
        place = Place(place_path)
        place.importTable('./databases/place.tab')
    else:
        place = Place(place_path)
    meme_motifs = iupac2meme(place.getMotifSeq())
    meme_slv.update(split_meme(meme_motifs))
    united.importDB('PLACE', place_path)

if options.denovo_path:
    meme_motifs = open(options.denovo_path).read()
    meme_slv.update(split_meme(meme_motifs))
    united.importDenovo(options.denovo_path, 'MEME')
