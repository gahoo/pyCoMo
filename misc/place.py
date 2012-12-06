#! /usr/bin/env python
#coding=utf-8
from sqlitedb import sqlitedb
from collections import defaultdict
import sqlite3
import os
import re


class Place(sqlitedb):
    """docstring for Place"""

    tables = [
    """CREATE TABLE Motif (
            id INTEGER PRIMARY KEY,
            Motif VARCHAR(25) NOT NULL UNIQUE,
            Accession CHAR(7),
            Description VARCHAR(520),
            Seq VARCHAR(35) NOT NULL UNIQUE)
    """,
    """CREATE TABLE Organism (
            id INTEGER NOT NULL REFERENCES Motif(id),
            Organism VARCHAR(220) NOT NULL)
    """,
    """CREATE TABLE Keyword (
            id INTEGER NOT NULL REFERENCES Motif(id),
            Keyword VARCHAR(40) NOT NULL)
    """,
    """CREATE TABLE Pubmed (
            id INTEGER NOT NULL REFERENCES Motif(id),
            Pubmed INTEGER NOT NULL)
    """]

    views = [
    """CREATE VIEW Main AS
        SELECT Motif.id AS id,Motif,Accession,Description,
               group_concat(Distinct Organism) AS Organism,
               group_concat(Distinct Keyword) AS Keyword,
               group_concat(Distinct Pubmed) AS Pubmed
        FROM Motif,Organism,Keyword,Pubmed
        WHERE Motif.id=Organism.id
        AND Motif.id=Keyword.id
        AND Motif.id=Pubmed.id
        GROUP BY Motif.id
    """]

    (Motif, Keyword, Description, Accession, Seq, Organism) = \
        (defaultdict(str), defaultdict(str), defaultdict(str),
        defaultdict(str), defaultdict(str), defaultdict(str))

    def __init__(self, dbfile):
        if not os.path.exists(dbfile):
            #若不存在数据库则新建之
            sqlitedb.__init__(self, dbfile)
            self.newDB(dbfile)
        else:
            sqlitedb.__init__(self, dbfile)
            self.loadDB()

    def newDB(self, dbfile):
        '''
        根据tables和views属性建立表和视图
        '''
        for SQL in self.tables + self.views:
            try:
                self.cur.execute(SQL)
            except sqlite3.OperationalError, Error:
                print "ERROR:", Error

    def loadDB2dict(self, tname, kcol, vcol, dic, where=None, group=None):
        '''
        查询数据库中的表，指定两列分别作为字典的键和值
        tname：表名
        kcol：作为键的列
        vcol：作为值的列
        return 字典
        '''
        for rs in self.select(tname, [kcol, vcol], where, group, distinct=True):
            #dic[rs[kcol].encode('gbk')] = rs[vcol]
            dic[rs[0]] = rs[1]

    def loadDB(self):
        self.loadDB2dict('Motif', "id", 'Motif', self.Motif)
        self.loadDB2dict('Motif', "id", 'Accession', self.Accession)
        self.loadDB2dict('Motif', "id", 'Description', self.Description)
        self.loadDB2dict('Motif', "id", 'Seq', self.Seq)
        self.loadDB2dict('Keyword', "Keyword", 'group_concat(Distinct id)', self.Keyword, group="Keyword")
        for k in self.Keyword.keys():
            self.Keyword[k] = map(int, self.Keyword[k].split(','))
        self.loadDB2dict('Organism', "Organism", 'group_concat(Distinct id)', self.Organism, group="Organism")
        for k in self.Organism.keys():
            self.Organism[k] = map(int, self.Organism[k].split(','))

    def importTable(self, place_tab_file):
        regex = re.compile(r'PubMed: *(\d+)')
        #get the cur of table
        i = len(self.select('Motif'))
        for line in open(place_tab_file).readlines():
            line = line.strip().split('\t')
            (motif, accession, date, description, keyword, Organism, seq) = line[:7]
            #i add 1 when insert successfully
            if self.insert('Motif', ["Motif", "Accession", "Description", "Seq"], \
                        [motif, accession, description.replace('"', "'"), seq.strip(' ')]):
                i += 1
            else:
                continue
            pubmed = regex.findall("\t".join(line[7:]))
            self.addOtherTable('Pubmed', i, pubmed)
            keyword = keyword.strip('; ').replace('"', "'").split('; ')
            self.addOtherTable('Keyword', i, keyword)
            Organism = Organism.strip('; ').split('; ')
            self.addOtherTable('Organism', i, Organism)
            self.commit()
        self.loadDB()

    def addOtherTable(self, tname, cur, vals):
        for v in vals:
            self.insert(tname, ['id', tname], [cur, v])

    def getMotifSeq(self, where={}):
        #print self.select('Organism', None, where={'Organism': 'rice (Oryza sativa)'}, like=True)
        #print self.select('Organism', None, where={'Organism': ['rice (Oryza sativa)']}, like=False)
        IDs = set(self.Motif.keys())
        Motif2id = dict([(v, k) for (k, v) in self.Motif.items()])
        Accession2id = dict([(v, k) for (k, v) in self.Accession.items()])
        for (col, vals) in where.items():
            if col == 'Motif':
                IDs = IDs & set([Motif2id[motif] for motif in vals])
            elif col == 'Accession':
                IDs = IDs & set([Accession2id[accession] for accession in vals])
            elif col == 'Organism':
                tmp_ids = []
                for val in vals:
                    tmp_ids.extend(self.Organism[val])
                IDs = IDs & set(tmp_ids)
            elif col == 'Keyword':
                tmp_ids = []
                for val in vals:
                    tmp_ids.extend(self.Keyword[val])
                IDs = IDs & set(tmp_ids)
        return [self.Seq[ID] for ID in IDs]

if __name__ == '__main__':
    pl = Place('place.db')
    #pl.importTable('place.tab')
    pl.getMotifSeq()
    print pl.Keyword
    print pl.getMotifSeq({'Keyword': ['phyA', 'chloroplast']})
    print pl.getMotifSeq()
