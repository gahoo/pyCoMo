#! /usr/bin/env python
#coding=utf-8
from sqlitedb import sqlitedb
from collections import defaultdict
import sqlite3
import os


class PlantCARE(sqlitedb):
    '''
    具体的PlantCARE数据库类，封装了对数据库的基本操作
    方法：
    loadDB2dict：加载数据库到字典LOC_ID/Organism：Motif：Sequence：中
    newDB：初始化数据库（建表和索引）
    buildDB：导入指定目录下的FIMO结果文件进入数据库
    importMotifs：导入PlantCARE那里得到的Motif信息
    addEntry：Scan结果的一个条目（一行）进入数据库（查询数据库保证无重复）
    readCARE：读取PlantCARE文件中的内容至列表中
    属性：
    tables：表的结构
    views：视图的结构
    LOC_ID：LOC_ID表的内容（字典）
    Organism：Organism表的内容（字典）
    Motif：Motif表的内容（字典）
    Seq：Seq表的内容（字典）
    '''
    tables = [
        """CREATE TABLE Organism (
            id INTEGER PRIMARY KEY,
            Organism VARCHAR(50) NOT NULL UNIQUE)
        """,
        """CREATE TABLE Motif (
            id INTEGER PRIMARY KEY,
            Motif VARCHAR(20) NOT NULL UNIQUE,
            Description VARCHAR(120),
            Note VARCHAR(100))
        """,
        """CREATE TABLE Accession (
            id INTEGER PRIMARY KEY,
            Accession VARCHAR(50) NOT NULL UNIQUE)
        """,
        """CREATE TABLE Instance (
            id INTEGER PRIMARY KEY,
            REF_Motif INTEGER REFERENCES Motif(id),
            Type VARCHAR(12) NOT NULL,
            REF_Organism INTEGER REFERENCES Organism(id),
            Description VARCHAR(120),
            REF_Accession INTEGER REFERENCES Accession(id),
            REF_MotifSeq  INTEGER REFERENCES MotifSeq(id))
        """,
        """CREATE TABLE MotifSeq (
            id INTEGER PRIMARY KEY,
            Seq VARCHAR(120) NOT NULL UNIQUE)
        """
        ]

    views = [
    """CREATE VIEW Main AS
        SELECT Motif,Type,Organism,Motif.Description,Accession,Seq
        FROM Instance,Motif,Organism,Accession,MotifSeq
        WHERE Motif.id=Instance.REF_Motif
            AND Organism.id=Instance.REF_Organism
            AND Accession.id=Instance.REF_Accession
            AND MotifSeq.id=Instance.REF_MotifSeq
    """,
    """CREATE VIEW Motif_Seq AS
        SELECT DISTINCT Motif,Seq
        FROM Instance,Motif,MotifSeq
        WHERE Motif.id=Instance.REF_Motif
            AND MotifSeq.id=Instance.REF_MotifSeq
    """
    ]

    # Accession = {}
    # Organism = {}
    # Motif = {}
    # MotifSeq = {}
    # SeqName = {}
    (Accession, Organism, Motif, MotifSeq, Description, Motif_Seq) = \
     (defaultdict(str), defaultdict(str), defaultdict(str), defaultdict(str), defaultdict(str), defaultdict(str))
    dbfile = ''

    def __init__(self, dbfile):
        '''
        dbfile数据库文件路径
        '''
        self.dbfile = dbfile
        if not os.path.exists(dbfile):
            #若不存在数据库则新建之
            sqlitedb.__init__(self, dbfile)
            self.newDB(dbfile)
        else:
            sqlitedb.__init__(self, dbfile)
            self.loadDB()

    def loadDB(self):
        self.loadDB2dict('Motif', "Motif", 'id', self.Motif)
        self.loadDB2dict('MotifSeq', "Seq", 'id', self.MotifSeq)
        self.loadDB2dict('Organism', "Organism", 'id', self.Organism)
        self.loadDB2dict('Accession', "Accession", 'id', self.Accession)
        self.loadDB2dict('Motif_Seq', "Seq", 'Motif', self.Motif_Seq)

    def loadDB2dict(self, tname, kcol, vcol, dic, where=None):
        '''
        查询数据库中的表，指定两列分别作为字典的键和值
        tname：表名
        kcol：作为键的列
        vcol：作为值的列
        return 字典
        '''
        for rs in self.select(tname, [kcol, vcol], where, distinct=True):
            #dic[rs[kcol].encode('gbk')] = rs[vcol]
            dic[rs[0]] = rs[1]

    def newDB(self, dbfile):
        '''
        根据tables和view   s属性建立表和视图
        '''
        for SQL in self.tables + self.views:
            try:
                self.cur.execute(SQL)
            except sqlite3.OperationalError, Error:
                print "ERROR:", Error

    def importMotifs(self, care_file):
        '''
        导入PlantCARE那里得到的Motif信息
        '''

        def multi2iupac(Sequence):
            '''
            将形如(A/T)的序列转成IUPAC字符
            '''
            iupac = {
                'M': ['A', 'C'],
                'R': ['A', 'G'],
                'W': ['A', 'T'],
                'S': ['C', 'G'],
                'Y': ['C', 'T'],
                'K': ['G', 'T']}
            Sequence = Sequence.upper()
            for key in iupac.keys():
                Sequence = Sequence.replace("(" + "/".join(iupac[key]) + ")", key)
                Sequence = Sequence.replace("(" + "/".join(iupac[key][::-1]) + ")", key)
            return Sequence

        for CARE in self.readCARE(care_file):
            if CARE[0] == '':
                #is Instance
                if CARE[4] == '':
                    #has no Sequence
                    continue
                REF_Accession = self.addEntry('Accession',
                ['Accession'],
                [CARE[3]],
                self.Accession, 0)
                REF_MotifSeq = self.addEntry('MotifSeq',
                    ['Seq'],
                    [multi2iupac(CARE[4])],
                    self.MotifSeq, 0)
                self.addEntry('Instance',
                ['REF_Motif', 'Type', 'REF_Organism', 'Description', 'REF_Accession', 'REF_MotifSeq'],
                [REF_Motif, CARE[1], REF_Organism, CARE[2], REF_Accession, REF_MotifSeq])
            else:
                #is Motif
                REF_Organism = self.addEntry('Organism',
                    ['Organism'],
                    [CARE[1]],
                    self.Organism, 0)
                REF_Motif = self.addEntry('Motif',
                    ['Motif', 'Description', 'Note'],
                    [CARE[0], CARE[2], CARE[3]],
                    self.Motif, 0)
        self.commit()
        self.loadDB()
        return

    def addEntry(self, tname, colums, values, dic=None, kcol=None):
        '''
        向指定表格插入条目，用字典控制重复
        效率较高
        tname：表名
        colums：要插入的列
        values：要插入的值
        dic：要插入的表所对应的属性，如LOC_ID、Motif等
        kcol：根据values[kcol]来检查是否重复
        return 插入条目的id
        '''
        if kcol == None and dic == None:
            #即导入CARE的情况下
            self.insert(tname, colums, values, False)
        elif values[kcol] in dic.keys():
            #有key的情况下
            return dic[values[kcol]]
        else:
            #没key就插入
            self.insert(tname, colums, values, False)
            id = len(dic) + 1
            dic[values[kcol]] = id
            return id

    def readCARE(self, care_file):
        '''
        读取PlantCARE的内容至列表care中
        care_file：care文件路径
        return care（列表）
        '''
        file = open(care_file)
        care = []
        for line in file.readlines():
            care.append(line.strip('\n').split('\t'))
        return care

    def getMotifSeq(self, where=None):
        if where:
            for key in where.keys():
                if key == 'Motif':
                    where['REF_Motif'] = [str(self.Motif[motif]) for motif in where.pop('Motif')]
                elif key == 'Organism':
                    where['REF_Organism'] = [str(self.Organism[organism]) for organism in where.pop('Organism')]
                elif key == 'Accession':
                    where['REF_Accession'] = [str(self.Accession[accession]) for accession in where.pop('Accession')]
        REF_MotifSeq = [str(ref[0]) for ref in self.select('Instance', ['REF_MotifSeq'], where=where, distinct=True)]
        Seqs = [seq[0] for seq in self.select('MotifSeq', ['Seq'], {'id': REF_MotifSeq})]
        return Seqs


if __name__ == '__main__':
    care = PlantCARE('plantcare.db')
    care.importMotifs('CARE.txt')
    print care.Accession
    print care.Organism.keys()
    print care.getMotifSeq({'Organism': [u'Antirrhinum majus', u'Flaveria trinervia', u'Nicotiana plumbaginifolia']})
    print care.getMotifSeq()
