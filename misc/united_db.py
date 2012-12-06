#! /usr/bin/env python
#coding=utf-8
from sqlitedb import sqlitedb
from collections import defaultdict
from hashlib import md5
from meme_motif_parsing import join_meme
import sqlite3
#import shelve
import zshelve
import os
import sys


class UnitiedDB(sqlitedb):
    '''
    combined database
    '''

    tables = [
    """CREATE TABLE Motif(
            uid INTEGER PRIMARY KEY,
            Motif VARCHAR (48) NOT NULL,
            Description VARCHAR(520),
            Database VARCHAR(9),
            REF_ID INTEGER NOT NULL,
            MEME_ID VARCHAR(120) NOT NULL)
    """,
    """CREATE TABLE Scan(
            id INTEGER PRIMARY KEY,
            REF_uid INTEGER NOT NULL REFERENCES Motif(uid),
            REF_SeqName INTEGER NOT NULL REFERENCES SeqName(id),
            start INTEGER NOT NULL,
            stop INTEGER NOT NULL,
            strand BOOLEAN NOT NULL,
            score REAL NOT NULL,
            pValue REAL NOT NULL,
            qValue REAL NOT NULL,
            REF_FoundMotif INTEGER REFERENCES FoundMotif(id))
    """,
    """CREATE TABLE SeqName(
            id INTEGER PRIMARY KEY,
            SeqName VARCHAR(30) NOT NULL UNIQUE)
    """,
    """CREATE TABLE FoundMotif (
            id INTEGER PRIMARY KEY,
            FoundMotif VARCHAR(120) NOT NULL UNIQUE)
    """,
    """CREATE TABLE MotifCounts (
            REF_uid INTEGER NOT NULL REFERENCES Motif(uid),
            Counts INTEGER NOT NULL,
            Range INTEGER,
            strand BOOLEAN NOT NULL,
            pValue REAL,
            qValue REAL,
            CONSTRAINT uidRangeCnt PRIMARY KEY (REF_uid,Counts,Range,strand,pValue,qValue))
    """,
    """CREATE TABLE CoMotifCounts (
            REF_uidA INTEGER NOT NULL REFERENCES Motif(uid),
            REF_uidB INTEGER NOT NULL REFERENCES Motif(uid),
            Counts INTEGER NOT NULL,
            Range INTEGER,
            strand BOOLEAN NOT NULL,
            pValue REAL,
            qValue REAL,
            dist_min INTEGER NOT NULL,
            dist_max INTEGER NOT NULL,
            CONSTRAINT uidCoRangeCnt PRIMARY KEY (REF_uidA,REF_uidB,Counts,Range,strand,pValue,qValue,dist_min,dist_max))
    """]

    views = []
    MEME_ID2uid = defaultdict(int)
    Motif2uid = defaultdict(int)
    uid2Motif = defaultdict(str)
    uid2Desc = defaultdict(str)
    uid2MEME_ID = defaultdict(str)
    Motif2MEME_ID = defaultdict(str)
    SeqName2REF_SeqName = defaultdict(int)
    FoundMotif2REF_FoundMotif = defaultdict(int)

    def __init__(self, dbfile):
        self.dbname = ".".join(dbfile.split(".")[:-1])
        self.dbpath = os.path.abspath(dbfile)
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

    def importDB(self, dbname, dbfile):
        '''
        import PLACE, PlantCARE, Jaspar.
        dbname = ["PLACE", "PlantCARE", "Jaspar"]
        dbfile = The seperated db file of place, plantcare, Jaspar
        '''
        if dbname == "PLACE":
            SQL = """INSERT INTO Motif (Motif,Description,Database,REF_ID,MEME_ID)
                    SELECT Motif,Description,"PLACE",id,Seq
                    FROM PLACE.Motif"""
        elif dbname == "PlantCARE":
            SQL = """INSERT INTO Motif (Motif,Description,Database,REF_ID,MEME_ID)
                    SELECT DISTINCT Motif,Instance.Description,"PlantCARE",Instance.id,Seq
                    FROM PlantCARE.Motif,PlantCARE.Instance,PlantCARE.MotifSeq
                    WHERE Instance.REF_Motif=Motif.id
                    AND Instance.REF_MotifSeq=MotifSeq.id
                    GROUP BY Seq
                    HAVING COUNT(*)>=1"""
        elif dbname == "Jaspar":
            SQL = """INSERT INTO Motif (Motif,Description,Database,REF_ID,MEME_ID)
                    SELECT NAME,group_concat(VAL),"Jaspar",MATRIX.ID,BASE_ID || '.' || VERSION
                    FROM Jaspar.MATRIX,Jaspar.MATRIX_ANNOTATION
                    WHERE MATRIX.ID=MATRIX_ANNOTATION.ID
                    GROUP BY MATRIX.ID"""
        self.cur.execute("ATTACH '%s' AS %s" % (dbfile, dbname))
        self.cur.execute(SQL)
        self.commit()

    def importDenovo(self, motiffile, format):
        '''
        import the MEME format Motif information
        motiffile = the file that contains motif
        '''
        if format.upper() == "MEME":
            i = 0
            for line in open(motiffile).readlines():
                if line.count("MOTIF"):
                    try:
                        (motif_name, motif_desc) = line.strip().split(' ')[1:]
                    except ValueError:
                        motif_name = line.strip().split(' ')[1]
                        motif_desc = ""
                    i += 1
                    self.insert('motif',
                            ['Motif', "Description", "Database", "REF_ID", "MEME_ID"],
                            [motif_name, motif_desc, "Denovo", i, motif_name])
        self.commit()

    def exportMotifs(self, MEME_IDs, meme_slv, meme_file=""):
        '''
        export Motifs in the meme_slv with given MEME_IDs
        when meme_file is set, result will write into the files.
        MEME_IDs: matches the Motif.MEME_ID
        meme_slv[MEME_ID]=meme_motifs: python shelve object contains MEME format motifs
        meme_file: output filename
        '''
        if meme_file:
            motif_file = open(meme_file, 'w')
            motif_file.write(join_meme(meme_slv, MEME_IDs))
            motif_file.close()
        else:
            return join_meme(meme_slv, MEME_IDs)

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

    def loadDB(self):
        '''
        load MEME_ID2uid,Motif2MEME_ID,SeqName2REF_SeqName,FoundMotif2REF_FoundMotif
        '''
        self.loadDB2dict('Motif', "MEME_ID", 'uid', self.MEME_ID2uid)
        self.loadDB2dict('Motif', "Motif", 'uid', self.Motif2uid)
        self.loadDB2dict('Motif', "uid", 'Motif', self.uid2Motif)
        self.loadDB2dict('Motif', "uid", 'Description', self.uid2Desc)
        self.loadDB2dict('Motif', "Motif", 'MEME_ID', self.Motif2MEME_ID)
        self.loadDB2dict('Motif', "uid", 'MEME_ID', self.uid2MEME_ID)
        self.loadDB2dict('SeqName', "SeqName", 'id', self.SeqName2REF_SeqName)

        #somtimes too big to load
        #self.loadDB2dict('FoundMotif', "FoundMotif", 'id', self.FoundMotif2REF_FoundMotif)

    def getMEME_ID(self, where, like=False):
        '''
        Load MEME_IDs from database
        where={"Database":["Jaspar","PlantCARE"],"Motif":[...]}: SQL that subset the data.
        like: whether using the fuzzy like in SQL
        '''
        return [rs[0] for rs in self.select('Motif', ['MEME_ID'], where=where, like=like, distinct=True)]

    def importScan(self, fimo_result):
        '''
        import the results of fimo
        fimo_result=[[MEME_ID, SeqName, start, stop, score, pValue, qValue, FoundMotif],...]
        '''
        #print "importing scanned result"
        for (MEME_ID, SeqName, start, stop, score, pValue, qValue, FoundMotif) in fimo_result:
            #(start, stop) = map(int, [start, stop])
            REF_uid = self.MEME_ID2uid[MEME_ID]
            #SeqName table
            REF_SeqName = self.addEntry('SeqName', ['SeqName'], [SeqName],
                        self.SeqName2REF_SeqName, 0)
            #FoundMotif table
            REF_FoundMotif = self.addEntry('FoundMotif', ['FoundMotif'], [FoundMotif],
                        self.FoundMotif2REF_FoundMotif, 0)
            #Scan table
            if int(stop) > int(start):
                #plus strand
                strand = 1
            else:
                strand = 0
            self.addEntry('Scan',
            ['REF_uid', 'REF_SeqName', 'start', 'stop', 'strand', 'score', 'pValue', 'qValue', 'REF_FoundMotif'],
            [REF_uid, REF_SeqName, start, stop, strand, score, pValue, qValue, REF_FoundMotif])
        self.commit()

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
        elif dic.has_key(values[kcol]):
            #有key的情况下
            return dic[values[kcol]]
        else:
            #没key就插入
            self.insert(tname, colums, values, False)
            id = len(dic) + 1
            dic[values[kcol]] = id
            return id

    def countMotifs(self, up_bp=None, strand=None, pValue=None, qValue=None, pair=False, dist_min=None, dist_max=None, REF_SeqName=None, REF_uid=None):
        '''
        统计数据库内全部Motif对应的SeqName数目
        #TODO 必须注意，如果选取的序列长度不同应该重新进行计算，否则背景是不同的，会导致Enrichment的计算有误，可以考虑增加一列记录长度的字段
        up_bp is the upstream base num if not set, all seq will count.
        statistic the SeqName num of a motif or motif pair
        up_bp: x bp upstream
        strand, pValue, qValue: subset the data
        pair: motif pair or not
        dist_min, dist_max: the distanct constraintion
        '''
        (colums, where) = self.__buildScanColumsWhere(up_bp, strand, pValue, qValue, REF_SeqName, REF_uid)
        colums = ",".join(map(str, colums))
        if where:
            if pair:
                where = "AND " + " AND ".join(["S1.%s AND S2.%s" % (w, w) for w in where])
                where = where + " AND ABS(S1.start - S2.start) > %d" % dist_min
                where = where + " AND ABS(S1.start - S2.start) < %d" % dist_max
                colums = colums + ",%d,%d" % (dist_min, dist_max)
            else:
                where = "WHERE " + " AND ".join(where)
        else:
            where = ""
        if pair and not REF_SeqName:
            SQL = """INSERT INTO CoMotifCounts(REF_uidA,REF_uidB,Counts,Range,strand,pValue,qValue,dist_min,dist_max)
                SELECT S1.REF_uid,S2.REF_uid,COUNT(DISTINCT S1.REF_SeqName),%s
                FROM Scan S1,Scan S2
                WHERE S1.REF_SeqName=S2.REF_SeqName
                AND S1.REF_uid!=S2.REF_uid
                %s
                GROUP BY S1.REF_uid,S2.REF_uid"""
        elif pair and REF_SeqName:
            SQL = """SELECT S1.REF_uid,S2.REF_uid,COUNT(DISTINCT S1.REF_SeqName)
                FROM Scan S1,Scan S2
                WHERE S1.REF_SeqName in (%s)
                AND S2.REF_SeqName in (%s)
                AND S1.REF_SeqName=S2.REF_SeqName
                AND S1.REF_uid!=S2.REF_uid
                %s
                GROUP BY S1.REF_uid,S2.REF_uid"""
            REF_SeqName = ",".join(map(str, REF_SeqName))
            #print SQL % (REF_SeqName, REF_SeqName, where)
            self.cur.execute(SQL % (REF_SeqName, REF_SeqName, where))
            return self.cur.fetchall()
        else:
            SQL = """INSERT INTO MotifCounts(REF_uid,Counts,Range,strand,pValue,qValue)
                SELECT REF_uid,count(distinct REF_SeqName),%s
                FROM Scan
                %s
                GROUP BY REF_uid"""
        #print SQL % (colums, where)
        self.cur.execute(SQL % (colums, where))
        self.commit()

    def __buildScanColumsWhere(self, up_bp=None, strand=None, pValue=None, qValue=None, REF_SeqName=None, REF_uid=None):
        '''
        build where statments for SQL
        up_bp: upstream x bp
        strand: [0:-,1:+,None:both]
        pValue: <0.05
        qValue: <1.0
        REF_SeqName: REF_SeqName in [1,2,3...]
        REF_uid: REF_uid in [1,2,3...]
        '''
        where = []
        if strand == None or strand == "both":
            strand = 2
        elif strand == '+' or strand == '1':
            strand = 1
            where.append("strand=%s" % strand)
        elif strand == '-' or strand == '0':
            strand = 0
            where.append("strand=%s" % strand)
        if pValue:
            where.append("pValue<=%s" % pValue)
        else:
            pValue = "NULL"
        if qValue:
            where.append("qValue<=%s" % qValue)
        else:
            qValue = "NULL"
        if up_bp:
            where.append("start>%d" % (3000 - up_bp))
        else:
            up_bp = "NULL"
        if REF_uid:
            where.append("REF_uid IN ('%s')" % "','".join(map(str, REF_uid)))
        if REF_SeqName:
            where.append("REF_SeqName IN ('%s')" % "','".join(map(str, REF_SeqName)))
        return ([up_bp, strand, pValue, qValue], where)

    def dumpScan2Dict_(self, up_bp=None, strand=None, pValue=None, qValue=None, dist_min=None, dist_max=None, REF_SeqName=None, REF_uid=None):
        '''
        return Motif[uid][REF_SeqNames]=[(start,stop),]
        '''
        motifs = defaultdict(dict)
        (colums, where) = self.__buildScanColumsWhere(up_bp, strand, pValue, qValue, REF_SeqName)
        SQL = """SELECT REF_SeqName,start,stop FROM Scan WHERE REF_uid=%s %s"""
        if not REF_uid:
            REF_uid = self.MEME_ID2uid.values()
        for REF_uid in REF_uid:
            self.cur.execute(SQL % (REF_uid, "AND " + " AND ".join(where)))
            motifs[str(REF_uid)] = defaultdict(list)
            #print REF_uid
            for (REF_SeqName, start, stop) in self.cur.fetchall():
                motifs[str(REF_uid)][REF_SeqName].append((start, stop))
            if not motifs[str(REF_uid)]:
                motifs.pop(str(REF_uid))
        return motifs

    def dumpScan2Dict(self, up_bp=None, strand=None, pValue=None, qValue=None, dist_min=None, dist_max=None, REF_SeqName=None, REF_uid=None):
        '''
        return Motif[uid][REF_SeqNames]=[(start,stop),]
        low ram version which store obj in shelve can not be used because shelve['1'].append(1) will not change the list
        '''
        motifs = defaultdict(dict)
        if not REF_uid:
            REF_uid = self.MEME_ID2uid.values()
        (colums, where) = self.__buildScanColumsWhere(up_bp, strand, pValue, qValue, REF_SeqName, REF_uid)
        SQL = """SELECT REF_uid,REF_SeqName,start,stop FROM Scan %s"""
        #print SQL % "WHERE " + " AND ".join(where)
        self.cur.execute(SQL % "WHERE " + " AND ".join(where))
        for (REF_uid, REF_SeqName, start, stop) in self.cur.fetchall():
            if not motifs[str(REF_uid)]:
                motifs[str(REF_uid)] = defaultdict(list)
            motifs[str(REF_uid)][REF_SeqName].append((start, stop))
        return motifs

    def countMotifs2(self, up_bp=None, strand=None, pValue=None, qValue=None, dist_min=None, dist_max=None, REF_SeqName=None, REF_uid=None):
        '''
        Faster than count by sqlite
        '''
        from collections import Counter
        if not REF_uid:
            REF_uid = self.MEME_ID2uid.values()
        (colums, where) = self.__buildScanColumsWhere(up_bp, strand, pValue, qValue, REF_SeqName, REF_uid)
        SQL = """SELECT DISTINCT REF_uid,REF_SeqName FROM Scan %s"""
        print SQL % "WHERE " + " AND ".join(where)
        self.cur.execute(SQL % "WHERE " + " AND ".join(where))
        return Counter([REF_uid for (REF_uid, REF_SeqName) in self.cur.fetchall()])

    def countMotifs3(self, up_bp=None, strand=None, pValue=None, qValue=None, dist_min=None, dist_max=None, REF_SeqName=None, REF_uid=None):
        '''
        Faster than count by countMotifs2 after run once
        '''
        SeqNames = self.dumpScan2Dict(up_bp, strand, pValue, qValue, dist_min, dist_max, REF_SeqName, REF_uid)
        return dict([(uid, len(SeqNames[uid])) for uid in SeqNames.keys()])

    def cache(self, table='Scan', colums=None, REF_SeqName=None, REF_uid=None, exact=False):
        '''
        将Scanned表读入内存，带来约三倍的速度提升。
        exact only support plantcare and place
        '''
        if REF_SeqName or REF_uid or (exact and table == 'Scan'):
            where = {}
            if REF_SeqName:
                where['REF_SeqName'] = [str(sn_id) for sn_id in REF_SeqName]
            if REF_uid:
                where['REF_MotifSeq'] = [str(uid) for uid in REF_uid]
            if exact and table == 'Scan':
                SQL = "REF_FoundMotif in (select id from FoundMotif where FoundMotif in (select MEME_ID from Motif))"
                where['SQL'] = SQL
        else:
            where = None
        self.select(table, colums, where, cache=True)

    def loadScan(self, colums, where, cache=False):
        try:
            self.cur.execute("ATTACH ':memory:' AS cache")
        except sqlite3.OperationalError, Error:
            print "ERROR:", Error
        SQL = """CREATE TABLE cache.Scan AS
                    SELECT %s
                    FROM Scan
                    %s"""
        (colums, where) = self.__buildScanColumsWhere()
        self.cur.execute(SQL % (colums, where))

if __name__ == '__main__':
    united = UnitiedDB('../RAP3k.db')
    #print united.MEME_ID2uid.values()
    #united.importDB('PlantCARE', 'plantcare.db')
    #united.importDB('PLACE', 'place.db')
    #united.importDB('Jaspar', 'Jaspar.db')
    #united.importDenovo('k.meme', 'MEME')
    #meme_slv = shelve.open(meme_slv)
    #united.exportMotifs(['ATGTACGTGGAGG', 'CACACATGGAA', 'MA0346.1', 'TGACY', 'AAGATTGATTGAG', 'AATCTAATCC', 'ACGT', 'ATAATGGGCCACACTGTGGGGCAT'], meme_slv, 'motif4fimo.meme')
    #united.countMotifs(None,None,0.05,1)
    united.dumpScan2Dict(2000, 1, 0.01, 0.1)
    #print len(united.getMotif_SeqName(2000, 1, 0.05, 0.5).keys())
    # print united.countMotifs2(None,None,None,None)
    # print united.countMotifs3(None,None,None,None)
    #united.countMotifs(1000,"both",0.5,0.5,True,10,250)
    #united.countMotifs(1000,None,0.5,0.5,False)
    #print united.countMotifs(1000,"both",0.01,0.7,True,10,250,range(200,800))
    #united.cache(exact=True)
    #print united.select('cache.Scan', ["count(*)"])
