#! /usr/bin/env python
#coding=utf-8
from sqlitedb import sqlitedb
from collections import defaultdict
import sqlite3
import os


class JasparDB(sqlitedb):
    '''
    Exact map of Jaspar DataBase
    '''
    tables = [
    """CREATE TABLE MATRIX(
            ID INT NOT NULL,
            COLLECTION VARCHAR (16) DEFAULT '',
            BASE_ID VARCHAR (16)DEFAULT '' NOT NULL,
            VERSION TINYINT DEFAULT 1 NOT NULL,
            NAME VARCHAR (255) DEFAULT '' NOT NULL,
            PRIMARY KEY (ID))
    """,
    """CREATE TABLE MATRIX_DATA(
           ID INT NOT NULL,
           row VARCHAR(1) NOT NULL,
           col TINYINT(3) NOT NULL,
           val float(10,3),
           PRIMARY KEY (ID, row, col))
    """,
    """CREATE TABLE MATRIX_ANNOTATION(
           ID INT NOT NULL,
           TAG VARCHAR(255)DEFAULT '' NOT NULL,
           VAL VARCHAR(255) DEFAULT '',
           PRIMARY KEY (ID, TAG))
    """,
    """CREATE TABLE MATRIX_SPECIES(
           ID INT NOT NULL,
           TAX_ID VARCHAR(255)DEFAULT '' NOT NULL)
    """,
    """ CREATE TABLE MATRIX_PROTEIN(
           ID INT NOT NULL,
           ACC VARCHAR(255)DEFAULT '' NOT NULL)
    """]
    views = [
    """CREATE VIEW Main AS
        SELECT COLLECTION,BASE_ID,VERSION,NAME,group_concat(TAG) AS ANNO_COL,group_concat(VAL) AS ANNO,ACC,TAX_ID
        FROM MATRIX,MATRIX_ANNOTATION,MATRIX_PROTEIN,MATRIX_SPECIES
        WHERE MATRIX.ID=MATRIX_ANNOTATION.ID
            AND MATRIX.ID=MATRIX_PROTEIN.ID
            AND MATRIX.ID=MATRIX_SPECIES.ID
        GROUP BY BASE_ID
    """,
    """CREATE VIEW MATRIX_NON_REDUNDANT AS
        SELECT *
        FROM MATRIX
        GROUP BY BASE_ID
        HAVING COUNT(*)>=1
    """,
    """CREATE VIEW MATRIX_ANNOTATION_NON_REDUNDANT AS
        SELECT *
        FROM MATRIX_ANNOTATION
        WHERE ID IN (SELECT ID
                    FROM MATRIX
                    GROUP BY BASE_ID
                    HAVING COUNT(*)>=1)
    """]

    Names = defaultdict(str)
    Base_ID = defaultdict(str)
    #Annotation = {}
    Class = defaultdict(str)
    Description = defaultdict(str)
    Family = defaultdict(str)
    Species = defaultdict(str)
    Accession = defaultdict(str)
    Pazar_TF = defaultdict(str)

    def __init__(self, dbfile):
        if not os.path.exists(dbfile):
            #若不存在数据库则新建之
            sqlitedb.__init__(self, dbfile)
            self.newDB(dbfile)
        else:
            sqlitedb.__init__(self, dbfile)
            self.loadDB()

    def loadDB(self):
        self.loadDB2dict('MATRIX', "ID", 'NAME', self.Names)
        self.loadDB2dict('MATRIX', "ID", "(BASE_ID || '.' || VERSION)", self.Base_ID)
        self.loadDB2dict('MATRIX_ANNOTATION', 'ID', 'VAL', self.Class, {'TAG': ['class']})
        self.loadDB2dict('MATRIX_ANNOTATION', 'ID', 'VAL', self.Description, {'TAG': ['Description', 'description']})
        self.loadDB2dict('MATRIX_ANNOTATION', 'ID', 'VAL', self.Family, {'TAG': ['family']})
        self.loadDB2dict('MATRIX_ANNOTATION', 'ID', 'VAL', self.Pazar_TF, {'TAG': ['pazar_TF']})
        self.loadDB2dict('MATRIX_SPECIES', 'ID', 'TAX_ID', self.Species)
        self.loadDB2dict('MATRIX_PROTEIN', 'ID', 'ACC', self.Accession)

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

    def getBASE_ID(self, where, redundant=False):
        '''
        where={col:[val]}
        '''
        IDs = set(self.Base_ID.keys())
        for (col, vals) in where.items():
            if col in ["NAME", "COLLECTION"]:
                if redundant:
                    tname = "MATRIX"
                else:
                    tname = "MATRIX_NON_REDUNDANT"
                IDs = IDs & set([ID[0] for ID in self.select(tname, ["ID"], where={col: vals}, distinct=True)])
            else:
                tname = "MATRIX_ANNOTATION"
                IDs = IDs & set([ID[0] for ID in self.select(tname, ["ID"], where={"TAG": [col], "VAL": vals}, distinct=True)])
        return [self.Base_ID[ID] for ID in IDs]

    def newDB(self, dbfile):
        '''
        根据tables和views属性建立表和视图
        '''
        for SQL in self.tables + self.views:
            try:
                self.cur.execute(SQL)
            except sqlite3.OperationalError, Error:
                print "ERROR:", Error

    def importTables(self, sql_tables_path):
        tables = {'MATRIX': ['ID', 'COLLECTION', 'BASE_ID', 'VERSION', 'NAME'],
                'MATRIX_ANNOTATION': ['ID', 'TAG', 'VAL'],
                'MATRIX_DATA': ['ID', 'col', 'row', 'val'],
                'MATRIX_PROTEIN': ['ID', 'ACC'],
                'MATRIX_SPECIES': ['ID', 'TAX_ID']}
        for tname in tables.keys():
            print tname
            for values in self.loadTable(sql_tables_path + tname + '.txt'):
                try:
                    self.insert(tname, tables[tname], values, False)
                except AttributeError, err:
                    print err, values
                except sqlite3.OperationalError, err:
                    print err, values
            self.commit()
        self.loadDB()

    def loadTable(self, filepath):
        table = []
        for line in open(filepath).readlines():
            table.append(line.strip().split('\t'))
        return table


if __name__ == '__main__':
    j = JasparDB('Jaspar.db')
    #j.importTables('/home/gahoo/Project/jaspar/jaspar_CORE/non_redundant/all_species/sql_tables/')
    j.importTables('/home/gahoo/Project/jaspar/all_data/sql_tables/')
    #print j.Description
    print j.getBASE_ID({"COLLECTION": ["CORE"]})
    print j.getBASE_ID({"COLLECTION": ["CORE"], "tax_group": ["vertebrates", "insect"], "family": ["GATA"]})
