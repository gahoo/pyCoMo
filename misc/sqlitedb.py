#! /usr/bin/env python
#coding=utf-8
import sqlite3


class sqlitedb:
    '''
    连接sqlite数据库，对其操作的类
    方法：
    insert：插入
    createTable：建表
    createIndex：建立索引
    createView：建立视图
    select：查询
    commit：提交事务
    属性：
    dbfile：数据库文件路径
    conn：数据库连接对象
    cur：数据库游标对象
    '''
    dbfile = None
    conn = None
    cur = None

    def __init__(self, dbfile):
        '''
        初始化连接数据库
        '''
        self.dbfile = dbfile
        self.conn = sqlite3.connect(self.dbfile)
        self.cur = self.conn.cursor()
        return

    def insert(self, table, colums, values, commit=False):
        '''
        SQL insert语句的简单封装
        table：表名
        colums：列（列表）
        values：值（列表）
        commit：是否马上提交事务
        '''
        colums = ",".join(colums)
        values = [str(v) for v in values]
        values = '"' + '","'.join(values) + '"'
        SQL = 'INSERT INTO %s (%s) VALUES (%s)' % (table, colums, values)
        #print SQL
        try:
            self.cur.execute(SQL)
            if commit:
                self.conn.commit()
            return True
        except sqlite3.IntegrityError, Error:
            print SQL
            print "ERROR:", Error
            log = open('log.csv', 'a')
            log.write(values + "\n")
            return False

    def createTable(self, table, colums):
        '''
        SQL 建表的简单封装
        tname：表名
        colums：表的结构（字典）
        '''
        #Not safe
        colums = ",\n".join([k + " " + v for (k, v) in colums.items()])
        SQL = 'CREATE TABLE %s (\n%s\n)' % (table, colums)
        #print SQL
        try:
            self.cur.execute(SQL)
        except sqlite3.OperationalError, Error:
            print "ERROR:", Error
        return

    def createIndex(self, index_name, table, colums):
        '''
        建立索引的简单封装
        index_name：索引名
        table：表名
        colums：要索引的列（列表）
        '''
        #Not safe
        colums = ",".join(colums)
        SQL = 'CREATE INDEX %s ON  %s (%s)' % (index_name, table, colums)
        try:
            self.cur.execute(SQL)
        except sqlite3.OperationalError, Error:
            print "ERROR:", Error
        return

    def createView(self, view, tables, colums, equals):
        '''
        建立视图的简单封装
        view：视图名
        tables：表名（列表）
        colums：列名（列表）
        equals：连接的列[(a,b),(b,c)]
        '''
        colums = ",".join(colums)
        tables = ",".join(tables)
        equals = " AND ".join([" = ".join(equal) for equal in equals])
        SQL = 'CREATE VIEW %s AS SELECT DISTINCT %s FROM %s WHERE %s' % \
            (view, colums, tables, equals)
        #print SQL
        try:
            self.cur.execute(SQL)
        except sqlite3.OperationalError, Error:
            print "ERROR:", Error
        return

    def select(self, table, colums=None, where=None, group=None, cache=False, distinct=False, like=False):
        '''
        select语句的简单封装
        table：表名
        colums：列名（列表）
        where：条件（字典）{colum:value}
        cache：是否载入内存
        '''
        if cache:
            try:
                self.cur.execute("ATTACH ':memory:' AS cache")
            except sqlite3.OperationalError, Error:
                print "ERROR:", Error
            cache = 'CREATE TABLE cache.%s AS' % table
        else:
            cache = ''
        if colums:
            colums = ",".join(colums)
        else:
            colums = "*"
        if distinct:
            distinct = "DISTINCT"
        else:
            distinct = ""
        if group:
            group = "GROUP BY %s" % group
        else:
            group = ""
        if where:
            raw_sql = None
            if where.has_key('SQL'):
                raw_sql = where.pop('SQL')
            if like:
                where = [k + " LIKE '%" + v + "%'" for (k, v) in where.items()]
            else:
                where = [k + " IN ('" + "','".join(v) + "')" for (k, v) in where.items()]
            if raw_sql:
                where.append(raw_sql)
            where = " AND ".join(where)
            SQL = '%s SELECT %s %s FROM %s \nWHERE %s %s' % \
                (cache, distinct, colums, table, where, group)
        else:
            SQL = '%s SELECT %s %s FROM %s %s' % (cache, distinct, colums, table, group)
        #print SQL
        try:
            self.cur.execute(SQL)
        except sqlite3.OperationalError, Error:
            print SQL
            raise Exception, Error
        if cache:
            return
        else:
            return self.cur.fetchall()

    def commit(self):
        '''
        提交事务
        '''
        #不commit，insert的结果不会保存
        self.conn.commit()
