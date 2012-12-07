#! /usr/bin/env python
#coding=utf-8
from optparse import OptionParser
from collections import defaultdict
import sys
import os
import time
import operator
import xmlrpclib


def loadEdges(edgefile):
    '''
    加载Edge文件
    返回字典Edges
    '''
    edges = open(edgefile, 'r')
    Edges = defaultdict(list)
    for edge in edges.readlines():
        (motifA, motifB, lst_cnt, lst_size, db_cnt, db_size, enrich_ratio, p, adjp) = edge.strip().split('\t')
        Edges["motifA"].append(motifA)
        Edges["motifB"].append(motifB)
        Edges["lst_cnts"].append(int(lst_cnt))
        Edges["lst_sizes"].append(int(lst_size))
        Edges["db_cnts"].append(int(db_cnt))
        Edges["db_sizes"].append(int(db_size))
        Edges["enrich_ratios"].append(float(enrich_ratio))
        Edges["Ps"].append(float(p))
        Edges["adjPs"].append(float(adjp))
    return Edges


def loadNodes(nodefile):
    '''
    加载Node文件
    返回字典Nodes
    '''
    nodes = open(nodefile, 'r')
    Nodes = defaultdict(list)
    for node in nodes.readlines():
        (motif, desc, motifseq, lst_cnt, lst_size, db_cnt, db_size, enrich_ratio, p, adjp) = node.strip().split('\t')
        Nodes['motifs'].append(motif)
        Nodes['descs'].append(desc)
        Nodes['motifSeq'].append(motifseq)
        Nodes['lst_cnts'].append(int(lst_cnt))
        Nodes['lst_sizes'].append(int(lst_size))
        Nodes['db_cnts'].append(int(db_cnt))
        Nodes['db_sizes'].append(int(db_size))
        Nodes['enrich_ratios'].append(float(enrich_ratio))
        Nodes['Ps'].append(float(p))
        Nodes['adjPs'].append(float(adjp))
    return Nodes


def createNetwork(server, Nodes, Edges, network_name, layout, nodePvalue, edgePvalue):
    '''
    根据Node和Edge构建网络
    '''
    networkid = server.Cytoscape.createNetwork(network_name)
    importNodes(server, networkid, Nodes)
    importEdges(server, networkid, Edges)
    server.Cytoscape.setDefaultBackgroundColor('default', "#FFFFFF")
    server.Cytoscape.performLayout(networkid, layout)
    setMappers(server, Nodes, Edges)
    saveNet(server, networkid, network_name)

    if nodePvalue:
        rich_nodes = subNodes(server, networkid, network_name, Nodes, layout,
                    attribute='adjPs', op='<=', bound=nodePvalue)
        saveNet(server, rich_nodes, "%s_NodePs%s" % (network_name, nodePvalue))
    if edgePvalue:
        sig_edges = subEdges(server, networkid, network_name, Nodes, Edges, layout,
                    attribute='adjPs', op='<=', bound=edgePvalue)
        saveNet(server, sig_edges, "%s_EdgePs%s" % (network_name, edgePvalue))
        server.Cytoscape.destroyNetwork(sig_edges)
    if nodePvalue and edgePvalue:
        rich_sig = subEdges(server, rich_nodes, "%s_Node%sEdge%s" % (network_name, nodePvalue, edgePvalue),
                    Nodes, Edges, layout, attribute='adjPs', op='<=', bound=edgePvalue)
        saveNet(server, rich_sig, "%s_Node%sEdgePs%s" % (network_name, nodePvalue, edgePvalue))
        server.Cytoscape.destroyNetwork(rich_sig)
    if nodePvalue:
        server.Cytoscape.destroyNetwork(rich_nodes)
    saveCys(server, network_name)
    server.Cytoscape.destroyNetwork(networkid)


def importNodes(server, networkid, Nodes):
    '''
    导入节点属性信息
    '''
    attributes = ['descs', 'motifSeq', 'lst_cnts', 'lst_sizes', 'db_cnts', 'db_sizes', 'enrich_ratios', 'Ps', 'adjPs']
    attr_types = ['String', 'String', 'Integer', 'Integer', 'Integer', 'Integer', 'Double', 'Double', 'Double']
    addNodes = {'String': server.Cytoscape.addStringNodeAttributes,
                'Integer': server.Cytoscape.addIntegerNodeAttributes,
                'Double': server.Cytoscape.addDoubleNodeAttributes}
    for attribute, attr_types in zip(attributes, attr_types):
        addNodes[attr_types](attribute, Nodes['motifs'], Nodes[attribute])
    #节点类型，是否富集，用于设置节点形状
    nodes_types = map(setNodeType, Nodes['adjPs'])
    server.Cytoscape.addStringNodeAttributes('nodes_types', Nodes['motifs'], nodes_types)


def setNodeType(x):
    '''
    set node type!
    '''
    if x <= 0.05:
        if x <= 0.001:
            return 'sig_en'
        else:
            return 'enrich'
    else:
        return 'not'


def importEdges(server, networkid, Edges):
    '''
    导入边属性信息
    '''
    server.Cytoscape.createEdges(Edges['motifA'], Edges['motifB'])
    edgenames = ["%s (directed) %s" % (Edges['motifA'][i], Edges['motifB'][i]) for i in range(len(Edges['motifB']))]
    attributes = ['lst_cnts', 'lst_sizes', 'db_cnts', 'db_sizes', 'enrich_ratios', 'Ps', 'adjPs']
    attr_types = ['Integer', 'Integer', 'Integer', 'Integer', 'Double', 'Double', 'Double']
    addEdges = {'Integer': server.Cytoscape.addIntegerEdgeAttributes,
                'Double': server.Cytoscape.addDoubleEdgeAttributes}
    for attribute, attr_types in zip(attributes, attr_types):
        addEdges[attr_types](attribute, edgenames, Edges[attribute])
    # server.Cytoscape.addIntegerEdgeAttributes('lst_cnts', edgenames, Edges['lst_cnts'])
    # server.Cytoscape.addIntegerEdgeAttributes('lst_sizes', edgenames, Edges['lst_sizes'])
    # server.Cytoscape.addIntegerEdgeAttributes('db_cnts', edgenames, Edges['db_cnts'])
    # server.Cytoscape.addIntegerEdgeAttributes('db_sizes', edgenames, Edges['db_sizes'])
    # server.Cytoscape.addDoubleEdgeAttributes('enrich_ratios', edgenames, Edges['enrich_ratios'])
    # server.Cytoscape.addDoubleEdgeAttributes('Ps', edgenames, Edges['Ps'])
    # server.Cytoscape.addDoubleEdgeAttributes('adjPs', edgenames, Edges['adjPs'])
    edges_types = map(lambda x: 'highly_sign' if x <= 0.01 else 'sign', Edges['adjPs'])
    server.Cytoscape.addStringEdgeAttributes('edges_types', edgenames, edges_types)


def setMappers(server, Nodes, Edges):
    '''
    将属性映射至图像中，包括颜色，大小，形状等
    节点：三角形：富集、大小：出现次数、颜色：富集率
    边：实线：显著相关P<0.01、宽度：P显著程度、颜色：共现次数
    '''
    #Nodes Attributes
    Mapper = {'Discrete': server.Cytoscape.createDiscreteMapper,
            'Continuous': server.Cytoscape.createContinuousMapper}

    m_type = {'Discrete': 'default',
            'Continuous': 'range'}
    #this should place in an config file
    vizer = [{ #节点大小与lst_cnts，即列表中模体出现次数成正比
            'attr': 'nodes_types',
            'map': 'Node Shape',
            'type': 'Discrete',
            'default': "ellipse",
            #富集的节点为三角形，否则为原型
            'values': {'sig_en': 'diamond', 'enrich': "triangle", 'not': "ellipse"}},
            { #节点大小与lst_cnts，即列表中模体出现次数成正比
            'attr': 'lst_cnts',
            'map': 'Node Size',
            'type': 'Continuous',
            'range': [float(min(Nodes['lst_cnts'])), float(max(Nodes['lst_cnts']))],
            'values': [5.0, 10.0, 50.0, 60.0]},
            { #节点颜色深浅与enrich_ratios，即富集率成正比
            'attr': "enrich_ratios",
            'map': "Node Color",
            'type': 'Continuous',
            'range': [float(min(Nodes['enrich_ratios'])), 1.0, float(max(Nodes['enrich_ratios']))],
            'values': ['#FFFFFF', '#FFFFFF', '#FFFFFF', '#3030FF', '#3030FF']},
            { #显著相关的边为实线，一般相关为虚线
            'attr': 'edges_types',
            'map': 'Edge Line Style',
            'type': 'Discrete',
            'default': 'EQUAL_DASH',
            'values': {'highly_sign': "SOLID", 'sign': "EQUAL_DASH"}},
            { #边的宽度与lst_cnts，即列表中模体出现次数成正比
            'attr': 'lst_cnts',
            'map': 'Edge Line Width',
            'type': 'Continuous',
            'range': [float(min(Edges['lst_cnts'])), float(max(Edges['lst_cnts']))],
            'values': [0.5, 1.0, 7.0, 8.0]},
            { #边颜色深浅与enrich_ratios，即富集率成正比
            'attr': 'enrich_ratios',
            'map': 'Edge Color',
            'type': 'Continuous',
            'range': [float(min(Edges['enrich_ratios'])), float(max(Edges['enrich_ratios']))],
            'values': ['#F0F0FF', '#F0F0FF', '#0000FF', '#0000FF']}]

    for viz in vizer:
        Mapper[viz['type']]('default', viz['attr'], viz['map'],
                viz[m_type[viz['type']]], viz['values'])

    # server.Cytoscape.createDiscreteMapper('default', 'nodes_types', 'Node Shape', \
    #                                       "ellipse", \
    #                                       {'sig_en': 'diamond', 'enrich': "triangle", 'not': "ellipse"})
    # #节点大小与lst_cnts，即列表中模体出现次数成正比
    # server.Cytoscape.createContinuousMapper('default', 'lst_cnts', 'Node Size', \
    #                                         [float(min(Nodes['lst_cnts'])), float(max(Nodes['lst_cnts']))], \
    #                                         [5.0, 10.0, 50.0, 60.0])
    # #节点颜色深浅与enrich_ratios，即富集率成正比
    # server.Cytoscape.createContinuousMapper('default', 'enrich_ratios', 'Node Color', \
    #                                         [float(min(Nodes['enrich_ratios'])), 1.0, float(max(Nodes['enrich_ratios']))], \
    #                                         ['#FFFFFF', '#FFFFFF', '#FFFFFF', '#3030FF', '#3030FF'])
    # #Edges Attributes
    # #显著相关的边为实线，一般相关为虚线
    # server.Cytoscape.createDiscreteMapper('default', 'edges_types', 'Edge Line Style', \
    #                                       "EQUAL_DASH", \
    #                                       {'highly_sign': "SOLID", 'sign': "EQUAL_DASH"})
    # #边的宽度与P值成反比，P越小线越粗
    # # server.Cytoscape.createContinuousMapper('default', 'adjPs', 'Edge Line Width', \
    # #                                         [float(min(Edges['adjPs'])), float(max(Edges['adjPs']))], \
    # #                                         [8.0, 7.0, 1.0, 0.5])
    # server.Cytoscape.createContinuousMapper('default', 'lst_cnts', 'Edge Line Width', \
    #                                         [float(min(Edges['lst_cnts'])), float(max(Edges['lst_cnts']))], \
    #                                         [0.5, 1.0, 7.0, 8.0])
    # #边的颜色和共现次数co_nums成正比，颜色越深共现次数越多
    # server.Cytoscape.createContinuousMapper('default', 'enrich_ratios', 'Edge Color', \
    #                                         [float(min(Edges['enrich_ratios'])), float(max(Edges['enrich_ratios']))], \
    #                                         ['#F0F0FF', '#F0F0FF', '#0000FF', '#0000FF'])


def saveNet(server, networkid, network_name):
    '''
    保存网络到文件，导出为图片
    '''
    #子网络无节点
    if not server.Cytoscape.countNodes(networkid):
        return
    #计算合适的scale
    node_cnt = server.Cytoscape.countNodes(networkid)
    scale = 2 + node_cnt / 15.0
    #避免motif名称出边界
    zoom = 0.8 * server.Cytoscape.getZoom(networkid)
    server.Cytoscape.setZoom(networkid, zoom)
    netfile = "%s/%s_Cyto" % (os.getcwd(), network_name)
    print network_name, scale, zoom#,netfile
    #server.Cytoscape.saveNetwork(networkid,"%s.gml" % netfile)
    server.Cytoscape.executeCommand('network', 'export', {'file': "%s.xgmml" % netfile, 'type': "xgmml"})
    server.Cytoscape.exportView(networkid, "%s.png" % netfile, "png", scale)


def subNodes(server, networkid, network_name, Nodes, layout, attribute='adjPs', op='<=', bound=0.05):
    '''
    指定属性在制定范围内的节点子网络
    默认为P小于0.05
    '''
    oper = {'>=': operator.ge,
            '<=': operator.le,
            '>': operator.gt,
            '<': operator.lt,
            '==': operator.eq,
            '!=': operator.ne}

    enriched = [n for n, v in zip(Nodes['motifs'], Nodes[attribute]) if oper[op](v, bound)]
    #去掉那些虽然enrich但不在网络中的节点
    all_nodes = server.Cytoscape.getNodes(networkid)
    enriched = list(set(enriched) & set(all_nodes))
    server.Cytoscape.selectNodes(networkid, enriched)
    sub_networkid = server.Cytoscape.createNetworkFromSelection(networkid,
                    "%s Nodes %s%s%s" % (network_name, attribute, op, bound))
    server.Cytoscape.performLayout(sub_networkid, layout)
    #if server.Cytoscape.countNodes(sub_networkid):
    #    saveNet(server,sub_networkid,"%s_Node%s%s" % (network_name,attribute,str(bound)))
    #server.Cytoscape.destroyNetwork(sub_networkid)
    return sub_networkid


def subEdges(server, networkid, network_name, Nodes, Edges, layout, attribute='adjPs', op='<=', bound=0.05):
    '''
    指定属性在制定范围内的边子网络
    默认为P小于0.05
    '''
    oper = {'>=': operator.ge,
            '<=': operator.le,
            '>': operator.gt,
            '<': operator.lt,
            '==': operator.eq,
            '!=': operator.ne}

    rev_op = {'>=': '<', '<=': '>', '<': '>=', '>': '<=', '==': '!=', '!=': '=='}
    #由于不能直接选择边创建网络，只能迂回
    #先根据显著的边所连的所有节点创建一个临时网络
    significant_edges = ["%s (directed) %s" % (Edges['motifA'][i], Edges['motifB'][i]) \
                        for i in range(len(Edges[attribute])) \
                        if oper[op](Edges[attribute][i], bound)]
    server.Cytoscape.selectEdges(networkid, significant_edges)
    selectNodesBYEdges(server, networkid, significant_edges)
    tmp_networkid = server.Cytoscape.createNetworkFromSelection(networkid, "tmp")
    #再删除临时网络中不显著的边
    sub_edges = server.Cytoscape.getEdges(tmp_networkid)
    sub_edges_p = server.Cytoscape.getEdgesAttributes(attribute, sub_edges)
    not_sign_edges = [sub_edges[i] for i in range(len(sub_edges_p)) if oper[rev_op[op]](sub_edges_p[i], bound)]
    #not_sign_edges=map(lambda x: True if x<=0.001 else False, sub_edges_p)
    for not_sign_edge in not_sign_edges:
        server.Cytoscape.removeEdge(tmp_networkid, not_sign_edge)
    #最后选择剩余所有边所连接的节点，据此创建网络
    remaining_edges = server.Cytoscape.getEdges(tmp_networkid)
    selectNodesBYEdges(server, tmp_networkid, remaining_edges)
    sub_networkid = server.Cytoscape.createNetworkFromSelection(tmp_networkid,
                    "%s Edges %s%s%s" % (network_name, attribute, op, str(bound)))
    server.Cytoscape.performLayout(sub_networkid, layout)
    #if server.Cytoscape.countNodes(sub_networkid):
    #    saveNet(server,sub_networkid,"%s_Edge%s%s" % (network_name,attribute,str(bound)))
    server.Cytoscape.destroyNetwork(tmp_networkid)
    #server.Cytoscape.destroyNetwork(sub_networkid)
    return sub_networkid


def selectNodesBYEdges(server, networkid, edges):
    '''
    根据所选的边选择对应的节点
    '''
    targer_nodes = server.Cytoscape.getEdgeTargetNodes(networkid, edges)
    source_nodes = server.Cytoscape.getEdgeSourceNodes(networkid, edges)
    server.Cytoscape.selectNodes(networkid, targer_nodes)
    server.Cytoscape.selectNodes(networkid, source_nodes)


def saveCys(server, network_name):
    cysfile = "%s/%s_Cyto.cys" % (os.getcwd(), network_name)
    server.Cytoscape.saveSessionAsCys(cysfile)


def getOptions():
    parser = OptionParser(usage="Usage: drawCyto [options] -e <edgefile> -n <nodefile>")
    parser.add_option("-e", "--edgefile", dest="edgefile",
                help="Edges file generated by Stats", metavar="FILE")
    parser.add_option("-n", "--nodefile", dest="nodefile",
                help="Nodes file generated by Stats", metavar="FILE")
    parser.add_option("-m", "--network_name", dest="network_name",
                default="",
                type="string",
                help="The network name", metavar="FILE")
    parser.add_option("-p", "--nodepvalue", dest="nodePvalue",
                default=None,
                type="float",
                help="Pvalue Threshold to output Nodes",
                metavar="FLOAT")
    parser.add_option("-q", "--edgepvalue", dest="edgePvalue",
                default=None,
                type="float",
                help="Pvalue Threshold to output Edges",
                metavar="FLOAT")
    parser.add_option("-l", "--layout", dest="layout",
                default="Kamada-Kawai",
                help="Network Layouts:jgraph-circle,attribute-circle,jgraph-annealing,jgraph-radial-tree"
                ",Kamada-Kawai-Noweight,Fruchterman-Rheingold,Kamada-Kawai,jgraph-gem,hierarchical,circular"
                ",isom,jgraph-moen,jgraph-sugiyama,attributes-layout,grid,jgraph-tree,force-directed,"
                "degree-circle,jgraph-spring",
                metavar="STRING")
    parser.add_option("-s", "--server", dest="server",
                default="http://localhost:9000",
                help="Cytoscape Server Address",
                metavar="STRING")
    parser.add_option("-r", "--enrich_ratio", dest="enrich_ratio",
                default=None,
                type="float",
                help="The enrich ratio threshold",
                metavar="FLOAT")
    return parser


if __name__ == '__main__':
    parser = getOptions()
    (options, args) = parser.parse_args()

    server = xmlrpclib.ServerProxy(options.server)

    try:
        server.Cytoscape.test()
    except:
        print "CytoscapeRPC is not running"
        sys.exit(2)

    if not options.edgefile or not options.nodefile:
        parser.print_help()
        sys.exit(2)

    if not os.path.exists(options.edgefile):
        print "Edges file %s not exists." % options.edgefile
        sys.exit(2)
    if not os.path.exists(options.nodefile):
        print "Nodes file %s not exists." % options.nodefile
        sys.exit(2)

    if not options.network_name:
        options.network_name = ".".join(os.path.basename(options.edgefile).split('.')[:-1])

    #print network_name,layout,enrichratioThresh,pThresh
    timestamp = time.time()
    createNetwork(server, \
                loadNodes(options.nodefile), \
                loadEdges(options.edgefile), \
                options.network_name, \
                options.layout, \
                options.nodePvalue, \
                options.edgePvalue)
    print time.time() - timestamp
