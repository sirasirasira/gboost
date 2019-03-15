# -*- coding: utf-8 -*-

import sys


argvs = sys.argv  # コマンドライン引数を格納したリストの取得
argc = len(argvs) # 引数の個数
if (argc != 2):   # 引数が足りない場合は、その旨を表示
    print 'Usage: # python %s filename' % argvs[0]
    quit() 


g=-1
dict = {}
for line in open(argvs[1],"r"):
    #print line,
    #print line.split(' ')
    if line[0:5]=="GRAPH":
        g += 1
    if len(line.split(' ')) == 5:
        if line.split(' ')[3][0:3]=="alp":
            i = line.split(' ')[4][:-1]
            #print g,i
            if i.split("=")[1] in dict:
                dict[i.split("=")[1]] += ","+str(g)
            else:
                dict[i.split("=")[1]] = str(g)


                
#print dict
for k,v in sorted(dict.iteritems()):
    print str(k)+"\t"+v+"\n",
