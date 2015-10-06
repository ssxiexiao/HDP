# -*- coding:utf-8-*-
import math
import random
import numpy as np
import csv

def sample(p):
    for i in range(len(p)):
        if i > 0:
            p[i] += p[i-1]
    p[len(p)-1] = 1.0
    t = random.random()
    for i in range(len(p)):
        if p[i] > t:
            return i

class UserInfo:

    def __init__(self):
        self.id2Attr = {}
        self.id2Indice = {}
        self.indice2Id = []

    def loadData(self, data):
        self.id2Attr = data
        self.indice2Id = []
        self.indice2Id.append([])
        for id in self.id2Attr:
            self.indice2Id[0].append[id]

    def getAttrById(self, id):
        attr = self.id2Attr[id][:]
        return attr
    
    def getGroup(self):
        return len(self.indice2Id)
    
    def getNJK(self, j, k, clusterinfo):
        njk = 0
        for id in self.indice2Id[j]:
            z = clusterinfo.getLabel(id)
            if z == k:
                njk += 1
        return njk

    def getNJKAttr(self, j, k, clusterinfo):
        result = []
        for id in self.indice2Id[j]:
            z = clusterinfo.getLabel(id)
            if z == k:
                attr = self.getAttrById(id)
                if len(result == 0):
                    for d in range(len(attr)):
                        result.append(0)
                        result[d] = attr[d]
                else:
                    for d in range(len(attr)):
                        result[d] += attr[d]

        return result

    def getNJ(self):
        nj = []
        for j in range(len(self.indice2Id)):
            nj.push(len(self.indice2Id[j]))
        return nj

    def getJ(self, j):
        return self.indice2Id[j]


class ClusterInfo:
    
    def __init__(self):
        self._cluster = {}
        self.maxcluster = 0

    def addToCluster(self, id, z):
        if not z in self._cluster.keys():
            self._cluster[z] = {}
        self._cluster[z][id] = 1

    def removeFromCluster(self, id, z):
        self._cluster[z].pop(id) 
        if len(self._cluster[z].keys()) == 0:
            self._cluster.pop(z)
        return
    
    def getClusters(self):
        return self._cluster.keys()

    def getLabel(self, id):
        for z in self._cluster.keys():
            if id in self._cluster[z].keys():
                return z

class Operator:

    def __init__(self):
        self.stirling = []

    def getLogValue(self, a, b):
        maxv = max(a,b)
        minv = min(a,b)
        if maxv == minv:
            return 1
        result = 0
        for i in range(minv, maxv-1):
            result += math.log(i)
        if a < b:
            return -result
        else:
            return result

    def getMarginal(self, tf, sumtf, tf0, sumtf0):
        result = 0
        result += self.getLogValue(sumtf0, sumtf0+sumtf)
        for i in range(len(tf)):
            result += self.getLogValue(tf0[i]+tf[i], tf0[i])
        return result

    def getNewMarginal(self, tf, sumtf):
        result = 0
        result += self.getLogValue(sumtf0, sumtf0+sumtf)
        for i in range(len(tf)):
            result += self.getLogValue(tf0[i]+tf[i], tf0[i])
        return result

    def samplingT(self, nd, alpha):
        if nd < 1:
            return 0
        t = 1
        for i in range(1, nd):
            if random.random() < (alpha / (nd + alpha)):
                t+=1
        return t

class HDP:

    def __init__(self, clusterinfo, userinfo):
        self.clusterinfo = clusterinfo
        self.userinfo = userinfo

    def sampling(self):
        maxIter = 100
        for i in range(maxIter):
            self.samplingZ()
            self.samplingBeta()
            self.samplingAlpha()
            self.samplingGama()

    def samplingZ(self):
        for j in self.userinfo.getGroup():
            for id in self.userinfo.getJ(j):
                zold = self.clusterinfo.getLabel(id)
                self.clusterinfo.removeFromCluster(id, zold)
                gk = self.clusterinfo.getClusters()
                p = []
                tf = self.userinfo.getAttrById(id)
                sumtf = sum(tf)
                for i in range(len(gk)):
                    njk = self.userinfo.getNJK(j, gk[i], self.clusterinfo)
                    betak = self.beta[gk[i]]
                    operator = Operator()
                    tf0 = self.userinfo.getNJKAttr(j, gk[i], self.clusterinfo)
                    sumtf0 = sum(tf0)
                    logmarg = operator.getMarginal(tf, sumtf, tf0, sumtf0)
                    result = math.log(njk + self.alpha * betak) + logmarg
                    p.append(result)
                logmarg = operator.getNewMarginal(tf, sumtf)
                result = math.log(self.alpha * self.beta_u) + logmarg
                p.append(result)
                s = sampling(p)
                if s == len(p):
                    znew = self.clusterinfo.maxcluster
                    self.clusterinfo.maxcluster += 1
                    nu = np.random.beta(self.gama, 1.0)
                    self.beta[znew] = self.beta_u * (1.0 - nu)
                    self.beta_u = self.beta_u * nu
                else:
                    znew = gk[s]
                self.clusterinfo.addToCluster(id, znew)

        return

    def samplingBeta(self):
        mk = []
        operator = Operator()
        for k in self.clusterinfo.getClusters():
            t = 0
            for j in self.userinfo.getGroup():
                njk = self.userinfo.getNJK(j, k, self.clusterinfo)
                t += operator.samplingT(njk, self.alpha * self.beta[k])
            mk.append(t)
        self.totalTable = sum(mk)
        mk.append(self.gama)
        _beta = np.random.dirichlet(mk, 1)
        _beta = _beta[0]
        count = 0
        for k in self.clusterinfo.getClusters():
            self.beta[k] = _beta[count]
            count+=1
        self.beta_u = _beta[len(_beta)-1]
        return

    def samplingAlpha(self):
        iter = 20
        J = len(self.userinfo.getGroup())
        nj = self.userinfo.getNJ()
        m = self.totalTable
        w = 0
        s = 0
        t = 0
        al = self.alpha
        for i in range(iter):
            a = self.alpha_a + m
            b = self.alpha_b
            #sampling wj and sj
            for j in range(J):
                w = np.random.beta(al + 1, nj[j])
                t = nj[j] / al
                if random.random() < t / (t + 1):
                    s = 1
                else:
                    s = 0
                a -= s
                b -= math.log(w)
            
            #sampling alpha
            al = np.random.gamma(a, 1 / b)
        
        self.alpha = al
        return

    def samplingGama(self):
        iter = 20
        ga = 0 
        eta = 0
        m = self.totalTable        
        K = len(self.clusterinfo.getClusters())
        kai = 0
        beta = 0
        ga = self.gama
        for i in range(iter):
            eta = np.random.beta(ga + 1, m)
            kai = self.gama_a + K - 1
            beta = self.gama_b - math.log(eta)
            if random.random() < (kai / (kai + m * beta)):
                kai += 1

            ga = np.random.gamma(kai, 1/beta)       
        self.gama = ga
        return
 
#LOAD DATA
rawdata = {}
reader = csv.reader(file('E:\\DATA\\sever230final(fitered).csv', 'r'))
count = 0
title = []
for line in reader:
    if count == 0:
        title = line
        count += 1
    else:
        time = line[0]
        id = line[1]
        line = line[2:]
        for i in range(len(line)):
            line[i] = float(line[i])
        if time not in rawdata.keys():
            rawdata[time] = {}
        rawdata[time][id] = line[:]
idlist = []
count = 0
timeIndex = rawdata.keys()
timeIndex.sort()
timeIndex = timeIndex[:8]
print timeIndex
for time in timeIndex:
    if count == 0:
        count += 1
        idlist = rawdata[time].keys()
    else:
        newidlist = []
        for id in idlist:
            if id in rawdata[time].keys():
                newidlist.append(id)
        idlist = newidlist
print len(idlist) 
data = {}
for time in timeIndex:
    data[time] = {}
for time in timeIndex:
    for id in idlist:
        data[time][id] = rawdata[time][id][:]


userinfo = UserInfo()
clusterinfo = ClusterInfo()
hdp = HDP(clusterinfo, userinfo)
count = 0
for i in timeIndex:
    #CHANGE ATTRIBUTE
    hdp.userinfo.loadData(data[i])
    if count == 0:
        count += 1
        #INITIAL CLUSTERING
        INITIAL_K = 5
        p = []
        for l in range(INITIAL_K):
            p.append(1.0/INITIAL_K)
        for id in data[i]:
            z = sample(p) + 1
            hdp.clusterinfo.addToCluster(id, z)
        hdp.clusterinfo.maxcluster = INITIAL_K + 1

    #RE-CLUSTERING
    hdp.sampling()

    #OUTPUT
