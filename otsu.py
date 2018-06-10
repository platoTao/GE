from random import randint
from random import uniform
import random
import numpy as np
from PIL import Image
from skimage import data,filters
import matplotlib.pyplot as plt
image = data.camera()
scored={}#存储分数，不用重复计算

def creat_chromosome(n,chromlength):
    '''生成n条染色体'''
    chromosomes=[]
    for i1 in range(n):
        g=[]
        for j in range(chromlength):
            x1=randint(0,1)
            g.append(x1)
        chromosomes.append(g)
    return chromosomes

def vairation(chromosome,pa):
    '''突变'''
    for i2 in range(len(chromosome)):
        j2=uniform(0,1)
        if j2<=pa:
            chromosome[i2]=1-chromosome[i2]
        else:
            pass
    return chromosome

def crossove(a,b,pc):
    '''两条染色体重组,采用一点交叉'''
    x2=uniform(0,1)
    if x2<=pc:
        x3=randint(0,5)
        a[x3:x3+2],b[x3:x3+2]=b[x3:x3+2],a[x3:x3+2]
    return a,b

import math
def decodechrom(chromosomes):
    '''将二进制转化为十进制'''
    decodes=[]
    for k1 in chromosomes:
        x4=0
        for i3 in range(len(k1)-1,-1,-1):
            x4+=math.pow(2,len(k1)-1-i3)*k1[i3]
        decodes.append(x4)
    return decodes

def d(best_chromosome):
    '''将二进制转化为十进制'''
    x4=0
    for r1 in range(len(best_chromosome)-1,-1,-1):
        x4+=math.pow(2,len(best_chromosome)-1-r1)*best_chromosome[r1]
    print(x4)

def transform(image):
    image_1=[]
    for value1 in image:
        for value2 in value1:
            image_1.append(value2)
    image_1.sort()
    return image_1
def cal_scores(chromosomes,image_1,decodes,scored):
    '''计算每条染色体的得分'''
    length=len(image_1)
    mean=np.mean(image_1)
    scores=[]
    for i4 in decodes:
        if i4 in scored:
            scores.append(scored[i4])
        else:
            t1=0
            while i4 >=image_1[t1] and t1<length-1:
                t1+=1
            keys1=image_1[0:t1]
            keys2=image_1[t1:]
            w1=t1/length
            w2=1-w1
            if len(keys1):
                mean_1=np.mean(keys1)
            else:
                mean_1=0
            if len(keys2):
                mean_2=np.mean(keys2)
            else:
                mean_2=0
            score=w1*(mean_1-mean)**2+w2*(mean_2-mean)**2
            scored[i4]=score
            scores.append(score)
    return scores


def cal_fitvals(scores):
    '''计算每条染色体的适应度'''
    total_score=sum(scores)
    return [i5/total_score for i5 in scores]

def cal_fitvalsum(fitvals):
    '''计算累计概率'''
    fitvalsum=[]
    for i6 in range(len(fitvals)):
        t2=0
        j6=0
        while(j6<=i6):
            t2+=fitvals[j6]
            j6+=1
        fitvalsum.append(t2)
    return fitvalsum

def selection(n,fitvalsum,chromosomes,scores,high_score,best_chromosome,enough):
    '''选择——复制'''
    m=max_fitval(scores)
    chromosomes_1=chromosomes[:]
    parents=[]
    if scores[m] >high_score:
        high_score=scores[m]
        best_chromosome=chromosomes_1[m]
        enough=0
        #print("改变后最佳",best_chromosome,high_score)
    for i7 in range(n-1):
        t3=uniform(0,1)
        j7=0
        while fitvalsum[j7]<=t3:
            j7+=1
        parents.append(chromosomes[j7])
    return parents,high_score,best_chromosome,enough


def generate(n,parents):
    '''生成下一代'''
    chromosomes=[]
    parents=random.sample(parents,len(parents))

    for i8 in range(0,n-1,2):
        if i8==n-2:
            a=parents[i8]
            a=vairation(a,pa)#突变
            chromosomes.append(a)
        else:
            a=parents[i8]
            b=parents[i8+1]
            a,b=crossove(a,b,pc)#交叉
            vairation(a,pa)#突变
            chromosomes.append(a)
            vairation(b,pa)#突变
            chromosomes.append(b)
    return chromosomes

def max_fitval(scores):
    '''得分最高的'''
    t=-1
    for i9,s in enumerate(scores):
        if s > t:
            m=i9
            t=s
        else:
            pass
    return m

def ge(n,chromlength,pa,pc,times,image,scored):
    chromosomes=creat_chromosome(n,chromlength)#随机生成第一代
    image_1=transform(image)
    decodes=decodechrom(chromosomes)
    Scores=cal_scores(chromosomes,image_1,decodes,scored)
    M=max_fitval(Scores)
    high_score=Scores[M]
    best_chromosome=chromosomes[M]
    print('结果如下所示：\n')
    enough=0
    for i in range(times):
        decodes=decodechrom(chromosomes)
        scores=cal_scores(chromosomes,image_1,decodes,scored)#计算得分
        max_fitvalue=decodes[max_fitval(scores)]
        contents='第{0}代适应度最大时灰度为{1}\n'.format(i,max_fitvalue)
        with open('result.txt','a') as f:
            f.write(contents)
        print(contents)
        fitvals=cal_fitvals(scores)#计算适应度
        fitvalsum=cal_fitvalsum(fitvals)#计算累计适应度
        parents,high_score,best_chromosome,enough=selection(n,fitvalsum,chromosomes,scores,high_score,best_chromosome,enough)#选择父代
        d(best_chromosome)
        best_chromosome_vice=best_chromosome[:]#不懂为什么best_cheomosome经过generate之后会改变，这里保存一下值
        chromosomes=generate(n,parents)#生成子代
        chromosomes.append(best_chromosome_vice)
        best_chromosome=best_chromosome_vice[:]#校正best_chromosome的值
        enough+=1
        print(enough)
        if enough==10:
            return max_fitvalue
    return max_fitvalue


n=8#族群大小
chromlength=8#染色体长度
times=50#迭代次数
pa=0.01#突变几率
pc=0.6#交叉概率



thresh=ge(n,chromlength,pa,pc,times,image,scored)
