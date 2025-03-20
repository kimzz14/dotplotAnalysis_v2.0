from multiprocessing import Pool, shared_memory
from dataclasses import dataclass, field
from lib.KJH_SVG.KJH_SVG import element
from itertools import groupby
from scipy.stats import chi2
import math
import numpy as np
###########################################################################################
def chi_squared_test(variance, number):
    sigma0_squared = 200000 * 200000

    chi_squared_stat = (number - 1) * variance / sigma0_squared

    alpha = 0.05

    df = number - 1
    critical_value = chi2.ppf(1 - alpha, df)

    return chi_squared_stat < critical_value

def median_with_outliers(array):
    median = np.median(array)
    variance = np.var(array)

    return median, variance

def median_without_outliers(array):
    sorted_array = np.sort(array)
    # 사분위수 계산
    Q1 = np.percentile(sorted_array, 25)
    Q3 = np.percentile(sorted_array, 75)
    IQR = Q3 - Q1
    
    # 아웃라이어 범위 설정
    lower_bound = Q1 - 1.5 * IQR
    upper_bound = Q3 + 1.5 * IQR

    # 아웃라이어 제거
    filtered_array = sorted_array[(sorted_array >= lower_bound) & (sorted_array <= upper_bound)]

    # 중앙값 계산
    median = np.median(filtered_array)
    variance = np.var(filtered_array)

    #print(len(array), len(filtered_array))

    #print("Filtered Median:", median)

    return median, variance
###########################################################################################
class FAIDX_READER:
    def __init__(self, fileName):
        self.fileName = fileName
        self.seqLen_DICT = {}
        self.seqName_LIST = []
        self.totalN = 0

        fin = open(self.fileName)
        for line in fin:
            seqName, seqLen = line.rstrip('\n').split('\t')[0:2]
            self.seqName_LIST += [seqName]
            seqLen = int(seqLen)
            self.seqLen_DICT[seqName] = seqLen
            self.totalN += seqLen
        fin.close()
    
    def get(self, seqName):
        return self.seqLen_DICT[seqName]

@dataclass
class Match:
    qname: str
    rname: str
    strand: str
    intercept: int
    number: int
    variance: float
    std_deviation: float
    qsPos: int
    qePos: int
    rsPos: int
    rePos: int
    coverage: float

class CONTIG:
    def __init__(self, qname):
        self.qname = qname

        self.block_DICT = {}
        #self.blockN_DICT = {}
        #self.block_LIST = []

        self.bestMatch = None

    def destory(self):
        del self.block_DICT
    
    def add(self, rname, strand, qpos, rpos):
        if strand == 0:
            intercept = -qpos + rpos
        elif strand == 1:
            intercept = qpos + rpos

        if not rname in self.block_DICT:
            self.block_DICT[rname] = {0:{}, 1:{}}
        
        if not intercept in self.block_DICT[rname][strand]:
            self.block_DICT[rname][strand][intercept]  = [[qpos, rpos, False]]
        else:
            self.block_DICT[rname][strand][intercept] += [[qpos, rpos, False]]

    def removeRepeat(self, maxRepeat):
        _qblock_DICT = {}
        _rblock_DICT = {}
        for rname, sub1_DICT in self.block_DICT.items():
            _qblock_DICT[rname] = {}
            _rblock_DICT[rname] = {}
            for strand, sub2_DICT in sub1_DICT.items():
                for intercept, block_LIST in sub2_DICT.items():
                    for block in block_LIST:
                        #_qblock_DICT
                        if not block[0] in _qblock_DICT[rname]:
                            _qblock_DICT[rname][block[0]]  = [block]
                        else:
                            _qblock_DICT[rname][block[0]] += [block]

                        #_rblock_DICT
                        if not block[1] in _rblock_DICT[rname]:
                            _rblock_DICT[rname][block[1]]  = [block]
                        else:
                            _rblock_DICT[rname][block[1]] += [block]

        for rname, sub1_DICT in _qblock_DICT.items():
            for qpos, block_LIST in sub1_DICT.items():
                if len(block_LIST) < maxRepeat: continue
                for block in block_LIST:
                    block[2] = True

        for rname, sub1_DICT in _rblock_DICT.items():
            for rpos, block_LIST in sub1_DICT.items():
                if len(block_LIST) < maxRepeat: continue
                for block in block_LIST:
                    block[2] = True

        del _qblock_DICT
        del _rblock_DICT

        _contig = CONTIG(self.qname)
        for rname, sub1_DICT in self.block_DICT.items():
            for strand, sub2_DICT in sub1_DICT.items():
                for intercept, block_LIST in sub2_DICT.items():
                    for block in block_LIST:
                        if block[2] == True: continue
                        _contig.add(rname, strand, block[0], block[1])

        return _contig
    
    def extract(self, minBlockN):
        _contig = CONTIG(self.qname)
        for rname, sub1_DICT in self.block_DICT.items():
            for strand, sub2_DICT in sub1_DICT.items():
                for intercept, block_LIST in sub2_DICT.items():
                    if len(block_LIST) < minBlockN: continue
                    for block in block_LIST:
                        _contig.add(rname, strand, block[0], block[1])

        return _contig

    def calc_best(self, coverage):
        bestMatch = None

        rname_LIST = list(self.block_DICT.keys())
        for rname in rname_LIST:
            match = self.calc(rname, 1)

            if bestMatch == None or bestMatch.number < match.number:
                bestMatch = match

        for rname in rname_LIST:
            if bestMatch.rname == rname: continue

            del self.block_DICT[rname]

        self.bestMatch = bestMatch
        return self.bestMatch

    def calc(self, rname, coverage):
        target_LIST = []
        for strand, sub1_DICT in self.block_DICT[rname].items():
            for intercept, block_LIST in sub1_DICT.items():
                target_LIST += [(strand, intercept, block_LIST)]

        target_LIST.sort(key=lambda x : len(x[2]), reverse=True)

        totalN = sum([len(target[2]) for target in target_LIST])
        subN = 0

        qsPos = 1000000000
        qePos = 0
        rsPos = 1000000000
        rePos = 0

        interceptP_LIST = []
        interceptM_LIST = []
        for [strand, intercept, block_LIST] in target_LIST:
            blockN = len(block_LIST)
            subN += blockN

            if strand == 0:
                interceptP_LIST += [intercept]* blockN

            if strand == 1:
                interceptM_LIST += [intercept]* blockN

            if subN/totalN > coverage: break
        
        #print(rname, len(interceptP_LIST), len(interceptM_LIST))
        if len(interceptP_LIST) > len(interceptM_LIST):
            STRAND = 0
            meanIntercept, variance = median_with_outliers(interceptP_LIST)
            std_deviation = variance ** 0.5
        else:
            STRAND = 1
            meanIntercept, variance = median_with_outliers(interceptM_LIST)
            std_deviation = variance ** 0.5

        
        #calculate boundery
        for [strand, intercept, block_LIST] in target_LIST:
            if STRAND != strand: continue
            if abs(intercept - meanIntercept) >  std_deviation * 3: continue
            for block in block_LIST:
                block[2] = True

                qsPos = min(qsPos, block[0])
                qePos = max(qePos, block[0])

                rsPos = min(rsPos, block[1])
                rePos = max(rePos, block[1])

        number = len(interceptP_LIST) + len(interceptM_LIST)

        return Match(self.qname, rname, STRAND, meanIntercept, number, variance, std_deviation, qsPos, qePos, rsPos, rePos, float(qePos - qsPos) / (rePos - rsPos))
            
class Engine:
    def __init__(self, prefix, batchN, batchIDX):
        self.batchN = batchN
        self.batchIDX = batchIDX

        self.fin = open(prefix + '.sam')
    
    def close(self):
        self.fin.close()
    
    def run(self):
        for line in self.fin:
            if line.startswith('@PG') == True:
                break
        
        contigIDX = -1
        for qname, group1 in groupby(self.fin, lambda line: line.split('$')[0]):
            contigIDX += 1
            if contigIDX%self.batchN != self.batchIDX: continue
            #print('----------------------------------------------------------', QNAME)
            contig = CONTIG(qname)
            for key2, group2 in groupby(group1, lambda line: line.split('\t')[0]):
                #print('----------------------------------------------------------', key2)
                qpos = int(key2.split('$')[1])

                for data in group2:
                    flag, rname, rpos = data.split('\t')[1:4]
                    if rname == '*': continue

                    rpos = int(rpos)
                    strand = 1 if int(flag)&16 == 16 else 0

                    #block = BLOCK(strand, qpos, rpos)

                    contig.add(rname, strand, qpos, rpos)
            
            yield contig

class IMAGE:
    def __init__(self, seqName, seqLen, posRate):
        margin = 100

        self.seqName = seqName
        self.seqLen = seqLen
        self.html = element('html', None)

        self.posRate = posRate


        self.svg = element('svg', self.html)
        self.svg.attr('viewBox', '0 0 {0} {1}'.format(self.cal_pos(self.seqLen) + margin * 2, self.cal_pos(self.seqLen) + margin * 2))
        self.svg.attr('height', self.cal_pos(self.seqLen) + margin * 2)
        self.svg.attr('width',  self.cal_pos(self.seqLen) + margin * 2)
        self.svg.style('background', 'white')
        self.group = element('g', self.svg)
        self.group.attr('transform', 'translate({0},{1}) scale(1)'.format(margin, margin)) #scale(100)
        self.set_ruler()

        self.dotplot_DICT = {}

    def cal_pos(self, pos):
        return self.posRate * pos

    def set_ruler(self):
        ref_ruler = element('rect', self.group)
        ref_ruler.attr('x', -2)
        ref_ruler.attr('y', 0)
        ref_ruler.attr('height', self.cal_pos(self.seqLen))
        ref_ruler.attr('width', 1)
        ref_ruler.attr('fill', 'black')
        query_ruler = element('rect', self.group)
        query_ruler.attr('x', 0)
        query_ruler.attr('y', self.cal_pos(self.seqLen))
        query_ruler.attr('height', 1)
        query_ruler.attr('width', self.cal_pos(self.seqLen))
        query_ruler.attr('fill', 'black')

    def set(self, seqName, seqLen):
        if not seqName in self.dotplot_DICT:
            self.dotplot_DICT[seqName] = DOTPLOT(self.group, seqName, seqLen, self.posRate)
        return self.dotplot_DICT[seqName]
    
    def get(self, seqName):
        return self.dotplot_DICT[seqName]

class DOTPLOT:
    def __init__(self, container, seqName, seqLen, posRate):
        self.seqName = seqName
        self.seqLen = seqLen
        self.posRate = posRate

        self.group = element('g', container)
        self.text = element('text', self.group)
        self.text.add(self.seqName)
        self.text.attr('font-size', 10)
        
        self.strand = 'N'

    def cal_pos(self, pos):
        return self.posRate * pos
    
    def set_position(self, strand, intercept):
        if strand == 0:
            self.strand = '+'
        else:
            self.strand = '-'
        
        x = self.cal_pos(intercept)
        y = 0

        if self.strand == '+':
            self.group.attr('transform', 'translate({0:.3f},{1:.3f}) '.format(x, y))
        else:
            self.group.attr('transform', 'translate({0:.3f},{1:.3f}) scale(-1, 1) '.format(x, y))
    
    def set_border(self, qsPos, qePos, rsPos, rePos):
        self.text.add('   ' + self.strand*5)

        rect = element('rect ', self.group)
        
        x1 = self.cal_pos(qsPos)
        x2 = self.cal_pos(qePos)
        y1 = self.cal_pos(rsPos)
        y2 = self.cal_pos(rePos)

        rect.attr('x', "{0:.3f}".format(x1))
        rect.attr('y', "{0:.3f}".format(y1))
        rect.attr('width',  "{0:.3f}".format(x2 - x1))
        rect.attr('height', "{0:.3f}".format(y2 - y1))
        rect.attr('fill', 'rgba(0,0,0,0.02)')
        rect.attr('stroke', 'gray')
        rect.attr('stroke-width', '0.3')

        if self.strand == '+':
            self.text.attr('x', x1 + 4)
            self.text.attr('y', y1 -  2)
        else:
            self.text.attr('x', x1 + 4 - (x2 - x1))
            self.text.attr('y', y1 - 2)
            self.text.attr('transform', 'scale(-1, 1)')

    def add_block(self, block):
        dot = element('circle ', self.group)

        r = 0.05

        cx = self.cal_pos(block[0]) - r/2
        cy = self.cal_pos(block[1]) - r/2
 
        dot.attr('cx', "{0:.4f}".format(cx))
        dot.attr('cy', "{0:.4f}".format(cy))
        dot.attr('r', r)

        if block[2] == True:
            dot.attr('fill', 'red')
        else:
            dot.attr('fill', 'gray')


###########################################################################################
from optparse import OptionParser
import sys, gzip
#option parser
parser = OptionParser(usage="""Run annotation.py \n Usage: %prog [options]""")
parser.add_option("-n","--name",action = 'store',type = 'string',dest = 'target',help = "")
parser.add_option("-t","--threadN",action = 'store',type = 'int',dest = 'threadN',help = "")
(opt, args) = parser.parse_args()
if opt.threadN == None:
    print('Basic usage')
    print('')
    print('     python 03.draw_dotplot_mp.py -t 24 -n Chr01 (optional)')
    print('')
    sys.exit()

prefix = 'query_100'
qFAR = FAIDX_READER('query.fa.fai')
rFAR = FAIDX_READER('ref/ref.fa.fai')

image_DICT = {}
imageLow_DICT = {}

batchN = opt.threadN
targetName = opt.target
###########################################################################################
def make_contig(batchIDX):
    result = []

    engine = Engine(prefix, batchN, batchIDX)
    for contig in engine.run():
        rContig = contig.removeRepeat(5)
        contig.destory()
        del contig

        fContig = rContig.extract(10)
        rContig.destory()
        del rContig

        fContig.calc_best(1)

        if targetName != None and fContig.bestMatch.rname != targetName:
            fContig.destory()
            del fContig
            continue

        result += [fContig]
    engine.close()

    return result

with Pool(processes=batchN) as pool:
    result_LIST = pool.map(make_contig, range(batchN))

###########################################################################################
fout = open(prefix + '.contig', 'w')
contig_DICT = {}
for contig_LIST in result_LIST:
    for contig in contig_LIST:
        if contig.bestMatch == None: continue

        context  = [contig.bestMatch.qname]
        context += [contig.bestMatch.rname]
        context += ['+' if contig.bestMatch.strand == 0 else '-']
        context += [contig.bestMatch.intercept]
        context += [contig.bestMatch.number]
        context += [contig.bestMatch.variance]
        context += [contig.bestMatch.std_deviation]
        context += [contig.bestMatch.qsPos]
        context += [contig.bestMatch.qePos]
        context += [contig.bestMatch.rsPos]
        context += [contig.bestMatch.rePos]
        context += [contig.bestMatch.coverage]
        fout.write('\t'.join(map(str,context)) + '\n')

        rname = contig.bestMatch.rname
        if not rname in contig_DICT: contig_DICT[rname] = []

        contig_DICT[rname] += [contig]

fout.close()
###########################################################################################
def draw_image(rname):
    rsize = rFAR.get(rname)

    posRate_S = 0.00001
    posRate_L = 0.0003

    image = IMAGE(rname, rsize, posRate_S)
    for contig in contig_DICT[rname]:
        bestMatch = contig.bestMatch
        qsize = qFAR.get(bestMatch.qname)

        dotplot = image.set(bestMatch.qname, qsize)
        dotplot.set_position(bestMatch.strand, bestMatch.intercept)
        dotplot.set_border(1, qsize, bestMatch.rsPos, bestMatch.rePos)

        for _strand, sub_DICT in contig.block_DICT[rname].items():
            for _intercept, block_LIST in sub_DICT.items():
                for block in block_LIST:
                    dotplot.add_block(block)

    fout = open(f'dotplot_S/{rname}.html', 'w')
    fout.write(str(image.html))
    fout.close()
    del image

    image = IMAGE(rname, rsize, posRate_S)
    for contig in contig_DICT[rname]:
        bestMatch = contig.bestMatch
        qsize = qFAR.get(bestMatch.qname)

        dotplot = image.set(bestMatch.qname, qsize)
        dotplot.set_position(bestMatch.strand, bestMatch.intercept)
        dotplot.set_border(1, qsize, bestMatch.rsPos, bestMatch.rePos)

        for _strand, sub_DICT in contig.block_DICT[rname].items():
            for _intercept, block_LIST in sub_DICT.items():
                for blockIDX, block in enumerate(block_LIST):
                    if blockIDX%10 == 0:
                        dotplot.add_block(block)

    fout = open(f'dotplot_S/{rname}_Low.html', 'w')
    fout.write(str(image.html))
    fout.close()
    del image

    image = IMAGE(rname, rsize, posRate_L)
    for contig in contig_DICT[rname]:
        bestMatch = contig.bestMatch
        qsize = qFAR.get(bestMatch.qname)

        dotplot = image.set(bestMatch.qname, qsize)
        dotplot.set_position(bestMatch.strand, bestMatch.intercept)
        dotplot.set_border(1, qsize, bestMatch.rsPos, bestMatch.rePos)

        for _strand, sub_DICT in contig.block_DICT[rname].items():
            for _intercept, block_LIST in sub_DICT.items():
                for block in block_LIST:
                    dotplot.add_block(block)

    fout = open(f'dotplot_L/{rname}.html', 'w')
    fout.write(str(image.html))
    fout.close()
    del image

    image = IMAGE(rname, rsize, posRate_L)
    for contig in contig_DICT[rname]:
        bestMatch = contig.bestMatch
        qsize = qFAR.get(bestMatch.qname)

        dotplot = image.set(bestMatch.qname, qsize)
        dotplot.set_position(bestMatch.strand, bestMatch.intercept)
        dotplot.set_border(1, qsize, bestMatch.rsPos, bestMatch.rePos)

        for _strand, sub_DICT in contig.block_DICT[rname].items():
            for _intercept, block_LIST in sub_DICT.items():
                for blockIDX, block in enumerate(block_LIST):
                    if blockIDX%10 == 0:
                        dotplot.add_block(block)

    fout = open(f'dotplot_L/{rname}_Low.html', 'w')
    fout.write(str(image.html))
    fout.close()
    del image

with Pool(processes=batchN) as pool:
    pool.map(draw_image, contig_DICT.keys())