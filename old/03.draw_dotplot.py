from lib.KJH_SVG.KJH_SVG import element
from lib.JobTimer.JobTimer import JobTimer
from itertools import groupby
from scipy.stats import chi2
import math
import numpy as np

###########################################################################################
def cal_pos(pos):
    POS_RATE = 0.00003
    return pos * POS_RATE
###########################################################################################
def chi_squared_test(variance, number):
    sigma0_squared = 200000 * 200000

    chi_squared_stat = (number - 1) * variance / sigma0_squared

    alpha = 0.05

    df = number - 1
    critical_value = chi2.ppf(1 - alpha, df)

    return chi_squared_stat < critical_value
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
        self.totalN = 0

        fin = open(self.fileName)
        for line in fin:
            seqName, seqLen = line.rstrip('\n').split('\t')[0:2]
            seqLen = int(seqLen)
            self.seqLen_DICT[seqName] = seqLen
            self.totalN += seqLen
        fin.close()
    
    def get(self, seqName):
        return self.seqLen_DICT[seqName]

class BLOCK:
    def __init__(self, rname, strand, qpos, rpos):
        self.rname = rname
        self.strand = strand
        self.qpos = qpos
        self.rpos = rpos

        self.isIn = False
        self.isRepeat = False

        if strand == '+':
            self.intercept = -qpos + rpos
        elif strand == '-':
            self.intercept =  qpos + rpos

class CONTIG:
    def __init__(self, qname):
        self.qname = qname
        self.block_DICT = {}
        self.blockN_DICT = {}
        self.block_LIST = []
    
    def add(self, block):
        self.block_LIST += [block]
        if not block.rname in self.block_DICT:
            self.block_DICT[block.rname] = {'+':{}, '-':{}}
            self.blockN_DICT[block.rname] = 0
        
        self.blockN_DICT[block.rname] += 1
        
        if not block.intercept in self.block_DICT[block.rname][block.strand]:
            self.block_DICT[block.rname][block.strand][block.intercept]  = [block]
        else:
            self.block_DICT[block.rname][block.strand][block.intercept] += [block]

    def removeRepeat(self, maxRepeat):
        _qblock_DICT = {}
        for rname, sub1_DICT in self.block_DICT.items():
            _qblock_DICT[rname] = {}
            for strand, sub2_DICT in sub1_DICT.items():
                for intercept, block_LIST in sub2_DICT.items():
                    for block in block_LIST:
                        if not block.qpos in _qblock_DICT[rname]:
                            _qblock_DICT[rname][block.qpos]  = [block]
                        else:
                            _qblock_DICT[rname][block.qpos] += [block]

        for rname, sub1_DICT in _qblock_DICT.items():
            for qpos, block_LIST in sub1_DICT.items():
                if len(block_LIST) < maxRepeat: continue
                for block in block_LIST:
                    block.isRepeat = True

        _rblock_DICT = {}
        for rname, sub1_DICT in self.block_DICT.items():
            _rblock_DICT[rname] = {}
            for strand, sub2_DICT in sub1_DICT.items():
                for intercept, block_LIST in sub2_DICT.items():
                    for block in block_LIST:
                        if not block.rpos in _rblock_DICT[rname]:
                            _rblock_DICT[rname][block.rpos]  = [block]
                        else:
                            _rblock_DICT[rname][block.rpos] += [block]

        for rname, sub1_DICT in _rblock_DICT.items():
            for rpos, block_LIST in sub1_DICT.items():
                if len(block_LIST) < maxRepeat: continue
                for block in block_LIST:
                    block.isRepeat = True

        _contig = CONTIG(self.qname)
        for rname, sub1_DICT in self.block_DICT.items():
            for strand, sub2_DICT in sub1_DICT.items():
                for intercept, block_LIST in sub2_DICT.items():
                    for block in block_LIST:
                        if block.isRepeat == True: continue
                        _contig.add(block)

        return _contig
    
    def extract(self, minBlockN):
        _contig = CONTIG(self.qname)
        for rname, sub1_DICT in self.block_DICT.items():
            for strand, sub2_DICT in sub1_DICT.items():
                for intercept, block_LIST in sub2_DICT.items():
                    blockN = len(block_LIST)
                    if blockN < minBlockN: continue

                    for block in block_LIST:
                        _contig.add(block)

        return _contig

        #self.block_DICT = _contig.block_DICT
        #self.blockN_DICT = _contig.blockN_DICT

    def get_most(self):
        rname = None
        nubmer = 1
        for _rname in self.block_DICT.keys():
            qname, strand, intercept, _number, variance, std_deviation, qsPos, qePos, rsPos, rePos, coverage = self.calc(_rname, 1)
            if _number > nubmer:
                nubmer = _number
                rname = _rname
        
        return rname

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

            if strand == '+':
                interceptP_LIST += [intercept]* blockN

            if strand == '-':
                interceptM_LIST += [intercept]* blockN

            if subN/totalN > coverage: break
        
        #print(rname, len(interceptP_LIST), len(interceptM_LIST))
        if len(interceptP_LIST) > len(interceptM_LIST):
            STRAND = '+'
            meanIntercept, variance = median_without_outliers(interceptP_LIST)
            std_deviation = variance ** 0.5
        else:
            STRAND = '-'
            meanIntercept, variance = median_without_outliers(interceptM_LIST)
            std_deviation = variance ** 0.5

        #calculate boundery
        for [strand, intercept, block_LIST] in target_LIST:
            if STRAND != strand: continue
            if abs(intercept - meanIntercept) >  std_deviation * 3: continue
            for block in block_LIST:
                block.isIn = True

                qsPos = min(qsPos, block.qpos)
                qePos = max(qePos, block.qpos)

                rsPos = min(rsPos, block.rpos)
                rePos = max(rePos, block.rpos)

        number = len(interceptP_LIST) + len(interceptM_LIST)
        return self.qname, STRAND, meanIntercept, number, variance, std_deviation, qsPos, qePos, rsPos, rePos, float(qePos - qsPos) / (rePos - rsPos)
            
class Engine:
    def __init__(self, prefix):
        self.READ_SIZE = 100

        self.fin = open(prefix + '.sam')
    
    def run(self):
        #skip head
        for line in self.fin:
            if line.startswith('@PG') == True:
                break

        for qname, group1 in groupby(self.fin, lambda line: line.split('$')[0]):
            #print('----------------------------------------------------------', QNAME)
            contig = CONTIG(qname)
            for key2, group2 in groupby(group1, lambda line: line.split('\t')[0]):
                #print('----------------------------------------------------------', key2)
                qpos = int(key2.split('$')[1])

                for data in group2:
                    flag, rname, rpos = data.split('\t')[1:4]
                    if rname == '*': continue

                    rpos = int(rpos)
                    strand = '-' if int(flag)&16 == 16 else '+'

                    block = BLOCK(rname, strand, qpos, rpos)

                    contig.add(block)
            
            yield contig


class IMAGE:
    def __init__(self, seqName, seqLen):
        margin = 100

        self.seqName = seqName
        self.seqLen = seqLen
        self.html = element('html', None)


        self.svg = element('svg', self.html)
        self.svg.attr('viewBox', '0 0 {0} {1}'.format(cal_pos(self.seqLen) + margin * 2, cal_pos(self.seqLen) + margin * 2))
        self.svg.attr('height', cal_pos(self.seqLen) + margin * 2)
        self.svg.attr('width',  cal_pos(self.seqLen) + margin * 2)
        self.svg.style('background', 'white')
        self.group = element('g', self.svg)
        self.group.attr('transform', 'translate({0},{1}) scale(1)'.format(margin, margin)) #scale(100)
        self.set_ruler()

        self.dotplot_DICT = {}

    def set_ruler(self):
        ref_ruler = element('rect', self.group)
        ref_ruler.attr('x', -2)
        ref_ruler.attr('y', 0)
        ref_ruler.attr('height', cal_pos(self.seqLen))
        ref_ruler.attr('width', 1)
        ref_ruler.attr('fill', 'black')
        query_ruler = element('rect', self.group)
        query_ruler.attr('x', 0)
        query_ruler.attr('y', cal_pos(self.seqLen))
        query_ruler.attr('height', 1)
        query_ruler.attr('width', cal_pos(self.seqLen))
        query_ruler.attr('fill', 'black')

    def set(self, seqName, seqLen):
        if not seqName in self.dotplot_DICT:
            self.dotplot_DICT[seqName] = DOTPLOT(self.group, seqName, seqLen)
        return self.dotplot_DICT[seqName]
    
    def get(self, seqName):
        return self.dotplot_DICT[seqName]

class DOTPLOT:
    def __init__(self, container, seqName, seqLen):
        self.seqName = seqName
        self.seqLen = seqLen

        self.group = element('g', container)
        self.text = element('text', self.group)
        self.text.add(self.seqName)
        self.text.attr('font-size', 10)
        
        self.strand = '+'

    def set_position(self, strand, intercept):
        self.strand = strand
        
        x = cal_pos(intercept)
        y = 0

        if strand == '+':
            self.group.attr('transform', 'translate({0:.3f},{1:.3f}) '.format(x, y))
        else:
            self.group.attr('transform', 'translate({0:.3f},{1:.3f}) scale(-1, 1) '.format(x, y))
    
    def set_border(self, qsPos, qePos, rsPos, rePos):
        self.text.add('   ' + strand*5)

        rect = element('rect ', self.group)
        
        x1 = cal_pos(qsPos)
        x2 = cal_pos(qePos)
        y1 = cal_pos(rsPos)
        y2 = cal_pos(rePos)

        rect.attr('x', "{0:.3f}".format(x1))
        rect.attr('y', "{0:.3f}".format(y1))
        rect.attr('width',  "{0:.3f}".format(x2 - x1))
        rect.attr('height', "{0:.3f}".format(y2 - y1))
        rect.attr('fill', 'rgba(0,0,0,0.02)')
        rect.attr('stroke', 'gray')
        rect.attr('stroke-width', '0.3')

        if strand == '+':
            self.text.attr('x', x1 + 4)
            self.text.attr('y', y1 -  2)
        else:
            self.text.attr('x', x1 + 4 - (x2 - x1))
            self.text.attr('y', y1 - 2)
            self.text.attr('transform', 'scale(-1, 1)')

    def add_block(self, block):
        dot = element('circle ', self.group)

        r = 0.05

        cx = cal_pos(block.qpos) - r/2
        cy = cal_pos(block.rpos) - r/2
 
        dot.attr('cx', "{0:.4f}".format(cx))
        dot.attr('cy', "{0:.4f}".format(cy))
        dot.attr('r', r)

        if block.isIn == True:
            dot.attr('fill', 'red')
        else:
            dot.attr('fill', 'gray')


###########################################################################################
jobTimer = JobTimer()

prefix = 'query_100'
qFAR = FAIDX_READER('query.fa.fai')
rFAR = FAIDX_READER('ref/ref.fa.fai')
engine = Engine(prefix)

image_DICT = {}
imageLow_DICT = {}

jobTimer.reset()
processedN = 0

fout = open(prefix + '.contig', 'w')
for contig in engine.run():
    rContig = contig.removeRepeat(5)
    fContig = rContig.extract(10)

    print(len(contig.block_LIST), len(rContig.block_LIST), len(fContig.block_LIST))

    qname = contig.qname
    rname = fContig.get_most()
    if rname == None: continue

    rsize = rFAR.get(rname)
    qsize = qFAR.get(qname)

    qname, strand, intercept, number, variance, std_deviation, qsPos, qePos, rsPos, rePos, coverage = fContig.calc(rname, 1)
    
    if not rname in image_DICT: 
        image_DICT[rname] = IMAGE(rname, rFAR.get(rname))

    #High Resolution
    if not rname in image_DICT: 
        image_DICT[rname] = IMAGE(rname, rFAR.get(rname))

    dotplot = image_DICT[rname].set(qname, qsize)
    dotplot.set_position(strand, intercept)
    dotplot.set_border(1, qsize, rsPos, rePos)

    for _strand, sub_DICT in fContig.block_DICT[rname].items():
        for _intercept, block_LIST in sub_DICT.items():
            for block in block_LIST:
                dotplot.add_block(block)
    
    #Low Resolution
    if not rname in imageLow_DICT: 
        imageLow_DICT[rname] = IMAGE(rname, rFAR.get(rname))
        
    dotplotLow = imageLow_DICT[rname].set(qname, qsize)
    dotplotLow.set_position(strand, intercept)
    dotplotLow.set_border(1, qsize, rsPos, rePos)

    for _strand, sub_DICT in fContig.block_DICT[rname].items():
        for _intercept, block_LIST in sub_DICT.items():
            for blockIDX, block in enumerate(block_LIST):
                if blockIDX%10 == 0:
                    dotplotLow.add_block(block)

    processedN += qsize
    jobTimer.check()
    percentage = float(processedN)/qFAR.totalN
    fout.write('\t'.join(map(str,[qname, rname, strand, intercept, number, variance, std_deviation, qsPos, qePos, rsPos, rePos, coverage])) + '\n')
    fout.flush()
    print("Process... [{0:6.2f}%] remainTime: {1}  ----   {2:>10}{3:>20}{4:>20.5f}{5:>20.5f}{6:>20}".format(percentage*100, jobTimer.remainTime(percentage), qname, rname, std_deviation, coverage, number))

fout.close()
for rname, image in image_DICT.items():
    fout = open('dotplot_S/{0}.html'.format(rname), 'w')
    fout.write(str(image.html))
    fout.close()

for rname, image in imageLow_DICT.items():
    fout = open('dotplot_S/{0}_Low.html'.format(rname), 'w')
    fout.write(str(image.html))
    fout.close()