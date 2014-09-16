import csv
from collections import defaultdict
import codecs
import os
from operator import attrgetter

class Segment:
    def __init__(self, l, prefix=''):
        self.chromosome = ('0'+l['CHROMOSOME'])[-2:]
        if l['CHROMOSOME'] == 'X':
            self.chromosomeN = 23
        else:
            self.chromosomeN = int(l['CHROMOSOME']) 
        self.name = prefix + codecs.decode(l['MATCHNAME'], 'utf-8').replace('  ', ' ')
        self.L1 = int(l['START LOCATION'])
        self.L2 = int(l['END LOCATION'])
        self.cm = float(l['CENTIMORGANS'])
        self.full_sort_key = '%02d/%09d-%09d' % (self.chromosomeN, self.L1, self.L2)

    def contains(self, chromosome, loc):
        return chromosome == self.chromosome and self.L1 <= loc and loc <= self.L2

    def overlaps(self, segment2):
        if self.chromosome != segment2.chromosome:
            return None
    
        overlapL1 = max(self.L1, segment2.L1)
        overlapL2 = min(self.L2, segment2.L2)
        if overlapL1 >= overlapL2:
            return None

        # ok, so they overlap. Can we estimate the cM value of that overlap?
        # we'll do it using both as the basis for pro-rating, and take the arithmetic mean

        return (0.5 * self.cm * (overlapL2 - overlapL1) / (self.L2 - self.L1) +
                0.5 * segment2.cm * (overlapL2 - overlapL1) / (segment2.L2 - segment2.L1))

class FTDNA:
    def __init__(self, filename = None, chromosome = None, minCM = 0.0, prefix = '', excludes = [], kit = ''):
        self.segments = []
        self.load(filename, chromosome, minCM, prefix, excludes, kit)

    def load(self, filename = None, chromosome = None, minCM = 0.0, prefix = '', excludes = [], kit=''):
        if filename == None:
            filename = max([fn for fn in os.listdir('.') if '_chromosome_browser_results' in fn and fn.startswith(kit)])
        self.filename = filename

        lastL = ""
        for l in csv.DictReader(open(filename, 'rb')):
            if l == lastL:
                continue
            lastL = l
            segment = Segment(l, prefix)
            if chromosome and chromosome != segment.chromosomeN:
                continue
            if segment.cm < minCM:
                continue
            if segment.name in excludes:
                continue
            self.segments.append(segment)

    def segments_sorted_by(self, criterion, reverse=False):
        return sorted(self.segments, key=attrgetter(criterion), reverse=reverse)

    def segments_grouped_by(self, criterion, reverse=False):
        group = defaultdict(list)
        ag = attrgetter(criterion)
        for s in self.segments:
            k = ag(s)
            group[k].append(s)
        return group
