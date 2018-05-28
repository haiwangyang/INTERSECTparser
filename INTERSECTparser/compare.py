#!/usr/bin/env python

from __future__ import print_function, division

def get_id_from_name(name):
    """ get gene id from bed name """
    return(name.split("ID=")[1].split("-")[0])

def get_start_end_from_string(start_end):
    """ 111.222 => 111, 222 """
    start, end = start_end.split(".")
    start = int(start)
    end = int(end)
    return(start, end)

def write_dct_table(dct, path):
    with open(path, "w") as f:
        for i in sorted(dct.keys()):
            f.write(str(i) + "\t" + str(dct[i]) + "\n")

class Bed:
    """ Bed object (now only use for maker) """
    def __init__(self, feature):
        self.feature = feature
        self.get_feature_ids()
 
    def get_feature_ids(self):
        st = set()
        dct = dict()
        with open("bed/" + self.feature + ".bed", "r") as f:
            for line in f.readlines():
                scaffold, start, end, name, one, strand = line.rstrip().split("\t")
                id = get_id_from_name(name)
                st.add(id)
                if id not in dct.keys():
                    dct[id] = list()
                dct[id].append(start + "." + end)
        self.ids = st
        self.id2ses = dct # id to start_ends


class Intersect:
    """ Intersect object between feature1 (now only maker) and feature2 """
    def __init__(self, feature1, feature2):
        self.feature1 = feature1
        self.feature2 = feature2
        self.get_intersect_dct()        
        self.get_percentage_dct()

    def get_intersect_dct(self):
        """ get gene-level intersect dict """
        dct = dict()
        with open("intersect/" + self.feature1 + "." + self.feature2 + ".intersect", "r") as f:
            for line in f.readlines():
                scaffold1, start1, end1, name1, one1, strand1, scaffold2, start2, end2, name2, one2, strand2, overlap = line.rstrip().split("\t")
                id1 = get_id_from_name(name1)
                if not id1 in dct.keys():
                    dct[id1] = dict()
                if not start1 + "." + end1 in dct[id1].keys():
                    dct[id1][start1 + "." + end1] = list()
                dct[id1][start1 + "." + end1].append(start2 + "." + end2)
        self.intersect_dct = dct

    def get_percentage_dct(self):
        """ get percertage of coverage on gene-level """
        b = Bed("maker")
        dct = dict()
        for id1 in b.ids:
            if id1 in self.intersect_dct.keys():
                dct_overlap = dict()
                for start1_end1 in self.intersect_dct[id1].keys():
                    start1, end1 = get_start_end_from_string(start1_end1)
                    #print("START1: " + str(start1) + "\t" + "END1: " + str(end1))
                    for i in range(start1, end1, 1):
                        dct_overlap[i] = 0
                    for start2_end2 in self.intersect_dct[id1][start1_end1]:
                        start2, end2 = get_start_end_from_string(start2_end2)
                        #print("START2: " + str(start2) + "\t" + "END2: " + str(end2))
                        for ii in range(start2, end2, 1):
                            if ii in dct_overlap.keys():
                                dct_overlap[ii] = 1

                count_total = 0
                count_overlap = 0
                for iii in dct_overlap.keys():
                    if dct_overlap[iii] == 0:
                        count_total += 1
                    if dct_overlap[iii] == 1:
                        count_total += 1
                        count_overlap += 1
                dct[id1] = count_overlap / count_total
                #print(count_overlap / count_total)
                #print("\n\n")
            else:
                dct[id1] = 0
        self.percentage_dct = dct 


def main():
    b = Bed("maker")
    for feature1 in ["maker", ]:
        for feature2 in ["repeatmasker", "protein2genome", "genemark", "ori_snap", "ori_augustus"]:
            i = Intersect(feature1, feature2)
            write_dct_table(i.percentage_dct, "output/" + feature1 + "." + feature2 + ".percentage.txt")

if __name__ == '__main__':
    main()
