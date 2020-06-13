#!/usr/bin/env python3

LIBVERSION = "1.0"

g_verbose = 0

class PasttrecDefaults:
    c_cable = [0x00 << 19, 0x00 << 19, 0x00 << 19]
    c_asic = [0x2000, 0x4000]

#                Bg_int,K,Tp      TC1      TC2      Vth
    c_config_reg = [0x00000, 0x00100, 0x00200, 0x00300]
    c_bl_reg = [0x00400, 0x00500, 0x00600, 0x00700,
                0x00800, 0x00900, 0x00a00, 0x00b00]

    c_base_w = 0x0050000
    c_base_r = 0x0051000


class PasttrecRegs(PasttrecDefaults):
    bg_int = 1
    gain = 0
    peaking = 0
    tc1c = 0
    tc1r = 0
    tc2c = 0
    tc2r = 0
    vth = 0
    bl = [0] * 8

    def __init__(self, bg_int=1, gain=0, peaking=0,
                 tc1c=0, tc1r=0, tc2c=0, tc2r=0,
                 vth=0, bl=[0]*8):
        self.bg_int = bg_int
        self.gain = gain
        self.peaking = peaking
        self.tc1c = tc1c
        self.tc1r = tc1r
        self.tc2c = tc2c
        self.tc2r = tc2r
        self.vth = vth
        self.bl = [i for i in bl]
	

    @staticmethod
    def load_asic_from_dict(d, test_version=None):
        if (test_version is not None) and (test_version != LIBVERSION):
            return False
        p = PasttrecRegs()
        for k, v in d.items():
            setattr(p, k, v)
        return p

    def dump_config(self):
        r_all = [0] * 12
        t = (self.bg_int << 4) | (self.gain << 2) | self.peaking
        r_all[0] = self.c_config_reg[0] | t
        t = (self.tc1c << 3) | self.tc1r
        r_all[1] = self.c_config_reg[1] | t
        t = (self.tc2c << 3) | self.tc2r
        r_all[2] = self.c_config_reg[2] | t
        r_all[3] = self.c_config_reg[3] | self.vth

        for i in range(8):
            r_all[4+i] = self.c_bl_reg[i] | self.bl[i]

        return r_all

    def dump_config_hex(self):
        return [hex(i) for i in self.dump_config()]

    def dump_bl_hex(self):
        return [hex(i) for i in self.dump_config()[4:]]


class PasttrecCard():
    name = None
    asic1 = None
    asic2 = None

    def __init__(self, name, asic1=None, asic2=None):
        self.name = name
        self.asic1 = asic1
        self.asic2 = asic2

    def set_asic(self, pos, asic):
        if pos == 0:
            self.asic1 = asic
        elif pos == 1:
            self.asic2 = asic

    def export(self):
        return {
            'name': self.name,
            'asic1': self.asic1.__dict__ if self.asic1 is not None else None,
            'asic2': self.asic2.__dict__ if self.asic2 is not None else None
        }

    def export_script(self, cable):
        regs = []
        if self.asic1:
            regs.extend(self.asic1.dump_config(cable, 0))
        if self.asic2:
            regs.extend(self.asic2.dump_config(cable, 1))
        return regs

    @staticmethod
    def load_card_from_dict(d, test_version=None):
        if (test_version is not None) and (test_version != LIBVERSION):
            return False, LIBVERSION

        if d is None:
            return False, None

        pc = PasttrecCard(d['name'],
                          PasttrecRegs().load_asic_from_dict(d['asic1']),
                          PasttrecRegs().load_asic_from_dict(d['asic2']))

        return True, pc


class TdcConnection():
    id = 0
    cable1 = None
    cable2 = None
    cable3 = None

    def __init__(self, id, cable1=None, cable2=None, cable3=None):
        self.id = hex(id) if isinstance(id, int) else id
        self.cable1 = cable1
        self.cable2 = cable2
        self.cable3 = cable3

    def set_card(self, pos, card):
        if pos == 0:
            self.cable1 = card
        elif pos == 1:
            self.cable2 = card
        elif pos == 2:
            self.cable3 = card

    def export(self):
        c1 = self.cable1.export() if isinstance(self.cable1, PasttrecCard) \
            else None
        c2 = self.cable2.export() if isinstance(self.cable2, PasttrecCard) \
            else None
        c3 = self.cable3.export() if isinstance(self.cable3, PasttrecCard) \
            else None

        return self.id, {
            'cable1': c1,
            'cable2': c2,
            'cable3': c3
        }

    def export_script(self):
        c1 = self.cable1.export_script(0) \
            if isinstance(self.cable1, PasttrecCard) else None
        c2 = self.cable2.export_script(1) \
            if isinstance(self.cable2, PasttrecCard) else None
        c3 = self.cable3.export_script(2) \
            if isinstance(self.cable3, PasttrecCard) else None

        c = []
        if c1:
            c.extend(c1)
        if c2:
            c.extend(c2)
        if c3:
            c.extend(c3)
        return self.id, c


def dump(tdcs):
    d = {'version': LIBVERSION}
    if isinstance(tdcs, list):
        for t in tdcs:
            k, v = t.export()
            d[k] = v
    elif isinstance(tdcs, TdcConnection):
        k, v = tdcs.export()
        d[k] = v

    return d


def dump_script(tdcs):
    d = []
    if isinstance(tdcs, list):
        for t in tdcs:
            k, v = t.export_script()
            for _v in v:
                d.append("trbcmd w {:s} 0xa000 {:s}".format(k, hex(_v)))
    elif isinstance(tdcs, TdcConnection):
        k, v = tdcs.export_script()
        for _v in v:
            d.append("trbcmd w {:s} 0xa000 {:s}".format(k, hex(_v)))

    return d


def load(d, test_version=True):
    if test_version:
        if 'version' in d:
            if d['version'] != LIBVERSION:
                return False, d['version']
        else:
            return False, '0.0.0'

    connections = []
    for k, v in d.items():
        if k == 'version':
            continue

        id = int(k, 16)
        r1, _c1 = PasttrecCard.load_card_from_dict(v['cable1'])
        r2, _c2 = PasttrecCard.load_card_from_dict(v['cable2'])
        r3, _c3 = PasttrecCard.load_card_from_dict(v['cable3'])

        c1 = _c1 if r1 else None
        c2 = _c2 if r2 else None
        c3 = _c3 if r3 else None

        connections.append(TdcConnection(id, cable1=c1, cable2=c2, cable3=c3))

    return True, connections


def_max_bl_registers = 32
def_pastrec_channel_range = 8


class Baselines:
    baselines = None
    config = None

    def __init__(self):
        self.baselines = {}

    def add_trb(self, trb):
        if trb not in self.baselines:
            w = def_max_bl_registers
            h = def_pastrec_channel_range
            a = len(PasttrecDefaults.c_asic)
            c = len(PasttrecDefaults.c_cable)
            self.baselines[trb] = [[[[0 for x in range(w)] for y in range(h)]
                                    for _a in range(a)] for _c in range(c)]

import argparse
from colorama import Fore, Style
import copy
import json
import sys
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit
from scipy.interpolate import splrep, sproot, splev
import os
import math
import sys,glob
import ROOT
#from array import array
#from ROOT import *
#f = ROOT.TFile("outfile.root", "recreate")
#c_counts = ROOT.TCanvas("c_counts","c_counts")

def bl_list_with_marker(l, pos):
    s = ""
    for i in range(len(l)):
        if i == pos:
            s += Fore.YELLOW + "{:d}".format(l[i]) + Style.RESET_ALL + ", "
        else:
            s += "{:d}, ".format(l[i])
    return s

if __name__=="__main__":
    parser=argparse.ArgumentParser(description='Calculate baseline width')
    parser.add_argument('json_file', help='list of arguments', type=str)
    args=parser.parse_args()

    print(args)
    
          
    base=os.path.basename(args.json_file)
    out_txt = os.path.splitext(base)[0] + ".txt"
    file = open(out_txt,"w")

    counts_txt = os.path.splitext(base)[0] + "cnts.txt"
    file2 = open(counts_txt,"w")
    
    Fullwidth_list = []
    Peak_norm = []


    tdc_cal0 = [0,1,2]
    tdc_cal1 = [0,1,2]
    tdc_cal2 = [0]

    plt.figure(1)

    with open(args.json_file) as json_data:
        d = json.load(json_data)
        json_data.close()

    tdc_cal = []
    file_name = args.json_file
    print("just_checking ",file_name)
   

    print("TDC CAL : ",tdc_cal)

    bls = d['baselines']
    cfg = d['config']

    cable_counter = 0

    x_ind =[]
    tlist = []
    p = PasttrecRegs()


    for k,v in cfg.items():
        setattr(p, k, v)

    print(cfg)

    x = list(range(0,32))
    idx = 1
    
    for k,v in bls.items():
        for reg in list(range(32)):
            for c in [0,1,2]:
                for a in [0,1]:
                    for ch in list(range(8)): 
                        file2.write("{:d}\t".format(v[c][a][ch][reg]))
            file2.write("\n")
    
    
    for k,v in bls.items():
        t = TdcConnection(k)
        plt.figure(idx)
        print("K : ",k)
        if (k == "0x6464" or k == "0x6465"):
            tdc_cal = tdc_cal1
        else:
            tdc_cal = tdc_cal0
        for c in tdc_cal:
            card = PasttrecCard("noname")
            cable_counter = cable_counter + 1
            for a in [0,1]:
                print(Fore.YELLOW + "Scanning {:s}  CARD: {:d}  ASIC: {:d}".format(k, c, a) + Style.RESET_ALL)
                bl = [0] * 8
                for ch in list(range(8)):
                    b = v[c][a][ch]
                    s = 0
                    w = 0                   
                    
                    for i in range(32):
                        s = s + (i+1) * b[i]
                        w += b[i]
                    if w == 0:
                        b = 0
                    else:
                        b = s/w - 1
                    bl[ch] = int(round(b))

                    print(ch,
                        " bl:", Fore.GREEN, "{:2d}".format(bl[ch]), Style.RESET_ALL,
                        "(0x{:s})".format(hex(bl[ch])[2:].zfill(2)),
                        Fore.RED, "{:>+3d} mV".format(-31 + 2 * bl[ch]), Style.RESET_ALL,
                        " [ ", bl_list_with_marker(v[c][a][ch], bl[ch]), "]")
                    
                    dd = v[c][a][ch]

                    x_valu =[]
                    
                    count_sum = 0
                    bl_sum = 0
                    bl_count_sum =0
                    bl_mean =0
                    bl_blmean_sqr =0;
                    sigma = 0
                    sqr_wave =0
                    
                    #Calc mean baseline value in mV
                    for b_pos,b_cnt in enumerate(dd):
                        bl_pos_mV = ( b_pos * 2 ) - 32  # conver lsb to mV
                        if (b_cnt > 0 ):
                            count_sum += b_cnt
                            bl_sum += bl_pos_mV 
                            bl_count_sum += (b_cnt * bl_pos_mV) #summation (count_i * Bl_Value_i)
                            sqr_wave = sqr_wave+1

                    if (bl_count_sum != 0) :

                        bl_mean = bl_count_sum/count_sum
                        
                    # Calculate Standard Deviation    
                    for b_poss ,b_cntt in enumerate(dd):
                        bl_poss_mV = ( b_poss * 2 ) - 32  # conver lsb to mV
                        if (b_cntt > 0) :
                           bl_blmean_sqr += (((bl_poss_mV - bl_mean )**2 ) * b_cntt )                           
                           

                    if (bl_count_sum != 0) :
                        if(sqr_wave ==1) :
                            sigma = 0.2886 * 2 #mV Vp* 1/srt(12)
                        else :
                            sigma = math.sqrt(bl_blmean_sqr / count_sum) #standard deviation
                    
                    #coutstopeak = {}
                    #postopeak = array('d')
                    graphs = ROOT.TGraph()
                    for i, j in enumerate(dd):
                        if j > 0: 
                            x_valu.append(i)
                        #if (cable_counter == 1):
                            #if j > max(dd):
                                #coutstopeak.append(j)
                                #postopeak.append(int(i))  
                            #graphs = ROOT.TGraph(len(postopeak), postopeak, coutstopeak)
                            #graphs.Draw("same")
                            #graphs.Write()
                    ch_status = 1
                    


                    if dd[0]==0 and dd[31]==0 and max(dd)>0:
                        Full_width = len(x_valu) * 2 # conver lsb to mV 
                        file.write("{:s} {:d} {:d} {:d} {:d} {:d} {:.3f} {:.3f} {:d} {:s}\n".format(k,c,a,(8*a)+ch,Full_width,np.argmax(dd),bl_mean,sigma,ch_status,file_name))

                    elif dd[0]!=0 and dd[31]!=0:
                        Full_width = len(x_valu) * 2
                        ch_status = 4;

                        file.write("{:s} {:d} {:d} {:d} {:d} {:d} {:.3f} {:.3f} {:d} {:s}\n".format(k,c,a,(8*a)+ch,Full_width,np.argmax(dd),bl_mean,sigma,ch_status,file_name))

                    elif dd[0]==0 and dd[31]!=0:
                        Full_width = len(x_valu) * 2
                        ch_status = 2;

                        file.write("{:s} {:d} {:d} {:d} {:d} {:d} {:.3f} {:.3f} {:d} {:s}\n".format(k,c,a,(8*a)+ch,Full_width,np.argmax(dd),bl_mean,sigma,ch_status,file_name))

                    elif dd[0]!=0 and dd[31]==0:
                        Full_width = len(x_valu) * 2
                        ch_status = 3;

                        file.write("{:s} {:d} {:d} {:d} {:d} {:d} {:.3f} {:.3f} {:d} {:s}\n".format(k,c,a,(8*a)+ch,Full_width,np.argmax(dd),bl_mean,sigma,ch_status,file_name))

                    else:
                        Full_width = len(x_valu) * 2
                        ch_status = 0;

                        file.write("{:s} {:d} {:d} {:d} {:d} {:d} {:.3f} {:.3f} {:d} {:s}\n".format(k,c,a,(8*a)+ch,Full_width,np.argmax(dd),bl_mean,sigma,ch_status,file_name))

                    Ch_sq_no = (16*c)+(8*a)+ch
                    Fullwidth_list.insert(Ch_sq_no,Full_width)

            t.set_card(c, card)
            #file2.writelines("\n")
        tlist.append(t)
        idx += 1
        #all_peaks =[]
        #max_of_peak = max(Peak_norm)
        #for q in Peak_norm:
            #all_peaks.append((q*5)/max_of_peak)

        #xx = [a for a in range((len(filelist)*cable_counter*16))] # -16 if/since one of the files contains only the results for 1 tdc
        xx = [a for a in range((cable_counter*16))] # -16 if/since one of the files contains only the results for 1 tdc
   #     print(Fullwidth_list)
   #     print(len(xx))
        print(len(Fullwidth_list))
        plt.plot(xx, Fullwidth_list, label="some_label1")
        plt.xticks([i*16.0 for i in range(0,2)])
        plt.yticks([0.,0.5,1.,1.5,2,2.5,3,3.5,4,4.5,5,5.5,6,7,8,9,10,11,12,13,14,15])
        plt.grid(True)

        file.close()
        file2.close()
      #  plt.show()

