# ----------------------------------------------------------------------------
# Authors: Seung-been "Steven" Lee, Aparna Radhakrishnan, Sean McGee
# The license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import copy
import datetime
import operator
import os
import statistics
import subprocess
import sys
import timeit
import types
import shutil
import logging
import gzip
import pkgutil
import time

import pandas as pd
import numpy as np
from phenotyper import phenotyper
from bam2gdf import bam2gdf

from version import __version__
from _get_parser import get_parser
from _call_sv import call_sv
from _plot import plot
from _get_loci import get_loci
from _report import create_report

LINEBREAK = '-' * 70

class VCF:
    def __init__(self):
        self.meta = []
        self.header = []
        self.data = []

class GDF:
    def __init__(self):
        self.header = []
        self.data = []

class SNPAllele:
    def __init__(self):
        self.pos = '' # reference genome position
        self.wt = '' # wild type (*1) allele
        self.var = '' # variant allele
        self.rs = '' # rs ID
        self.het = False # heterozygous
        self.ad = 0 # allelic depth
        self.td = 0 # total depth
        self.n = '' # SNP table number
        self.hg = '' # reference genome allele
        self.so = '' # sequence ontology
        self.effect = '' # coding effect
        self.impact = '' # variant impact
        self.rev = False # reverting variant

    @property
    def key(self):
        return (self.pos, self.wt, self.var)

    @property
    def af(self):
        return 0 if self.td == 0 else self.ad / self.td

    def __eq__(self, other):
        return self.key == other.key

    def __hash__(self):
        return hash(self.key)

    def summary(self):
        return '<{}:{}>{}:{}/{}:{:.2f}:{}:{}:{}>'.format(self.pos, self.wt, self.var, self.ad, self.td, self.af, self.so, self.impact, self.effect)

class StarAllele:
    def __init__(self):
        self.name = ''
        self.score = -100.0
        self.core = []
        self.tag = []
        self.sv = ''

    @property
    def ranked_as(self):
        """
        Unknown function alleles should be broken ties with normal function
        alleles using attributes other than activity score. Increased
        function alleles should come before normal function alleles.
        """

        if self.score < 0:
            return 1.0
        elif self.score > 1:
            return 0.99
        else:
            return self.score

    @property
    def rank(self):
        return (self.ranked_as, -1 * int(bool(self.sv)), -1 * len(self.core))

    def __str__(self):
        return self.name

    def __repr__(self):
        return str(self)

    def __eq__(self, other):
        return self.name == other.name

    def __hash__(self):
        return hash(self.name)

class Haplotype:
    def __init__(self):
        self.cand = []
        self.obs = []

    @property
    def sv(self):
        sv = 'no_sv'
        sv_list = []
        for star_allele in self.cand:
            if star_allele.sv and star_allele.sv not in sv_list:
                sv_list.append(star_allele.sv)
        if len(sv_list) > 1:
            raise ValueError('haplotype contains multiple structural variant calls')
        if len(sv_list) == 1:
            sv = sv_list[0]
        return sv

    @property
    def af(self):
        return [0] if not self.obs else [x.af for x in self.obs]

    @property
    def af_mean_main(self):
        filtered = [x.af for x in self.obs if x.td > 10 and x.het and x in [y for y in self.cand[0].core]]
        return -1 if not filtered else statistics.mean(filtered)

    def af_mean_gene(self, start, end):
        filtered = [x.af for x in self.obs if x.td > 10 and x.het and start <= x.pos <= end]
        return -1 if not filtered else statistics.mean(filtered)

    def fit_data(self, total_cn, start, end):
        """Return the fit MAF and CN."""
        maf_choices = []
        for i in range(1, total_cn):
            maf_choices.append(i / total_cn)
        fit_maf = min(maf_choices, key = lambda x: abs(x - self.af_mean_gene(start, end)))
        fit_cn = maf_choices.index(fit_maf) + 1
        return fit_maf, fit_cn

    def remove_star(self, sX):
        """Remove the given star allele from the candidates list."""
        for i, star in enumerate(self.cand):
            if star.name == sX.name:
                del self.cand[i]
                break

    def add_dup(self, cn):
        """Duplicate the main star allele by the given CN."""
        if cn == 1:
            return
        if cn > 10:
            cn = 10
        sX = self.cand[0]
        score = sX.score * cn
        name = sX.name + 'x' + str(cn)
        sY = StarAllele()

        sY.name = name; sY.score = score; sY.core = copy.deepcopy(sX.core); sY.sv = 'cnv{}'.format(cn)

        self.cand.insert(0, sY)
        self.remove_star(sX)

class Person:
    def __init__(self):
        self.name = '' # sample ID
        self.gt = False # true if genotyped
        self.sv = ['', ''] # SV calls
        self.pt = '' # predicted phenotype
        self.ssr = '' # sum of squared residuals
        self.dip_cand = [] # candidate stars
        self.hap = [Haplotype(), Haplotype()]
        self.bad = False # true if QC failed
        self.af_df = None

        

    def get_status(self):
        return 'g' if self.gt else 'ng'

    def get_hap1_main(self):
        if self.gt:
            return self.hap[0].cand[0].name
        else:
            return '.'

    def get_hap2_main(self):
        if self.gt:
            return self.hap[1].cand[0].name
        else:
            return '.'

    def get_hap1_cand(self):
        return ','.join([x.name for x in self.hap[0].cand])

    def get_hap2_cand(self):
        return ','.join([x.name for x in self.hap[1].cand])

    def get_hap1_score(self):
        if self.gt:
            return str(self.hap[0].cand[0].score)
        else:
            return '.'

    def get_hap2_score(self):
        if self.gt:
            return str(self.hap[1].cand[0].score)
        else:
            return '.'

    def get_dip_score(self):
        if self.gt:
            return str(self.hap[0].cand[0].score + self.hap[1].cand[0].score)
        else:
            return '.'

    def get_phenotype(self):
        if self.gt:
            return self.pt
        else:
            return '.'

    def get_dip_sv(self):
        return ','.join(self.sv)

    def get_hap1_sv(self):
        if self.gt:
            return self.hap[0].sv
        else:
            return '.'

    def get_hap2_sv(self):
        if self.gt:
            return self.hap[1].sv
        else:
            return '.'

    def get_ssr(self):
        return self.ssr

    def get_dip_cand(self):
        return ','.join([x.name for x in self.dip_cand])

    def get_hap1_main_core(self):
        if self.gt:
            return [x.summary() for x in self.hap[0].obs if x in self.hap[0].cand[0].core]
        else:
            return '.'

    def get_hap2_main_core(self):
        if self.gt:
            return [x.summary() for x in self.hap[1].obs if x in self.hap[1].cand[0].core]
        else:
            return '.'

    def get_hap1_main_tag(self):
        if self.gt:
            return [x.summary() for x in self.hap[0].obs if x in self.hap[0].cand[0].tag]
        else:
            return '.'

    def get_hap2_main_tag(self):
        if self.gt:
            return [x.summary() for x in self.hap[1].obs if x in self.hap[1].cand[0].tag]
        else:
            return '.'


    def get_hap1_af_mean_gene(self, target_locus):
        if self.gt:
            return "{:.2f}".format(self.hap[0].af_mean_gene(target_locus.gene_start, target_locus.gene_end))
        else:
            return '.'

    def get_hap2_af_mean_gene(self, target_locus):
        if self.gt:
            return "{:.2f}".format(self.hap[1].af_mean_gene(target_locus.gene_start, target_locus.gene_end))
        else:
            return '.'

    def get_hap1_af_mean_main(self):
        if self.gt:
            return "{:.2f}".format(self.hap[0].af_mean_main)
        else:
            return '.'

    def get_hap2_af_mean_main(self):
        if self.gt:
            return "{:.2f}".format(self.hap[1].af_mean_main)
        else:
            return '.'



    def summarize(self, target_locus):
        list2str = lambda x: '.' if not x else ','.join([str(x) for x in x])

        results = [
            self.name,
            self.get_status(),
            self.get_hap1_main(),
            self.get_hap2_main(),
            self.get_hap1_cand(),
            self.get_hap2_cand(),
            self.get_hap1_score(),
            self.get_hap2_score(),
            self.get_dip_score(),
            self.get_phenotype(),
            self.get_dip_sv(),
            self.get_hap1_sv(),
            self.get_hap2_sv(),
            self.get_ssr(),
            self.get_dip_cand(),
            list2str(self.get_hap1_main_core()),
            list2str(self.get_hap2_main_core()),
            list2str(self.get_hap1_main_tag()),
            list2str(self.get_hap2_main_tag()),
            self.get_hap1_af_mean_gene(target_locus),
            self.get_hap2_af_mean_gene(target_locus),
            self.get_hap1_af_mean_main(),
            self.get_hap2_af_mean_main(),
        ]
        
        return results


def copy_vcf(original_vcf, items):
    copied_vcf = VCF()
    if 'meta' in items:
        copied_vcf.meta = copy.deepcopy(original_vcf.meta)
    if 'header' in items:
        copied_vcf.header = copy.deepcopy(original_vcf.header)
    if 'data' in items:
        copied_vcf.data = copy.deepcopy(original_vcf.data)
    return copied_vcf

def parse_region(region):
    return {'chr': region.split(':')[0].replace('chr', ''), 'start': int(region.split(':')[1].split('-')[0]), 'end': int(region.split(':')[1].split('-')[1])}

def read_vcf_simple(file):
    f = gzip.open(file, 'rt') if '.gz' in file else open(file)
    vcf = VCF()
    for line in f:
        if '##' in line:
            vcf.meta.append(line)
            continue
        fields = line.strip().split('\t')
        if fields[0] == '#CHROM':
            vcf.header = fields
            continue
        chr = fields[0].replace('chr', '')
        vcf.data.append([chr] + fields[1:])
    f.close()
    return vcf

def read_vcf_region(file, region):
    vcf = VCF()
    region_dict = parse_region(region)
    f = gzip.open(file, 'rt') if '.gz' in file else open(file)
    for line in f:
        if '##' in line:
            vcf.meta.append(line)
            continue
        fields = line.strip().split('\t')
        if fields[0] == '#CHROM':
            vcf.header = fields
            continue
        chr, pos = fields[0].replace('chr', ''), int(fields[1])
        if chr != region_dict['chr'] or pos < region_dict['start']:
            continue
        if pos > region_dict['end']:
            break
        vcf.data.append([chr] + fields[1:])
    f.close()
    return vcf

def write_vcf(vcf, file):
    with open(file, 'w') as f:
        for line in vcf.meta:
            f.write(line)
        f.write('\t'.join(vcf.header) + '\n')
        for fields in vcf.data:
            f.write('\t'.join([str(x) for x in fields]) + '\n')

def get_snp_list(target_gene, program_dir, genome_build):
    df = pd.read_table(f"{program_dir}/snp_table.tsv")
    df = df[df['gene'] == target_gene]
    def func(row):
        snp_allele = SNPAllele()
        snp_allele.n = row['sg_id']
        snp_allele.effect = row['functional_effect']
        snp_allele.pos = row[f'{genome_build}_pos']
        snp_allele.id = row['rs_id']
        snp_allele.hg = row[f'{genome_build}_allele']
        snp_allele.var = row['var_allele']
        snp_allele.wt = row['wt_allele']
        snp_allele.so = row['sequence_ontology']
        snp_allele.impact = row['variant_impact']
        snp_allele.rev = row[f'{genome_build}_revertant'] == 'revertant'
        return snp_allele
    snp_alleles = df.apply(func, axis=1).to_list()

    return snp_alleles

def get_star_dict(target_gene, snp_list, program_dir, genome_build, star_table):
    #df = pd.read_table(f"{program_dir}/star_table.tsv")
    df = pd.read_table(star_table)
    df = df[df["gene"] == target_gene]
    def func(row):
        star_allele = StarAllele()
        star_allele.name = row["name"]
        star_allele.score = row["score"]
        if row[f"{genome_build}_core"] in ['.', 'ref']:
            star_allele.core = []
        else:
            star_allele.core = copy.deepcopy([x for x in snp_list if '{}:{}>{}'.format(x.pos, x.wt, x.var) in row[f"{genome_build}_core"].split(',')])

        if row[f"{genome_build}_tag"] == '.':
            star_allele.tag = []
        else:
            star_allele.tag = copy.deepcopy([x for x in snp_list if '{}:{}>{}'.format(x.pos, x.wt, x.var) in row[f"{genome_build}_tag"].split(',')])

        if row["sv"] == '.':
            star_allele.sv = ''
        else:
            star_allele.sv = row["sv"]

        return star_allele

    star_alleles = df.apply(func, axis=1).to_list()

    star_dict = {}

    for star_allele in star_alleles:
        star_dict[star_allele.name] = star_allele


    return star_dict

def vcf2samples(vcf):
    samples = []
    for name in vcf.header[9:]:
        sample = Person()
        sample.name = name
        i = vcf.header.index(name)
        for fields in vcf.data:
            pos, rs, ref, alt, inf, fmt = int(fields[1]), fields[2], fields[3], fields[4].split(','), fields[7].split(';'), fields[8]

            if not any(['PS=D' in x for x in inf]):
                continue

            #gt = [int(x) for x in fields[i].split(":")[0].split("|")]
            gt = [x for x in fields[i].split(':')[0].split('|')]
            gt = [int(y) if y is not "." else y for y in gt]

            al = [ref] + alt
            vi_list = ['no_change'] + [x for x in inf if 'VI=' in x][0].replace('VI=', '').split(',')
            so_list = ['no_change'] + [x for x in inf if 'SO=' in x][0].replace('SO=', '').split(',')
            fe_list = ['no_change'] + [x for x in inf if 'FE=' in x][0].replace('FE=', '').split(',')

            for j in [0, 1]:

                snp = SNPAllele()
                if gt[j] != ".":
                    #print(pos,gt[j])
                    #print(gt)
                    snp.pos, snp.wt, snp.var, snp.rs, snp.het, snp.so, snp.impact, snp.effect = pos, ref, al[gt[j]], rs, gt[0] != gt[1], so_list[gt[j]], vi_list[gt[j]], fe_list[gt[j]]

                if 'AD' in fmt:
                    ad_ind=fmt.split(":").index("AD")
                    ad_list=[x for x in fields[i].split(':')[ad_ind].split(',')]
                    ad_list = [int(y) if y is not "." else 0 for y in ad_list]
                    #ad_list = [int(x) for x in fields[i].split(':')[1].split(',')]
                    #print(ad_list)

                    if not ad_list:
                         snp.ad=0; snp.td=0
                    elif gt[j] == ".":
                         snp.ad=0; snp.td=0
                    else:
                     #    print(gt[j],ad_list)
                         snp.ad = ad_list[gt[j]]; snp.td = sum(ad_list)

                sample.hap[j].obs.append(snp)
        samples.append(sample)
    return samples



def sort_star_names(names):
    def f(x):
        cn = 1
        if '*' not in x or x == '*DEL':
            n = 999
        else:
            n = int(''.join([y for y in x.split('+')[0].split('x')[0] if y.isdigit()]))
            if 'x' in x.split('+')[0]:
                cn = int(x.split('+')[0].split('x')[1])
        return (n, cn, len(x))

    return sorted(names, key = f)

def read_sv_table(fn):
    result = {}

    with open(fn) as f:
        header = next(f).strip().split()
        for line in f:
            #print(line)
            fields = line.strip().split()
            gene = fields[0]
            name = fields[2]

            if gene not in result:
                result[gene] = {}

            result[gene][name] = dict(zip(header, fields))

    return result

def cyp2a6(sample, star_dict):
    call_sv1(sample, "gc_e1e4", "*34", star_dict)
    call_sv1(sample, "gc_e1e2", "*12", star_dict)
    call_tandem(sample, "dup2", "*2", "*S1", star_dict)
    call_tandem(sample, "dup2", "*1", "*S1", star_dict)
    call_tandem(sample, "gc_e9", "*1", "*S2", star_dict)
    call_tandem(sample, "dup7", "*1", "*S3", star_dict)
    call_tandem(sample, "dup7b", "*1", "*S6", star_dict)
    cyp2a6_svcomb(sample, star_dict)

def cyp2b6(sample, star_dict):
    call_sv1(sample, "gc_i4e9", "*29", star_dict)

def cyp2d6(sample, star_dict):
    call_tandem(sample, 'del1', '*S2', '*1', star_dict, ordered = True)
    call_tandem(sample, 'gc_i1e9_5', '*68x5', '*4', star_dict)
    call_tandem(sample, 'gc_i1e9', '*68', '*4', star_dict)
    call_tandem(sample, 'gc_i1e9', '*S1', '*1', star_dict)
    call_tandem(sample, 'gc_e9', '*4N', '*4', star_dict)
    call_tandem(sample, 'gc_e9_3', '*36x3', '*10', star_dict)
    call_tandem(sample, 'gc_e9_7', '*36x7', '*10', star_dict)
    call_tandem(sample, 'gc_e9', '*36', '*10', star_dict)
    call_tandem(sample, 'gc_e9', '*83', '*2', star_dict)
    call_tandem(sample, 'gc_7to6_i4', '*13A', '*2', star_dict)
    call_sv1(sample, 'gc_e1e7', '*13C', star_dict)
    call_sv1(sample, 'gc_7to6_i1', '*13B', star_dict)
    cyp2d6_svcomb(sample, star_dict)

def cyp2e1(sample, star_dict):
    call_sv1(sample, "dup_e7e9", "*S1", star_dict)

def gstm1(sample, star_dict):
    gstm1_svcomb(sample, star_dict)

def gstt1(sample, star_dict):
    gstt1_svcomb(sample, star_dict)

def slc22a2(sample, star_dict):
    call_sv1(sample, "del_i9", "*S1", star_dict)
    call_sv1(sample, "del_e11", "*S2", star_dict)

def slco1b1(sample, star_dict):
    call_sv1(sample, "dup1", "*S3", star_dict)

def ugt1a4(sample, star_dict):
    call_sv1(sample, 'del_i1', '*S1', star_dict)
    call_sv1(sample, 'del2', '*S2', star_dict)

def ugt2b15(sample, star_dict):
    call_sv1(sample, "del_i3e6", "*S1", star_dict)

def ugt2b17(sample, star_dict):
    ugt2b17_svcomb(sample, star_dict)

def new_tandem(sv, star_names, star_dict):
    stars = [star_dict[x] for x in star_names]
    score = sum([x.score for x in stars])
    core = list(set([x for y in [z.core for z in stars] for x in y]))
    tandem = StarAllele()
    tandem.name, tandem.score, tandem.core, tandem.sv = '+'.join(star_names), score, copy.deepcopy(core), sv
    return tandem

def new_dup(sX, cnv):
    times = int(cnv.replace('cnv', ''))
    score = sX.score * times
    dup = StarAllele()
    dup.name, dup.score, dup.core, dup.sv = sX.name + 'x' + str(times), score, copy.deepcopy(sX.core), cnv
    return dup

def remove_select(hap, stars):
    '''
    args:
        hap (list of Stars)
        stars (list of Stars)
    '''
    for i in reversed(range(len(hap))):
        if hap[i].name in [x.name for x in stars]:
            del hap[i]

def remove_sv(hap, l = []):
    for i in reversed(range(len(hap))):
        if hap[i].sv:
            if hap[i].name in l:
                continue
            else:
                del hap[i]

def which_has(sample, stars):
    '''
    args:
        sample (Person)
        stars (list of str)
    Returns:
        i (int)
    '''
    h1 = set(stars).issubset([x.name for x in sample.hap[0].cand])
    h2 = set(stars).issubset([x.name for x in sample.hap[1].cand])
    if h1 and h2:
        i = 3
    elif h1 and not h2:
        i = 1
    elif not h1 and h2:
        i = 2
    else:
        i = 0
    return i

def which_severe(sample):
    '''
    args:
        sample (Person)
    Returns:
        i (int)
    '''
    h1 = copy.deepcopy(sample.hap[0].cand)
    h2 = copy.deepcopy(sample.hap[1].cand)
    remove_sv(h1)
    remove_sv(h2)
    if h1[0].name == h2[0].name:
        i = 3
    else:
        h3 = sorted([h1[0], h2[0]], key = lambda x: x.rank)
        if h3[0].name == h1[0].name:
            i = 1
        else:
            i = 2
    return i

def call_sv1(sample, sv, x, star_dict):
    '''
    This function calls a sample's final genotype if the sample has only one SV.
    If a SV-carrying allele is in LD with other alleles, the funtion takes those other alleles as input.
    x = the name of the star allele with the SV
    '''
    if sample.gt or sample.sv != ["no_sv", sv]:
        return
    if star_dict[x].core:
        h1 = set(star_dict[x].core).issubset(sample.hap[0].obs)
        h2 = set(star_dict[x].core).issubset(sample.hap[1].obs)
    else:
        h1 = True
        h2 = True
    if not h1 and not h2:
        return
    elif h1 and not h2:
        i, j = 1, 0
    elif not h1 and h2:
        i, j = 0, 1
    else:
        l1 = copy.deepcopy(sample.hap[0].cand)
        l2 = copy.deepcopy(sample.hap[1].cand)
        remove_sv(l1)
        remove_sv(l2)
        if l1[0].name == l2[0].name:
            i, j = 0, 1
        else:
            l3 = sorted([l1[0], l2[0]], key = lambda x: x.rank)
            if l3[0].name == l1[0].name:
                i, j = 0, 1
            else:
                i, j = 1, 0
    remove_sv(sample.hap[i].cand)
    remove_sv(sample.hap[j].cand, [x])
    sample.gt = True

def call_tandem(sample, sv, x, y, star_dict, ordered=False):
    """
    Calls a tandem duplication allele containing two gene copies (e.g., CYP2D6*36+*10)
    x = the name of the 1st star allele in the tandem (e.g., '*36')
    y = the name of the 2nd star allele in the tandem (e.g., '*10')
    """
    if sample.gt or sample.sv != ["no_sv", sv]:
        return
    h1 = set([x, y]).issubset([_.name for _ in sample.hap[0].cand])
    h2 = set([x, y]).issubset([_.name for _ in sample.hap[1].cand])
    if not h1 and not h2:
        return
    elif h1 and not h2:
        i, j = 0, 1
    elif not h1 and h2:
        i, j = 1, 0
    else:
        l1 = copy.deepcopy(sample.hap[0].cand)
        l2 = copy.deepcopy(sample.hap[1].cand)
        remove_sv(l1)
        remove_sv(l2)
        if l1[0] == l2[0]:
            i, j = 1, 0
        else:
            l3 = sorted([l1[0], l2[0]], key = lambda x: x.rank)
            if l3[0] == l1[0]:
                i, j = 1, 0
            else:
                i, j = 0, 1
    sX = star_dict[x]
    sY = star_dict[y]

    # find SNPs shared by both star alleles
    overlap = []
    for snp in sX.core:
        if snp in sY.core:
            overlap.append(snp)

    # return if allele fraction in any of the shared SNPs is less than 0.4
    for snp in overlap:
        if [x for x in sample.hap[i].obs if x == snp][0].af < 0.4:
            return

    tandem = new_tandem(sv, [sX.name, sY.name], star_dict)
    sample.hap[i].cand.insert(0, tandem)
    remove_select(sample.hap[i].cand, [sX, sY])
    remove_sv(sample.hap[i].cand, [tandem.name])
    remove_sv(sample.hap[j].cand)
    if ordered:
        sample.hap[i].cand.sort(key = lambda x: x.rank)
    sample.gt = True

def call_cnv3(sample, target_locus):
    if sample.gt or sample.sv != ["no_sv", 'cnv2']:
        return
    sX = sample.hap[0].cand[0]
    sY = sample.hap[1].cand[0]

    # Simplest case (e.g., *2/*2x2)
    if sX.name == sY.name:
        sample.hap[0].add_dup(2)
        sample.gt = True
        return

    sX_gene = sample.hap[0].af_mean_gene(target_locus.gene_start, target_locus.gene_end)
    sY_gene = sample.hap[1].af_mean_gene(target_locus.gene_start, target_locus.gene_end)
    if sX_gene == -1:
        sX_gene = 1 - sY_gene
    if sY_gene == -1:
        sY_gene = 1 - sX_gene
    diff_gene = sY_gene - sX_gene

    sX_main = sample.hap[0].af_mean_main
    sY_main = sample.hap[1].af_mean_main
    if sX_main == -1:
         sX_main = 1 - sY_main
    if sY_main == -1:
        sY_main = 1 - sX_main
    diff_main = sY_main - sX_main


    fit_maf1, _ = sample.hap[0].fit_data(3, target_locus.gene_start, target_locus.gene_end)
    fit_maf2, _ = sample.hap[1].fit_data(3, target_locus.gene_start, target_locus.gene_end)


    means = [round(fit_maf1, 2), round(fit_maf2, 2)]


    f = lambda a, b: ((a == b) & (a == 0)) | (a * b > 0)

    if f(diff_gene, diff_main):
        if means == [0.33, 0.67]:
            sample.hap[1].add_dup(2)
            sample.gt = True
        elif means == [0.67, 0.33]:
            sample.hap[0].add_dup(2)
            sample.gt = True
    else:

        if abs(diff_main) > abs(diff_gene):

            if sY_main > sX_main:
                sample.hap[1].add_dup(2)
                sample.gt = True

            else:
                sample.hap[0].add_dup(2)
                sample.gt = True

        else:
            if means == [0.33, 0.67]:
                sample.hap[1].add_dup(2)
                sample.gt = True
            elif means == [0.67, 0.33]:
                sample.hap[0].add_dup(2)
                sample.gt = True

def call_cnv_plus(sample, target_locus):
    '''This function calls a final genotype with CN > 3 gene copies.'''
    if sample.gt:
        return
    if 'cnv' not in sample.sv[1]:
        return
    if sample.sv[0] == 'no_sv' and (sample.sv[1] == 'cnv0' or sample.sv[1] == 'cnv2'):
        return
    if sample.sv[0] == 'no_sv' and 'cnv' in sample.sv[1]:
        total_cn = int(sample.sv[1].replace('cnv', '')) + 1
    elif 'cnv' in sample.sv[0] and 'cnv' in sample.sv[1]:
        total_cn = int(sample.sv[0].replace('cnv', '')) + int(sample.sv[1].replace('cnv', ''))
        if total_cn < 4:
            return

    # allele fraction profile is not informative -- i.e. it's empty
    if sample.hap[0].af_mean_gene == -1:
        sample.hap[0].add_dup(total_cn - 1)
        sample.gt = True
        return

    fit_maf, fit_cn = sample.hap[0].fit_data(total_cn, target_locus.gene_start, target_locus.gene_end)
    sample.hap[0].add_dup(fit_cn)
    sample.hap[1].add_dup(total_cn - fit_cn)
    sample.gt = True

def cyp2a6_svcomb(sample, star_dict):
    if sample.gt:
        return
    gt = []
    for sv in sample.sv:
        if sv == 'cnv0':
            gt.append(star_dict['*4'])
        elif sv == 'cnv2':
            gt.append(new_dup(star_dict['*1'], sv))
        elif sv == 'gc_e1e2':
            gt.append(star_dict['*12'])
        elif sv == 'gc_e1e4':
            gt.append(star_dict['*34'])
        elif sv == 'dup2':
            gt.append(new_tandem(sv, ['*1', '*S1'], star_dict))
        elif sv == 'gc_e9':
            gt.append(new_tandem(sv, ['*1', '*S2'], star_dict))
        elif sv == 'dup7':
            gt.append(new_tandem(sv, ['*1', '*S3'], star_dict))
        elif sv == 'dup7x2':
            gt.append(new_tandem(sv, ['*1', '*S3', '*S3'], star_dict))
        elif sv == 'dup7b':
            gt.append(new_tandem(sv, ['*1', '*S6'], star_dict))
    if len(gt) == 2:
        sample.hap[0].cand = [gt[0]]
        sample.hap[1].cand = [gt[1]]
        sample.gt = True

def    cyp2d6_svcomb(sample, star_dict):
    if sample.gt:
        return
    gt = []
    for sv in sample.sv:
        if sv == 'cnv0':
            gt.append(star_dict['*5'])
        elif sv == 'gc_i1e9' and which_has(sample, ['*68', '*4']):
            gt.append(new_tandem(sv, ['*68', '*4'], star_dict))
        elif sv == 'gc_i1e9' and which_has(sample, ['*S1', '*1']):
            gt.append(new_tandem(sv, ['*S1', '*1'], star_dict))
        elif sv == 'gc_e9' and which_has(sample, ['*4N', '*4']):
            gt.append(new_tandem(sv, ['*4N', '*4'], star_dict))
        elif sv == 'gc_e9' and which_has(sample, ['*36', '*10']):
            gt.append(new_tandem(sv, ['*36', '*10'], star_dict))
        elif sv == 'gc_e9' and which_has(sample, ['*83', '*2']):
            gt.append(new_tandem(sv, ['*83', '*2'], star_dict))
        elif sv == 'gc_7to6_i4' and which_has(sample, ['13A', '*2']):
            gt.append(new_tandem(sv, ['13A', '*2'], star_dict))
        elif sv == 'gc_7to6_i1':
            gt.append(star_dict['*13B'])
        elif sv == 'gc_e1e7':
            gt.append(star_dict['*13C'])
    cnv = None
    for sv in sample.sv:
        if 'cnv' in sv and sv != 'cnv0':
            cnv = sv
    if cnv:
        if '*68+*4' in [x.name for x in gt]:
            svcomb_tandem_cnv(gt, sample, ['*68', '*4'], cnv)
        elif '*S1+*1' in [x.name for x in gt]:
            svcomb_tandem_cnv(gt, sample, ['*S1', '*1'], cnv)
        elif '*4N+*4' in [x.name for x in gt]:
            svcomb_tandem_cnv(gt, sample, ['*4N', '*4'], cnv)
        elif '*36+*10' in [x.name for x in gt]:
            svcomb_tandem_cnv(gt, sample, ['*36', '*10'], cnv)
        elif '*83+*2' in [x.name for x in gt]:
            svcomb_tandem_cnv(gt, sample, ['*83', '*2'], cnv)
        elif '*13A+*2' in [x.name for x in gt]:
            svcomb_tandem_cnv(gt, sample, ['*13A', '*2'], cnv)
        elif '*13B' in [x.name for x in gt]:
            svcomb_sv1_cnv(gt, sample, '*13B', cnv)
    if len(gt) == 2:
        sample.hap[0].cand = [gt[0]]
        sample.hap[1].cand = [gt[1]]
        sample.gt = True

def svcomb_sv1_cnv(gt, sample, sX_name, cnv):
    i = which_has(sample, [sX_name])
    if not i:
        return
    if i != 3:
        j = {0: 1, 1: 0}[i - 1]
    elif i == 3:
        j = 0
    l = copy.deepcopy(sample.hap[j].cand)
    remove_sv(l)
    sY = new_dup(l[0], cnv)
    gt.insert(j, sY)

def svcomb_tandem_cnv(gt, sample, tandem, cnv):
    i = which_has(sample, [tandem[0], tandem[1]])
    if not i:
        return
    if i != 3:
        j = {0: 1, 1: 0}[i - 1]
        l = copy.deepcopy(sample.hap[j].cand)
        remove_sv(l)
        sX = new_dup(l[0], cnv)
        gt.insert(j, sX)
    elif i == 3:
        for x in sample.hap[0].cand:
            if x.name == tandem[1]:
                gt.append(new_dup(x, cnv))
                break

def gstt1_svcomb(sample, star_dict):
    if sample.gt:
        return
    gt = []
    for sv in sample.sv:
        if sv == 'cnv0':
            gt.append(star_dictstar_dict["*2"])
    if len(gt) == 2:
        sample.hap[0].cand = [gt[0]]
        sample.hap[1].cand = [gt[1]]
        sample.gt = True

def gstm1_svcomb(sample, star_dict):
    if sample.gt:
        return
    gt = []
    for sv in sample.sv:
        if sv == 'cnv0':
            gt.append(star_dict["*2"])
    if len(gt) == 2:
        sample.hap[0].cand = [gt[0]]
        sample.hap[1].cand = [gt[1]]
        sample.gt = True

def ugt2b17_svcomb(sample, star_dict):
    if sample.gt:
        return
    gt = []
    for sv in sample.sv:
        if sv == 'cnv0':
            gt.append(star_dict["*2"])
    if len(gt) == 2:
        sample.hap[0].cand = [gt[0]]
        sample.hap[1].cand = [gt[1]]
        sample.gt = True

def partial(lst, query):
    for i,s in enumerate(lst):
         if query in s:
              return({'sample':i, 'gt':s})
         else:
              pass

def check_impute_calls(phaseme_vcf, phased_vcf,genome_build,program_dir):
   starTableName=program_dir+"/"+genome_build+"_star.txt"
   starTable=pd.read_table(starTableName)
   #print("starTable",starTable)
   check=[]
   for i,line in enumerate(phaseme_vcf.data):
      x=partial(line,"./.")
      if x is not None:
          x['pos']=line[1]
          x['ref']=line[3]
          x['alt']=line[4]
          check.append(x)      

   check2=[]
   #print("check:",check)
   phased_df=pd.DataFrame(phased_vcf.data, columns=phased_vcf.header)
   for l in check:
    #   print(l)
    #  print(starTable.head())
       tmp=starTable.loc[starTable['pos']==int(l['pos']),].reset_index()
       if tmp.shape[0]==0:
           next
       else:
           gt_phased=phased_df.loc[phased_df['POS']==l['pos'],].values.flatten().tolist()[l['sample']]
           if "0|0" in gt_phased:
               next
           elif tmp.shape[0]==1 and tmp.loc[0,'ref']==l['ref'] and tmp.loc[0,'alt']==l['alt']:
               l['sample']=phased_df.columns[l['sample']]
               l['gt_beagle']=gt_phased
               l['star']=tmp['name'].to_list()
               check2.append(l)
           else:
               tmp1=tmp[['pos','ref','alt']].drop_duplicates()
               for i,row in tmp1.iterrows():
                  if row['ref']==l['ref'] and row['alt']==l['alt']:
                     l['sample']=phased_df.columns[l['sample']]
                     l['gt_beagle']=gt_phased
                     stars=tmp.loc[((tmp['ref']==row['ref']) & (tmp['alt']==row['alt'])),'name'].to_list()
                     l['star']=stars
                     check2.append(l)
   #print("check2:",check2)
   return(check2)

def assess_vcf(input_vcf,
               gdf_file,
               data_type,
               logger):
    # Start logging.
    logger.info(LINEBREAK)
    logger.info("Step 2/10: Assessing input VCF file...")

    if gdf_file == None:
        pass
     #   print("No GDF file")
    elif not gdf_file:
        # check whether the sample list is identical between VCF and GDF
        with open(gdf_file) as f:
            gdf_samples = [x.replace('Depth_for_', '') for x in f.readline().strip().split('\t')[3:]]
        vcf_samples = input_vcf.header[9:]
        if len(gdf_samples) != len(vcf_samples):
            raise TypeError(f'The sample size differs between the VCF file (N={len(vcf_samples)}) and the GDF file (N={len(gdf_samples)})')
        if len(set(gdf_samples + vcf_samples)) != len(vcf_samples):
            raise TypeError(f'Two different sets of samples were detected from the VCF file and the GDF file')
        for i in range(len(vcf_samples)):
            if vcf_samples[i] != gdf_samples[i]:
                raise TypeError(f"The order of samples differs between the VCF file ('{vcf_samples[i]}') and the GDF file ('{gdf_samples[i]}') at sample index {i}")

        # make sure the sample size > 1 when using TS data
        if len(vcf_samples) < 5 and data_type == 'ts':
            raise TypeError(f"Genotyping with TS data requires at least five samples (the current sample size is {len(vcf_samples)})")

    log_dict = {'row': 0, 'AD': 0, 'phased': 0, 'unphased': 0, 'both': 0}

    for fields in input_vcf.data:
        chrom = fields[0].replace('chr', '')
        pos = int(fields[1])
        format = fields[8].split(':')

        # Check GT field
        if 'GT' not in format:
            ValueError('GT field not found [{}]'.format(pos))

        # Check AD field
        if 'AD' in format:
            log_dict['AD'] += 1

        # Check phasing status
        def f(x):
            gt = x.split(':')[format.index('GT')]
            if '/' in gt:
                return '/'
            elif '|' in gt:
                return '|'
            else:
                if chrom == 'X' or chrom == 'Y':
                    return
                else:
                    raise ValueError('Genotype separator not found for autosomal chromosome chr{}:{} GT=[{}]'.format(chrom, pos, gt))

        separators = set([f(x) for x in fields[9:] if f(x)])

        log_dict['row'] += 1
        if len(separators) == 1:
            if '|' in separators:
                log_dict['phased'] += 1
            else:
                log_dict['unphased'] += 1
        else:
            log_dict['both'] += 1

    # Check if input VCF is empty
    if log_dict['row'] == 0:
        vcf_empty = True
    else:
        vcf_empty = False

    # Determine AD mode
    if log_dict['row'] > 0 and log_dict['AD'] / log_dict['row'] > 0.8:
        vcf_ad = True
    else:
        vcf_ad = False

    # Determine phasing mode
    if log_dict['phased'] == log_dict['row']:
        vcf_sep = '|'
    elif log_dict['unphased'] == log_dict['row']:
        vcf_sep = '/'
    else:
        vcf_sep = 'b'

    logger.info(f"Samples total: {len(input_vcf.header[9:])}")
    logger.info(f"Markers total: {log_dict['row']}")
    logger.info(f"Markers with allelic depth: {log_dict['AD']}")
    logger.info(f"Markers unphased: {log_dict['unphased']}")
    logger.info(f"Markers phased: {log_dict['phased']}")
    logger.info(f"Markers partially phased: {log_dict['both']}")
    logger.info("Finished assessing input VCF file")
    logger.info(LINEBREAK)

    return vcf_ad, vcf_sep, vcf_empty

def process_vcf(input_vcf, vcf_ad, vcf_sep, logger,gb):
    start=time.process_time()
    # Start logging.
    logger.info(LINEBREAK)
    logger.info("Step 3/10: Processing input VCF...")

    log_dict = {'IA': 0, 'allelic_imbalance': 0, 's50': 0}
    processed_vcf = copy_vcf(input_vcf, ['header'])
    processed_vcf.meta = [
        '##fileformat=VCFv4.2\n',
        '##reference={}\n'.format(gb),
        '##INFO=<ID=FE,Number=A,Type=String,Description="Functional Effect">\n',
        '##INFO=<ID=PS,Number=1,Type=String,Description="Phasing Status (A, in preparation; B1, ready for phasing as is; B2, ready for phasing after conformation to reference VCF; C1, excluded from phasing because marker is absent in reference VCF; C2, excluded from phasing because marker has different REF allele; C3, excluded from phasing because marker has no overlapping ALT alleles; D1, statistically phased; D2, manually phased with certainty; D3, manually phased without certainty; D4, already phased; D5, manually phased by extension; E, omitted during statistical phasing)">\n',
        '##INFO=<ID=RV,Number=A,Type=String,Description="Reverting Variation">\n',
        '##INFO=<ID=SO,Number=A,Type=String,Description="=Sequence Ontology">\n',
        '##INFO=<ID=VI,Number=A,Type=String,Description="Variant Impact">\n',
        '##FILTER=<ID=IA,Description="Invalid Allele">\n',
        '##FILTER=<ID=s50,Description="Less than 50% of samples have data">\n',
        '##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths for the ref and alt alleles in the order listed">\n',
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
        '##FORMAT=<ID=HE,Number=4,Type=Integer,Description="Matching scores computed by the phase-by-extension algorithm">\n'
    ]

    for fields in input_vcf.data:
        chr, ref, alt, fmt, flt, inf = fields[0].replace('chr', ''), fields[3], fields[4].split(','), fields[8].split(':'), [], ['PS=D4'] if vcf_sep == '|' else ['PS=A']

        def f(x):
            gt_field = x.split(':')[fmt.index('GT')]

            # Determine genotype separator
            if '/' in gt_field:
                gt_sep = '/'
            elif '|' in gt_field:
                gt_sep = '|'
            else:
                gt_sep = ''

            # Unphase genotype if input VCF is partially phased
            if vcf_sep == 'b' and gt_sep == '|':
                if gt_field == '.|.':
                    gt_field = './.'
                else:
                    gt_field = '/'.join(sorted(gt_field.split('|'), key = lambda x: int(x)))

            # Conform genotype for sex chromosomes if necessary
            if not gt_sep and (chr == 'X' or chr == 'Y'):
                gt_field = '0|' + gt_field if vcf_sep == '|' else '0/' + gt_field

            # Get AD field information
            if vcf_ad:
                ad_field = x.split(':')[fmt.index('AD')]

                if gt_field == './.':
                    ad_field = ','.join(['0'] * (len(alt) + 1))

                # Some tools such as Genalice produce variable-length AD fields (i.e., single AD value if sample is homozygous for ALT allele)
                if len(ad_field.split(',')) == 1 and gt_sep and gt_field.split(gt_sep)[0] == gt_field.split(gt_sep)[1] and gt_field.split(gt_sep)[0] != '0' and len(alt) == 1:
                    if 'DP' in fmt:
                        dp_field = x.split(':')[fmt.index('DP')]
                        ad_field = str(int(dp_field) - int(ad_field)) + ',' + ad_field
                    else:
                        ad_field = '0,' + ad_field
                ad_field = ':' + ad_field

            else:
                ad_field = ''

            return gt_field + ad_field

        fields[9:] = [f(x) for x in fields[9:]]

        # Check invalid allele
        if ref == 'I' or ref == 'D' or '.' in alt or 'D' in alt:
            flt.append('IA')
            log_dict['IA'] += 1

        # Check high missingness
        if ['.' in x.split(':')[0] for x in fields[9:]].count(True) / len(fields[9:]) > 0.5:
            flt.append('s50')
            log_dict['s50'] += 1

        # Define quick function 3
        def qf3(x):
            gt = x.split(':')[0]
            if '.' in gt:
                return False
            if '|' in gt:
                gt = gt.split('|')
            else:
                gt = gt.split('/')
            if gt[0] == gt[1]:
                return False
            if "." in [y for y in x.split(':')[1].split(',')]:
                ad = [y for y in x.split(':')[1].split(',')]
                ad=list(filter(('.').__ne__, ad))
                ad = [int(y) for y in ad]
            else:
                ad = [int(y) for y in x.split(':')[1].split(',')]
            if sum(ad) == 0:
                return False
            return max(ad) / sum(ad)

        # Check allelic imbalance
        if vcf_ad:
            ratios = [qf3(x) for x in fields[9:] if qf3(x)]
            if ratios:
                median = statistics.median(ratios)
                if median > 0.8 or median < 0.2:
                    log_dict['allelic_imbalance'] += 1

        fields[5] = '.'
        fields[6] = ';'.join(flt) if flt else 'PASS'
        fields[7] = ';'.join(inf)
        fields[8] = 'GT:AD' if vcf_ad else 'GT'
        processed_vcf.data.append(fields)

    logger.info(f"Markers with allelic imbalance: {log_dict['allelic_imbalance']}")
    logger.info(f"Markers with high missingness: {log_dict['s50']}")
    logger.info(f"Markers with invalid allele: {log_dict['IA']}")
    logger.info("Finished processing input VCF")
    #logger.info(f"Processed VCF: {time.process_time() - start}")
    logger.info(LINEBREAK)

    #print("Processed VCF",time.process_time() - start)

    return processed_vcf

def adjust_vcf(processed_vcf,star_dict):
    start=time.process_time()
    '''
    Conform multiallelic loci containing more than one indels to the star allele table
    For example, the UGT1A1 gene has such locus which defines three star alleles: *28 (234668879:CAT>CATAT), *36 (234668879:CAT>C), and *37 (234668879:CAT>CATATAT)
    This locus can be represented in a VCF file in many different ways as shown below

    File        Star Alleles    POS            REF        ALT
    1st VCF        *28                234668879    C        CAT
    2nd VCF        *36                234668879    CAT        C
    3rd VCF        *37                234668879    C        CATAT
    4th VCF        *28,*36            234668879    CAT        CATAT,C
    5th VCF        *28,*37            234668879    C        CAT,CATAT
    6th VCF     *28,*36,*37        234668879    CAT        CATAT,C,CATATAT
    Ref VCF        *28                234668879    C        CAT    

    The following step would conform all these records to 234668879:CAT>CATAT,C,CATATAT
    '''

    adjusted_vcf = copy_vcf(processed_vcf, ['meta', 'header', 'data'])

    for i in range(len(adjusted_vcf.data)):
        fields = adjusted_vcf.data[i]
        pos = fields[1]
        ref = fields[3]
        alt = fields[4].split(',')

        # skip if the locus does not have an indel
        if len(ref) == 1 and all([len(x) == 1 for x in alt]):
            continue

        # skip if the locus is not used to define a star allele
        star_list = []
        for name, star in star_dict.items():
            pos=int(pos)
            if pos in [x.pos for x in star.core]:
                star_list.append(star)
        if not star_list:
            continue

        # skip if the locus is already formatted properly
        bool_list = []
        for x in alt:
            is_found = False
            for star in star_list:
                if f'{pos}:{ref}>{x}' in [f'{x.pos}:{x.hg}>{x.var}' for x in star.core]:
                    is_found = True
                    break
            bool_list.append(is_found)
        if all(bool_list):
            continue

        # try to find new sequences for REF and ALT
        snp_list = []
        for star in star_list:
            snp_list += star.core
        snp_list = list(set(snp_list))
        new_ref = []
        new_alt = []
        for x in alt:
            result = [y for y in snp_list if y.pos == pos and (len(y.hg) - len(ref) == len(y.var) - len(x))]
            if result:
                new_ref.append(result[0].hg)
                new_alt.append(result[0].var)
            else:
                new_alt.append(x)
        new_ref = list(set(new_ref))
        
        # skip if there are more than one sequence for REF
        if len(new_ref) > 1:
            continue
                    
        
        # skip if the number of sequences does not match for ALT -- WHY? 
        if len(new_alt) != len(alt):
            continue

        # update the adjusted record
        if len(new_ref)>0:
            fields[3] = new_ref[0]
        fields[4] = ','.join(new_alt)
        

    #print("Adjusted VCF:",time.process_time()-start)
    return adjusted_vcf

def conform_vcf(adjusted_vcf, ref_vcf, snp_list,  vcf_ad, vcf_empty, vcf_sep, vcf_only, logger):
    start=time.process_time()
    # Start logging.
    logger.info(LINEBREAK)
    logger.info("Step 4/10: Conforming input VCF...")

    log_dict = {'status': 'Skipped because input VCF is empty ', 'row': 0, 'filtered': 0, 'PS=B1': 0, 'PS=B2': 0, 'PS=C1': 0, 'PS=C2': 0, 'PS=C3': 0}

    def write_log():
        logger.info(f"Markers total: {len(adjusted_vcf.data)}")
        logger.info(f"Markers filtered: {log_dict['filtered']}")
        logger.info(f"Markers remaining: {log_dict['row']}")
        logger.info(f"Markers total: {len(adjusted_vcf.data)}")
        logger.info(f"Markers filtered: {log_dict['filtered']}")
        logger.info(f"Markers remaining: {log_dict['row']}")
        logger.info(f"Markers phasable: {log_dict['PS=B1'] + log_dict['PS=B2']}")
        logger.info(f"Markers ready: {log_dict['PS=B1']}")
        logger.info(f"Markers conformed: {log_dict['PS=B2']}")
        logger.info(f"Markers unphasable: {log_dict['PS=C1'] + log_dict['PS=C2'] + log_dict['PS=C3']}")
        logger.info(f"Markers absent in reference VCF: {log_dict['PS=C1']}")
        logger.info(f"Markers with different REF allele: {log_dict['PS=C2']}")
        logger.info(f"Markers with no overlapping ALT alleles: {log_dict['PS=C3']}")
        logger.info("Finished conforming input VCF")
        logger.info(LINEBREAK)

    if vcf_empty:
        write_log()
        return adjusted_vcf

    if vcf_sep == '|':
        log_dict['status'] = 'Skipped because input VCF is already fully phased'
        write_log()
        return adjusted_vcf

    log_dict['status'] = 'Completed'
    conformed_vcf = copy_vcf(adjusted_vcf, ['meta', 'header'])
    ref_pos=[int(x[1]) for x in ref_vcf.data]

    for fields1 in adjusted_vcf.data:
        chr1, pos1, ref1, alt1, flt1, inf1 = fields1[0], int(fields1[1]), fields1[3], fields1[4].split(','), fields1[6], fields1[7].split(';')
        if vcf_only:
            if flt1 !='PASS' and flt1 != 's50':
                log_dict['filtered'] += 1
                continue
        else:
            if flt1 != 'PASS':
                log_dict['filtered'] += 1
                continue

        log_dict['row'] += 1

        is_found = False
        if (pos1 not in [x.pos for x in snp_list]) and (pos1 not in ref_pos):
            #print(fields1[0:8])
            continue

        if (pos1 in [x.pos for x in snp_list]) or (pos1 in ref_pos):
            for i in range(len(ref_vcf.data)):
                fields2 = ref_vcf.data[i]
                chr2, pos2, ref2, alt2 = fields2[0], int(fields2[1]), fields2[3], fields2[4].split(',')

                # Keep looking if not found
                if chr1 != chr2 or pos1 != pos2:
                    continue

                # Although found, will not be phased
                if ref1 != ref2:
                    fields3 = ref_vcf.data[i + 1]
                    pos3 = int(fields3[1])
                    ref3 = fields3[3]

                    # Check if the next line matches
                    if pos1 == pos3 and ref1 == ref3:
                        continue
                    else:
                        log_dict['PS=C2'] += 1
                        inf1 = ['PS=C2' if x == 'PS=A' else x for x in inf1]
                        is_found = True
                        break

                # There are no overlapping ALT alleles
                if not list(set(alt1) & set(alt2)):
                    is_found = True
                    log_dict['PS=C3'] += 1
                    inf1 = ['PS=C3' if x == 'PS=A' else x for x in inf1]
                    break

                # Found and perfectly matched, no need to conform
                if alt1 == alt2 or len(set(alt1) & set(alt2)) == 0:
                    is_found = True
                    log_dict['PS=B1'] += 1
                    inf1 = ['PS=B1' if x == 'PS=A' else x for x in inf1]
                    break

                # Although found, missing one or more ALT alleles
                if set(alt1).issubset(alt2) and len(alt2) > len(alt1):
                    diff = len(alt2) - len(alt1)
                    for allele in alt2:
                        if allele not in alt1:
                            alt1.append(allele)
                    if vcf_ad:
                        fields1[9:] = [x + ',0' * diff for x in fields1[9:]]

                # Although same ALT alleles, wrong order
                if set(alt1) == set(alt2):
                    is_found = True
                    mapping = {0: 0}
                    for i in range(len(alt1)):
                        mapping[i + 1] = alt2.index(alt1[i]) + 1
                    fields1[4] = ','.join(alt2) # update ALT alleles
                    log_dict['PS=B2'] += 1
                    inf1 = ['PS=B2' if x == 'PS=A' else x for x in inf1]

                    def f(x):
                        gt = x.split(':')[0].split('/')
                        for i in [0, 1]:
                            if gt[i] != '0' and gt[i] != '.':
                                gt[i] = str(mapping[int(gt[i])])
                        if not vcf_ad:
                            return '/'.join(gt)
                        ad1 = x.split(':')[1].split(',')
                        ad2 = [0 for y in ad1]
                        for i in range(len(ad2)):
                            ad2[mapping[i]] = ad1[i]
                        return '/'.join(gt) + ':' + ','.join(ad2)

                    fields1[9:] = [f(x) for x in fields1[9:]]
                    break
            
                # Completely missing data with ad=0,0

            if not is_found:
                inf1 = ['PS=C1' if x == 'PS=A' else x for x in inf1]
                log_dict['PS=C1'] += 1
        fields1[7] = ';'.join(inf1)
        conformed_vcf.data.append(fields1)

    logger.info(f"Markers total: {len(adjusted_vcf.data)}")
    logger.info(f"Markers filtered: {log_dict['filtered']}")
    logger.info(f"Markers remaining: {log_dict['row']}")
    logger.info(f"Markers phasable: {log_dict['PS=B1'] + log_dict['PS=B2']}")
    logger.info(f"Markers ready: {log_dict['PS=B1']}")
    logger.info(f"Markers conformed: {log_dict['PS=B2']}")
    logger.info(f"Markers unphasable: {log_dict['PS=C1'] + log_dict['PS=C2'] + log_dict['PS=C3']}")
    logger.info(f"Markers absent in reference VCF: {log_dict['PS=C1']}")
    logger.info(f"Markers with different REF allele: {log_dict['PS=C2']}")
    logger.info(f"Markers with no overlapping ALT alleles: {log_dict['PS=C3']}")
    logger.info("Finished conforming input VCF")
    #logger.info(f"Conformed VCF: {time.process_time() - start}")
    logger.info(LINEBREAK)

    return conformed_vcf

def phase_vcf(conformed_vcf, program_dir, target_region, snp_list, impute_ungenotyped_markers, output_dir, hap_panel, vcf_ad, vcf_empty, vcf_sep, logger,genome_build):
    start=time.process_time()
    # Start logging.
    logger.info(LINEBREAK)
    logger.info("Step 5/10: Statistically phasing input VCF...")

    def write_log():
        logger.info(f"Markers attempted: {log_dict['attempted']}")
        logger.info(f"Markers phased: {log_dict['PS=D1']}")
        logger.info(f"Markers omitted: {log_dict['PS=E']}")
        logger.info("Finished statistically phasing input VCF")
        logger.info(LINEBREAK)

    log_dict = {'status': 'Skipped because input VCF is empty', 'attempted': 0, 'PS=D1': 0, 'PS=E': 0}
    #write_log()

    if vcf_empty:
        write_log()
        return({'combined_vcf':conformed_vcf,'impute_dict':[]})
        #return conformed_vcf

    if vcf_sep == '|':
        log_dict['status'] = 'Skipped because input VCF is already fully phased'
        write_log()
        return({'combined_vcf':conformed_vcf,'impute_dict':[]})
        #return conformed_vcf

    phaseme_vcf = copy_vcf(conformed_vcf, ['meta', 'header', 'data'])
    #print("Phaseme:",phaseme_vcf.data)        

    for i in reversed(range(len(phaseme_vcf.data))):
        fields = phaseme_vcf.data[i]
        inf = fields[7].split(';')
        if any(['PS=C' in x for x in inf]):
            del phaseme_vcf.data[i]

    log_dict['attempted'] = len(phaseme_vcf.data)

    if log_dict['attempted'] == 0:
        log_dict['status'] = 'Skipped because none of the markers are eligible for phasing'
        write_log()
        return({'combined_vcf':conformed_vcf,'impute_dict':[]})
        #return conformed_vcf

    if log_dict['attempted'] == 1:
        log_dict['status'] = 'Skipped because there is only single marker eligible for phasing (manually phased)'
        combined_vcf = copy_vcf(conformed_vcf, ['meta', 'header', 'data'])

        for i in range(len(combined_vcf.data)):
            fields = combined_vcf.data[i]
            inf = fields[7].split(';')
            if not any(['PS=C' in x for x in inf]):
                fields[7] = 'PS=D2'
                combined_vcf.data[i] = fields[:9] + [x.replace('/', '|') for x in fields[9:]]
                impute_dict=[]
                break
        write_log()
        return({'combined_vcf':combined_vcf,'impute_dict':impute_dict})
    

    log_dict['status'] = 'Completed'
    write_vcf(phaseme_vcf, f'{output_dir}/phaseme.vcf')
    subprocess.call(['java', '-Xmx2g', '-jar', program_dir + '/beagle.05May22.33a.jar', f'gt={output_dir}/phaseme.vcf', 'chrom=' + target_region, 'ref=' + hap_panel, f'out={output_dir}/phased', 'impute=' + str(impute_ungenotyped_markers).lower()], stdout = open(f'{output_dir}/beagle.log', 'a'), stderr = open(f'{output_dir}/beagle.log', 'a'))
    subprocess.call(['gunzip', f'{output_dir}/phased.vcf.gz'])
    phased_vcf = read_vcf_simple(f'{output_dir}/phased.vcf')
    #print("Phasing done",time.process_time() - start)
    impute_dict=check_impute_calls(phaseme_vcf,phased_vcf,genome_build,program_dir)
    #print("Impute dict:",impute_dict)
    phasing_time=time.process_time() - start
    start=time.process_time()
    combined_vcf = copy_vcf(conformed_vcf, ['meta', 'header'])
    
    phased_pos=[int(x[1]) for x in phased_vcf.data]
    for fields1 in conformed_vcf.data:
        pos1, ref1, alt1, inf1 = int(fields1[1]), fields1[3], fields1[4], fields1[7].split(';')

        if any(['PS=C' in x for x in inf1]):
            fields1[9:] = [x.replace('./.', '0/0') for x in fields1[9:]]

            if (pos1 in [x.pos for x in snp_list]) and all([x.split(':')[0].split('/')[0] == x.split(':')[0].split('/')[1] for x in fields1[9:]]):
                fields1[9:] = [x.replace('/', '|') for x in fields1[9:]]
                fields1[7] = 'PS=D2'

            combined_vcf.data.append(fields1)
            continue

        is_found = False

        if pos1 not in phased_pos:
            continue

        fields2=phased_vcf.data[phased_pos.index(pos1)]
        pos2, ref2, alt2 = int(fields2[1]), fields2[3], fields2[4]
        gt=lambda x: fields2[9:][x].split(':')[0]
        if "AD" in fields[8].split(":"):
            ad_ind=fields[8].split(":").index("AD")
            ad = lambda x: fields1[9:][x].split(':')[ad_ind]
        else:
            for i in list(range(len(fields2[9:]))):
                if gt(i)=="0|1" or gt(i)=="1|0":
                    ad=lambda x: "50,50"
                elif gt(i)=="1|1":
                    ad=lambda x: "0,50"
                else:
                    ad=lambda x: "0,0"
        idx=list(range(len(fields1[9:])))
        def redo_gt(gt,ad,vcf_ad):
            if vcf_ad!=True:
                s=gt+":"+ad
                return(s)
            if (gt == '0|1' or gt == '1|1' or gt == '1|0') and (ad =='0,0' or ad =='0,0,0' or ad =='0,0,0,0'):
                gt="0|0"
                s=gt+":"+ad
                return(s)
            else:
                s=gt+":"+ad
                return(s)
        
        if pos1 == pos2 and ref1 == ref2 and alt1 == alt2:
            is_found = True
            #fields1[9:] = [gt(i) + ':' + ad(i) for i in idx] if vcf_ad else [gt(i) for i in idx]
            fields1[9:] = [gt(i) + ':' + ad(i) for i in idx]
            fields1[9:] = [redo_gt(gt(i),ad(i),vcf_ad) for i in idx]
            log_dict['PS=D1'] += 1
            inf1 = ['PS=D1' if 'PS=B' in x else x for x in inf1]
            #break
            #continue
        if not is_found:
            log_dict['PS=E'] += 1
            inf1 = ['PS=E' if 'PS=B' in x else x for x in inf1]

        fields1[7] = ';'.join(inf1)
        combined_vcf.data.append(fields1)


    logger.info(f"Markers attempted: {log_dict['attempted']}")
    logger.info(f"Markers phased: {log_dict['PS=D1']}")
    logger.info(f"Markers omitted: {log_dict['PS=E']}")
    logger.info("Finished statistically phasing input VCF")
    #logger.info(f"Phased VCF: {phasing_time}")
    #logger.info(f"Combined VCF: {time.process_time() - start}")
    logger.info(LINEBREAK)
   
    return({'combined_vcf':combined_vcf,'impute_dict':impute_dict})

def annotate_vcf(combined_vcf, target_locus, snp_list, data_type, vcf_ad, logger):
    # Start logging.
    logger.info(LINEBREAK)
    logger.info("Step 6/10: Annotating input VCF...")

    log_dict = {'stargazer_membership': 0, 'low_impact': 0, 'moderate_impact': 0, 'high_impact': 0, 'reverting_variation': 0}
    annotated_vcf = copy_vcf(combined_vcf, ['meta', 'header', 'data'])
    undetected_revertants = [x for x in snp_list if x.rev]

    for fields in annotated_vcf.data:
        pos, ref, alt, inf = int(fields[1]), fields[3], fields[4].split(','), fields[7].split(';')
        vi_list = []; rv_list = []; so_list = []; fe_list = []

        for var in alt:
            filtered = [x for x in snp_list if pos == x.pos and ref == x.hg and (var == x.var or var == x.wt)]
            #print(filtered,pos,var)
            if filtered:
                #print(filtered[0].impact,filtered[0].n,filtered[0].effect,filtered[0].pos,filtered[0].id,filtered[0].hg,filtered[0].var,filtered[0].wt,filtered[0].so,filtered[0].impact,filtered[0].rev)
                vi, rv, so, fe = filtered[0].impact, 'revertant_true' if filtered[0].rev else 'revertant_false', filtered[0].so, filtered[0].effect
                #print (vi,rv,so,fe)
                log_dict[vi] += 1
                undetected_revertants = [x for x in undetected_revertants if x != filtered[0]]
                if filtered[0].rev:
                    log_dict['reverting_variation'] += 1
            else:
                vi, rv, so, fe = 'unknown_impact', 'revertant_unknown', 'unknown_variant', 'unknown_effect'
            vi_list.append(vi); rv_list.append(rv); so_list.append(so); fe_list.append(fe)

        if any([x != 'unknown_impact' for x in vi_list]):
            log_dict['stargazer_membership'] += 1
       
        inf.append('VI=' + ','.join(vi_list)); inf.append('RV=' + ','.join(rv_list)); inf.append('SO=' + ','.join(so_list)); inf.append('FE=' + ','.join(fe_list))
        fields[7] = ';'.join(inf)

    # Manually add variants that are part of the genome assembly
    if data_type != 'chip':
        chr = target_locus.chr
        dat, fmt = (['0|0:0,0' for x in annotated_vcf.header[9:]], 'GT:AD') if vcf_ad else (['0|0' for x in annotated_vcf.header[9:]], 'GT')
        for snp in undetected_revertants:
            inf = ';'.join(['PS=D2', 'VI=' + snp.impact, 'RV=revertant_true', 'SO=' + snp.so, 'FE=' + snp.effect])
            fields = [chr, snp.pos, snp.id, snp.hg, snp.wt, '.', 'PASS', inf, fmt] + dat
            annotated_vcf.data.append(fields)
        annotated_vcf.data.sort(key = lambda x: int(x[1]))

    # Manually phase, without certainty, any unphased revertants
    for fields in annotated_vcf.data:
        inf = fields[7].split(';')
        if any(['PS=C' in x for x in inf]) and 'RV=revertant_true' in inf:
            fields[7] = ';'.join(['PS=D3' if 'PS=C' in x else x for x in inf])
            fields[9:] = [x.replace('/', '|') for x in fields[9:]]


    logger.info(f"Markers with Stargazer membership: {log_dict['stargazer_membership']}")
    logger.info(f"Variants with low impact: {log_dict['low_impact']}")
    logger.info(f"Variants with moderate impact: {log_dict['moderate_impact']}")
    logger.info(f"Variants with high impact: {log_dict['high_impact']}")
    logger.info(f"Variants reverted to wild type: {log_dict['reverting_variation']}/{len([x for x in snp_list if x.rev])}")
    logger.info("Finished annotating input VCF")
    logger.info(LINEBREAK)

    return annotated_vcf

def account_vcf(annotated_vcf, target_gene, vcf_ad, gb):
    accounted_vcf = copy_vcf(annotated_vcf, ['meta', 'header', 'data'])
    accounted_vcf.meta = ['##reference={}-{}*1\n'.format(gb,target_gene.upper()) if '##reference=' in x else x for x in accounted_vcf.meta]

    for fields in accounted_vcf.data:
        ref, alt, inf = fields[3], fields[4].split(','), fields[7].split(';')
        rv_list = [x for x in inf if 'RV=' in x][0].replace('RV=', '').split(',')
        if 'revertant_true' not in rv_list:
            continue
        i = rv_list.index('revertant_true')
        fields[3] = alt[i]
        fields[4] = ','.join([ref if x == alt[i] else x for x in alt])

        def f(x):
            gt = x.split(':')[0].split('|')
            field = '|'.join([str(i + 1) if y == '0' else '0' if y == str(i + 1) else y for y in gt])
            if vcf_ad:
                ad = x.split(':')[1].split(',')
                field += ':' + ','.join([ad[i + 1] if y == 0 else ad[0] if y == i + 1 else ad[y] for y in range(len(ad))])
            return field

        fields[9:] = [f(x) for x in fields[9:]]

    return accounted_vcf

def extend_vcf(accounted_vcf, star_dict, logger):
    # Start logging.
    logger.info(LINEBREAK)
    logger.info("Step 7/10: Phasing input VCF by haplotype extension...")

    log_dict = {'attempted': 0, 'phased': 0, 'omitted': 0}
    finalized_vcf = copy_vcf(accounted_vcf, ['meta', 'header', 'data'])
    pseudo_samples = vcf2samples(finalized_vcf)

    for fields in finalized_vcf.data:

        pos, ref, alt, inf, fmt = int(fields[1]), fields[3], fields[4].split(','), fields[7].split(';'), fields[8]
        ps, vi = inf[0], inf[1]
        if 'PS=D' in ps or ('high_impact' not in vi and 'low_impact' not in vi and 'moderate_impact' not in vi and 'stargazer_membership' not in vi):
            continue
        log_dict['attempted'] += 1

        def pbe(i):
            x = fields[i]
            gt = x.split(':')[0].split('/')
            if gt[0] == gt[1]:
                return x.replace('/', '|')
            scores = [[0, 0], [0, 0]]
            for j in [0, 1]:
                if gt[j] == '0':
                    continue
                idx = int(gt[j])
                target_snp = SNPAllele()
                target_snp.pos = pos
                target_snp.wt = ref
                target_snp.var = alt[idx - 1]
                relevant_stars = [v for k, v in star_dict.items() if target_snp in v.core]
                name = finalized_vcf.header[i]
                for k in [0, 1]:
                    for star in relevant_stars:
                        score = 0
                        for snp in pseudo_samples[i-9].hap[k].obs:
                            if snp in star.core + star.tag:
                                score += 1
                        if score > scores[j][k]:
                            scores[j][k] = score
            a = scores[0]
            b = scores[1]
            flip = False
            if max(a) == max(b):
                if a[0] < a[1] and b[0] > b[1]:
                    flip = True
                elif a[0] == a[1] and b[0] > b[1]:
                    flip = True
                elif a[0] < a[1] and b[0] == b[1]:
                    flip = True
                else:
                    pass
            else:
                if max(a) > max(b):
                    if a[0] > a[1]:
                        pass
                    else:
                        flip = True
                else:
                    if b[0] > b[1]:
                        flip = True
                    else:
                        pass
            if flip:
                result = f'{gt[1]}|{gt[0]}'
            else:
                result = f'{gt[0]}|{gt[1]}'
            if 'AD' in fmt:
                ad_ind=fmt.split(":").index("AD")
                result = result + ':' + x.split(':')[ad_ind]
            result = result + f':{",".join([str(x) for x in a + b])}'
            return result

        new_fields = [pbe(i) for i in range(9, len(fields))]

        if not all(new_fields):
            log_dict['omitted'] += 1
            continue

        log_dict['phased'] += 1

        inf[0] = 'PS=D5'
        fields[7] = ';'.join(inf)
        fields[8] += ':HE'
        fields[9:] = new_fields

    logger.info(f"Markers attempted: {log_dict['attempted']}")
    logger.info(f"Markers phased: {log_dict['phased']}")
    logger.info(f"Markers omitted: {log_dict['omitted']}")
    logger.info("Finished phasing input VCF by haplotype extension")
    logger.info(LINEBREAK)

    return finalized_vcf


def pretty_vcf(processed_vcf, finalized_vcf, logger):
    # Start logging.
    logger.info(LINEBREAK)
    logger.info("Step 8/10: Creating a pretty VCF as a final product...")

    log_dict = {'saved': 0, 'filled': 0}
    pretty_vcf = copy_vcf(finalized_vcf, ['meta', 'header', 'data'])
    
    pos_all={}
    for fields in processed_vcf.data:
        pos,ref,alt=int(fields[1]), fields[3], fields[4]
        pos_all[pos]=[ref,alt]
    
    for fields in pretty_vcf.data:

        pos= int(fields[1])
        log_dict['saved'] += 1
        if pos in pos_all.keys():
            fields[3]=pos_all[pos][0]
            fields[4]=pos_all[pos][1]
            
    return pretty_vcf
            

    


def main():
    # Start the timer. Used later to compute elapsed run time.
    start_time = timeit.default_timer()

    # Get absolute path to the program directory.
    program_dir = os.path.dirname(os.path.realpath(__file__))
    

    # Import the gene table.
    gene_df = pd.read_table(f"{program_dir}/gene_table.tsv")
    target_genes = gene_df[gene_df["type"] == "target"]["name"].to_list()
    control_genes = gene_df[gene_df["control"] == "yes"]["name"].to_list()

    # Define input data types.
    data_types = {'wgs': 'Whole genome sequencing',
                  'ts': 'Targeted sequencing',
                  'chip': 'Single nucleotide polymorphism array'}

    # Parse user-provided arguments.
    parser = get_parser(target_genes, control_genes, list(data_types), program_dir)
    args = parser.parse_args()
    print(args)

    # Create the output directory. Will overwrite if already existing.
    try:
        shutil.rmtree(args.output_dir)
    except OSError:
        pass
    os.mkdir(args.output_dir)

    # Start logging.
    logging.basicConfig(
        level=logging.INFO,
        format="[%(levelname)s] %(message)s",
        handlers=[
            logging.FileHandler(f"{args.output_dir}/stargazer.log"),
            logging.StreamHandler()
        ]
    )
    logger = logging.getLogger()

    # Give the intro.
    logger.info(LINEBREAK)
    logger.info(f"Stargazer {__version__}")
    logger.info('Author: Seung-been "Steven" Lee')
    logger.info("Enter 'python3 Stargazer -h' to view command line arguments")
    logger.info("For more details, please visit "
                "https://stargazer.gs.washington.edu/stargazerweb")
    logger.info(f"Date: {datetime.datetime.now().strftime('%Y-%m-%d')}")
    logger.info(f"Start time: {datetime.datetime.now().strftime('%H:%M:%S')}")
    logger.info(f"Stargazer path: {sys.argv[0]}")
    logger.info(f"Command line: stargazer {' '.join(sys.argv[1:])}")
    logger.info(LINEBREAK)

    # Get the Locus instances.
    loci = get_loci(gene_df, args.target_gene,
        args.control_gene, control_genes, args.genome_build)
    target_locus, control_locus, paralog_locus = loci

    # Log the genotype mode.
    logger.info(LINEBREAK)
    if args.create_gdf_file is None:
        logger.info("Step 1/9: Determining genotype mode...")
    else:
        logger.info('Entering GDF file generation mode.')
        
    logger.info(f"Reference genome assembly: {args.genome_build}")
    logger.info(f"Target gene: {target_locus.name.upper()}")
    if paralog_locus is None:
        logger.info(f"Target paralog: None")
    else:
        logger.info(f"Target paralog: {paralog_locus.name.upper()}")
    logger.info(f"Target region: chr{target_locus.region}")
    if control_locus is None:
        logger.info(f"Control gene: None")
    else:
        logger.info(f"Control gene: {control_locus.name.upper()}")
        logger.info(f"Control region: chr{control_locus.region}")
    logger.info(f"Input data source: {data_types[args.data_type]}")

    # run as a GDF creator (then leave)
    if args.create_gdf_file is not None:
                       
        bam2gdf(target_locus.region,control_locus.region,**vars(args))
        
        logger.info(f"Finished creating GDF file: {args.output_dir}/{args.create_gdf_file}")
        logger.info(LINEBREAK)
        exit() 

    # Import the data files.
    snp_list = get_snp_list(args.target_gene, program_dir, args.genome_build)
    star_dict = get_star_dict(args.target_gene, snp_list,
        program_dir, args.genome_build, args.star_table)
    sv_table = read_sv_table(f"{program_dir}/sv_table.tsv")
    sv_dict = sv_table[args.target_gene]


    # Determine whether VCF-only mode should be turned on.
    vcf_only = not all([args.control_gene, args.gdf_file])


    del_allele = [v.name for k, v in star_dict.items() if v.sv == 'cnv0'][0]
    if args.hap_panel is None:
        panel_file = "{}/1kgp_vcf/{}/{}.vcf.gz".format(program_dir, args.genome_build, args.target_gene)
    else:
        panel_file = args.hap_panel


    logger.info(f"VCF-only mode is turned on: {vcf_only}")
    logger.info(f"Impute ungenotyped markers: {args.impute_ungenotyped_markers}")
    logger.info("Finished determining genotype mode")
    logger.info(LINEBREAK)


    ref_vcf = read_vcf_simple(panel_file)
    input_vcf = read_vcf_region(args.vcf_file, target_locus.region)
    vcf_ad, vcf_sep, vcf_empty = assess_vcf(input_vcf, args.gdf_file, args.data_type, logger) # step 2

    processed_vcf = process_vcf(input_vcf, vcf_ad, vcf_sep, logger,args.genome_build) # step 3
    adjusted_vcf = adjust_vcf(processed_vcf,star_dict)
    #write_vcf(adjusted_vcf, f'{args.output_dir}/adjusted.vcf')
    conformed_vcf = conform_vcf(adjusted_vcf, ref_vcf, snp_list, vcf_ad, vcf_empty, vcf_sep, vcf_only, logger) # step 4
    #write_vcf(conformed_vcf, f'{args.output_dir}/conformed.vcf')
    result_dict = phase_vcf(conformed_vcf, program_dir, target_locus.region, snp_list, args.impute_ungenotyped_markers, args.output_dir, panel_file, vcf_ad, vcf_empty, vcf_sep, logger, args.genome_build) # step 5
    #print(result_dict)
    combined_vcf = result_dict['combined_vcf']
    impute_dict=result_dict['impute_dict']
    #write_vcf(combined_vcf, f'{args.output_dir}/combined.vcf')
    annotated_vcf = annotate_vcf(combined_vcf, target_locus, snp_list, args.data_type, vcf_ad, logger) # step 6
    accounted_vcf = account_vcf(annotated_vcf, args.target_gene, vcf_ad, args.genome_build)
    finalized_vcf = extend_vcf(accounted_vcf, star_dict, logger) # step 7
    #write_vcf(finalized_vcf, f'{args.output_dir}/pre_finalized.vcf')
    #
    # step 8: create a VCF for the public that has original REF and ALT alleles
    prettied_vcf = pretty_vcf(processed_vcf,finalized_vcf,logger)
    write_vcf(prettied_vcf, f'{args.output_dir}/finalized.vcf')

    persons = vcf2samples(finalized_vcf)


    for sample in persons:
        if vcf_only:
            sample.sv = ['no_sv', 'no_sv']
            sample.ssr = '.'

    # Call structural variation.
    logger.info(LINEBREAK)
    logger.info("Step 9/10: Detecting structural variants...")
    if vcf_only:
        logger.info("Skipped because VCF-only mode is turned on")
    else:
        dataframes = call_sv(persons,
                             args.gdf_file,
                             target_locus,
                             control_locus,
                             args.data_type,
                             sv_dict,
                             logger,
                             args.genome_build,
                             args.sample_list)
        target_df, control_df, profile_df, unmasked_df = dataframes
        logger.info("Finished detecting structural variants")
    logger.info(LINEBREAK)

    for sample in persons:
        f = lambda x: sorted([v for k, v in star_dict.items() if set(v.core).issubset(x) and not (v.sv and v.sv not in sample.sv)], key = lambda x: x.rank)
        hap1_snp = [x for x in sample.hap[0].obs if x.wt != x.var]
        hap2_snp = [x for x in sample.hap[1].obs if x.wt != x.var]
        sample.hap[0].cand = f(hap1_snp)
        sample.hap[1].cand = f(hap2_snp)
        sample.dip_cand = f(list(set(hap1_snp + hap2_snp)))

    genotype_callers = {'cyp2a6': cyp2a6,
                        'cyp2b6': cyp2b6,
                        'cyp2d6': cyp2d6,
                        'cyp2e1': cyp2e1,
                        'gstm1': gstm1,
                        'gstt1': gstt1,
                        'slc22a2': slc22a2,
                        'slco1b1': slco1b1,
                        'ugt1a4': ugt1a4,
                        'ugt2b15': ugt2b15,
                        'ugt2b17': ugt2b17}

    for sample in persons:
        
        if sample.sv == ['no_sv', 'no_sv']: sample.gt = True
        call_sv1(sample, 'cnv0', del_allele, star_dict)
        call_cnv3(sample, target_locus)
        call_cnv_plus(sample, target_locus)

        if target_locus.name in genotype_callers:
            genotype_callers[args.target_gene](sample, star_dict)

    # Remove extra *1 alleles.
    def f(l,target_locus):
        if target_locus.name=="cyp2c19":
            wt_star="*38"
        else:
            wt_star="*1"
        if len(l) == 1:
            return
        for i in reversed(range(len(l))):
            if l[i].name == wt_star:
                del l[i]

    for sample in persons:
        f(sample.hap[0].cand,target_locus)
        f(sample.hap[1].cand,target_locus)
        f(sample.dip_cand,target_locus)

    # Order the haplotypes.
    for sample in persons:
        if not sample.gt:
            continue
        if sort_star_names([sample.hap[0].cand[0].name, sample.hap[1].cand[0].name])[0] == sample.hap[1].cand[0].name:
            sample.hap[0], sample.hap[1] = sample.hap[1], sample.hap[0]

    # Predict the phenotype.
    for person in persons:
        #print(target_locus.name,person.hap[0].cand[0].name,person.hap[1].cand[0].name)
        person.pt = phenotyper(target_locus.name,
                               person.hap[0].cand[0].name,
                               person.hap[1].cand[0].name)

    # Write the result file.
    result_df = pd.DataFrame.from_records(
        [x.summarize(target_locus) for x in persons],
        columns=['name', 'status', 'hap1_main', 'hap2_main',
                 'hap1_cand', 'hap2_cand', 'hap1_score', 'hap2_score',
                 'dip_score', 'phenotype', 'dip_sv', 'hap1_sv',
                 'hap2_sv','ssr', 'dip_cand', 'hap1_main_core',
                 'hap2_main_core', 'hap1_main_tag', 'hap2_main_tag',
                 'hap1_af_mean_gene', 'hap2_af_mean_gene',
                 'hap1_af_mean_main', 'hap2_af_mean_main']
    )
    result_df.to_csv(f"{args.output_dir}/genotype-calls.tsv", sep='\t', index=False)

    # Create and write report file.
    report_df=pd.DataFrame(create_report(result_df,f"{program_dir}/star_table.tsv",args.target_gene,impute_dict))
    report_df.to_csv(f"{args.output_dir}/report.tsv", sep='\t', index=False)


    # Plot read detph, copy number, and allele fraction.
    logger.info(LINEBREAK)
    logger.info("Step 10/10: Plotting various profiles...")
    if vcf_only:
        logger.info("Skipped because VCF-only mode is turned on")
    elif args.include_profiles:
        plot(target_df,
             control_df,
             profile_df,
             unmasked_df,
             target_locus,
             control_locus,
             args.output_dir,
             persons,
             logger,
             paralog_locus=paralog_locus)
        logger.info("Finished plotting various profiles")
    else:
        logger.info("Skipped due to user request")
    logger.info(LINEBREAK)

    # Stop the timer and calculate elasped run time.
    stop_time=timeit.default_timer()
    elapsed_time=str(datetime.timedelta(
        seconds=(stop_time - start_time))).split('.')[0]

    # Give the outro.
    logger.info(LINEBREAK)
    logger.info(f"Elapsed time: {elapsed_time}")
    logger.info("Stargazer finished")
    logger.info(LINEBREAK)

    import sysconfig
    print(sysconfig.get_paths())
if __name__ == '__main__':
    main()
