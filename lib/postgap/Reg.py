#! /usr/bin/env python

import Globals
from postgap.DataModel import *
import postgap.BedTools
from postgap.Utils import *
import logging
import requests
import subprocess
import tempfile


class Reg_source(object):
    def run(self, ld_snps):
        assert False, "This stub should be defined"


class S_DHS(Reg_source):
    display_name = "S_DHS"

    def run(self, ld_snps):
        snp_hash = dict((snp.rsID, snp) for snp in ld_snps)
        intersection = postgap.BedTools.overlap_snps_to_bed(
            ld_snps, postgap.Globals.DATABASES_DIR + "/S_DHS.bed")
        res = filter(lambda X: X.score, (self.get_regulome_evidence(
            feature, snp_hash) for feature in intersection))
        return res

    def get_regulome_evidence(self, feature, snp_hash):
        return Regulatory_Evidence(
            snp=snp_hash[feature[6]],
            source=self.display_name,
            score=1,
            study=None,
            tissue=None,
            info=None
        )


class S_TSS(Reg_source):
    display_name = "S_TSS"

    def run(self, ld_snps):
        snp_hash = dict((snp.rsID, snp) for snp in ld_snps)
        intersection = postgap.BedTools.overlap_snps_to_bed(
            ld_snps, postgap.Globals.DATABASES_DIR + "/S_TSS.bed")
        res = filter(lambda X: X.score, (self.get_regulome_evidence(
            feature, snp_hash) for feature in intersection))
        return res

    def get_regulome_evidence(self, feature, snp_hash):
        return Regulatory_Evidence(
            snp=snp_hash[feature[6]],
            source=self.display_name,
            score=1,
            study=None,
            tissue=None,
            info=None
        )


class S_coding(Reg_source):
    display_name = "S_coding"

    def run(self, ld_snps):
        snp_hash = dict((snp.rsID, snp) for snp in ld_snps)
        intersection = postgap.BedTools.overlap_snps_to_bed(
            ld_snps, postgap.Globals.DATABASES_DIR + "/S_coding.bed")
        res = filter(lambda X: X.score, (self.get_regulome_evidence(
            feature, snp_hash) for feature in intersection))
        return res

    def get_regulome_evidence(self, feature, snp_hash):
        return Regulatory_Evidence(
            snp=snp_hash[feature[6]],
            source=self.display_name,
            score=1,
            study=None,
            tissue=None,
            info=None
        )


def get_filtered_subclasses(subclasses_filter):
    return [subclass for subclass in sources if subclass.display_name in subclasses_filter]


sources = Reg_source.__subclasses__()
