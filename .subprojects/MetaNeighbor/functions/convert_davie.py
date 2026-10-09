import gzip
import numpy as np
import pandas as pd
import h5py
import anndata as ad
import scipy.sparse as sp

def load_annotation(path):
    sym2fbgn, annid2fbgn = {}, {}
    with gzip.open(path, "rt") as f:
        for line in f:
            if line.startswith("#"):
                continue
            p = line.rstrip("\n").split("\t")
            if len(p) < 5:
                continue
            sym, fbgn, annid = p[0], p[2], p[4]
            if sym:
                sym2fbgn.setdefault(sym, set()).add(fbgn)
            if annid:
                annid2fbgn.setdefault(annid, set()).add(fbgn)
    return sym2fbgn, annid2fbgn

def load_synonyms(path):
    syn = {}
    with gzip.open(path, "rt") as f:
        for line in f:
            if line.startswith("#"):
                continue
            p = line.rstrip("\n").split("\t")
            if len(p) < 6:
                continue
            fbgn = p[0]
            for s in p[5].split("|"):
                s = s.strip()
                if s:
                    syn.setdefault(s, set()).add(fbgn)
    return syn