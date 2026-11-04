"""Training-only operating-point selection with explicit unsupported status."""
import numpy as np


def select_operating_point(rows,n_gene_words,n_blank_words,blank_ceiling=.1,min_gene_calls=50):
    """Maximize training gene yield under a declared blank proxy ceiling.

    This does not match true molecular errors or supply a confidence bound.
    A zero-blank, low-count sample cannot qualify by itself. The chosen state
    must remain fixed across geometric arms and independent evaluation cells.
    """
    if not isinstance(n_gene_words,int) or not isinstance(n_blank_words,int) or min(n_gene_words,n_blank_words)<1:
        raise ValueError('positive integer codebook class sizes required')
    if not 0<blank_ceiling<1 or not isinstance(min_gene_calls,int) or min_gene_calls<1:
        raise ValueError('invalid operating point qualification')
    audited=[]
    for index,row in enumerate(rows):
        g=row['gene_calls'];b=row['blank_calls']
        if not isinstance(g,int) or not isinstance(b,int) or min(g,b)<0:raise ValueError('counts must be nonnegative integers')
        proxy=b/g*n_gene_words/n_blank_words if g else None
        eligible=g>=min_gene_calls and proxy is not None and proxy<=blank_ceiling
        audited.append(dict(index=index,**row,computed_blank_proxy=proxy,eligible=eligible))
    eligible=[r for r in audited if r['eligible']]
    if not eligible:return dict(qualified=False,abstained=True,reason='no_training_setting_meets_count_and_proxy_criteria',selected=None,rows=audited)
    chosen=min(eligible,key=lambda r:(-r['gene_calls'],r['computed_blank_proxy'],r['index']))
    return dict(qualified=True,abstained=False,reason='training_only_count_and_proxy_selection',selected=chosen,rows=audited)
