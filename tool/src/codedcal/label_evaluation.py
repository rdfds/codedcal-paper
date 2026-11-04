"""Label-consistent set recovery, including co-located different words."""
import numpy as np
from .evaluation import validate_calls
from scipy.spatial import cKDTree
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import maximum_bipartite_matching

def label_consistent_score(reference_xy,reference_words,call_xy,call_words,is_blank,radius=1.5,max_edges=2000000):
    """Maximum same-label spatial matches; all calls remain in the denominator.

    Unlike a position-first assignment this cannot swap identities arbitrarily
    at a demixed, co-located pair. Independent truth is still required for an
    accuracy interpretation; a decoded reference gives consistency only.
    Unsolved matching components yield unavailable F1, never apparent zeros.
    """
    blank=np.asarray(is_blank)
    if blank.ndim!=1 or blank.dtype!=bool or not blank.any() or blank.all():
        raise ValueError('boolean gene/blank book mask required')
    a,aw=validate_calls(reference_xy,reference_words,len(blank))
    b,bw=validate_calls(call_xy,call_words,len(blank))
    if blank[aw].any():raise ValueError('reference contains blanks')
    if not np.isfinite(radius) or radius<=0 or not isinstance(max_edges,int) or max_edges<1:
        raise ValueError('positive radius and integer edge budget required')
    correct=0;skipped=0;ambiguous=0;edge_count=0
    for word in np.unique(aw):
        aa=a[aw==word];bb=b[bw==word]
        if not len(bb):continue
        tree=cKDTree(bb)
        counts=tree.query_ball_point(aa,radius,return_length=True)
        edges=int(counts.sum());edge_count+=edges;ambiguous+=int((counts>1).sum())
        if edges>max_edges:
            skipped+=1;continue
        neighbors=tree.query_ball_point(aa,radius)
        indices=np.concatenate(neighbors).astype(int) if edges else np.array([],int)
        graph=csr_matrix((np.ones(edges,np.int8),indices,np.r_[0,np.cumsum(counts)]),shape=(len(aa),len(bb)))
        correct+=int((maximum_bipartite_matching(graph,perm_type='column')>=0).sum())
    supported=skipped==0;r=len(a);c=len(b)
    return dict(n_reference=r,calls=c,correct_label_matches=correct,
        supported=supported,skipped_label_graphs=skipped,candidate_edges=edge_count,
        ambiguous_same_label_references=ambiguous,
        identity_f1=2*correct/(r+c) if supported and r+c else None,
        identity_precision=correct/c if supported and c else None,
        identity_recall=correct/r if supported and r else None,
        gene_calls=int((~blank[bw]).sum()),blank_calls=int(blank[bw].sum()))

def training_support(words,local_traces):
    """Diagnose the active-code incidence graph for relative gain fitting.

    Inactive bits and disconnected groups have no shared amplitude scale from
    these selected words. This is structural support, not statistical precision
    or independent validation of inferred training labels.
    """
    from scipy.sparse.csgraph import connected_components
    from scipy.sparse import csr_matrix
    words=np.asarray(words,float)
    if len(local_traces)!=len(words):raise ValueError('training word mismatch')
    active=np.array([len(x)>0 for x in local_traces])
    on=words[active]>0;covered=on.any(0) if len(on) else np.zeros(words.shape[1],bool)
    adjacency=on.astype(int).T@on.astype(int)
    n,labels=connected_components(csr_matrix(adjacency),directed=False)
    return dict(active_words=int(active.sum()),covered_readouts=covered.tolist(),
        uncovered_readouts=np.flatnonzero(~covered).tolist(),incidence_components=int(n),
        connected_complete_support=bool(covered.all() and n==1))
