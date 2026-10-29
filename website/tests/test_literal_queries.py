"""Queries rewritten with checked literals (#146 speed-up) still answer, and bad ids are refused instead of reaching the SQL."""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
os.environ.setdefault('SCBRAIN_DATA', 'D:/超算/scBRAIN/db')
from fastapi.testclient import TestClient
from app.main import app, lit
from fastapi import HTTPException
import pytest

c = TestClient(app)


def test_meta_cached_and_same():
    a = c.get('/api/meta').json(); b = c.get('/api/meta').json()
    assert a == b and len(a['contexts']) == 110


def test_gene_page():
    d = c.get('/api/gene/GPNMB').json()
    assert len(d['eqtl']) == 110 and d['n_esnp_pairs'] > 0
    assert [x['ctx'] for x in d['eqtl']] == sorted(x['ctx'] for x in d['eqtl'])


def test_groups_has_both_sexes():
    d = c.get('/api/gene/GPNMB/groups?ctx=FCTX.ExN').json()
    assert {g['group'] for g in d['sex']} == {'Female', 'Male'} and sum(g['n'] for g in d['sex']) > 100


def test_esnps_genotype_variant():
    e = c.get('/api/gene/GPNMB/esnps?ctx=FCTX.ExN').json()['rows']
    assert e and e[0]['pval'] <= e[-1]['pval']
    g = c.get('/api/genotype?gene=GPNMB&rsid=rs156429').json()
    assert g['boxes'] and g['stats']
    v = c.get('/api/variant/rs156429').json()
    assert v['pairs'] and 'GPNMB' in {p['symbol'] for p in v['pairs']}


def test_bad_ids_refused():
    assert c.get('/api/gene/GPNMB/groups?ctx=FCTX.ExN%27%20OR%201=1').status_code == 404
    assert c.get("/api/variant/rs1'--").status_code == 404
    with pytest.raises(HTTPException):
        lit("ENSG1' OR '1'='1", 'gene')
