"""/api/celltype/{cell}: one cell type's estimated fraction in the 10 regions (PI 10-20「细胞比例就放首页，也像这个图可以点进去」).
python -m pytest db_site/tests -q   (SCBRAIN_DATA must point at the db tables)"""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
os.environ.setdefault('SCBRAIN_DATA', 'D:/超算/scBRAIN/db')
from fastapi.testclient import TestClient
from app.main import app, REG

c = TestClient(app)


def test_celltype_regions_in_order_with_stats():
    r = c.get('/api/celltype/Purk')
    assert r.status_code == 200
    d = r.json()
    assert d['cell'] == 'Purk' and d['cell_name'] == 'Purkinje cells'
    assert [x['region'] for x in d['regions']] == REG
    for x in d['regions']:
        assert x['lo'] <= x['q1'] <= x['median'] <= x['q3'] <= x['hi']
        assert x['n'] > 0 and x['ctx'] == x['region'] + '.Purk'


def test_purkinje_highest_in_cerebellum():
    d = c.get('/api/celltype/Purk').json()
    top = max(d['regions'], key=lambda x: x['median'])
    assert top['region'] == 'CRBL'


def test_every_cell_is_strict_json():
    import json
    for cell in ['ExN', 'OPC', 'Gran', 'Purk']:
        r = c.get('/api/celltype/' + cell)
        json.loads(r.text, parse_constant=lambda s: (_ for _ in ()).throw(ValueError(s)))   # NaN / Infinity break the browser


def test_unknown_cell_is_404():
    assert c.get('/api/celltype/Foo').status_code == 404
