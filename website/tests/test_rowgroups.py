"""RowGroups (row-group lookup in a key-sorted parquet, #146 speed-up) returns exactly the rows plain SQL returns, including keys that
span row-group boundaries, keys absent from the file, and refuses a file that is not sorted by the key."""
import os, sys, tempfile
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
os.environ.setdefault('SCBRAIN_DATA', 'D:/超算/scBRAIN/db')
import duckdb, pytest
from app.main import RowGroups

SRC = 'D:/超算/scBRAIN/db/tables/esnp_pairs.parquet'
TMP = tempfile.mkdtemp(prefix='rg_')
SORTED = TMP + '/esnp_sorted.parquet'
duckdb.execute(f"COPY (SELECT * FROM read_parquet('{SRC}') ORDER BY gene, ctx, pval) TO '{SORTED}' (FORMAT parquet, ROW_GROUP_SIZE 20000)")
IX = RowGroups(SORTED, 'gene')


def sql_rows(g):
    return duckdb.execute(f"SELECT * FROM read_parquet('{SORTED}') WHERE gene = '{g}' ORDER BY ctx, pval, rsid").df()


def canon(d):
    return d.sort_values(['ctx', 'pval', 'rsid']).reset_index(drop=True)


def test_boundary_and_random_keys_match_sql():
    genes = set(IX.lo[:40]) | set(IX.hi[:40]) | {'ENSG00000136235'}
    for g in genes:
        a, b = canon(IX.rows(g)), sql_rows(g)
        assert len(a) == len(b) > 0, g
        assert (a[['ctx', 'rsid']].values == b[['ctx', 'rsid']].values).all()
        assert (a.pval.values == b.pval.values).all()


def test_absent_key_is_empty():
    assert len(IX.rows('ENSG99999999999')) == 0
    assert len(IX.rows('ENSG00000000000')) == 0


def test_unsorted_file_refused():
    with pytest.raises(ValueError):
        RowGroups(SRC, 'rsid')
