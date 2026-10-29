"""疾病风险位点富集的 2×2 计数（照贺奕捷 desease.sh 的口径），三套并列：
  orig_asis  原 eSNP 名单 + 原脚本的数法（名单末尾两行空行会让交集多数 1、总数多数 2）——用来对回他的结果
  orig_clean 原 eSNP 名单，去掉空行后按真实交集数
  eg         只留 eGene 的新 eSNP 名单（韩教授 10-09 定），按真实交集数
non-eSNP 名单三套相同。输出长表 fisher_redo/counts.tsv。
"""
import os, sys, time

P = '/public/home/yucong/scbrain_in/paper2'
D = f'{P}/single-cell-eqtl/data'
DIS = '/public/home/yucong/scbrain_in/DISEASE_new'
EG = '/public/home/yucong/scbrain/esnp_eg/esnp'
OUT = os.environ.get('FISHER_OUT', '/public/home/yucong/scbrain/fisher_redo')
os.makedirs(OUT, exist_ok=True)
# 敏感性分析：去掉 6 号染色体 25–35 Mb（hg38，扩展 MHC 区）里的全部变异，只出正式口径 eg
NO_MHC = os.environ.get('EXCLUDE_MHC', '') == '1'
MHC = set()
if NO_MHC:
    with open(f'{D}/genotype/out.bim') as f:
        for line in f:
            a = line.split('\t')
            if a[0] == 'chr6' and 25_000_000 <= int(a[3]) <= 35_000_000:
                MHC.add(a[1])
    print('MHC variants excluded:', len(MHC), flush=True)
REGIONS = ["MEDU", "PUTM", "TCTX", "WHMT", "CRBL", "OCTX", "SNIG", "HIPP", "THAL", "FCTX"]
CELLS = ["Micro", "InN", "ExN", "Gran", "NPC", "Olig", "Endo", "OPC", "Perc", "Astro", "Purk"]
t0 = time.time()


def read_lines(path):
    if not os.path.exists(path):
        return None
    with open(path) as f:
        return [l.rstrip('\n') for l in f]


def load_diseases(kind, namefile):
    names = [l.strip() for l in open(namefile) if l.strip()]
    out = {}
    miss = []
    for n in names:
        p = f'{DIS}/{kind}-disease-snp/{n}2.txt'
        ls = read_lines(p)
        if ls is None:
            miss.append(n); continue
        out[n] = set(x for x in ls if x) - MHC
    print(kind, 'diseases', len(out), 'missing files', len(miss), miss[:5], flush=True)
    return out


dis = {'nervous': load_diseases('nervous', f'{P}/eqtl-gwas/fisher/nervous/nervous-disease-fisher/name.txt'),
       'non-nervous': load_diseases('non-nervous', f'{P}/eqtl-gwas/fisher/non-nervous/non-nervous-disease-fisher/name.txt')}

with open(f'{OUT}/counts.tsv', 'w') as w:
    w.write('variant\tclass\tdisease\tctx\tE.N\tE.NN\tNE.N\tNE.NN\n')
    for r in REGIONS:
        for c in CELLS:
            ctx = f'{r}.{c}'
            e0 = read_lines(f'{D}/esnp/{ctx}-esnp3.txt')
            n0 = read_lines(f'{D}/non-esnp/{ctx}-non-esnp4.txt')
            e1 = read_lines(f'{EG}/{ctx}-esnp3.txt')
            e0 = e0 or []; n0 = n0 or []; e1 = e1 or []
            E0 = set(x for x in e0 if x) - MHC; N0 = set(x for x in n0 if x) - MHC; E1 = set(x for x in e1 if x) - MHC
            # 原脚本的数法：cat 两个文件 | sort | uniq -c | 取计数为 2 的行；名单里恰有 2 行空行时空行也被算成 1 个交集
            e0_blank = sum(1 for x in e0 if not x); n0_blank = sum(1 for x in n0 if not x)
            e_bonus = 1 if e0_blank == 2 else 0; n_bonus = 1 if n0_blank == 2 else 0
            for cls, dd in dis.items():
                for name, S in dd.items():
                    ie0 = len(S & E0) if len(S) < len(E0) else sum(1 for x in E0 if x in S)
                    in0 = sum(1 for x in S if x in N0)
                    ie1 = len(S & E1) if len(S) < len(E1) else sum(1 for x in E1 if x in S)
                    a = ie0 + e_bonus; cc = in0 + n_bonus
                    if not NO_MHC:
                        w.write(f'orig_asis\t{cls}\t{name}\t{ctx}\t{a}\t{len(e0) - a}\t{cc}\t{len(n0) - cc}\n')
                        w.write(f'orig_clean\t{cls}\t{name}\t{ctx}\t{ie0}\t{len(E0) - ie0}\t{in0}\t{len(N0) - in0}\n')
                    w.write(f'eg\t{cls}\t{name}\t{ctx}\t{ie1}\t{len(E1) - ie1}\t{in0}\t{len(N0) - in0}\n')
            print(f'[{time.time() - t0:6.0f}s] {ctx} E0={len(E0)} E1={len(E1)} N={len(N0)}', flush=True)
print('counts done')
