#!/bin/bash
# 登录节点上并行下载新 GWAS（下载不计费）。先停掉顺序下载的旧进程，避免两个进程写同一个文件。
S=/public/home/yucong/scbrain
PY=/public/home/yucong/.conda/envs/scvi/bin/python
cd $S
pkill -u yucong -f "scripts/gwas_fetch.py" 2>/dev/null
sleep 2
pkill -u yucong -f "wget -q -c --tries=8 --timeout=90 -O external/gwas10" 2>/dev/null
sleep 1
for g in GCST90432156:MG GCST90558092:MS GCST90726344:MDD GCST009325:PD GCST90018844:FebrileSeizures; do
  nohup $PY scripts/gwas_fetch.py external/gwas10 $g > logs/gwas_fetch_${g#*:}.log 2>&1 < /dev/null &
done
mkdir -p external/gwas10/PGC3_SCZ external/gwas10/PGC_CUD
(
  cd external/gwas10/PGC3_SCZ
  wget -q -c --tries=8 --timeout=90 -O PGC3_SCZ_sumstats_release_v3.2_README.pdf https://ndownloader.figshare.com/files/34865091
  wget -q -c --tries=8 --timeout=90 -O PGC3_SCZ_wave3.european.autosome.public.v3.vcf.tsv.gz https://ndownloader.figshare.com/files/34517828
  echo "SCZ wget exit $? size $(stat -c %s PGC3_SCZ_wave3.european.autosome.public.v3.vcf.tsv.gz 2>/dev/null) (expect 239710564)"
) > logs/gwas_fetch_SCZ.log 2>&1 < /dev/null &
(
  cd external/gwas10/PGC_CUD
  wget -q -c --tries=8 --timeout=90 -O CUD.20210518.readme.docx https://ndownloader.figshare.com/files/28570821
  wget -q -c --tries=8 --timeout=90 -O CUD_EUR_full_public_11.14.2020.gz https://ndownloader.figshare.com/files/28570824
  echo "CUD wget exit $? size $(stat -c %s CUD_EUR_full_public_11.14.2020.gz 2>/dev/null) (expect 223085683)"
) > logs/gwas_fetch_CUD.log 2>&1 < /dev/null &
disown -a
echo launched
