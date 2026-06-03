#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
fetch_16S.py  （Biopython 版）
=============================
从 NCBI 批量获取菌种的 16S rRNA 基因序列，并以 FASTA 格式保存为 data.fasta。

策略
----
对每个物种名，检索 NCBI nucleotide 库中“最具代表性”的 16S rRNA 记录：
优先 RefSeq 的 16S rRNA 记录（NR_ 系列，长度 1200-1700 bp），
找不到时退回普通核酸记录。

输入文件每行形如：
    Acetivibrio thermocellus DSM 1313 (NC_017304)
脚本解析出物种名（属+种，必要时含 Candidatus），向 NCBI 查询其 16S rRNA。

运行环境
--------
* Python 3.6+
* Biopython： pip install biopython
* 运行机器需能访问 https://eutils.ncbi.nlm.nih.gov （正常联网即可）

使用方法
--------
NCBI 要求所有 E-utilities 请求带上邮箱，因此 --email 为必填项。

1) 小批量测试（推荐，先取前 20 条）：
       python3 fetch_16S.py --input 菌种名称.txt --email you@example.com --limit 20
2) 确认无误后跑全部：
       python3 fetch_16S.py --input 菌种名称.txt --email you@example.com
3) 若有 NCBI API key（NCBI 账号设置中申请，限速 3→10 次/秒）：
       python3 fetch_16S.py --input 菌种名称.txt --email you@example.com --api-key 你的KEY

输出
----
* data.fasta     —— 所有成功获取的序列（FASTA 格式）
* fetch_16S.log  —— 运行日志（成功 / 失败 / 跳过）
脚本支持断点续跑：若 data.fasta 已存在，已完成物种会被自动跳过。
"""

import argparse
import os
import re
import sys
import time

try:
    from Bio import Entrez, SeqIO
except ImportError:
    sys.exit("缺少 Biopython，请先运行:  pip install biopython")


# --------------------------------------------------------------------------
# 输入解析
# --------------------------------------------------------------------------
def parse_input(path):
    """读取菌种名称文件，返回 [(原始行, 物种全名, 基因组accession), ...]。"""
    entries = []
    with open(path, "r", encoding="utf-8") as fh:
        for raw in fh:
            line = raw.strip()
            if not line:
                continue
            acc = ""
            m = re.search(r"\(([A-Za-z0-9_]+)\)\s*$", line)
            if m:
                acc = m.group(1)
            # 物种名 = 去掉所有括号内容后的部分
            name = re.sub(r"\([^)]*\)", "", line).strip()
            entries.append((line, name, acc))
    return entries


def species_query(name):
    """由物种全名生成适合 NCBI 检索的物种名（属+种；Candidatus 取前三词）。"""
    n = name.replace("[", "").replace("]", "").strip()
    parts = n.split()
    if not parts:
        return n
    if parts[0].lower() == "candidatus" and len(parts) >= 3:
        return " ".join(parts[:3])
    if len(parts) >= 2:
        return " ".join(parts[:2])
    return n


# --------------------------------------------------------------------------
# NCBI 访问（Biopython Entrez）
# --------------------------------------------------------------------------
def esearch_16s(species, retries=3):
    """
    检索某物种最具代表性的 16S rRNA 记录，返回 UID。
    依次尝试：RefSeq 16S(1200-1700bp) -> 任意 16S(1200-1700bp) -> 任意 16S。
    """
    base = (f'"{species}"[Organism] AND '
            f'(16S ribosomal RNA[Title] OR 16S rRNA[Title])')
    terms = [
        base + ' AND refseq[Filter] AND ("1200"[SLEN] : "1700"[SLEN])',
        base + ' AND ("1200"[SLEN] : "1700"[SLEN])',
        base,
    ]
    for term in terms:
        for attempt in range(retries):
            try:
                handle = Entrez.esearch(
                    db="nuccore", term=term, retmax=1, sort="relevance")
                rec = Entrez.read(handle)
                handle.close()
                ids = rec.get("IdList", [])
                if ids:
                    return ids[0]
                break  # 该 term 无结果，换下一个 term
            except Exception:
                time.sleep(2.0 * (attempt + 1))
    return None


def efetch_fasta_record(uid, retries=3):
    """根据 UID 取一条 SeqRecord。"""
    for attempt in range(retries):
        try:
            handle = Entrez.efetch(
                db="nuccore", id=uid, rettype="fasta", retmode="text")
            record = SeqIO.read(handle, "fasta")
            handle.close()
            return record
        except Exception:
            time.sleep(2.0 * (attempt + 1))
    return None


# --------------------------------------------------------------------------
# 断点续跑
# --------------------------------------------------------------------------
def load_done(fasta_path):
    """从已存在的 data.fasta 读出已完成的物种名集合。"""
    done = set()
    if not os.path.exists(fasta_path):
        return done
    with open(fasta_path, "r", encoding="utf-8") as fh:
        for line in fh:
            if line.startswith(">"):
                m = re.match(r">([^|]+)\|", line)
                if m:
                    done.add(m.group(1).strip())
    return done


# --------------------------------------------------------------------------
# 主流程
# --------------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser(
        description="基于 Biopython 批量获取 16S rRNA 序列并存为 data.fasta")
    ap.add_argument("--input", required=True, help="菌种名称 txt 文件路径")
    ap.add_argument("--email", required=True,
                    help="联系邮箱（NCBI E-utilities 强制要求）")
    ap.add_argument("--output", default="data.fasta", help="输出 FASTA 文件名")
    ap.add_argument("--log", default="fetch_16S.log", help="日志文件名")
    ap.add_argument("--limit", type=int, default=0,
                    help="只处理前 N 条（0=全部，用于测试）")
    ap.add_argument("--api-key", default="", help="NCBI API key（可选，提速）")
    args = ap.parse_args()

    # 配置 Entrez
    Entrez.email = args.email
    if args.api_key:
        Entrez.api_key = args.api_key
    # 有 key 时 10 次/秒，否则 3 次/秒，留余量
    delay = 0.12 if args.api_key else 0.40

    entries = parse_input(args.input)
    if args.limit > 0:
        entries = entries[:args.limit]

    # 按物种名去重，保留首次出现顺序
    seen, unique = set(), []
    for raw, name, acc in entries:
        sp = species_query(name)
        if sp in seen:
            continue
        seen.add(sp)
        unique.append((sp, acc))

    done = load_done(args.output)
    print(f"输入条目: {len(entries)}  去重后物种: {len(unique)}  "
          f"已完成: {len(done)}")

    ok = fail = skip = 0
    with open(args.output, "a", encoding="utf-8") as out, \
         open(args.log, "a", encoding="utf-8") as log:
        log.write(f"\n=== run {time.strftime('%Y-%m-%d %H:%M:%S')} ===\n")
        for i, (sp, acc) in enumerate(unique, 1):
            if sp in done:
                skip += 1
                continue
            tag = f"[{i}/{len(unique)}] {sp}"
            try:
                uid = esearch_16s(sp)
                time.sleep(delay)
                if not uid:
                    fail += 1
                    msg = f"FAIL  {tag}  -- 未找到 16S 记录"
                    print(msg); log.write(msg + "\n")
                    continue
                record = efetch_fasta_record(uid)
                time.sleep(delay)
                if record is None or len(record.seq) == 0:
                    fail += 1
                    msg = f"FAIL  {tag}  -- 序列获取失败"
                    print(msg); log.write(msg + "\n")
                    continue
                # 写出 FASTA：header 含物种名，便于后续识别
                out.write(f">{sp} | {record.description}\n")
                seq = str(record.seq)
                for j in range(0, len(seq), 70):       # 每行 70 个碱基
                    out.write(seq[j:j + 70] + "\n")
                out.flush()
                ok += 1
                msg = f"OK    {tag}  -- {record.id} ({len(seq)} bp)"
                print(msg); log.write(msg + "\n")
            except Exception as e:
                fail += 1
                msg = f"ERROR {tag}  -- {e}"
                print(msg); log.write(msg + "\n")
                time.sleep(2.0)

        summary = (f"\n完成: 成功 {ok}  失败 {fail}  跳过(已存在) {skip}\n"
                   f"输出文件: {os.path.abspath(args.output)}\n")
        print(summary)
        log.write(summary)


if __name__ == "__main__":
    main()
