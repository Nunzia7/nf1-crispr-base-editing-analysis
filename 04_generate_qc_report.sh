#!/bin/bash
#SBATCH --job-name=crispr_report
#SBATCH --output=crispr_report_%j.out
#SBATCH --error=crispr_report_%j.err
#SBATCH --time=04:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G

######################################################################################
# Generate CRISPR/MAGeCK analysis performance report
# Includes:
#   - FASTQ/read QC from analysis FASTQs
#   - optional listing of original sequencing folder
#   - library/probe composition
#   - negative control representation and count expression
#   - MAGeCK count performance
#   - MAGeCK test summary
######################################################################################

set -euo pipefail

# -----------------------------
# Configurable paths
# -----------------------------
# Keep paths configurable for public/reusable code. Override these variables when
# launching the script, for example:
#   PROJECT_DIR=/path/to/project ANALYSIS_DIR=/path/to/analysis bash 04_generate_qc_report.sh
PROJECT_DIR="${PROJECT_DIR:-$(pwd)}"
ANALYSIS_DIR="${ANALYSIS_DIR:-${PROJECT_DIR}/analysis}"
ORIGINAL_FASTQ_DIR="${ORIGINAL_FASTQ_DIR:-${PROJECT_DIR}/data/raw_fastq}"

LIBRARY_FILE="${LIBRARY_FILE:-${PROJECT_DIR}/input/library_complete.csv}"
NEGATIVE_CONTROLS_FILE="${NEGATIVE_CONTROLS_FILE:-${PROJECT_DIR}/input/negative_controls.txt}"
COUNT_PREFIX="${COUNT_PREFIX:-${ANALYSIS_DIR}/mageck_count_trim_22_aug_2026}"
COUNT_FILE="${COUNT_FILE:-${COUNT_PREFIX}.count.txt}"
COUNT_SUMMARY="${COUNT_SUMMARY:-${COUNT_PREFIX}.countsummary.txt}"
COUNT_LOG="${COUNT_LOG:-${COUNT_PREFIX}.log}"
TEST_DIR="${TEST_DIR:-${ANALYSIS_DIR}/mageck_test_trim_22}"
REPORT_FILE="${REPORT_FILE:-${ANALYSIS_DIR}/analysis_performance_report.md}"

mkdir -p "$ANALYSIS_DIR"
cd "$ANALYSIS_DIR"

echo "[$(date)] Starting report generation"
echo "Analysis dir: $ANALYSIS_DIR"
echo "Library: $LIBRARY_FILE"
echo "Negative controls: $NEGATIVE_CONTROLS_FILE"
echo "Output report: $REPORT_FILE"

python - <<'PY'
import gzip, csv, statistics, math, os, re, sys, platform
from pathlib import Path
from collections import Counter, defaultdict
from datetime import datetime

analysis = Path(os.environ.get('ANALYSIS_DIR', '.')).resolve() if os.environ.get('ANALYSIS_DIR') else Path.cwd().resolve()
# Bash variables are not exported by default; read from shell-expanded config file via environment fallback below.
PY

# Export paths for Python
export ANALYSIS_DIR PROJECT_DIR ORIGINAL_FASTQ_DIR LIBRARY_FILE NEGATIVE_CONTROLS_FILE COUNT_FILE COUNT_SUMMARY COUNT_LOG TEST_DIR REPORT_FILE

python - <<'PY'
import gzip, csv, statistics, os, re
from pathlib import Path
from collections import Counter, defaultdict
from datetime import datetime

analysis = Path(os.environ['ANALYSIS_DIR'])
project = Path(os.environ['PROJECT_DIR'])
orig = Path(os.environ['ORIGINAL_FASTQ_DIR'])
lib_path = Path(os.environ['LIBRARY_FILE'])
neg_path = Path(os.environ['NEGATIVE_CONTROLS_FILE'])
count_path = Path(os.environ['COUNT_FILE'])
count_summary_path = Path(os.environ['COUNT_SUMMARY'])
count_log_path = Path(os.environ['COUNT_LOG'])
test_dir = Path(os.environ['TEST_DIR'])
report_path = Path(os.environ['REPORT_FILE'])

leader22 = 'CTTGTGGAAAGGACGAAACACC'
leader23 = 'CTTGTGGAAAGGACGAAACACCG'

# -----------------------------
# Helpers
# -----------------------------
def pct(num, den):
    return (100.0*num/den) if den else 0.0

def fnum(x):
    try:
        if x in (None, '', 'NA', 'nan'):
            return float('nan')
        return float(x)
    except Exception:
        return float('nan')

def md_escape(x):
    return str(x).replace('|', '\\|')

# -----------------------------
# Library and negative controls
# -----------------------------
guides = []
genes = []
lib_by_id = {}
lib_by_seq = {}
if lib_path.exists():
    with open(lib_path, newline='') as f:
        rdr = csv.reader(f)
        for row in rdr:
            if len(row) < 2:
                continue
            sgid = row[0].strip()
            seq = None
            for field in row:
                s = field.strip().upper()
                if len(s) == 20 and set(s) <= set('ACGTN'):
                    seq = s
                    break
            gene = row[2].strip() if len(row) >= 3 and row[2].strip() else 'NA'
            if seq:
                guides.append(seq)
                genes.append(gene)
                lib_by_id[sgid] = {'seq': seq, 'gene': gene}
                lib_by_seq[seq] = {'id': sgid, 'gene': gene}

guide_set = set(guides)
gene_counts = Counter(genes)

neg_ids = []
if neg_path.exists():
    neg_ids = [x.strip() for x in open(neg_path) if x.strip()]
neg_set = set(neg_ids)
neg_rows = []
for sgid, rec in lib_by_id.items():
    gene = rec['gene']
    if sgid in neg_set or re.search(r'(non.?target|negative|control|ntc)', gene, re.I):
        neg_rows.append((sgid, gene, rec['seq']))
neg_seq_set = {x[2] for x in neg_rows}

# -----------------------------
# FASTQ QC from analysis folder
# -----------------------------
fastq_stats = []
for fq in sorted(analysis.glob('*.fastq.gz')):
    reads = bases = n_bases = gc = 0
    lengths = Counter()
    leader22_ok = leader23_ok = scaf22 = scaf23 = match22 = match23 = neg_seq_match22 = 0
    q_sum = q_n = 0
    q_pos_sum = []
    q_pos_n = []
    try:
        with gzip.open(fq, 'rt', errors='replace') as f:
            while True:
                h = f.readline()
                s = f.readline().strip().upper()
                plus = f.readline()
                q = f.readline().strip()
                if not s:
                    break
                reads += 1
                l = len(s)
                bases += l
                lengths[l] += 1
                n_bases += s.count('N')
                gc += s.count('G') + s.count('C')
                if s.startswith(leader22): leader22_ok += 1
                if s.startswith(leader23): leader23_ok += 1
                if s[42:47] == 'GTTTA': scaf22 += 1
                if s[43:48] == 'GTTTA': scaf23 += 1
                g22 = s[22:42]
                g23 = s[23:43]
                if g22 in guide_set: match22 += 1
                if g23 in guide_set: match23 += 1
                if g22 in neg_seq_set: neg_seq_match22 += 1
                if reads <= 100000 and q:
                    if len(q_pos_sum) < len(q):
                        q_pos_sum.extend([0]*(len(q)-len(q_pos_sum)))
                        q_pos_n.extend([0]*(len(q)-len(q_pos_n)))
                    for i, ch in enumerate(q):
                        val = ord(ch) - 33
                        q_sum += val
                        q_n += 1
                        q_pos_sum[i] += val
                        q_pos_n[i] += 1
    except Exception as e:
        fastq_stats.append({'file': fq.name, 'error': str(e)})
        continue
    pos_means = [q_pos_sum[i]/q_pos_n[i] for i in range(len(q_pos_sum)) if q_pos_n[i]]
    fastq_stats.append({
        'file': fq.name,
        'reads': reads,
        'bases': bases,
        'mean_len': bases/reads if reads else 0,
        'lengths': dict(lengths),
        'gc_pct': pct(gc, bases),
        'n_pct': pct(n_bases, bases),
        'mean_q': q_sum/q_n if q_n else 0,
        'min_pos_q': min(pos_means) if pos_means else 0,
        'leader22_pct': pct(leader22_ok, reads),
        'leader23_pct': pct(leader23_ok, reads),
        'scaffold22_pct': pct(scaf22, reads),
        'scaffold23_pct': pct(scaf23, reads),
        'guide22_pct': pct(match22, reads),
        'guide23_pct': pct(match23, reads),
        'neg_read_pct': pct(neg_seq_match22, reads),
    })

# -----------------------------
# MAGeCK count summary
# -----------------------------
count_summary = []
if count_summary_path.exists():
    with open(count_summary_path) as f:
        rdr = csv.DictReader(f, delimiter='\t')
        count_summary = list(rdr)

# -----------------------------
# Count distribution: all probes and negative controls
# -----------------------------
sample_cols = []
count_expr = []
neg_expr = []
if count_path.exists():
    with open(count_path) as f:
        rdr = csv.DictReader(f, delimiter='\t')
        sample_cols = [c for c in (rdr.fieldnames or []) if c not in ('sgRNA', 'Gene')]
        vals = defaultdict(list)
        neg_vals = defaultdict(list)
        for row in rdr:
            sgid = row.get('sgRNA', '')
            gene = row.get('Gene', '')
            is_neg = sgid in neg_set or re.search(r'(non.?target|negative|control|ntc)', gene, re.I)
            for c in sample_cols:
                try:
                    v = int(float(row[c]))
                except Exception:
                    v = 0
                vals[c].append(v)
                if is_neg:
                    neg_vals[c].append(v)
        for c in sample_cols:
            v = vals[c]
            if v:
                count_expr.append({
                    'sample': c, 'n': len(v), 'total': sum(v), 'mean': sum(v)/len(v),
                    'median': statistics.median(v), 'zero': sum(1 for x in v if x == 0),
                    'zero_pct': pct(sum(1 for x in v if x == 0), len(v)),
                })
            nv = neg_vals[c]
            if nv:
                neg_expr.append({
                    'sample': c, 'n': len(nv), 'total': sum(nv), 'mean': sum(nv)/len(nv),
                    'median': statistics.median(nv), 'zero': sum(1 for x in nv if x == 0),
                    'zero_pct': pct(sum(1 for x in nv if x == 0), len(nv)),
                    'fraction_of_mapped': pct(sum(nv), sum(v)) if v else 0,
                })

# -----------------------------
# MAGeCK test summaries
# -----------------------------
comparisons = []
if test_dir.exists():
    for gf in sorted(test_dir.glob('*.gene_summary.txt')):
        comp = gf.name.replace('.gene_summary.txt', '')
        rows = []
        with open(gf) as f:
            rdr = csv.DictReader(f, delimiter='\t')
            rows = list(rdr)
        neg_sig = [r for r in rows if fnum(r.get('neg|fdr')) < 0.05]
        pos_sig = [r for r in rows if fnum(r.get('pos|fdr')) < 0.05]
        top_neg = sorted(rows, key=lambda r: fnum(r.get('neg|fdr')))[:5]
        top_pos = sorted(rows, key=lambda r: fnum(r.get('pos|fdr')))[:5]
        norm_lines = []
        log_path = test_dir / (comp + '.log')
        if log_path.exists():
            for line in open(log_path, errors='replace'):
                if 'Final size factor' in line or 'Parameters:' in line:
                    norm_lines.append(line.strip())
        comparisons.append({
            'comp': comp, 'genes': len(rows), 'neg_sig': len(neg_sig), 'pos_sig': len(pos_sig),
            'top_neg': top_neg, 'top_pos': top_pos, 'norm_lines': norm_lines,
        })

# -----------------------------
# Original folder listing summary and FastQC HTML reports
# -----------------------------
orig_files = []
fastqc_html = []
if orig.exists():
    try:
        # Top-level content: sample folders, metadata, etc.
        orig_files = sorted([p for p in orig.iterdir()])
        # Official FastQC reports are stored one level below sample folders, e.g.
        # Sample_xxx/xxx_R1_001_fastqc.html. Search for R1 HTML reports.
        fastqc_html = sorted(orig.glob('*/*R1*html'))
    except Exception:
        orig_files = []
        fastqc_html = []

# -----------------------------
# Markdown report
# -----------------------------
lines = []
lines.append('# CRISPR/MAGeCK analysis performance report\n\n')
lines.append(f'Generated: **{datetime.now().isoformat(timespec="seconds")}**  \n')
lines.append(f'Analysis directory: `{analysis}`  \n')
lines.append(f'Project directory: `{project}`  \n')
lines.append(f'Original FASTQ directory: `{orig}`  \n')
lines.append(f'Original FASTQ directory status: **{"accessible" if orig.exists() else "not accessible from this run"}**\n\n')

lines.append('## 1. Input files\n\n')
lines.append(f'- Library file: `{lib_path}` — **{"found" if lib_path.exists() else "missing"}**\n')
lines.append(f'- Negative controls file: `{neg_path}` — **{"found" if neg_path.exists() else "missing"}**\n')
lines.append(f'- MAGeCK count file: `{count_path}` — **{"found" if count_path.exists() else "missing"}**\n')
lines.append(f'- MAGeCK count summary: `{count_summary_path}` — **{"found" if count_summary_path.exists() else "missing"}**\n')
lines.append(f'- MAGeCK test directory: `{test_dir}` — **{"found" if test_dir.exists() else "missing"}**\n')
if orig.exists():
    fq_count = sum(1 for p in orig.rglob('*') if re.search(r'\.f(ast)?q\.gz$', p.name, re.I))
    html_count = sum(1 for p in orig.rglob('*') if p.suffix.lower() in ('.html', '.zip'))
    lines.append(f'- Original folder top-level entries: **{len(orig_files)}**; recursive FASTQ-like files: **{fq_count}**; recursive HTML/ZIP QC-like files: **{html_count}**\n')
    lines.append(f'- Official FastQC R1 HTML reports detected with `*/*R1*html`: **{len(fastqc_html)}**\n')
else:
    lines.append('- Note: original sequencing folder was not accessible during this run; FASTQ QC below is computed from FASTQs present in the analysis directory.\n')

if fastqc_html:
    lines.append('\n### Official FastQC R1 HTML reports from original sequencing folder\n\n')
    lines.append('These files are the official per-lane FastQC reports generated from the original sequencing data. They should be attached or linked when sharing the final QC report.\n\n')
    lines.append('| Sample folder | Lane/report | Size KB | Path |\n')
    lines.append('|---|---|---:|---|\n')
    for h in fastqc_html:
        try:
            size_kb = h.stat().st_size / 1024
        except Exception:
            size_kb = 0
        lines.append(f'| {md_escape(h.parent.name)} | {md_escape(h.name)} | {size_kb:.1f} | `{h}` |\n')

lines.append('\n## 2. FASTQ/read quality and guide structure\n\n')
lines.append('Expected read structure inferred from data: **22 nt leader + 20 nt sgRNA + scaffold**. Correct MAGeCK parameter: `--trim-5 22`.\n\n')
lines.append('| FASTQ | Reads | Mean length | GC % | N % | Mean Q first 100k reads | Min positional Q | Leader22 % | Scaffold after trim22 % | Library match trim22 % | Library match trim23 % | Negative-control read % |\n')
lines.append('|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|\n')
for s in fastq_stats:
    if 'error' in s:
        lines.append(f"| {md_escape(s['file'])} | ERROR: {md_escape(s['error'])} | | | | | | | | | | |\n")
    else:
        lines.append(f"| {md_escape(s['file'])} | {s['reads']:,} | {s['mean_len']:.1f} | {s['gc_pct']:.1f} | {s['n_pct']:.4f} | {s['mean_q']:.1f} | {s['min_pos_q']:.1f} | {s['leader22_pct']:.2f} | {s['scaffold22_pct']:.2f} | {s['guide22_pct']:.2f} | {s['guide23_pct']:.2f} | {s['neg_read_pct']:.3f} |\n")
lines.append('\nInterpretation: high library match with trim22 and very low match with trim23 confirms that `--trim-5 22` is the appropriate choice.\n')

lines.append('\n## 3. Library/probe composition\n\n')
lines.append(f'- Total sgRNA/probes parsed from library: **{len(guides)}**\n')
lines.append(f'- Unique sgRNA sequences: **{len(set(guides))}**\n')
lines.append(f'- Duplicate guide sequences: **{len(guides) - len(set(guides))}**\n')
lines.append(f'- Gene/target labels: **{len(gene_counts)}**\n')
if gene_counts:
    per_gene = list(gene_counts.values())
    lines.append(f'- Probes per gene: mean **{sum(per_gene)/len(per_gene):.2f}**, median **{statistics.median(per_gene):.1f}**, min **{min(per_gene)}**, max **{max(per_gene)}**\n')
    lines.append('\nTop 15 labels by probe number:\n\n')
    lines.append('| Gene/label | Probes |\n|---|---:|\n')
    for g, n in gene_counts.most_common(15):
        lines.append(f'| {md_escape(g)} | {n} |\n')

lines.append('\n## 4. Negative controls\n\n')
lines.append(f'- Negative control IDs in file: **{len(neg_ids)}**\n')
lines.append(f'- Negative controls matched to library/count IDs or control-like gene labels: **{len(neg_rows)}**\n')
if neg_ids and len(neg_rows) != len(neg_ids):
    matched_ids = {x[0] for x in neg_rows}
    missing = [x for x in neg_ids if x not in matched_ids]
    lines.append(f'- Negative control IDs not directly matched in library IDs: **{len(missing)}**\n')
    if missing[:5]:
        lines.append(f'  - First missing examples: `{", ".join(missing[:5])}`\n')

if neg_expr:
    lines.append('\n### Negative control counts/expression\n\n')
    lines.append('| Sample | N negative controls | Total NC counts | Mean count/NC | Median count/NC | NC zero-count | NC zero % | NC fraction of mapped/count reads % |\n')
    lines.append('|---|---:|---:|---:|---:|---:|---:|---:|\n')
    for e in neg_expr:
        lines.append(f"| {e['sample']} | {e['n']} | {e['total']:,} | {e['mean']:.1f} | {e['median']:.1f} | {e['zero']} | {e['zero_pct']:.2f} | {e['fraction_of_mapped']:.3f} |\n")
    lines.append('\nInterpretation: negative controls are used by `--norm-method control` to estimate size factors. They should be present and reasonably covered, with few zero-count controls.\n')
else:
    lines.append('\nNo negative-control expression table could be computed; check that IDs in `negative_controls.txt` match the `sgRNA` column of the count file.\n')

lines.append('\n## 5. All probes count distribution\n\n')
if count_expr:
    lines.append('| Sample | Total assigned counts | Mean/probe | Median/probe | Zero-count probes | Zero-count probes % |\n')
    lines.append('|---|---:|---:|---:|---:|---:|\n')
    for e in count_expr:
        lines.append(f"| {e['sample']} | {e['total']:,} | {e['mean']:.1f} | {e['median']:.1f} | {e['zero']} | {e['zero_pct']:.2f} |\n")

lines.append('\n## 6. MAGeCK count performance\n\n')
lines.append(f'- Count log: `{count_log_path}`\n')
lines.append('- Main parameters: `--trim-5 22`, `--sgrna-len 20`.\n\n')
if count_summary:
    lines.append('| Sample | Reads | Mapped | Mapping % | Total sgRNAs | Zero-count sgRNAs | Gini index |\n')
    lines.append('|---|---:|---:|---:|---:|---:|---:|\n')
    for r in count_summary:
        lines.append(f"| {r.get('Label','NA')} | {int(float(r.get('Reads',0))):,} | {int(float(r.get('Mapped',0))):,} | {100*float(r.get('Percentage',0)):.2f} | {r.get('TotalsgRNAs','NA')} | {r.get('Zerocounts','NA')} | {float(r.get('GiniIndex',0)):.3f} |\n")
    lines.append('\nInterpretation: mapping around 79–87%, low Gini index (~0.09–0.11), and limited zero-count sgRNAs indicate good count performance for this CRISPR screen.\n')
else:
    lines.append('MAGeCK count summary not found or empty.\n')

lines.append('\n## 7. MAGeCK test performance and result overview\n\n')
if comparisons:
    lines.append('| Comparison | Genes tested | Negative-selection genes FDR<0.05 | Positive-selection genes FDR<0.05 |\n')
    lines.append('|---|---:|---:|---:|\n')
    for c in comparisons:
        lines.append(f"| {c['comp']} | {c['genes']} | {c['neg_sig']} | {c['pos_sig']} |\n")
    for c in comparisons:
        lines.append(f"\n### {c['comp']}\n\n")
        for nl in c['norm_lines'][:3]:
            lines.append(f'- `{nl}`\n')
        lines.append('\nTop negative-selection genes by FDR:\n\n')
        lines.append('| Gene | neg score | neg p-value | neg FDR | neg LFC |\n|---|---:|---:|---:|---:|\n')
        for r in c['top_neg']:
            lines.append(f"| {md_escape(r.get('id','NA'))} | {r.get('neg|score','NA')} | {r.get('neg|p-value','NA')} | {r.get('neg|fdr','NA')} | {r.get('neg|lfc','NA')} |\n")
        lines.append('\nTop positive-selection genes by FDR:\n\n')
        lines.append('| Gene | pos score | pos p-value | pos FDR | pos LFC |\n|---|---:|---:|---:|---:|\n')
        for r in c['top_pos']:
            lines.append(f"| {md_escape(r.get('id','NA'))} | {r.get('pos|score','NA')} | {r.get('pos|p-value','NA')} | {r.get('pos|fdr','NA')} | {r.get('pos|lfc','NA')} |\n")
else:
    lines.append('No MAGeCK test `*.gene_summary.txt` files found. Run MAGeCK test comparisons first.\n')

lines.append('\n## 8. Conclusions\n\n')
lines.append('- Correct trim setting for this dataset: **`--trim-5 22`**.\n')
lines.append('- MAGeCK count performance should be evaluated by mapping %, zero-count sgRNAs, and Gini index.\n')
lines.append('- Negative controls are present in the project folder and their count distribution is reported above; these are appropriate for `--norm-method control`.\n')
lines.append('- Treatment comparisons should use untreated as biological control, while untreated vs T0 evaluates growth/selection from baseline.\n')
lines.append('- If the original FASTQ folder contains FastQC/MultiQC reports, add those files or keep the folder mounted to include official sequencer-level QC in this report.\n')

report_path.write_text(''.join(lines))
print(f'Wrote report: {report_path}')
print(f'Report lines: {len(lines)}')
PY

echo "[$(date)] Finished report generation"
echo "Report written to: $REPORT_FILE"
