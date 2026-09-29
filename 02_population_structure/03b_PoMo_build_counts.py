# 科学目的：按 46 个 Pop_L3 汇总 final 359 个体 all-sites VCF 为 PoMo COUNTSFILE。
# 输入：publication sample_metadata.csv、外部 populations.all.vcf；输出：populations.cf。
# 流程：匹配 VCF 样本 → 对每条 VCF record 逐群统计 A/C/G/T → 保存固定 NSITES header。

import csv
from pathlib import Path

# Load data
with Path("03_Public_Code/data/metadata/sample_metadata.csv").open() as stream:
    sample_to_pop = {}
    for row in csv.DictReader(stream):
        sample = row["SampleID"]
        assert sample not in sample_to_pop, f"重复 SampleID: {sample}"
        sample_to_pop[sample] = row["Pop_L3"]

vcf_path = Path("01_Data/vcf_260905_pomo/populations.all.vcf")
with vcf_path.open() as stream:
    for line in stream:
        if line.startswith("#CHROM"):
            vcf_samples = line.rstrip("\n").split("\t")[9:]
            break
    else:
        raise ValueError("VCF #CHROM header missing")

# Prepare data: sample order comes from the VCF, not the metadata table.
assert len(vcf_samples) == 359 and len(set(vcf_samples)) == 359
assert all(sample in sample_to_pop for sample in vcf_samples)
sample_pops = [sample_to_pop[sample] for sample in vcf_samples]
populations = sorted(set(sample_pops))
assert len(populations) == 46
base_index = {"A": 0, "C": 1, "G": 2, "T": 3}

# Analysis / Save results: preserve every record, including fixed sites; missing GT adds no counts.
output_path = Path("04_Results_Archive/03_ML_Tree/PoMo_population/populations.cf")
output_path.parent.mkdir(parents=True, exist_ok=True)
with vcf_path.open() as vcf, output_path.open("w") as out:
    out.write(f"COUNTSFILE NPOP {len(populations)} NSITES 751805\n")
    out.write("\t".join(["CHROM", "POS"] + populations) + "\n")
    for line in vcf:
        if line.startswith("#"):
            continue
        columns = line.rstrip("\n").split("\t")
        chrom, pos, ref, alt = columns[0], columns[1], columns[3], columns[4]
        counts = {pop: [0, 0, 0, 0] for pop in populations}
        for pop, field in zip(sample_pops, columns[9:]):
            gt = field.split(":")[0]
            if gt == "./.":
                continue
            elif gt == "0/0":
                alleles = (ref, ref)
            elif gt == "0/1":
                alleles = (ref, alt)
            elif gt == "1/1":
                alleles = (alt, alt)
            else:
                raise ValueError(f"Unexpected genotype: {gt}")
            for base in alleles:
                counts[pop][base_index[base]] += 1
        row = [chrom, pos] + [",".join(map(str, counts[pop])) for pop in populations]
        out.write("\t".join(row) + "\n")
