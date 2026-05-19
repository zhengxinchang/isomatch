# parse the TSS and TES in merged isomatch gtf
import gzip
import sys
data = []
count = 0
with gzip.open(sys.argv[1], 'rt') as f:
    for line in f:
        if line.startswith('#'):
            continue
        fields = line.strip().split('\t')
        chrom, source, feature, start, end, score, strand, frame, attributes = fields
        if feature != 'transcript':
            continue
        attr_dict = {}
        for attr in attributes.split(';'):
            attr = attr.strip()
            if not attr:
                continue
            key, value = attr.split(' ', 1)
            value = value.strip('"')
            attr_dict[key] = value
        isom_src = attr_dict.get('ISOM_SRC', '')
        isom_count = attr_dict.get('ISOM_COUNT', '')
        isom_exons = attr_dict.get('ISOM_EXONS', '')
        if isom_count == '1' or isom_count == '2':
            continue
        if isom_exons == "MONO":
            continue
        count += 1
        if count > 5000:
            break
        lefts = []
        rights = []
        for tx in isom_src.split('|'):
            tx_fields = tx.split(':')
            lefts.append(int(tx_fields[2]))
            rights.append(int(tx_fields[3]))
        min_left = min(lefts)
        max_right = max(rights)
        lefts_offsets = [x - min_left for x in lefts if x!=min_left]
        rights_offsets = [max_right -x for x in rights if x!=max_right]
        for (l, r) in zip(lefts_offsets, rights_offsets):

            data.append( (l, r) )

# plot lefts_offsets and rights_offsets
import matplotlib.pyplot as plt
import numpy as np
lefts_offsets = np.array([x[0] for x in data])
rights_offsets = np.array([x[1] for x in data])
plt.scatter(lefts_offsets, rights_offsets, alpha=0.5)
plt.xlabel('TSS offset')
plt.ylabel('TES offset')

plt.title('TSS vs TES offsets for merged transcripts')
plt.grid()
plt.savefig('tss_tes_offsets.png')

# make a 2d histogram
plt.figure(figsize=(10, 8))
plt.hist2d(lefts_offsets, rights_offsets, bins=2000, cmap='Blues')
plt.colorbar(label='Count in bin')

plt.xlabel('TSS offset')
plt.ylabel('TES offset')
plt.title('2D histogram of TSS vs TES offsets for merged transcripts')
plt.grid()
plt.savefig('tss_tes_offsets_hist2d.png')   