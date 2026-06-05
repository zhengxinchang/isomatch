isomatch largely follows SQANTI3 classification rules, assigning transcripts to FSM, ISM, NIC, NNC, and other categories with compatible subcategories. Key differences are noted below.

1. **Unstranded query transcripts.** SQANTI3 raises an error; isomatch instead assigns them to a new category `invalid_query_tx` (subcategory `unstranded_query_tx`), improving robustness.

2. **NNC and genic overlap.** SQANTI3 does not require splice-site overlap for NNC. isomatch requires at least one query splice site to overlap a reference splice site, ensuring the NNC transcript is more confidently associated with the assigned gene. This distinction is also discussed in [SQANTI3 issue #427](https://github.com/ConesaLab/SQANTI3/issues/427).

3. **Intron retention at the 5′ end.** Consider a transcript whose first and second exons are merged (intron retention) while all remaining junctions match a reference transcript. SQANTI3 classifies this as NIC (`intron_retention`) because no anchor exon is found during junction matching. isomatch classifies it as ISM (`intron_retention`), since all query junctions are a subset of the reference junctions.