# Merge Flowchart

This diagram summarizes the implemented merge logic in `src/merge`.
It includes only parameters that are actively used by the current merge code path.

```mermaid
flowchart TD
    A[Load indexed inputs .isomx] --> B[Collect chromosome names]
    B --> C[For each chromosome]
    C --> D[K-way merge PTIRs by genomic order]
    D --> E{Overlaps current super cluster?}
    E -- Yes --> F[Append PTIR to current super cluster<br/>Update cluster_max_end]
    E -- No --> G[Process previous super cluster]
    F --> E
    G --> H[Partition by strand and exon count<br/>key = strand + n_exons]

    H --> I{n_exons == 1?}

    I -- Yes --> J[Mono-exon merge]
    J --> K[Union-find by reciprocal overlap<br/>threshold: mono_ovlp]
    K --> L[Select representative boundaries<br/>policy: mono_policy]
    L --> M[Write merged transcript]

    I -- No --> N[Split transcripts by PTIR type]
    N --> O[Canonical: ALLC]
    N --> P[Non-canonical pool: PRTC + NOTC]

    O --> Q[Canonical splice merge]
    Q --> R[Union-find by junction similarity<br/>thresholds: wob_a, wob_d, wob_u]
    R --> S[Optional terminal refinement<br/>mode: terminal_merge<br/>thresholds: tss_wob, tes_wob]
    S --> T[Profile canonical group]
    T --> U[Select representative splice junctions<br/>policy: splice_policy]
    U --> V[Select representative TSS/TES]

    V --> W{Guide hits available?}
    W -- Both TSS and TES hit --> X[Choose among doubly-guided entries<br/>using tss_policy + tes_policy]
    W -- Partial guide hit --> Y[Choose highest guide-hit score<br/>Tie-break by longer transcript]
    W -- No guide hit --> Z[Choose by tss_policy + tes_policy]
    X --> AA[Canonical backbone ready]
    Y --> AA
    Z --> AA

    P --> AB[Try to absorb non-canonical transcripts<br/>into canonical backbones]
    AA --> AB
    AB --> AC{Matches a canonical backbone?}
    AC -- Yes --> AD[Attach to best backbone<br/>minimum donor+acceptor diff<br/>thresholds: wob_a_nc, wob_d_nc, wob_u_nc]
    AC -- No --> AE[Keep as residual non-canonical]

    AE --> AF[Merge residual non-canonical transcripts]
    AF --> AG[Union-find by junction similarity<br/>thresholds: wob_a_nc, wob_d_nc, wob_u_nc]
    AG --> AH[Optional terminal refinement<br/>mode: terminal_merge_nc<br/>thresholds: tss_wob_nc, tes_wob_nc]
    AH --> AI[Profile non-canonical group]
    AI --> AJ[Select representative splice junctions<br/>policy: splice_policy]
    AJ --> AK[Select representative TSS/TES<br/>guide-aware, same as canonical]
    AD --> AL[Write merged transcript]
    AK --> AL
    M --> AM[Continue until chromosome ends]
    AL --> AM
```

## Implemented Merge Parameters

- Guide-aware representative selection:
  - `guide_tss`
  - `guide_tes`
  - `guide_tss_flank`
  - `guide_tes_flank`
  - `chrmap`
- Canonical multi-exon merge:
  - `wob_a`
  - `wob_d`
  - `wob_u`
  - `tss_wob`
  - `tes_wob`
  - `terminal_merge`
- Non-canonical multi-exon merge:
  - `wob_a_nc`
  - `wob_d_nc`
  - `wob_u_nc`
  - `tss_wob_nc`
  - `tes_wob_nc`
  - `terminal_merge_nc`
- Representative selection:
  - `splice_policy`
  - `tss_policy`
  - `tes_policy`
  - `mono_policy`
- Mono-exon merge:
  - `mono_ovlp`

## Explicitly Excluded

These arguments exist in `MergeArgs` but are not currently used by the implemented merge path:

- `sx_max`
- `junc_diff`
- `shift_rescue`
