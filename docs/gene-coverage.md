# Gene coverage

Generated on 2026-11-14 from `docker/genes.tsv` (tool support and SV flags) and the `AUTHORITATIVE` map in
`docker/pgx-compare.py`. Regenerate it when either changes. ✓ = the tool runs for the gene;
★ = that method is the verdict authority when it makes a call (Figure 2 in the README).
PharmCAT also runs once per sample and is shown as a cross-check on other genes; ABO is typed
by the ABO read pileup and reported through the VCF-Check path. CYPOR is an alias of POR and is
not run separately.

| Gene | Region (GRCh38) | PyPGx | Stargazer | Aldy | StellarPGx | Cyrius | PharmCAT | VCF-Check | OptiType | arcasHLA | mutserve | SV-aware |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| CYP1A1 | chr15:74716541-74728528 | ✓ | ✓ | ✓ | ✓ |  |  |  |  |  |  |  |
| CYP1A2 | chr15:74745844-74759607 | ✓ | ✓ | ✓ | ✓ |  |  |  |  |  |  |  |
| CYP2A6 | chr19:40833540-40890447 | ✓ | ✓ | ✓ | ✓ |  |  |  |  |  |  | yes |
| CYP2B6 | chr19:40921281-41028398 | ✓ | ✓ | ✓ | ✓ |  | ★ |  |  |  |  | yes |
| CYP2C8 | chr10:95033771-95072497 | ✓ | ✓ | ✓ | ✓ |  |  |  |  |  |  |  |
| CYP2C9 | chr10:94935657-94993091 | ✓ | ✓ | ✓ | ✓ |  |  |  |  |  |  |  |
| CYP2C19 | chr10:94759680-94858547 | ✓ | ✓ | ✓ | ✓ |  |  |  |  |  |  |  |
| CYP2D6 | chr22:42116498-42155810 | ✓ | ✓ | ✓ | ✓ | ★ |  |  |  |  |  | yes |
| CYP2E1 | chr10:133517362-133549123 | ✓ | ✓ | ✓ | ✓ |  |  |  |  |  |  | yes |
| CYP3A4 | chr7:99753966-99787184 | ✓ | ✓ | ✓ | ✓ |  |  |  |  |  |  |  |
| CYP3A5 | chr7:99645193-99682996 | ✓ | ✓ | ✓ | ✓ |  |  |  |  |  |  |  |
| CYP4F2 | chr19:15863022-15913074 | ✓ | ✓ | ✓ | ✓ |  | ★ |  |  |  |  | yes |
| DPYD | chr1:97074742-97924034 | ✓ | ✓ | ✓ |  |  |  |  |  |  |  |  |
| G6PD | chrX:154528389-154550018 | ✓ | ✓ | ✓ |  |  |  | ★ |  |  |  | yes |
| GSTM1 | chr1:109684816-109696745 | ✓ | ✓ | ✓ | ✓ |  |  |  |  |  |  | yes |
| GSTT1 | ALT_CONTIG | ✓ |  |  | ✓ |  |  |  |  |  |  | yes |
| IFNL3 | chr19:39240552-39253525 | ✓ | ✓ | ✓ |  |  |  | ★ |  |  |  |  |
| NAT1 | chr8:18207108-18226689 | ✓ | ✓ | ✓ | ✓ |  |  |  |  |  |  |  |
| NAT2 | chr8:18388281-18404218 | ✓ | ✓ | ✓ | ✓ |  |  |  |  |  |  |  |
| NUDT15 | chr13:48034725-48050221 | ✓ | ✓ | ✓ | ✓ |  |  |  |  |  |  |  |
| POR | chr7:75912154-75989855 | ✓ | ✓ |  | ✓ |  |  |  |  |  |  |  |
| RYR1 | chr19:38430690-38590564 | ✓ | ✓ | ✓ |  |  |  | ★ |  |  |  |  |
| SLCO1B1 | chr12:21128193-21242796 | ✓ | ✓ | ✓ | ✓ |  |  |  |  |  |  |  |
| TPMT | chr6:18125310-18158169 | ✓ | ✓ | ✓ | ✓ |  |  |  |  |  |  |  |
| UGT1A1 | chr2:233754269-233779300 | ✓ | ✓ | ✓ | ✓ |  | ★ | ✓ |  |  |  |  |
| VKORC1 | chr16:31087853-31097797 | ✓ | ✓ | ✓ |  |  |  | ★ |  |  |  |  |
| ABCG2 | chr4:88085265-88236626 |  |  | ✓ |  |  |  |  |  |  |  |  |
| HLA-A | chr6:28510020-33480577 |  |  |  |  |  |  |  | ✓ |  |  |  |
| HLA-B | chr6:28510020-33480577 |  |  |  |  |  |  |  | ✓ |  |  |  |
| HLA-C | chr6:28510020-33480577 |  |  |  |  |  |  |  | ✓ |  |  |  |
| CACNA1S | chr1:201006956-201095000 | ✓ |  |  | ✓ |  |  | ★ |  |  |  |  |
| MT-RNR1 | chrM:648-1601 |  |  |  |  |  |  |  |  |  | ✓ |  |
| CYP2C-CLUSTER | chr10:94644745-94646745 |  |  |  |  |  |  | ★ |  |  |  |  |
| GSTP1 | chr11:67580811-67589653 | ✓ |  | ✓ |  |  |  |  |  |  |  |  |
| SLC15A2 | chr3:121891400-121947188 | ✓ |  |  |  |  |  |  |  |  |  |  |
| SLC22A2 | chr6:160206754-160268821 | ✓ |  |  |  |  |  |  |  |  |  | yes |
| SLCO2B1 | chr11:75148106-75209549 | ✓ |  |  |  |  |  |  |  |  |  |  |
| UGT2B15 | chr4:68640596-68676652 | ✓ |  |  |  |  |  |  |  |  |  | yes |
| HLA-DQA1 | chr6:28510020-33480577 |  |  |  |  |  |  |  |  | ✓ |  |  |
| HLA-DQB1 | chr6:28510020-33480577 |  |  |  |  |  |  |  |  | ✓ |  |  |
| HLA-DRB1 | chr6:28510020-33480577 |  |  |  |  |  |  |  |  | ✓ |  |  |
| ABO | chr9:133255000-133258000 |  |  |  |  |  |  | ★ |  |  |  |  |

42 genes and loci.
