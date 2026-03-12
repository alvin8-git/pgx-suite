# PGx Pipeline Documentation

This pipeline performs star-allele calling and diplotype assignment across 31 pharmacogenomically relevant genes using four complementary tools: PyPGx 0.26.0, Stargazer 2.0.3, Aldy 4.8.3, and StellarPGx 1.2.7, with OptiType for HLA typing and mutserve for mitochondrial variants. All analyses are performed against the GRCh38 reference genome. Results are aggregated into a unified HTML report with concordance scoring across tools.

---

## Table of Contents

1. [ABCG2](#abcg2)
2. [CACNA1S](#cacna1s)
3. [CYP1A1](#cyp1a1)
4. [CYP1A2](#cyp1a2)
5. [CYP2A6](#cyp2a6)
6. [CYP2B6](#cyp2b6)
7. [CYP2C8](#cyp2c8)
8. [CYP2C9](#cyp2c9)
9. [CYP2C19](#cyp2c19)
10. [CYP2D6](#cyp2d6)
11. [CYP2E1](#cyp2e1)
12. [CYP3A4](#cyp3a4)
13. [CYP3A5](#cyp3a5)
14. [CYP4F2](#cyp4f2)
15. [DPYD](#dpyd)
16. [G6PD](#g6pd)
17. [GSTM1](#gstm1)
18. [GSTT1](#gstt1)
19. [HLA-A](#hla-a)
20. [HLA-B](#hla-b)
21. [IFNL3](#ifnl3)
22. [MT-RNR1](#mt-rnr1)
23. [NAT1](#nat1)
24. [NAT2](#nat2)
25. [NUDT15](#nudt15)
26. [POR](#por)
27. [RYR1](#ryr1)
28. [SLCO1B1](#slco1b1)
29. [TPMT](#tpmt)
30. [UGT1A1](#ugt1a1)
31. [VKORC1](#vkorc1)

---

## ABCG2

### Overview

ABCG2 (ATP-binding cassette subfamily G member 2), also known as BCRP (Breast Cancer Resistance Protein), is a half-transporter expressed on the apical surface of intestinal epithelium, hepatocytes, renal proximal tubules, and the blood-brain barrier. It mediates efflux of a broad range of xenobiotics and endogenous substrates, including uric acid, into the intestinal lumen, bile, and urine. ABCG2 plays a central role in drug absorption, distribution, and elimination, and also contributes to cellular defence against toxic compounds.

### Pharmacogenomic Significance

The common missense variant rs2231142 (p.Q141K, c.421C>A) reduces ABCG2 protein surface expression and transport activity by approximately 50%, leading to increased systemic exposure of substrate drugs. For rosuvastatin specifically, the Q141K variant is associated with significantly higher plasma AUC, necessitating dose reduction to mitigate statin-induced myopathy risk. ABCG2 also influences exposure to sulfasalazine, topotecan, and several tyrosine kinase inhibitors.

### Key Variants

| Allele | Variant | Effect | Clinical Relevance |
|--------|---------|--------|--------------------|
| *2 | rs2231142 (Q141K) | ~50% reduced efflux activity | Increased rosuvastatin AUC; dose reduction advised |
| *5 | rs72552713 (Q126X) | Loss of function | Rare; full transporter loss |
| *1 (reference) | — | Normal function | Standard dosing |

### Population Frequencies

The Q141K variant (rs2231142) is strikingly common in East Asian populations (~30–35% allele frequency) compared to Europeans (~10%) and Africans (~1–2%). This difference has direct implications for statin dosing guidelines, as East Asian patients carry a substantially higher burden of reduced-function ABCG2 alleles. The Q126X variant is rare across all populations (<1%).

### CPIC Level

**CPIC Level B** — Rosuvastatin (statin-induced myopathy). ABCG2 is co-genotyped alongside SLCO1B1 in statin pharmacogenomics panels. Some CPIC warrants lower evidence because the rosuvastatin guideline is classified Level B (actionable but lower confidence than Level A).

### Tool Support in pgx-suite

| Tool | Supported |
|------|-----------|
| PyPGx | No |
| Stargazer | No |
| Aldy | Yes |
| StellarPGx | Yes |
| OptiType | No |
| mutserve | No |

---

## CACNA1S

### Overview

CACNA1S encodes the alpha-1S subunit of the L-type voltage-gated calcium channel (CaV1.1), which is the principal calcium channel of skeletal muscle transverse tubules. CaV1.1 is the voltage sensor that couples membrane depolarisation to ryanodine receptor (RYR1)-mediated calcium release from the sarcoplasmic reticulum, a process termed excitation-contraction (EC) coupling. Mutations in CACNA1S are a less common but established cause of malignant hyperthermia susceptibility (MHS) and periodic paralysis.

### Pharmacogenomic Significance

Pathogenic variants in CACNA1S confer susceptibility to malignant hyperthermia (MH), a potentially fatal pharmacogenetic syndrome triggered by volatile halogenated anaesthetics (e.g., halothane, sevoflurane, desflurane) and the depolarising neuromuscular blocking agent succinylcholine. During an MH episode, uncontrolled calcium release from the sarcoplasmic reticulum leads to muscle rigidity, hyperthermia, metabolic acidosis, and rhabdomyolysis. Dantrolene, which blocks RYR1-mediated calcium release, is the specific antidote.

### Key Variants

| Variant | cDNA | Protein | MH Susceptibility |
|---------|------|---------|-------------------|
| rs772226819 | c.520C>T | p.Arg174Trp | MHS (pathogenic) |
| rs28933981 | c.3257G>A | p.Arg1086His | MHS (pathogenic) |
| c.3668G>A | — | p.Arg1223His | MHS (pathogenic) |

CACNA1S variants account for roughly 1% of MH-susceptible families; the majority (~70%) carry RYR1 variants.

### Population Frequencies

CACNA1S pathogenic variants are rare in all populations, with a combined estimated prevalence of MH susceptibility of 1 in 2,000 to 1 in 10,000 individuals. No significant ethnic skew has been reported for CACNA1S-specific MH variants, unlike some RYR1 variants.

### CPIC Level

**CPIC Level A** — Volatile anaesthetic agents and succinylcholine. The CPIC guideline (jointly with CPNDS) recommends avoidance of triggering agents in genetically susceptible individuals and use of dantrolene prophylaxis protocols where indicated.

### Tool Support in pgx-suite

| Tool | Supported |
|------|-----------|
| PyPGx | **Yes** |
| Stargazer | No |
| Aldy | No |
| StellarPGx | **Yes** |
| OptiType | No |
| mutserve | No |

> Implemented in Phase 7 (config-only: PyPGx and StellarPGx both support CACNA1S natively with standard WGS BAM input).

---

## CYP1A1

### Overview

CYP1A1 (Cytochrome P450 Family 1 Subfamily A Member 1) is an extrahepatic phase I monooxygenase expressed predominantly in the lung, placenta, and lymphocytes. It is responsible for the metabolic activation of polycyclic aromatic hydrocarbons (PAHs), heterocyclic amines, and other procarcinogens found in cigarette smoke and environmental pollutants into reactive epoxide intermediates capable of forming DNA adducts. CYP1A1 expression is strongly inducible by aryl hydrocarbon receptor (AhR) ligands.

### Pharmacogenomic Significance

CYP1A1 is not a primary target for conventional drug dosing guidelines, but its polymorphisms modify individual susceptibility to carcinogen-induced cancers (lung, colorectal, breast) and influence the toxicity of certain chemotherapeutic agents that are CYP1A1 substrates or whose elimination depends on AhR pathway activity. The *2A (m1) and *2C (m2/Ile462Val) variants are the most studied in cancer risk epidemiology. CYP1A1 also contributes to the metabolism of erlotinib and other tyrosine kinase inhibitors in extrahepatic tissues.

### Key Variants

| Allele | Variant | Protein Change | Functional Effect |
|--------|---------|---------------|-------------------|
| *2A (m1) | rs4646903 (3'-UTR) | None | Increased inducibility |
| *2C (m2) | rs1048943 (Ile462Val) | Ile462Val | ~2-fold increased activity toward PAHs |
| *4 | rs1799814 (Thr461Asn) | Thr461Asn | Altered substrate specificity |

### Population Frequencies

The *2C (Ile462Val) variant is more common in East Asians (~30–40%) compared to Europeans (~5–10%) and Africans (~3%). The *2A variant shows similarly elevated frequency in East Asian populations. These differences partially explain population variation in carcinogen-related cancer susceptibility.

### CPIC Level

**Not assigned** — CYP1A1 does not have a standalone CPIC clinical pharmacogenomics guideline. Its role is primarily in cancer pharmacogenomics and environmental toxicology research contexts.

### Tool Support in pgx-suite

| Tool | Supported |
|------|-----------|
| PyPGx | No |
| Stargazer | No |
| Aldy | No |
| StellarPGx | No |
| OptiType | No |
| mutserve | No |

---

## CYP1A2

### Overview

CYP1A2 is the predominant hepatic CYP1A family member, accounting for approximately 13% of total hepatic CYP content. It metabolises caffeine (used as a CYP1A2 phenotyping probe), theophylline, clozapine, olanzapine, fluvoxamine, and a number of procarcinogens. CYP1A2 expression is highly variable between individuals, driven by both genetic polymorphisms and environmental inducers such as cigarette smoking and cruciferous vegetables, leading to up to 60-fold inter-individual variation in activity.

### Pharmacogenomic Significance

CYP1A2 is the primary metabolising enzyme for clozapine and olanzapine; poor metabolisers or those with low inducibility may experience elevated plasma concentrations and adverse effects (seizures, sedation, metabolic syndrome). Fluvoxamine is a potent CYP1A2 inhibitor and markedly elevates clozapine exposure when co-prescribed. CYP1A2 ultrarapid metabolisers (smokers and *1F carriers) require higher clozapine doses to achieve therapeutic levels.

### Key Variants

| Allele | Variant | Effect |
|--------|---------|--------|
| *1F | rs762551 (-163C>A, intron 1) | Increased inducibility by smoking; ultrarapid in smokers |
| *1C | rs2069514 (-3860G>A) | Decreased inducibility |
| *1K | rs12720461 (intron 1) | Decreased activity |
| *6 | rs56276455 (R431W) | Reduced/no function |
| *7 | rs72547513 (M348V, splice) | No function |

### Population Frequencies

CYP1A2 *1F (rs762551 -163C>A) is the most common functional variant, with minor allele frequency ~35–40% in Europeans and East Asians. Loss-of-function alleles (*6, *7) are rare globally (<2%). The *1C allele is more frequent in Africans and Asians than Europeans.

### CPIC Level

**CPIC Level B** — Clozapine, fluvoxamine, and olanzapine. The evidence supports genotype-guided dose adjustment but with moderate confidence. CPIC Level C for theophylline.

### Tool Support in pgx-suite

| Tool | Supported |
|------|-----------|
| PyPGx | Yes |
| Stargazer | Yes |
| Aldy | Yes |
| StellarPGx | Yes |
| OptiType | No |
| mutserve | No |

---

## CYP2A6

### Overview

CYP2A6 is the primary enzyme responsible for the hepatic metabolism of nicotine to cotinine and subsequently to trans-3'-hydroxycotinine (3HC), the primary urinary metabolite of nicotine. It also metabolises coumarin, tegafur (prodrug of 5-FU), and nitrosamines derived from tobacco. CYP2A6 is expressed mainly in the liver, with minor expression in lung and nasal mucosa, and its activity is the main determinant of nicotine clearance rate.

### Pharmacogenomic Significance

The 3HC/cotinine urinary ratio (the nicotine metabolite ratio, NMR) is a validated biomarker of CYP2A6 metabolic phenotype and predicts smoking cessation response to pharmacotherapy. Normal metabolisers clear nicotine rapidly and show better response to varenicline; slow metabolisers sustain higher nicotine levels from fewer cigarettes and may respond better to nicotine replacement therapy. CYP2A6 also activates tegafur, a fluoropyrimidine prodrug, so CYP2A6 poor metabolisers may have reduced efficacy with tegafur-based regimens.

### Key Variants

| Allele | Type | Effect | Frequency (Europeans) |
|--------|------|--------|-----------------------|
| *2 | SNP (L160H) | No function | ~2% |
| *4 | Whole-gene deletion | No function | ~1–2% (higher in East Asians ~15%) |
| *9 | Promoter TATA-box | Reduced expression | ~5–10% |
| *12 | Gene conversion | Reduced function | ~2% |
| *17 | SNP (I471T) | Reduced function | ~8% in Africans |

### Population Frequencies

CYP2A6 loss-of-function alleles are substantially more frequent in East Asians (combined ~35–40% reduced/no-function allele frequency) compared to Europeans (~10–15%) and Africans (~10–15% for different alleles). This explains why East Asian smokers tend to smoke fewer cigarettes per day on average.

### CPIC Level

**CPIC Level B** — Smoking cessation pharmacotherapy (nicotine replacement therapy, varenicline). Moderate evidence for NMR-guided pharmacotherapy selection.

### Tool Support in pgx-suite

| Tool | Supported |
|------|-----------|
| PyPGx | Yes |
| Stargazer | Yes |
| Aldy | Yes |
| StellarPGx | Yes |
| OptiType | No |
| mutserve | No |

---

## CYP2B6

### Overview

CYP2B6 accounts for approximately 2–10% of total hepatic CYP content and is notable for its high inter-individual variability (>100-fold). It is the primary enzyme responsible for the metabolism of several clinically important drugs, including efavirenz, nevirapine, methadone, bupropion, and cyclophosphamide. CYP2B6 is also expressed in the brain, intestine, and lung, and is highly inducible by rifampicin and other nuclear receptor (PXR/CAR) ligands.

### Pharmacogenomic Significance

CYP2B6 genotype is the strongest genetic predictor of efavirenz plasma concentrations, with slow metabolisers (carrying *6/*6 or *6/*18 diplotypes) reaching drug exposures 3–5 times higher than normal metabolisers. This leads to significantly elevated rates of central nervous system toxicity (nightmares, dizziness, depression) and, paradoxically, can improve virological suppression at the cost of tolerability. CPIC guidelines recommend dose reduction in poor metabolisers receiving efavirenz. CYP2B6 also influences methadone clearance, with slow metabolisers at increased risk of QTc prolongation.

### Key Variants

| Allele | Key Variants | Phenotype | Freq (Europeans) |
|--------|-------------|-----------|-----------------|
| *6 | Q172H + K262R (rs3745274 + rs2279343) | Reduced function | ~25–30% |
| *18 | I328T (rs28399499) | No function | ~5% (higher in Africans ~10%) |
| *22 | Promoter deletion | Reduced expression | ~5% in Africans |
| *4 | K262R only | Slightly reduced | ~15% |

### Population Frequencies

The *6 allele is common across all populations (Europeans ~28%, Africans ~33%, East Asians ~16%). The *18 (no-function) allele is predominantly observed in sub-Saharan Africans (~10%), making African patients disproportionately represented among CYP2B6 poor metabolisers treated with efavirenz-containing HIV regimens.

### CPIC Level

**CPIC Level A** — Efavirenz, nevirapine, and methadone. CPIC provides specific dose-adjustment recommendations by metaboliser phenotype for efavirenz-based antiretroviral therapy.

### Tool Support in pgx-suite

| Tool | Supported |
|------|-----------|
| PyPGx | Yes |
| Stargazer | Yes |
| Aldy | Yes |
| StellarPGx | Yes |
| OptiType | No |
| mutserve | No |

---

## CYP2C8

### Overview

CYP2C8 is a hepatic cytochrome P450 responsible for the oxidative metabolism of approximately 5% of clinically used drugs, particularly those with a high hepatic first-pass effect. Major substrates include paclitaxel (primary metabolic pathway), rosiglitazone, pioglitazone, repaglinide, amodiaquine, and chloroquine. CYP2C8 shares overlapping substrate specificity with CYP2C9 and is co-regulated by similar nuclear receptors (PXR, CAR).

### Pharmacogenomic Significance

The CYP2C8*3 allele, which carries two co-inherited amino acid substitutions (R139K and K399R), is associated with reduced enzyme activity toward paclitaxel and rosiglitazone, potentially affecting both efficacy and toxicity profiles. Paclitaxel neuropathy has been linked to CYP2C8 poor metaboliser phenotype, as slower clearance leads to prolonged drug exposure. Gemfibrozil glucuronide is a potent mechanism-based inhibitor of CYP2C8 and can cause dangerous drug interactions with repaglinide.

### Key Variants

| Allele | Variants | Protein Effect | Function |
|--------|---------|---------------|----------|
| *2 | I269F (rs11572103) | Structural alteration | Reduced |
| *3 | R139K + K399R (rs11572080 + rs10509681) | Altered active site | Reduced (~50%) |
| *4 | I264M (rs1058930) | Altered substrate access | Reduced |

### Population Frequencies

CYP2C8*3 is most prevalent in Europeans (~13% allele frequency) and much rarer in Africans (~1–2%) and East Asians (<1%). CYP2C8*2 is predominantly found in African populations (~18%). This population stratification affects which reduced-function allele dominates in different ethnic groups.

### CPIC Level

**CPIC Level B** — Paclitaxel (neuropathy risk) and rosiglitazone. Evidence is actionable but the guideline confidence is classified below Level A.

### Tool Support in pgx-suite

| Tool | Supported |
|------|-----------|
| PyPGx | Yes |
| Stargazer | Yes |
| Aldy | Yes |
| StellarPGx | Yes |
| OptiType | No |
| mutserve | No |

---

## CYP2C9

### Overview

CYP2C9 is a major hepatic cytochrome P450, accounting for ~20% of hepatic CYP content and metabolising approximately 15% of all clinically prescribed drugs. It is responsible for the hydroxylation of S-warfarin (the more potent enantiomer), phenytoin, celecoxib, ibuprofen, diclofenac, losartan, glipizide, and tolbutamide. CYP2C9 activity is tightly regulated by genetic polymorphism and can be inhibited by fluconazole and amiodarone, or induced by rifampicin.

### Pharmacogenomic Significance

CYP2C9 reduced-function alleles markedly decrease warfarin clearance, leading to elevated plasma concentrations and increased bleeding risk. CPIC and IWPC dosing algorithms incorporate CYP2C9 genotype alongside VKORC1 and CYP4F2 to predict therapeutic warfarin maintenance dose. For phenytoin, CYP2C9 poor metabolisers require substantially lower doses to avoid neurotoxicity. NSAIDs metabolised by CYP2C9 (celecoxib, ibuprofen) have prolonged exposure in poor metabolisers, raising gastrointestinal and cardiovascular risk.

### Key Variants

| Allele | Variant | Protein Change | Residual Activity | Freq (Europeans) |
|--------|---------|---------------|-------------------|-----------------|
| *2 | rs1799853 | R144C | ~30% of normal | ~12% |
| *3 | rs1057910 | I359L | ~5% of normal | ~6% |
| *5 | rs28371686 | D360E | Reduced | <1% |
| *6 | rs9332131 | 818delA | No function | <1% |
| *8 | rs7900194 | R150H | Reduced | ~1% (Africans ~4%) |
| *11 | rs28371685 | R335W | Reduced | <1% (Africans ~2%) |

### Population Frequencies

CYP2C9*2 and *3 (reduced-function) are common in Europeans but rare in East Asians and Africans. African populations carry different reduced-function alleles (*5, *6, *8, *11) that are often not included in commercial pharmacogenomic panels, creating health equity gaps in warfarin dosing accuracy.

### CPIC Level

**CPIC Level A** — Warfarin, acenocoumarol, phenprocoumon, phenytoin, and fosphenytoin. CYP2C9 is a cornerstone gene in clinical pharmacogenomics implementation.

### Tool Support in pgx-suite

| Tool | Supported |
|------|-----------|
| PyPGx | Yes |
| Stargazer | Yes |
| Aldy | Yes |
| StellarPGx | Yes |
| OptiType | No |
| mutserve | No |

---

## CYP2C19

### Overview

CYP2C19 is a polymorphic hepatic cytochrome P450 involved in the metabolism of approximately 10–15% of prescribed drugs, with particularly high clinical significance for antiplatelet therapy, proton pump inhibitors (PPIs), and psychiatric medications. Key substrates include clopidogrel (prodrug activation), omeprazole, pantoprazole, escitalopram, citalopram, amitriptyline, voriconazole, and diazepam. The enzyme exhibits a unique combination of both no-function (poor metaboliser) and gain-of-function (ultrarapid metaboliser) alleles of major clinical consequence.

### Pharmacogenomic Significance

CYP2C19 is essential for converting clopidogrel to its active thiol metabolite that irreversibly inhibits platelet P2Y12 receptors. Poor metabolisers (CYP2C19*2/*2 or *2/*3) fail to adequately activate clopidogrel, resulting in residual platelet reactivity and increased risk of stent thrombosis and major adverse cardiovascular events (MACE). CPIC recommends alternative antiplatelet therapy (prasugrel or ticagrelor) for poor and intermediate metabolisers following acute coronary syndrome. Ultrarapid metabolisers (*17/*17) show enhanced clopidogrel response but normal PPI and SSRI metabolism.

### Key Variants

| Allele | Variant | Effect | Freq (Europeans) |
|--------|---------|--------|-----------------|
| *2 | rs4244285 (splice defect) | No function | ~15% |
| *3 | rs4986893 (W212X) | No function | <1% (East Asians ~5%) |
| *17 | rs12248560 (promoter -806C>T) | Increased expression; ultrarapid | ~21% |
| *4 | rs28399504 (I331V) | No function | <1% |

### Population Frequencies

CYP2C19*2 poor-metaboliser allele is common across all populations (Europeans ~15%, East Asians ~29–35%, Africans ~17%). CYP2C19*3 is almost exclusively found in East and Southeast Asians (5–9%), where the combined poor-metaboliser frequency is ~15–20%. CYP2C19*17 is most prevalent in Africans (~16–28%) and Europeans (~21%), making ultrarapid metabolism particularly common in these groups.

### CPIC Level

**CPIC Level A** — Clopidogrel, PPIs (omeprazole, pantoprazole), tricyclic antidepressants, SSRIs (escitalopram, citalopram, sertraline), and voriconazole.

### Tool Support in pgx-suite

| Tool | Supported |
|------|-----------|
| PyPGx | Yes |
| Stargazer | Yes |
| Aldy | Yes |
| StellarPGx | Yes |
| OptiType | No |
| mutserve | No |

---

## CYP2D6

### Overview

CYP2D6 is one of the most intensively studied pharmacogenes, metabolising an estimated 20–25% of all clinically prescribed drugs despite representing only ~2–5% of total hepatic CYP content. Key substrate classes include opioid analgesics (codeine, tramadol, oxycodone), antidepressants (tricyclics, SSRIs, SNRIs), antipsychotics (haloperidol, risperidone), beta-blockers (metoprolol, propranolol), and the breast cancer drug tamoxifen. CYP2D6 is unusual among drug-metabolising enzymes in being non-inducible and in displaying extreme allelic diversity including gene deletion, duplication, and hybrid alleles.

### Pharmacogenomic Significance

CYP2D6 phenotype determines both the efficacy and toxicity of multiple drug classes. Codeine is a prodrug requiring CYP2D6-mediated conversion to morphine; poor metabolisers receive no analgesic benefit while ultrarapid metabolisers may generate potentially fatal morphine levels (CPIC recommends codeine avoidance in both extremes). Tamoxifen requires CYP2D6-mediated activation to endoxifen, its active metabolite; poor and intermediate metabolisers may have compromised breast cancer outcomes. Atomoxetine exposure is 10-fold higher in poor metabolisers, requiring dose adjustment.

### Key Variants

| Allele | Type | Function | Activity Score | Freq (Europeans) |
|--------|------|----------|---------------|-----------------|
| *3 | 2549A>del frameshift | No function | 0 | ~1–2% |
| *4 | Splice site defect | No function | 0 | ~20% |
| *5 | Gene deletion | No function | 0 | ~4% |
| *6 | T1707del frameshift | No function | 0 | ~1% |
| *10 | P34S (splice effect) | Reduced | 0.25 | ~2% (East Asians ~35–50%) |
| *17 | T107I + R296C | Reduced | 0.5 | ~20% in Africans |
| *41 | 2988G>A (splice) | Reduced | 0.5 | ~9% |
| *2xN | Gene duplication | Ultrarapid | 2.0 | ~2–5% (higher in East Africa) |

### Population Frequencies

CYP2D6*4 (no-function) is the most common loss-of-function allele in Europeans (~20%), making the European poor-metaboliser frequency ~7–10%. East Asians have lower poor-metaboliser rates (~1–2%) but very high *10 (reduced-function) frequencies. Africans have the highest prevalence of gene duplication (ultrarapid metabolism) and diverse reduced-function alleles (*17, *29), making phenotype prediction challenging without comprehensive genotyping.

### CPIC Level

**CPIC Level A** — Codeine, tramadol, oxycodone, tamoxifen, atomoxetine, tricyclic antidepressants, ondansetron, tropisetron, and multiple antidepressants/antipsychotics.

### Tool Support in pgx-suite

| Tool | Supported |
|------|-----------|
| PyPGx | Yes |
| Stargazer | Yes |
| Aldy | Yes |
| StellarPGx | Yes |
| OptiType | No |
| mutserve | No |

---

## CYP2E1

### Overview

CYP2E1 is a constitutively expressed cytochrome P450 located in the endoplasmic reticulum of hepatocytes that metabolises ethanol, acetaminophen (paracetamol), and a range of low-molecular-weight compounds including halogenated anaesthetics, benzene, and dimethylnitrosamine. CYP2E1 is unusual in generating a high rate of reactive oxygen species (ROS) during its catalytic cycle, contributing significantly to oxidative stress-mediated hepatotoxicity. Its expression is inducible by ethanol, obesity, and diabetes (via insulin signalling pathways).

### Pharmacogenomic Significance

CYP2E1 converts acetaminophen to the hepatotoxic metabolite N-acetyl-p-benzoquinone imine (NAPQI), which is normally detoxified by glutathione but causes hepatocellular necrosis at toxic doses, particularly in individuals with induced CYP2E1 (chronic alcohol users) or depleted glutathione reserves. CYP2E1 polymorphisms (notably the *5B/*5B and *6/*6 genotypes) may modulate susceptibility to acetaminophen-induced liver injury, anaesthetic-related hepatitis, and chemotherapy-induced organ toxicity, though the clinical evidence remains largely in the research domain.

### Key Variants

| Allele | Variant | Effect |
|--------|---------|--------|
| *5B | rs3813867 (5'-flanking) | Increased transcription in vitro |
| *6 | rs2031920 (PstI RFLP) | Altered transcription regulation |
| *7A | rs6413419 (exon 7 G7632A) | Altered function |

### Population Frequencies

The CYP2E1 *5B and *6 variant alleles are more common in East Asians (~25–35% frequency) than in Europeans (~5%) or Africans (~3%). This has been studied in the context of increased susceptibility to alcohol-related liver disease in some East Asian cohorts, though confounders are significant.

### CPIC Level

**Not assigned** — CYP2E1 does not have a dedicated CPIC clinical guideline. Its role is primarily in toxicology and research pharmacogenomics contexts.

### Tool Support in pgx-suite

| Tool | Supported |
|------|-----------|
| PyPGx | No |
| Stargazer | No |
| Aldy | No |
| StellarPGx | No |
| OptiType | No |
| mutserve | No |

---

## CYP3A4

### Overview

CYP3A4 is the most abundantly expressed cytochrome P450 in the human liver and small intestine, accounting for ~30–40% of total hepatic CYP protein and metabolising an estimated 40–50% of all clinically used drugs. Its unusually large and flexible active site accommodates substrates ranging from small molecules to large macrolide antibiotics. Key substrates include statins, immunosuppressants (tacrolimus, cyclosporine), benzodiazepines, calcium channel blockers, antiretrovirals, and many oncology agents. CYP3A4 activity is highly variable (up to 40-fold), driven more by environmental factors and drug interactions than by genetic variation.

### Pharmacogenomic Significance

Unlike most drug-metabolising CYPs, genetic variation in CYP3A4 explains only a fraction of its inter-individual variability; the majority is attributable to drug-drug interactions (DDIs), particularly with potent inhibitors (e.g., ketoconazole, ritonavir, grapefruit juice) and inducers (e.g., rifampicin, phenytoin, St John's Wort). The CYP3A4*22 allele (intronic SNP rs35599367) reduces hepatic mRNA expression by ~50% and is associated with modestly reduced metabolic capacity for tacrolimus, atorvastatin, and other substrates. CYP3A4 is therefore clinically important primarily in the context of DDI management and as a secondary modifier for CYP3A5 substrates.

### Key Variants

| Allele | Variant | Effect |
|--------|---------|--------|
| *22 | rs35599367 (intron 6) | ~50% reduced mRNA; reduced function |
| *1B | rs2740574 (5'-UTR) | Mildly increased expression |
| *20 | 14690delA (frameshift) | No function (rare) |

### Population Frequencies

CYP3A4*22 allele frequency is ~5–7% in Europeans and rare (<2%) in East Asians and Africans. The *1B variant is more common in individuals of African ancestry (~50–60% frequency) vs. Europeans (~5–10%), though its functional significance is debated.

### CPIC Level

**Not assigned as standalone** — CYP3A4 is included as a modifier in tacrolimus (via CYP3A5), simvastatin, and other guidelines but does not have a dedicated Level A CPIC guideline. CYP3A5 is the primary actionable gene in the CYP3A locus.

### Tool Support in pgx-suite

| Tool | Supported |
|------|-----------|
| PyPGx | Yes |
| Stargazer | No |
| Aldy | Yes |
| StellarPGx | Yes |
| OptiType | No |
| mutserve | No |

---

## CYP3A5

### Overview

CYP3A5 shares ~84% amino acid identity with CYP3A4 and overlapping substrate specificity, but unlike CYP3A4, its expression is polymorphically controlled: most individuals of European or East Asian descent express little or no CYP3A5 protein (*3/*3 genotype), while it is abundantly expressed in the majority of individuals of African descent (*1 carrier). CYP3A5 is a major contributor to extra-hepatic CYP3A activity in the kidney and small intestine, and is the primary determinant of tacrolimus clearance variability.

### Pharmacogenomic Significance

CYP3A5 expressors (*1 carriers) have substantially higher tacrolimus clearance than non-expressors (*3/*3), requiring approximately 1.5–2 times higher dose to achieve the same trough concentration. Failure to account for CYP3A5 genotype in tacrolimus dosing leads to prolonged time outside therapeutic range post-transplant, increasing risk of rejection or nephrotoxicity. CPIC recommends CYP3A5-guided tacrolimus starting doses, stratifying patients into expressors and non-expressors.

### Key Variants

| Allele | Variant | Effect | Freq (Europeans) | Freq (Africans) |
|--------|---------|--------|-----------------|-----------------|
| *3 | rs776746 (splice, intron 3) | No function (premature stop) | ~94% | ~27% |
| *6 | rs10264272 (exon 7 splice) | No function | ~1% | ~17% |
| *7 | rs41303343 (exon 11 frameshift) | No function | <1% | ~8% |

### Population Frequencies

CYP3A5 expressors (*1/*1 or *1/*3) are the minority in European (~10%) and East Asian (~25%) populations, but the majority in African populations (~70–75%). This difference has profound implications for tacrolimus dosing standards that were historically derived from predominantly European populations.

### CPIC Level

**CPIC Level A** — Tacrolimus. CPIC guideline recommends genotype-guided initial dose selection, with expressors requiring higher doses and more frequent monitoring.

### Tool Support in pgx-suite

| Tool | Supported |
|------|-----------|
| PyPGx | Yes |
| Stargazer | Yes |
| Aldy | Yes |
| StellarPGx | Yes |
| OptiType | No |
| mutserve | No |

---

## CYP4F2

### Overview

CYP4F2 is a member of the cytochrome P450 4F subfamily involved primarily in the omega-hydroxylation of arachidonic acid metabolites (leukotrienes), long-chain fatty acids, and vitamin K1 (phylloquinone) and vitamin K2 (menaquinone-4). Its role in the hepatic catabolism of vitamin K is pharmacogenomically critical because it modulates the size of the intrahepatic vitamin K pool available for vitamin K-dependent clotting factor gamma-carboxylation. CYP4F2 is expressed primarily in the liver and kidney.

### Pharmacogenomic Significance

The CYP4F2*3 variant (V433M, rs2108622) reduces vitamin K oxidation by approximately 50%, leading to higher hepatic vitamin K levels, more complete carboxylation of clotting factors, and a consequently lower sensitivity to vitamin K antagonist (VKA) anticoagulants such as warfarin. Patients carrying *3 require a modestly higher warfarin maintenance dose (~10–15% increase per *3 allele) to achieve target anticoagulation. CYP4F2 is therefore incorporated as a secondary modifier in multi-gene warfarin dosing algorithms alongside CYP2C9 and VKORC1.

### Key Variants

| Allele | Variant | Protein | Effect on Warfarin Dose |
|--------|---------|---------|------------------------|
| *3 | rs2108622 (V433M) | Val433Met | Higher dose required (~+10–15%/allele) |
| *1 | Reference | Wild type | Standard dosing |

### Population Frequencies

CYP4F2*3 allele frequency is approximately 26–28% in Europeans, 18–25% in East Asians, and 6–7% in Africans. The higher frequency in European populations is one reason it was initially identified in European warfarin cohort studies.

### CPIC Level

**Incorporated in CPIC Level A warfarin guideline** — CYP4F2 is included as a tertiary modifier in the CPIC/IWPC warfarin dosing algorithm. It is not classified independently but modifies dose predictions when included alongside CYP2C9 and VKORC1.

### Tool Support in pgx-suite

| Tool | Supported |
|------|-----------|
| PyPGx | Yes |
| Stargazer | Yes |
| Aldy | Yes |
| StellarPGx | Yes |
| OptiType | No |
| mutserve | No |

---

## DPYD

### Overview

DPYD encodes dihydropyrimidine dehydrogenase (DPD), the rate-limiting enzyme in the catabolism of uracil, thymine, and the fluoropyrimidine drugs 5-fluorouracil (5-FU) and capecitabine. DPD is expressed predominantly in hepatocytes and accounts for >80% of 5-FU elimination, converting it to the inactive metabolite dihydrofluorouracil (DHFU). Partial or complete DPD deficiency dramatically prolongs 5-FU half-life, resulting in prolonged exposure to the active drug and potentially fatal toxicity.

### Pharmacogenomic Significance

DPYD reduced-function variants cause severe fluoropyrimidine toxicity including mucositis, neutropaenia, hand-foot syndrome, diarrhoea, and treatment-related death. Four DPYD variants are designated as High Clinical Significance by CPIC and the European Medicines Agency (EMA): *2A, *13, rs67376798 (D949V), and HapB3 (rs56038477, c.1129-5923C>G splice). CPIC recommends dose reduction of 25–50% for heterozygous carriers and alternative therapy for homozygous null individuals. DPYD pre-treatment genotyping is now mandated by the EMA.

### Key Variants

| Allele / Variant | cDNA Change | Effect | Clinical Significance |
|-----------------|------------|--------|----------------------|
| *2A | IVS14+1G>A (rs3918290) | Exon 14 skipping; no function | CPIC High; 50% dose reduction |
| *13 | c.1679T>G (rs55886062, I560S) | No function | CPIC High; 50% dose reduction |
| rs67376798 | c.2846A>T (D949V) | Reduced function | CPIC High; 25–50% dose reduction |
| HapB3 (c.1129-5923C>G) | rs56038477 | Reduced function (splice) | CPIC High; 25% dose reduction |

### Population Frequencies

DPYD*2A has an allele frequency of ~1% in Europeans and is extremely rare in Africans and East Asians. HapB3 is predominantly observed in Europeans (~3–5%). Overall, ~3–5% of Europeans carry at least one DPYD reduced-function allele. African populations carry mostly different, less characterised variants, representing an area of active research and health equity concern.

### CPIC Level

**CPIC Level A** — Fluoropyrimidines: 5-fluorouracil (5-FU), capecitabine, and tegafur. Pre-treatment DPYD genotyping is now regulatory-mandated in the European Union.

### Tool Support in pgx-suite

| Tool | Supported |
|------|-----------|
| PyPGx | Yes |
| Stargazer | Yes |
| Aldy | Yes |
| StellarPGx | Yes |
| OptiType | No |
| mutserve | No |

---

## G6PD

### Overview

G6PD (glucose-6-phosphate dehydrogenase) is encoded on the X chromosome and is the rate-limiting enzyme of the pentose phosphate pathway, generating NADPH needed to maintain reduced glutathione in erythrocytes. Because mature red blood cells lack mitochondria and cannot generate NADPH by any other pathway, G6PD is essential for protecting erythrocytes from oxidative haemolysis. G6PD deficiency is the most common human enzyme disorder globally, affecting an estimated 400–500 million individuals worldwide.

### Pharmacogenomic Significance

G6PD-deficient individuals are at risk of acute haemolytic anaemia when exposed to oxidative drugs including rasburicase, primaquine, dapsone, nitrofurantoin, and certain sulphonamides. CPIC recommends rasburicase contraindication in G6PD-deficient patients due to the risk of severe, potentially life-threatening haemolysis. Primaquine and tafenoquine (anti-malarial radical cure agents) require G6PD testing prior to administration; deficient males should not receive these drugs. Notably, G6PD-deficient individuals may have reduced susceptibility to malaria, explaining the high allele frequency in malaria-endemic regions.

### Key Variants

| Variant | Common Name | Class | Residual Activity |
|---------|------------|-------|------------------|
| c.202G>A (Val68Met) | G6PD A- | III | ~10% |
| c.376A>G (Asn126Asp) | G6PD A+ (modifier) | III/IV | Variable |
| c.563C>T (Ser188Phe) | G6PD Mediterranean | II | <1% |
| c.1003G>A (Val335Met) | G6PD Canton | II | <5% |
| c.1376G>T (Arg459Leu) | G6PD Kaiping | II | ~10% |

WHO classifies G6PD variants into classes I–V by severity (I = most severe; V = increased activity).

### Population Frequencies

G6PD deficiency allele frequencies reach 10–35% in malaria-endemic regions of sub-Saharan Africa, the Middle East, South Asia, and Southeast Asia. G6PD A- (202A mutation combined with 376G) is the predominant variant in sub-Saharan Africans. G6PD Mediterranean is prevalent in Mediterranean populations, the Middle East, and South Asia. Because G6PD is X-linked, hemizygous males are fully deficient for a single variant allele, while females can be heterozygous with intermediate enzyme activity.

### CPIC Level

**CPIC Level A** — Rasburicase (contraindicated in deficiency), primaquine, tafenoquine, dapsone, and nitrofurantoin.

### Tool Support in pgx-suite

| Tool | Supported |
|------|-----------|
| PyPGx | Yes |
| Stargazer | No |
| Aldy | Yes |
| StellarPGx | No |
| OptiType | No |
| mutserve | No |

---

## GSTM1

### Overview

GSTM1 encodes glutathione S-transferase mu 1, a phase II conjugation enzyme that catalyses the transfer of glutathione to electrophilic substrates, enabling their subsequent elimination. GSTM1 is expressed in the liver, lung, kidney, and lymphocytes and participates in the detoxification of reactive oxygen species, carcinogens (including PAH epoxides and aflatoxin metabolites), and certain chemotherapeutic agents. The enzyme functions as a homodimer and belongs to the mu class of the glutathione transferase superfamily.

### Pharmacogenomic Significance

The most clinically significant GSTM1 variant is the *0 (null) allele, produced by a homozygous whole-gene deletion that results in complete absence of GSTM1 enzyme activity. GSTM1 null homozygotes have reduced capacity to detoxify electrophilic carcinogens, which is associated with modestly increased risk of lung, bladder, and colorectal cancers in smokers and individuals with occupational carcinogen exposure. In the context of chemotherapy, GSTM1 null genotype may influence response to platinum-based drugs and cyclophosphamide, as these drugs are partly inactivated by glutathione conjugation.

### Key Variants

| Genotype | Variant | Enzyme Activity |
|----------|---------|----------------|
| *0/*0 (null) | Homozygous whole-gene deletion | None (~0%) |
| *A/*A or *A/*0 | At least one copy present | Normal |
| *B/*B or heterozygous | Variant with one copy | Normal to reduced |

### Population Frequencies

GSTM1 null homozygosity frequency is approximately 47–65% in European populations, 35–55% in East Asians, and 26–35% in Africans. Given this high null frequency, GSTM1 null is the most common pharmacogenomically relevant homozygous loss-of-function variant in the human genome. Whole-gene copy number analysis is required for accurate GSTM1 genotyping.

### CPIC Level

**Not formally assigned (research/exploratory)** — GSTM1 null genotype has been studied extensively in pharmacogenomics but does not yet have a standalone CPIC clinical guideline. It is included in panels for cancer pharmacogenomics research.

### Tool Support in pgx-suite

| Tool | Supported |
|------|-----------|
| PyPGx | No |
| Stargazer | No |
| Aldy | No |
| StellarPGx | No |
| OptiType | No |
| mutserve | No |

> GSTM1 null genotyping requires copy number variant (CNV) analysis not currently implemented in this pipeline.

---

## GSTT1

### Overview

GSTT1 encodes glutathione S-transferase theta 1, a phase II enzyme involved in the conjugation of glutathione to small electrophilic substrates including halogenated compounds, epoxides, and reactive oxygen species. GSTT1 is expressed in liver, erythrocytes, and kidney and participates in the detoxification of environmental carcinogens (ethylene oxide, methyl bromide, dichloromethane) and is the primary glutathione conjugating enzyme for certain chemotherapy-related species. Like GSTM1, the predominant pharmacogenomically relevant variant is a null allele produced by whole-gene deletion.

### Pharmacogenomic Significance

GSTT1 null homozygotes have reduced capacity to conjugate and eliminate reactive metabolites of halogenated compounds and certain chemotherapeutic drugs (cyclophosphamide, busulfan, thiotepa). GSTT1 null genotype has been studied as a modifier of cancer risk (particularly bladder, brain, and haematological malignancies) and chemotherapy toxicity and efficacy in haematopoietic stem cell transplant settings. Clinical implementation remains primarily in research contexts.

### Key Variants

| Genotype | Enzyme Activity |
|----------|----------------|
| *0/*0 (null) | None |
| Present (at least one copy) | Normal |

> **GRCh38 Caveat**: GSTT1 is located on the chr22_KI270879v1_alt alternate contig in the GRCh38 reference assembly. Standard WGS aligners configured for primary chromosomes only may fail to align reads to GSTT1, causing apparent homozygous deletion calls (false-null genotype) in individuals who carry functional GSTT1 copies on the alternate contig. Careful aligner configuration including alt contigs, or orthogonal CNV methods, is required for accurate GSTT1 genotyping in GRCh38 pipelines.

### Population Frequencies

GSTT1 null homozygosity occurs in approximately 10–25% of Europeans, 55–65% of East Asians, and 25–35% of Africans. The high frequency in East Asians is notable. The substantial ethnic variation mirrors but is quantitatively distinct from GSTM1 null frequencies.

### CPIC Level

**Not formally assigned** — GSTT1 does not have a dedicated CPIC clinical guideline. Genotyping is performed in research and exploratory clinical pharmacogenomics panels.

### Tool Support in pgx-suite

| Tool | Supported |
|------|-----------|
| PyPGx | No |
| Stargazer | No |
| Aldy | No |
| StellarPGx | No |
| OptiType | No |
| mutserve | No |

> GRCh38 alt-contig placement is a specific technical barrier to GSTT1 implementation in this WGS pipeline.

---

## HLA-A

### Overview

HLA-A (Human Leukocyte Antigen A) encodes a major histocompatibility complex (MHC) class I cell-surface glycoprotein that presents intracellular peptide antigens (typically 8–10 amino acids) to CD8+ cytotoxic T lymphocytes, enabling immune surveillance of infected and malignant cells. The HLA-A locus is located within the MHC region on chromosome 6p21.3 and is the most polymorphic locus in the human genome, with thousands of distinct alleles defined at the 4-digit level. HLA type is co-dominantly expressed, with individuals typically expressing two distinct HLA-A alleles.

### Pharmacogenomic Significance

Specific HLA-A alleles are strongly associated with severe immune-mediated adverse drug reactions (ADRs) through a hapten or pharmacological interaction (p-i) mechanism, where the drug or its metabolite alters the MHC-peptide-TCR interaction. HLA-A*31:01 is a risk factor for carbamazepine-induced severe cutaneous adverse reactions (SCARs) including drug reaction with eosinophilia and systemic symptoms (DRESS), Stevens-Johnson syndrome (SJS), and toxic epidermal necrolysis (TEN), particularly in European and Japanese populations. HLA-A*68:01 has been linked to abacavir hypersensitivity in some cohorts, though HLA-B*57:01 remains the primary abacavir marker.

### Key Allele–Drug Associations

| HLA-A Allele | Drug | Reaction | Population |
|-------------|------|----------|------------|
| *31:01 | Carbamazepine | DRESS, SJS/TEN, MPE | European, Japanese, Han Chinese |
| *68:01 | Abacavir | Hypersensitivity (secondary association) | Reported in some cohorts |
| *02:01 | Clozapine | Agranulocytosis (weak association) | Exploratory |

### Population Frequencies

HLA-A*31:01 has allele frequency ~5% in Japanese, ~3–5% in Northern Europeans, and is generally rarer in other populations. The HLA region shows striking population stratification due to ancestral selective pressure from infectious pathogens, making population-specific allele frequency tables essential for pre-test probability calculations.

### CPIC Level

**CPIC Level A** — Carbamazepine (HLA-A*31:01; European and Japanese populations). CPIC recommends avoiding carbamazepine in *31:01 carriers when equally effective alternatives exist. HLA-B*15:02 screening is recommended in Asian populations for SJS/TEN risk.

### Tool Support in pgx-suite

| Tool | Supported |
|------|-----------|
| PyPGx | No |
| Stargazer | No |
| Aldy | No |
| StellarPGx | No |
| OptiType | Yes (4-digit typing) |
| mutserve | No |

> OptiType performs HLA-A 4-digit typing from WGS reads extracted from the MHC region (chr6:28510020–33480577).

---

## HLA-B

### Overview

HLA-B encodes the most polymorphic of the classical MHC class I molecules, with over 7,000 alleles defined at the 4-digit level. Like HLA-A, it presents intracellular peptides to CD8+ T cells and plays a central role in adaptive immune surveillance. HLA-B is located immediately centromeric to HLA-A on chromosome 6p21.3. Several HLA-B alleles have the strongest known pharmacogenomic associations of any gene in the human genome, linking specific alleles to life-threatening drug hypersensitivity reactions with positive predictive values approaching near-certainty for some allele-drug combinations.

### Pharmacogenomic Significance

HLA-B*57:01 is essentially the sole genetic determinant of abacavir (ABC) hypersensitivity syndrome (AHS), a systemic hypersensitivity reaction occurring in ~8% of *57:01-negative carriers but approaching 0% clinical events when *57:01-negative individuals receive abacavir (negative predictive value ~100% in prospective studies). CPIC and regulatory agencies mandate pre-treatment HLA-B*57:01 screening before abacavir initiation. HLA-B*58:01 confers strong risk of allopurinol-induced SJS/TEN, particularly in Han Chinese and Southeast Asian populations. HLA-B*15:02 confers risk of carbamazepine-induced SJS/TEN in Southeast Asian but not European populations.

### Key Allele–Drug Associations

| HLA-B Allele | Drug | Reaction | OR/RR | Population |
|-------------|------|----------|-------|------------|
| *57:01 | Abacavir | Hypersensitivity syndrome | ~900–1800 | All (freq ~8% Europeans) |
| *58:01 | Allopurinol | SJS/TEN, DRESS | >100 | Han Chinese, Thai, Korean |
| *15:02 | Carbamazepine | SJS/TEN | ~80 | Southeast Asian, Han Chinese |
| *15:02 | Phenytoin | SJS/TEN | ~10 | Southeast Asian |
| *57:03 | Flucloxacillin | DILI | ~80 | Europeans |

### Population Frequencies

HLA-B*57:01 allele frequency is ~5–8% in Europeans, ~6% in Africans, and ~2–3% in South Asians, but <1% in East Asians. HLA-B*58:01 frequency is approximately 6–8% in Han Chinese, ~8% in Koreans, ~8% in Thais, compared to <1% in Europeans — explaining why allopurinol-associated SJS/TEN is predominantly observed in Asian populations and motivating Asian-specific pre-treatment screening guidelines.

### CPIC Level

**CPIC Level A** — Abacavir (B*57:01), allopurinol (B*58:01), carbamazepine (B*15:02), and phenytoin/fosphenytoin (B*15:02).

### Tool Support in pgx-suite

| Tool | Supported |
|------|-----------|
| PyPGx | No |
| Stargazer | No |
| Aldy | No |
| StellarPGx | No |
| OptiType | Yes (4-digit typing) |
| mutserve | No |

> OptiType performs HLA-B 4-digit typing. HLA-B allele calling is the primary clinical application of HLA typing in this pipeline.

---

## IFNL3

### Overview

The IFNL3/IFNL4 locus on chromosome 19q13.13 encodes type III interferons (lambda interferons) involved in innate antiviral immunity, particularly in hepatocytes and intestinal epithelium. IFNL3 (IL28B) encodes interferon lambda-3. Historically, SNPs in the IFNL3 region (rs12979860 C>T and rs8099917 T>G) were identified as the strongest predictors of spontaneous and treatment-induced HCV clearance. Subsequent fine-mapping revealed that the causal gene is in fact IFNL4: the ancestral rs368234815 TT>dG variant (ss469415590) creates an open reading frame encoding a functional IFNL4 protein in non-responders, while the dG allele is in near-complete linkage disequilibrium with the rs12979860 T allele (poor-response genotype).

### Pharmacogenomic Significance

Patients carrying the favourable genotype (rs12979860 CC or IFNL4 TT/TT, non-dG) have significantly higher rates of sustained virological response (SVR) to pegylated interferon alfa-2a/2b and ribavirin therapy for chronic HCV infection, particularly genotype 1 HCV (~80% SVR in CC vs ~30–40% in TT). This distinction was crucial in the era of interferon-based therapy and remains relevant as a modifier of response and as a prognostic marker. With the advent of directly acting antivirals (DAAs) achieving >95% SVR regardless of IFNL3/IFNL4 genotype, the clinical utility of this test has declined substantially for HCV treatment but retains value for predicting spontaneous HCV clearance after acute infection.

### Key Variants

| Variant | Gene | Allele | Association |
|---------|------|--------|-------------|
| rs12979860 C>T | IFNL3 (proxy) | C = favourable | SVR predictor (interferon + ribavirin) |
| rs8099917 T>G | IFNL3 region | T = favourable | SVR predictor (interferon + ribavirin) |
| ss469415590 TT>dG | IFNL4 (causal) | TT = favourable | Creates functional IFNL4 protein if dG |

### Population Frequencies

The favourable rs12979860 C allele (non-IFNL4-expressing) is highly prevalent in East Asians (~93%), which correlates with their historically superior interferon-based HCV treatment outcomes compared to Europeans (~60–65% C allele) and Africans (~38%), where the T/dG allele predominates. This dramatic ethnic stratification makes IFNL3/IFNL4 a textbook example of pharmacogenomics-driven therapeutic disparities.

### CPIC Level

**CPIC Level A** — Peginterferon alfa-2a and alfa-2b (for HCV treatment). Clinical utility is now predominantly restricted to settings where direct-acting antivirals are unavailable.

### Tool Support in pgx-suite

| Tool | Supported |
|------|-----------|
| PyPGx | No |
| Stargazer | No |
| Aldy | No |
| StellarPGx | No |
| OptiType | No |
| mutserve | No |

> IFNL3/IFNL4 SNP genotyping (rs12979860, ss469415590) can be performed via targeted bcftools variant calling; pipeline integration is planned.

---

## MT-RNR1

### Overview

MT-RNR1 encodes the mitochondrial 12S ribosomal RNA (rRNA), an essential structural component of the small subunit of the mitochondrial ribosome. The 12S rRNA participates in protein synthesis within mitochondria and contains a functionally critical decoding site analogous to the prokaryotic 16S rRNA aminoglycoside-binding domain. Because aminoglycoside antibiotics structurally mimic aminoacyl-tRNA and target the prokaryotic ribosome decoding site, pathogenic 12S rRNA variants that increase structural similarity to the prokaryotic rRNA confer hypersensitivity to aminoglycosides. MT-RNR1 variants are maternally inherited and are present in all cells of an affected individual.

### Pharmacogenomic Significance

The variants m.1555A>G and m.1494C>T in MT-RNR1 dramatically increase susceptibility to aminoglycoside-induced sensorineural hearing loss (AIHL), a permanent and often severe adverse effect. In individuals carrying m.1555A>G, even a single standard course of aminoglycoside therapy (gentamicin, tobramycin, amikacin) can trigger irreversible bilateral deafness. CPIC recommends avoiding all aminoglycosides in confirmed m.1555A>G or m.1494C>T carriers unless no alternative exists, and performing urgent genetic counselling of maternal relatives. Because MT-RNR1 is mitochondrial, allele frequency is reported as heteroplasmy fraction rather than a diploid allele frequency.

### Key Variants

| Variant | mtDNA Position | Effect | Clinical Significance |
|---------|---------------|--------|----------------------|
| m.1555A>G | 12S rRNA | Mimics prokaryotic rRNA | CPIC High — aminoglycoside deafness |
| m.1494C>T | 12S rRNA decoding site | Same mechanism | CPIC High — aminoglycoside deafness |

### Population Frequencies

m.1555A>G is estimated to be present in approximately 1 in 500 (0.2%) individuals of European ancestry, making it one of the more common pathogenic mitochondrial variants. Prevalence is higher in some isolated populations (Spanish: ~0.7%). Because mitochondrial inheritance is strictly maternal, a single affected individual implies all maternal-line relatives are also carriers. Heteroplasmy (mixture of mutant and wild-type mtDNA) can occur but most clinically significant cases are homoplasmic.

### CPIC Level

**CPIC Level A** — Aminoglycosides (gentamicin, tobramycin, amikacin, streptomycin). Pre-treatment screening is recommended, particularly in neonates/ICU settings where aminoglycosides are frequently used empirically.

### Tool Support in pgx-suite

| Tool | Supported |
|------|-----------|
| PyPGx | No |
| Stargazer | No |
| Aldy | No |
| StellarPGx | No |
| OptiType | No |
| mutserve | **Yes** |

> Implemented in Phase 7 via `pgx-mt.sh` (mutserve v2.0.0 standalone JAR). chrM reads are extracted from the WGS BAM, mutserve calls variants at AF ≥1%, and the result JSON reports carrier/non-carrier status for m.1555A>G and m.1494C>T with allele fraction.

---

## NAT1

### Overview

NAT1 (N-acetyltransferase 1) is a phase II conjugation enzyme that catalyses the N-acetylation of arylamine and hydrazine compounds using acetyl-CoA as the acetyl donor. NAT1 is expressed ubiquitously across human tissues including erythrocytes, lymphocytes, liver, intestine, and bladder. While it shares considerable substrate overlap with NAT2, NAT1 has distinct substrate specificities for para-aminosalicylic acid, para-aminobenzoylglutamate (a folate metabolite), and the cancer drug amonafide. NAT1 activity is constitutive and much less variable than NAT2 at the population level.

### Pharmacogenomic Significance

NAT1 polymorphisms can influence the efficacy and toxicity of amonafide (a DNA topoisomerase II inhibitor) by modulating conversion of the drug to its N-acetyl metabolite. NAT1 also acetylates the hydralazine metabolite hydrazinohydralazine and contributes to heterocyclic amine bioactivation. NAT1 null alleles have been studied as cancer susceptibility modifiers (bladder cancer in smokers), though the clinical pharmacogenomics evidence base is considerably smaller than for NAT2.

### Key Variants

| Allele | Key Variants | Effect |
|--------|-------------|--------|
| *14 | rs4986782 (R187Q) | Reduced function |
| *15 | rs4986783 (R33Stop) | No function |
| *17 | rs4987076 | Reduced function |
| *22 | rs4986781 | Reduced activity |

### Population Frequencies

NAT1 reduced-function alleles are relatively uncommon (<5% individual allele frequencies) across all major population groups, with modest ethnic variation. NAT1 activity is substantially less variable at the population level compared to NAT2, with most individuals falling within the intermediate-to-extensive range.

### CPIC Level

**CPIC Level B** — Amonafide and hydralazine (in combination with NAT2 genotype). Evidence supports genotype-informed dosing consideration but at lower confidence than Level A.

### Tool Support in pgx-suite

| Tool | Supported |
|------|-----------|
| PyPGx | Yes |
| Stargazer | No |
| Aldy | Yes |
| StellarPGx | No |
| OptiType | No |
| mutserve | No |

---

## NAT2

### Overview

NAT2 (N-acetyltransferase 2) encodes a phase II enzyme that acetylates a clinically important set of arylamine and hydrazine drugs. Unlike NAT1, NAT2 is expressed primarily in the liver and intestine and is the main determinant of the classical "slow acetylator" pharmacogenetic trait, one of the earliest pharmacogenetic observations (identified in the 1950s during isoniazid therapy). NAT2 phenotype is classified as slow, intermediate, or rapid (fast) acetylator based on the haplotypic combination of SNPs across the gene.

### Pharmacogenomic Significance

Slow NAT2 acetylators have impaired hepatic acetylation of isoniazid (first-line anti-tuberculosis drug), resulting in higher plasma concentrations, prolonged exposure, and significantly elevated risk of isoniazid-induced peripheral neuropathy and hepatotoxicity. CPIC recommends pyridoxine (vitamin B6) supplementation for all patients on isoniazid, with enhanced monitoring or dose adjustment in slow acetylators, particularly those with renal impairment. Hydralazine (antihypertensive and heart failure drug) also exhibits NAT2-dependent clearance, with slow acetylators at higher risk of hydralazine-induced lupus syndrome. Sulfasalazine and procainamide are additional NAT2 substrates.

### Key Variants

| Allele | Key Variants | Phenotype Effect |
|--------|-------------|-----------------|
| *5 (*5A–*5G) | rs1801280 (T341C, I114T) | Slow acetylator component |
| *6 (*6A, *6B) | rs1799930 (G590A, R197Q) | Slow acetylator component |
| *7 | rs1799931 (G857A, G286E) | Slow acetylator component |
| *14 | rs1801279 (G191A, R64Q) | Slow acetylator (Africans) |

Slow acetylator phenotype requires homozygosity for two slow-allele haplotypes.

### Population Frequencies

Slow acetylator phenotype is highly prevalent across all populations: ~50–60% of Europeans and Africans, ~10–15% of East Asians (lower due to lower *5/*6/*7 frequencies), and ~40–50% of South Asians. Middle Eastern populations show some of the highest slow acetylator rates (~80–90% in some cohorts). NAT2 demonstrates one of the widest documented population-frequency variations of any pharmacogene.

### CPIC Level

**CPIC Level A** — Isoniazid and hydralazine. CPIC provides specific recommendations for NAT2-guided dosing and toxicity monitoring.

### Tool Support in pgx-suite

| Tool | Supported |
|------|-----------|
| PyPGx | Yes |
| Stargazer | No |
| Aldy | Yes |
| StellarPGx | No |
| OptiType | No |
| mutserve | No |

---

## NUDT15

### Overview

NUDT15 (nudix hydrolase 15) encodes a diphosphohydrolase that metabolises thiopurine nucleoside triphosphate metabolites, specifically 6-thioguanosine triphosphate (6-TGTP) and 6-thio-deoxyguanosine triphosphate (6-TdGTP), converting them to their corresponding monophosphates and thereby reducing intracellular levels of the cytotoxic thioguanine nucleotides. NUDT15 acts as a critical detoxifying enzyme that protects cells from thiopurine-induced DNA incorporation and cytotoxicity. It was first identified as a pharmacogene through genome-wide association studies in Asian populations with thiopurine toxicity.

### Pharmacogenomic Significance

NUDT15 loss-of-function variants lead to excessive accumulation of thioguanine nucleotides (TGN) incorporated into DNA, causing severe thiopurine-induced myelosuppression (leukopenia, neutropenia). The risk phenotype with NUDT15 is indistinguishable clinically from that caused by TPMT deficiency, and both genes must be tested together to capture the full range of at-risk patients. In East and Southeast Asian populations, NUDT15*3 (R139C) is a substantially more prevalent cause of thiopurine toxicity than TPMT variants, reversing the relative clinical importance of these two genes compared to European populations.

### Key Variants

| Allele | Variant | Protein | Function | Freq (East Asians) |
|--------|---------|---------|----------|--------------------|
| *2 | rs746071566 (p.Val18_Val19insGlyVal) | Insertion | No function | Rare |
| *3 | rs116855232 (R139C) | Arg139Cys | No function | ~4–9% |
| *4 | rs147390019 (p.Arg139His) | Arg139His | No function | Rare |
| *5 | rs186364861 (p.Val18Ile + p.Leu146=) | Combined | No function | Rare |

### Population Frequencies

NUDT15*3 (R139C) has allele frequency ~9–10% in East Asians (Koreans, Japanese, Han Chinese) and ~4% in South Asians, but <1% in Europeans and Africans. This dramatically different frequency distribution means that in East Asian patient populations, NUDT15 testing detects more thiopurine-risk individuals than TPMT testing alone, and combined NUDT15 + TPMT testing is essential for complete risk stratification.

### CPIC Level

**CPIC Level A** — Thiopurines: azathioprine, mercaptopurine (6-MP), and thioguanine (6-TG). CPIC provides dose-reduction guidance for heterozygous carriers (reduce by 30–50%) and alternative therapy recommendations for homozygous poor metabolisers.

### Tool Support in pgx-suite

| Tool | Supported |
|------|-----------|
| PyPGx | Yes |
| Stargazer | No |
| Aldy | Yes |
| StellarPGx | No |
| OptiType | No |
| mutserve | No |

---

## POR

### Overview

POR (P450 oxidoreductase) encodes the electron transfer protein that serves as the obligate electron donor for all microsomal cytochrome P450 enzymes in the endoplasmic reticulum, including CYP1A1, CYP1A2, CYP2C9, CYP2C19, CYP2D6, CYP3A4, and CYP3A5. POR transfers electrons from NADPH to the CYP enzyme active site, making it an essential co-factor without which no microsomal CYP activity can occur. Severe POR mutations cause the rare Antley-Bixler syndrome with disorders of sex development; milder variants modulate CYP activity across multiple substrates simultaneously.

### Pharmacogenomic Significance

POR variants modify the activity of multiple CYP enzymes simultaneously, distinguishing it from single-gene pharmacogenomic effects. The POR*28 variant (rs1057868, A503V) has been associated with modestly increased CYP2C9 and CYP3A4 activity, potentially affecting warfarin and tacrolimus metabolism independent of the primary metabolising gene genotype. Because POR influences all microsomal CYPs, POR variants can produce complex phenotypic interactions that complicate prediction from single-gene models. POR is therefore proposed as a pan-CYP modifier gene in multi-gene dosing algorithms.

### Key Variants

| Allele | Variant | Effect |
|--------|---------|--------|
| *28 | rs1057868 (A503V) | Modestly increased CYP2C9/CYP3A4 activity |
| Multiple | Various LOF mutations | Antley-Bixler syndrome (severe; rare) |

### Population Frequencies

POR*28 (A503V) allele frequency is approximately 27–28% in Europeans, ~25% in East Asians, and ~20% in Africans — making it a common variant across populations. Its moderate functional impact (rather than null effect) and ubiquity across populations makes it a candidate for inclusion in comprehensive warfarin and tacrolimus dosing algorithms.

### CPIC Level

**Not independently assigned** — POR is not a standalone CPIC guideline gene but is included as a modifier in pharmacogenomics research studies for CYP2C9, CYP2C19, and CYP3A4/3A5 substrates. Active area of research for multi-gene dosing model inclusion.

### Tool Support in pgx-suite

| Tool | Supported |
|------|-----------|
| PyPGx | Yes |
| Stargazer | No |
| Aldy | No |
| StellarPGx | No |
| OptiType | No |
| mutserve | No |

---

## RYR1

### Overview

RYR1 (ryanodine receptor 1) encodes the skeletal muscle sarcoplasmic reticulum calcium release channel, the largest known ion channel protein (5,037 amino acids). RYR1 is the principal calcium release channel of the sarcoplasmic reticulum and is the central effector of excitation-contraction coupling in skeletal muscle, responding to signals from CACNA1S (CaV1.1) voltage sensors in the transverse tubule membrane. Pathogenic RYR1 variants cause a spectrum of myopathies including malignant hyperthermia susceptibility (MHS), central core disease (CCD), multi-minicore disease, and nemaline myopathy.

### Pharmacogenomic Significance

RYR1 gain-of-function variants are the primary genetic cause of malignant hyperthermia susceptibility (MHS), accounting for ~70% of MH-susceptible pedigrees. Exposure of RYR1 MHS individuals to volatile halogenated anaesthetics (halothane, sevoflurane, isoflurane, desflurane) or succinylcholine triggers uncontrolled calcium release from the sarcoplasmic reticulum, resulting in muscle rigidity, life-threatening hyperthermia (>40°C), metabolic acidosis, rhabdomyolysis, and hyperkalemia. Without prompt treatment with dantrolene (IV, 2.5 mg/kg bolus), MH has a mortality approaching 70%. Pre-operative RYR1 genotyping allows identification of at-risk individuals and selection of total intravenous anaesthesia (TIVA) protocols.

### Key Variants

| Variant | cDNA | Protein | Class |
|---------|------|---------|-------|
| rs118192165 | c.7300G>A | p.Val2434Ile | MHS (European founder) |
| rs193922771 | c.1840C>T | p.Arg614Cys | MHS (Class I) |
| rs118192166 | c.14387A>G | p.Tyr4796Cys | MHS |
| rs193922844 | c.6617C>T | p.Thr2206Met | MHS |

The European Malignant Hyperthermia Group (EMHG) has catalogued >300 RYR1 variants associated with MHS, of which ~35 are classified as causative (class I/II evidence).

### Population Frequencies

RYR1 pathogenic MH variants are present across all ethnic populations, with an aggregate estimated MHS prevalence of 1 in 2,000 to 1 in 10,000 (genotype-based estimates vary). No single MH founder variant predominates globally, though some variants (e.g., p.Val2434Ile) are over-represented in European pedigrees. The clinical MH incidence is lower than genetic prevalence, as environmental triggers are required.

### CPIC Level

**CPIC Level A** — Volatile anaesthetic agents (halothane, sevoflurane, isoflurane, desflurane) and succinylcholine. CPIC/CPNDS guidelines recommend genetically informed anaesthetic planning and avoidance of triggering agents.

### Tool Support in pgx-suite

| Tool | Supported |
|------|-----------|
| PyPGx | **Yes** |
| Stargazer | **Yes** |
| Aldy | **Yes** |
| StellarPGx | No |
| OptiType | No |
| mutserve | No |

---

## SLCO1B1

### Overview

SLCO1B1 encodes organic anion transporting polypeptide 1B1 (OATP1B1), a hepatocyte-specific uptake transporter located on the sinusoidal (basolateral) membrane of hepatocytes. OATP1B1 mediates the sodium-independent uptake of a wide range of endogenous substrates (bilirubin, bile acids, thyroid hormones) and xenobiotics from portal blood into hepatocytes, serving as the critical first step in hepatic first-pass elimination. Major drug substrates include statins (simvastatin, atorvastatin, rosuvastatin, pravastatin), the anticoagulant repaglinide, the antibiotic rifampicin, and several HIV protease inhibitors.

### Pharmacogenomic Significance

The SLCO1B1*5 variant (rs4149056, Val174Ala) reduces OATP1B1 transport function, impairing hepatic uptake of statins and dramatically increasing their systemic plasma exposure. Homozygous *5/*5 carriers receiving simvastatin have an 8-fold increased risk of simvastatin-induced myopathy (myalgia, myositis, rhabdomyolysis) compared to wild-type individuals. The overall CPIC-estimated odds ratio for myopathy in *5 carriers on simvastatin is approximately 2.6 per allele. CPIC recommends lower simvastatin doses (≤20 mg daily) or alternative statins (pravastatin, fluvastatin) for *5 carriers.

### Key Variants

| Allele | Variant | Protein | Function | Myopathy Risk |
|--------|---------|---------|----------|--------------|
| *1A | Reference | Val174 | Normal | Baseline |
| *1B | rs2306283 (Asn130Asp) | Asn130Asp | Normal/slightly increased | Not elevated |
| *5 | rs4149056 (Val174Ala) | Val174Ala | ~50% reduced | ~2.6x per allele |
| *15 | *1B + *5 | Both changes | Reduced | ~2.6x per *5 allele |
| *17 | rs4149056 + rs2900478 | Multiple | Reduced | Elevated |

### Population Frequencies

The SLCO1B1*5 (rs4149056) allele frequency is approximately 14–17% in Europeans, 13–15% in East Asians, and 2–4% in Africans. The lower frequency in Africans, combined with the higher frequency of other SLCO1B1 variants in this population, again highlights the importance of comprehensive genotyping beyond the single-SNP assay approach that is predominantly validated in Europeans.

### CPIC Level

**CPIC Level A** — Simvastatin (statin-induced myopathy). CPIC provides specific dose-capping and statin-switching guidance by SLCO1B1 phenotype. Also incorporated in statin guidelines for atorvastatin, pravastatin, and rosuvastatin.

### Tool Support in pgx-suite

| Tool | Supported |
|------|-----------|
| PyPGx | Yes |
| Stargazer | Yes |
| Aldy | Yes |
| StellarPGx | Yes |
| OptiType | No |
| mutserve | No |

---

## TPMT

### Overview

TPMT (thiopurine S-methyltransferase) encodes a cytosolic enzyme that catalyses the S-methylation of aromatic and heterocyclic thiol compounds, including the thiopurine drugs azathioprine, mercaptopurine (6-MP), and thioguanine (6-TG). TPMT competes with the enzyme hypoxanthine phosphoribosyltransferase (HPRT) for thiopurine substrates: HPRT channels drugs toward active thioguanine nucleotide (TGN) formation (cytotoxic), while TPMT inactivates them via S-methylation. TPMT was one of the first pharmacogenes characterised biochemically, with its polymorphic nature described in the early 1980s.

### Pharmacogenomic Significance

TPMT poor metabolisers (homozygous for loss-of-function alleles, ~1 in 300 individuals) who receive standard thiopurine doses accumulate extremely high TGN concentrations, causing life-threatening myelosuppression. CPIC recommends drastically reduced doses (10-fold reduction or therapy avoidance) in TPMT poor metabolisers, with 30–50% dose reductions in heterozygous intermediate metabolisers. For complete thiopurine risk stratification, TPMT genotyping must be performed alongside NUDT15 genotyping, as both genes contribute independently and additively to myelosuppression risk (particularly relevant in Asian populations where NUDT15*3 predominates).

### Key Variants

| Allele | Variant | Protein | Function | Freq (Europeans) |
|--------|---------|---------|----------|-----------------|
| *2 | rs1800462 (G238C) | Ala80Pro | No function | ~0.5% |
| *3A | rs1800460 + rs1142345 | Ala154Thr + Tyr240Cys | No function | ~5% |
| *3B | rs1800460 (G460A) | Ala154Thr | No function | <0.5% |
| *3C | rs1142345 (A719G) | Tyr240Cys | No function | ~0.8% |
| *4 | rs1800584 (splice) | — | No function | Rare |

### Population Frequencies

TPMT poor metaboliser frequency is approximately 1 in 300 (0.3%) across European and African populations, with *3A the predominant no-function allele in Europeans. *3C is the predominant no-function allele in Africans. East Asians have lower TPMT poor metaboliser frequency (~1 in 1,500) but higher NUDT15 risk, again emphasising the need for combined testing in Asian populations.

### CPIC Level

**CPIC Level A** — Thiopurines: azathioprine, mercaptopurine (6-MP), and thioguanine (6-TG). CPIC guidelines recommend mandatory pre-treatment genotyping for all patients initiating thiopurine therapy.

### Tool Support in pgx-suite

| Tool | Supported |
|------|-----------|
| PyPGx | Yes |
| Stargazer | Yes |
| Aldy | Yes |
| StellarPGx | Yes |
| OptiType | No |
| mutserve | No |

---

## UGT1A1

### Overview

UGT1A1 (UDP-glucuronosyltransferase 1A1) encodes the primary hepatic enzyme responsible for the glucuronidation and biliary elimination of bilirubin (converting unconjugated to conjugated bilirubin). UGT1A1 also glucuronidates the active metabolite of irinotecan (SN-38), estrogens, and the HIV integrase inhibitor atazanavir. Complete UGT1A1 deficiency causes Crigler-Najjar syndrome; partial deficiency causes Gilbert syndrome, a benign condition of mild unconjugated hyperbilirubinaemia affecting approximately 5–10% of Western populations.

### Pharmacogenomic Significance

SN-38, the active metabolite of the chemotherapy drug irinotecan, is glucuronidated by UGT1A1 to the inactive SN-38G for biliary excretion. UGT1A1 poor metabolisers (particularly *28/*28 and *6/*6 homozygotes) have impaired SN-38 clearance, leading to prolonged exposure and severe irinotecan toxicity: grade 3/4 neutropenia and diarrhoea. FDA labelling of irinotecan recommends reduced starting doses in UGT1A1*28 homozygotes. Atazanavir competitively inhibits UGT1A1, causing predictable hyperbilirubinaemia that is more pronounced in *28 carriers, though this is generally a benign cosmetic effect rather than clinical toxicity.

### Key Variants

| Allele | Variant | Effect | Freq (Europeans) |
|--------|---------|--------|-----------------|
| *28 | rs8175347 (TA)7 in promoter TATA box | ~70% reduced expression | ~36% |
| *6 | rs4148323 (G211A, Gly71Arg) | Reduced activity | <1% (East Asians ~13%) |
| *27 | rs35350960 (C686A, Pro229Gln) | Reduced | Rare |
| *36 | (TA)5 repeat | Increased expression | Rare |
| *80 | rs887829 (C>T, intron) | Proxy for *28 | — |

### Population Frequencies

UGT1A1*28 (7-TA repeat) is common in Europeans (~36% allele frequency), Africans (~40%), and less frequent in East Asians (~9%). East Asian populations have instead a high prevalence of UGT1A1*6 (G71R, ~13% allele frequency), which is rare in Europeans. Irinotecan risk therefore must be assessed with *28 primarily in Europeans/Africans and *6 in East Asian populations to avoid missing at-risk patients.

### CPIC Level

**CPIC Level A** — Irinotecan (chemotherapy-induced toxicity). Also relevant for atazanavir (benign hyperbilirubinaemia). CPIC provides specific dose-reduction guidance for *28/*28 and *6/*6 homozygotes receiving irinotecan.

### Tool Support in pgx-suite

| Tool | Supported |
|------|-----------|
| PyPGx | Yes |
| Stargazer | Yes |
| Aldy | Yes |
| StellarPGx | Yes |
| OptiType | No |
| mutserve | No |

---

## VKORC1

### Overview

VKORC1 (vitamin K epoxide reductase complex subunit 1) encodes the catalytic subunit of the enzyme complex responsible for recycling vitamin K from its inactive epoxide form (vitamin K 2,3-epoxide) back to the reduced form (vitamin KH2) that is required as a cofactor for gamma-carboxylation of coagulation factors II, VII, IX, and X, as well as proteins C, S, and Z. VKORC1 is the direct molecular target of coumarin anticoagulants (warfarin, acenocoumarol, phenprocoumon), which competitively inhibit vitamin K recycling and thereby deplete functional clotting factor activity.

### Pharmacogenomic Significance

VKORC1 is the single most important pharmacogenomic determinant of warfarin dose requirement, explaining approximately 25–30% of warfarin dose variability. The promoter variant rs9923231 (-1639G>A) reduces VKORC1 mRNA transcription by approximately 30–45%, resulting in lower VKORC1 enzyme levels and consequently greater sensitivity to warfarin inhibition. Patients homozygous for the -1639A allele (haplotype A/A) require approximately 40% lower warfarin doses to achieve therapeutic anticoagulation compared to G/G individuals. VKORC1 genotyping is incorporated into all evidence-based warfarin dosing algorithms (IWPC, Gage).

### Key Variants

| Variant | Region | Haplotype | Effect | Freq (Europeans) |
|---------|--------|-----------|--------|-----------------|
| rs9923231 (-1639G>A) | Promoter | A = low expression | Lower warfarin dose | ~40% (A allele) |
| rs9934438 (1173C>T) | Intron 1 | Linked to promoter -1639 | Same effect | ~40% |
| rs2359612 (2255C>T) | Intron 2 | VKORC1 haplotyping | Modifier | ~40% |

The -1639G>A SNP (rs9923231) is the primary VKORC1 clinical reporting variant and is in strong linkage disequilibrium with several other VKORC1 intronic variants used for haplotype assignment.

### Population Frequencies

The VKORC1 -1639A (low-expression, warfarin-sensitive) allele frequency is approximately 38–42% in Europeans, 88–92% in East Asians (explaining why Asians require markedly lower warfarin doses), and ~5–10% in Africans (explaining why Africans generally require higher doses). This trimodal ethnic distribution of VKORC1 allele frequencies is a principal driver of population differences in warfarin therapeutic dose requirements and is one of the most striking examples of pharmacogenomics-driven health disparities.

### CPIC Level

**CPIC Level A** — Warfarin, acenocoumarol, and phenprocoumon. VKORC1 is a cornerstone component of all warfarin pharmacogenomics dosing algorithms, alongside CYP2C9 and CYP4F2.

### Tool Support in pgx-suite

| Tool | Supported |
|------|-----------|
| PyPGx | Yes |
| Stargazer | Yes |
| Aldy | Yes |
| StellarPGx | Yes |
| OptiType | No |
| mutserve | No |

---

*Document generated for pgx-suite pipeline — GRCh38 reference — 31 pharmacogenes covering CPIC Level A and selected Level B/research genes. For the most current CPIC guidelines, consult [cpicpgx.org](https://cpicpgx.org).*
