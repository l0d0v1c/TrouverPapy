# Résumé de l'analyse phénotypique et pathologique - AITI_43 & AITI_55

**Date**: 15 novembre 2025
**Fichier analysé**: AITI_43_55.vcf.gz
**Nombre de variants**: 1,233,013 SNPs

## Individus analysés

1. **Germany_Lech_EBA_AITI_43** - Bonne couverture génomique
2. **Germany_Lech_EBA_father.or.son.AITI_43_lc_AITI_55_d** - Couverture faible (peu de données)

---

## AITI_43 - Profil phénotypique complet

### 🎨 Pigmentation et apparence physique

| Trait | SNP | Génotype | Phénotype prédit |
|-------|-----|----------|------------------|
| **Couleur des yeux** | rs12913832 | A/A | **Marrons/bruns** (ancestral) |
| | rs12896399 | G/G | Yeux foncés |
| **Couleur des cheveux** | rs12203592 | C/C | **Cheveux foncés** |
| **Couleur de la peau** | rs1426654 | A/A | **Peau foncée** (génotype ancestral) |
| | rs16891982 | C/C | Peau normale/foncée |

**Interprétation**: AITI_43 possède le profil de pigmentation ancestral typique - yeux marrons, cheveux foncés et peau plus foncée que les populations européennes modernes. Ceci est cohérent avec un individu de l'Âge du Bronze ancien, avant la fixation complète des allèles de peau claire en Europe.

### 🥛 Métabolisme

| Trait | SNP | Génotype | Phénotype prédit |
|-------|-----|----------|------------------|
| **Tolérance au lactose** | rs4988235 | G/G | ✅ **TOLÉRANT** (mutation européenne) |
| **Métabolisme alcool** | rs671 | G/G | Normal |

**⭐ Découverte importante**: AITI_43 est **homozygote pour la persistance de la lactase** (G/G), ce qui indique qu'il était tolérant au lactose. Cette mutation est caractéristique des populations européennes qui ont développé l'élevage laitier. Sa présence à l'Âge du Bronze ancien en Allemagne du Sud (région de Lech) est cohérente avec l'adoption de l'agriculture et de l'élevage dans cette zone.

### 💪 Performance physique

| Trait | SNP | Génotype | Phénotype prédit |
|-------|-----|----------|------------------|
| **Type musculaire** | rs1815739 | T/T | Type **endurance** (perte ACTN3) |

**Interprétation**: Le génotype T/T au locus ACTN3 indique une perte de la protéine ACTN3 dans les fibres musculaires rapides, ce qui est associé à une meilleure performance en endurance plutôt qu'en sprint.

### 🏥 Santé et variants pathogènes

| Condition | SNP | Génotype | Risque |
|-----------|-----|----------|--------|
| **Hémochromatose (C282Y)** | rs1800562 | G/G | ✅ Aucun risque |
| **Hémochromatose (H63D)** | rs1799945 | C/C | ✅ Aucun risque |
| **Facteur V Leiden** | rs6025 | ./. | Pas de données |
| **MTHFR C677T** | rs1801133 | G/G | ✅ Normal |

**✅ Bilan**: Aucun variant pathogène majeur détecté. AITI_43 ne présente pas les mutations courantes associées à l'hémochromatose, condition relativement fréquente en Europe.

### 🧬 Autres traits

| Trait | SNP | Génotype | Interprétation |
|-------|-----|----------|----------------|
| **Type de cérumen** | rs17822931 | ./. | Pas de données |
| **Goût amer (PTC)** | rs713598 | C/C | Probablement capable de goûter l'amertume |

---

## AITI_55_d - Couverture limitée

**⚠️ Note importante**: AITI_55_d présente une très faible couverture génomique. La plupart des SNPs analysés retournent "./." (pas de données). Ceci est typique de l'ADN ancien de faible qualité ou de couverture de séquençage insuffisante.

**Recommandation**: Pour obtenir un profil phénotypique fiable pour AITI_55_d, il faudrait:
- Augmenter la profondeur de séquençage
- Utiliser des techniques d'enrichissement (capture ciblée)
- Ou se concentrer uniquement sur AITI_43 pour les analyses phénotypiques

---

## Contexte archéologique et génétique

### Période et lieu
- **Culture**: Âge du Bronze ancien (Early Bronze Age - EBA)
- **Région**: Lech, Allemagne du Sud
- **Période estimée**: ~2200-1800 BCE

### Cohérence des résultats avec le contexte

1. **Pigmentation foncée**: Cohérent avec les populations de l'EBA en Europe centrale, avant la fixation complète des allèles de peau claire

2. **Persistance de la lactase**: Indique que la mutation de tolérance au lactose était déjà présente et possiblement sélectionnée dans cette population pastorale

3. **Pas de variants nordiques récents**: L'absence de certains allèles (comme la peau très claire) est cohérente avec une population pré-migration massive des populations nordiques

---

## Outils créés pour cette analyse

1. **phenotype_analysis.py**: Script automatisé analysant ~20 SNPs clés
2. **search_custom_snps.py**: Outil de recherche flexible pour explorer des SNPs spécifiques
3. **snp_lists/**: Collection de listes de SNPs par catégorie:
   - cardiovascular.txt
   - metabolism.txt
   - immune_system.txt
   - neurological.txt
   - physical_traits.txt

### Utilisation

```bash
# Analyse complète automatique
python3 phenotype_analysis.py AITI_43_55.vcf.gz

# Recherche personnalisée (ex: traits cardiovasculaires)
python3 search_custom_snps.py AITI_43_55.vcf.gz --file snp_lists/cardiovascular.txt

# Recherche de SNPs individuels
python3 search_custom_snps.py AITI_43_55.vcf.gz rs12913832 rs1426654
```

---

## Limitations et avertissements

⚠️ **Important**:
1. Ces analyses sont basées sur des SNPs individuels et ne capturent pas la complexité polygénique de nombreux traits
2. Les interprétations phénotypiques sont probabilistes, non déterministes
3. L'ADN ancien peut avoir des dommages post-mortem affectant les génotypes
4. Ces résultats sont à but de recherche uniquement et ne constituent pas un diagnostic médical
5. De nombreux SNPs pathogènes ne sont pas présents dans ce panel 1240K (conçu pour la génétique des populations)

---

## Références et bases de données utilisées

- **dbSNP**: Base de données des variants génétiques
- **ClinVar**: Variants cliniquement significatifs
- **GWAS Catalog**: Associations génotype-phénotype
- **SNPedia**: Wiki communautaire des SNPs
- **Publications scientifiques** sur la pigmentation, le métabolisme et l'évolution des populations européennes

---

## Suggestions pour analyses futures

1. **Analyse de parenté**: Confirmer la relation père-fils entre AITI_43 et AITI_55
2. **Haplotypes**: Analyser les haplotypes du chromosome Y et de l'ADN mitochondrial
3. **Admixture**: Comparer avec les populations contemporaines et anciennes
4. **Sélection**: Rechercher des signaux de sélection positive
5. **Annotation fonctionnelle**: Annoter tous les variants avec VEP ou SnpEff

Pour plus d'informations, consultez **README_PHENOTYPE_ANALYSIS.md**
