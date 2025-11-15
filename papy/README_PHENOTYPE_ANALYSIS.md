# Analyse Phénotypique et Pathologique à partir de VCF

Ce répertoire contient des outils pour analyser les marqueurs génétiques (SNPs) dans un fichier VCF et les associer avec des traits phénotypiques et des pathologies.

## Fichiers

### 1. `phenotype_analysis.py`
Script principal qui analyse automatiquement ~20 SNPs importants couvrant:

**Traits de pigmentation:**
- Couleur des yeux (HERC2/OCA2, SLC24A4)
- Couleur des cheveux (MC1R, IRF4)
- Pigmentation de la peau (SLC24A5, SLC45A2)

**Métabolisme:**
- Tolérance au lactose (LCT/MCM6) - rs4988235
- Métabolisme de l'alcool (ALDH2)

**Maladies génétiques:**
- Hémochromatose (HFE) - surcharge en fer
- Facteur V Leiden (F5) - risque de thrombose
- Mucoviscidose (CFTR)

**Autres traits:**
- Performance musculaire (ACTN3)
- Type de cérumen (ABCC11)
- Groupe sanguin (ABO)

#### Utilisation:
```bash
python3 phenotype_analysis.py AITI_43_55.vcf.gz
```

Le rapport est sauvegardé dans `phenotype_report.txt`

### 2. `search_custom_snps.py`
Outil de recherche flexible pour explorer des SNPs spécifiques

#### Utilisation:

**Rechercher des SNPs individuels:**
```bash
python3 search_custom_snps.py AITI_43_55.vcf.gz rs12913832 rs1426654
```

**Rechercher à partir d'un fichier:**
```bash
python3 search_custom_snps.py AITI_43_55.vcf.gz --file ma_liste_snps.txt
```

**Exemple avec aide:**
```bash
python3 search_custom_snps.py --help
```

## Résultats pour AITI_43

### Caractéristiques physiques:
- **Yeux**: Marrons/bruns (A/A au rs12913832)
- **Cheveux**: Foncés (C/C au rs12203592)
- **Peau**: Foncée/normale - génotype ancestral (A/A au rs1426654)

### Métabolisme:
- **Lactose**: ✓ TOLÉRANT (G/G au rs4988235)
  - Ce résultat est intéressant car montre la mutation européenne de persistance de la lactase, déjà présente à l'Âge du Bronze ancien
- **Alcool**: Métabolisme normal (G/G au rs671)

### Performance physique:
- **Type musculaire**: Endurance (T/T au rs1815739 - perte ACTN3)

### Santé:
- ✓ Aucun variant pathogène majeur détecté
- Pas de risque d'hémochromatose
- Pas de mutation Facteur V Leiden

### Notes:
AITI_55_d a très peu de données disponibles (couverture faible)

## Ajouter de nouveaux SNPs d'intérêt

Pour analyser d'autres SNPs, vous pouvez:

1. **Modifier `phenotype_analysis.py`**: Ajouter des SNPs dans le dictionnaire `PHENOTYPE_SNPS`

2. **Utiliser `search_custom_snps.py`**: Pour une recherche rapide sans modifier le code

## SNPs supplémentaires intéressants

### Maladies cardiovasculaires:
- rs1333049 (9p21.3) - Risque de maladie coronarienne
- rs10757274 (9p21) - Infarctus du myocarde

### Métabolisme:
- rs9939609 (FTO) - Prédisposition à l'obésité
- rs7903146 (TCF7L2) - Diabète de type 2
- rs1801282 (PPARG) - Sensibilité à l'insuline

### Réponse immunitaire:
- rs6897932 (IL7R) - Sclérose en plaques
- rs2476601 (PTPN22) - Maladies auto-immunes

### Détoxification:
- rs1695 (GSTP1) - Métabolisme des toxines
- rs1051740 (CYP1A1) - Métabolisme des carcinogènes

### Santé osseuse:
- rs2228570 (VDR) - Récepteur de la vitamine D
- rs1800247 (COL1A1) - Ostéoporose

## Bases de données utiles

Pour trouver plus d'informations sur les SNPs:

- **dbSNP**: https://www.ncbi.nlm.nih.gov/snp/
- **ClinVar**: https://www.ncbi.nlm.nih.gov/clinvar/ (variants pathogènes)
- **GWAS Catalog**: https://www.ebi.ac.uk/gwas/ (associations phénotypiques)
- **SNPedia**: https://www.snpedia.com/ (wiki des SNPs)
- **Ensembl**: https://www.ensembl.org/ (annotation génomique)

## Avertissements

⚠️ **Important:**
- Ces analyses sont à but informatif et de recherche uniquement
- Ne constituent PAS un diagnostic médical
- Les interprétations sont basées sur les données scientifiques actuelles
- Pour l'ADN ancien, la couverture peut être faible (beaucoup de données manquantes)
- Les associations génotype-phénotype sont souvent probabilistes et influencées par de nombreux autres facteurs

## Contact et contributions

Pour ajouter de nouveaux SNPs ou améliorer les analyses, modifiez les scripts et testez avec vos données.
