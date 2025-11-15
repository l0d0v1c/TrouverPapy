#!/usr/bin/env python3
"""
Calcule des statistiques de couverture pour un fichier VCF
Utile pour évaluer la qualité des données d'ADN ancien
"""

import gzip
import sys


def calculate_vcf_stats(vcf_file):
    """
    Calcule les statistiques de couverture génotypique par échantillon
    """

    # Lire les noms d'échantillons
    samples = []
    with gzip.open(vcf_file, 'rt') as f:
        for line in f:
            if line.startswith('#CHROM'):
                fields = line.strip().split('\t')
                samples = fields[9:]
                break

    print(f"{'='*80}")
    print(f"STATISTIQUES VCF - {vcf_file}")
    print(f"{'='*80}\n")
    print(f"Nombre d'échantillons: {len(samples)}\n")

    # Initialiser les compteurs
    total_variants = 0
    genotype_counts = {sample: {'total': 0, 'missing': 0, 'hom_ref': 0, 'het': 0, 'hom_alt': 0}
                      for sample in samples}
    chromosome_counts = {}

    # Parcourir le VCF
    with gzip.open(vcf_file, 'rt') as f:
        for line in f:
            if line.startswith('#'):
                continue

            fields = line.strip().split('\t')
            chrom = fields[0]
            total_variants += 1

            # Compter par chromosome
            chromosome_counts[chrom] = chromosome_counts.get(chrom, 0) + 1

            # Analyser les génotypes
            sample_data = fields[9:]
            for i, sample in enumerate(samples):
                gt = sample_data[i].split(':')[0]
                genotype_counts[sample]['total'] += 1

                if gt == './.':
                    genotype_counts[sample]['missing'] += 1
                elif gt == '0/0':
                    genotype_counts[sample]['hom_ref'] += 1
                elif gt in ['0/1', '1/0']:
                    genotype_counts[sample]['het'] += 1
                elif gt == '1/1':
                    genotype_counts[sample]['hom_alt'] += 1

    # Afficher les statistiques générales
    print(f"STATISTIQUES GÉNÉRALES")
    print(f"{'-'*80}")
    print(f"Nombre total de variants: {total_variants:,}")
    print(f"Nombre de chromosomes: {len(chromosome_counts)}")
    print()

    # Statistiques par chromosome
    print(f"RÉPARTITION PAR CHROMOSOME")
    print(f"{'-'*80}")
    for chrom in sorted(chromosome_counts.keys(), key=lambda x: (len(x), x)):
        count = chromosome_counts[chrom]
        percent = (count / total_variants) * 100
        print(f"  chr{chrom:>2s}: {count:>8,} variants ({percent:>5.2f}%)")
    print()

    # Statistiques par échantillon
    print(f"COUVERTURE PAR ÉCHANTILLON")
    print(f"{'='*80}\n")

    for sample in samples:
        stats = genotype_counts[sample]
        total = stats['total']
        missing = stats['missing']
        called = total - missing
        coverage = (called / total) * 100 if total > 0 else 0

        print(f"Échantillon: {sample}")
        print(f"{'-'*80}")
        print(f"  Total de sites analysés: {total:,}")
        print(f"  Sites avec données:      {called:,} ({coverage:.2f}%)")
        print(f"  Sites manquants (./.):   {missing:,} ({(missing/total)*100:.2f}%)")
        print()
        print(f"  Répartition des génotypes appelés:")
        print(f"    Homozygote référence (0/0): {stats['hom_ref']:>8,} ({(stats['hom_ref']/total)*100:>5.2f}%)")
        print(f"    Hétérozygote (0/1):          {stats['het']:>8,} ({(stats['het']/total)*100:>5.2f}%)")
        print(f"    Homozygote alternatif (1/1): {stats['hom_alt']:>8,} ({(stats['hom_alt']/total)*100:>5.2f}%)")
        print()

        # Qualité pour ADN ancien
        if coverage < 30:
            quality = "⚠️  TRÈS FAIBLE - Données insuffisantes pour la plupart des analyses"
        elif coverage < 50:
            quality = "⚠️  FAIBLE - Utilisable pour quelques analyses ciblées"
        elif coverage < 70:
            quality = "🟡 MODÉRÉE - Acceptable pour certaines analyses"
        elif coverage < 90:
            quality = "🟢 BONNE - Utilisable pour la plupart des analyses"
        else:
            quality = "🟢 EXCELLENTE - Données de haute qualité"

        print(f"  Évaluation de qualité: {quality}")
        print(f"{'='*80}\n")

    # Calcul de l'hétérozygotie
    print(f"TAUX D'HÉTÉROZYGOTIE")
    print(f"{'-'*80}")
    for sample in samples:
        stats = genotype_counts[sample]
        called = stats['total'] - stats['missing']
        if called > 0:
            het_rate = (stats['het'] / called) * 100
            print(f"  {sample}: {het_rate:.3f}%")

            # Commentaire sur l'hétérozygotie
            if het_rate < 15:
                print(f"    ⚠️  Faible - Possible consanguinité ou erreurs de génotypage")
            elif het_rate < 25:
                print(f"    🟢 Normal pour populations européennes anciennes")
            else:
                print(f"    🟡 Élevé - Possiblement population admixée")
        print()


def main():
    if len(sys.argv) < 2:
        print("Usage: python vcf_stats.py <fichier.vcf.gz>")
        sys.exit(1)

    vcf_file = sys.argv[1]
    calculate_vcf_stats(vcf_file)


if __name__ == "__main__":
    main()
