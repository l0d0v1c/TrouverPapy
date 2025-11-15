#!/usr/bin/env python3
"""
Script pour rechercher des SNPs spécifiques dans un fichier VCF
Utile pour explorer de nouveaux marqueurs génétiques
"""

import gzip
import sys
import argparse


def search_snps(vcf_file, snp_list):
    """
    Recherche une liste de SNPs (rs IDs) dans le VCF

    Args:
        vcf_file: chemin vers le fichier VCF.gz
        snp_list: liste des rs IDs à rechercher
    """

    # Extraire les noms d'échantillons
    samples = []
    with gzip.open(vcf_file, 'rt') as f:
        for line in f:
            if line.startswith('#CHROM'):
                fields = line.strip().split('\t')
                samples = fields[9:]
                break

    print(f"Échantillons: {', '.join(samples)}\n")
    print("="*80)

    # Rechercher les SNPs
    found_snps = set()

    with gzip.open(vcf_file, 'rt') as f:
        for line in f:
            if line.startswith('#'):
                continue

            fields = line.strip().split('\t')
            chrom = fields[0]
            pos = fields[1]
            rs_id = fields[2]
            ref = fields[3]
            alt = fields[4]

            if rs_id in snp_list:
                found_snps.add(rs_id)

                print(f"\n✓ SNP TROUVÉ: {rs_id}")
                print(f"  Chromosome: {chrom}")
                print(f"  Position: {pos}")
                print(f"  Allèle de référence: {ref}")
                print(f"  Allèle alternatif: {alt}")
                print(f"  {'-'*76}")

                # Afficher les génotypes pour chaque échantillon
                sample_data = fields[9:]
                for i, sample in enumerate(samples):
                    gt = sample_data[i].split(':')[0]

                    # Convertir en format lisible
                    if gt == './.':
                        genotype = "Pas de données"
                    elif gt == '0/0':
                        genotype = f"{ref}/{ref} (homozygote référence)"
                    elif gt == '0/1' or gt == '1/0':
                        genotype = f"{ref}/{alt} (hétérozygote)"
                    elif gt == '1/1':
                        genotype = f"{alt}/{alt} (homozygote alternatif)"
                    else:
                        genotype = gt

                    print(f"  {sample}: {genotype}")

                print("="*80)

    # SNPs non trouvés
    not_found = set(snp_list) - found_snps
    if not_found:
        print(f"\n\n⚠ SNPs NON TROUVÉS dans le VCF ({len(not_found)}):")
        for snp in sorted(not_found):
            print(f"  - {snp}")

    print(f"\n\nRésumé: {len(found_snps)}/{len(snp_list)} SNPs trouvés")


def main():
    parser = argparse.ArgumentParser(
        description='Recherche de SNPs spécifiques dans un fichier VCF',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Exemples d'utilisation:

1. Rechercher des SNPs à partir d'une liste:
   python search_custom_snps.py AITI_43_55.vcf.gz rs12913832 rs1426654 rs4988235

2. Rechercher des SNPs à partir d'un fichier:
   python search_custom_snps.py AITI_43_55.vcf.gz --file snp_list.txt

SNPs d'intérêt courants:
  - rs12913832 : Couleur des yeux (HERC2/OCA2)
  - rs1426654  : Pigmentation de la peau (SLC24A5)
  - rs4988235  : Tolérance au lactose (LCT)
  - rs1815739  : Performance musculaire (ACTN3)
  - rs1800562  : Hémochromatose (HFE)
  - rs6025     : Facteur V Leiden (thrombose)
        """
    )

    parser.add_argument('vcf_file', help='Fichier VCF.gz à analyser')
    parser.add_argument('snps', nargs='*', help='Liste de SNPs (rs IDs) à rechercher')
    parser.add_argument('--file', '-f', help='Fichier texte contenant les rs IDs (un par ligne)')

    args = parser.parse_args()

    # Collecter les SNPs à rechercher
    snp_list = []

    if args.snps:
        snp_list.extend(args.snps)

    if args.file:
        try:
            with open(args.file, 'r') as f:
                for line in f:
                    line = line.strip()
                    # Ignorer les lignes vides et les commentaires
                    if not line or line.startswith('#'):
                        continue
                    # Extraire seulement le rs ID (avant tout commentaire sur la ligne)
                    rs_id = line.split()[0]  # Prend le premier mot
                    if rs_id.startswith('rs'):
                        snp_list.append(rs_id)
        except FileNotFoundError:
            print(f"Erreur: Fichier {args.file} non trouvé")
            sys.exit(1)

    if not snp_list:
        parser.print_help()
        print("\nErreur: Vous devez spécifier au moins un SNP à rechercher")
        sys.exit(1)

    # Supprimer les doublons
    snp_list = list(set(snp_list))

    print(f"Recherche de {len(snp_list)} SNPs dans {args.vcf_file}...\n")
    search_snps(args.vcf_file, snp_list)


if __name__ == "__main__":
    main()
