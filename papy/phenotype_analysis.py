#!/usr/bin/env python3
"""
Analyse phénotypique et pathologique à partir d'un fichier VCF
Extrait les génotypes pour des SNPs associés à des traits et pathologies connus
"""

import gzip
import sys
from collections import defaultdict

# Base de données des SNPs d'intérêt organisés par catégorie
PHENOTYPE_SNPS = {
    "Pigmentation - Yeux": {
        "rs12913832": {
            "gene": "HERC2/OCA2",
            "chr": "15",
            "pos": "28365618",
            "trait": "Couleur des yeux",
            "alleles": {
                "A/A": "Yeux marrons/bruns (ancestral)",
                "A/G": "Yeux verts/noisette (intermédiaire)",
                "G/G": "Yeux bleus"
            }
        },
        "rs12896399": {
            "gene": "SLC24A4",
            "chr": "14",
            "trait": "Couleur des yeux",
            "alleles": {
                "G/G": "Yeux foncés",
                "G/T": "Yeux intermédiaires",
                "T/T": "Contribution aux yeux clairs"
            }
        }
    },

    "Pigmentation - Cheveux": {
        "rs1805007": {
            "gene": "MC1R",
            "chr": "16",
            "trait": "Cheveux roux",
            "alleles": {
                "C/C": "Pas de cheveux roux",
                "C/T": "Porteur variant roux",
                "T/T": "Cheveux roux probable"
            }
        },
        "rs12203592": {
            "gene": "IRF4",
            "chr": "6",
            "trait": "Couleur des cheveux",
            "alleles": {
                "C/C": "Cheveux foncés",
                "C/T": "Cheveux intermédiaires",
                "T/T": "Cheveux clairs/blonds"
            }
        }
    },

    "Pigmentation - Peau": {
        "rs1426654": {
            "gene": "SLC24A5",
            "chr": "15",
            "trait": "Pigmentation de la peau",
            "alleles": {
                "A/A": "Peau foncée (ancestral)",
                "A/G": "Peau intermédiaire",
                "G/G": "Peau claire (dérivé européen)"
            }
        },
        "rs16891982": {
            "gene": "SLC45A2",
            "chr": "5",
            "trait": "Pigmentation de la peau",
            "alleles": {
                "C/C": "Peau normale/foncée",
                "C/G": "Peau intermédiaire",
                "G/G": "Peau très claire"
            }
        }
    },

    "Métabolisme - Lactose": {
        "rs4988235": {
            "gene": "MCM6/LCT",
            "chr": "2",
            "trait": "Persistance de la lactase (tolérance au lactose)",
            "alleles": {
                "A/A": "INTOLERANT au lactose (ancestral)",
                "A/G": "Tolérance partielle",
                "G/G": "TOLERANT au lactose (mutation européenne)"
            }
        },
        "rs182549": {
            "gene": "MCM6/LCT",
            "chr": "2",
            "trait": "Persistance de la lactase (variant secondaire)",
            "alleles": {
                "C/C": "Pas de persistance",
                "C/T": "Persistance partielle",
                "T/T": "Persistance de la lactase"
            }
        }
    },

    "Métabolisme - Alcool": {
        "rs671": {
            "gene": "ALDH2",
            "chr": "12",
            "trait": "Métabolisme de l'alcool",
            "alleles": {
                "G/G": "Métabolisme normal de l'alcool",
                "G/A": "Métabolisme réduit (flush asiatique)",
                "A/A": "Métabolisme très réduit, intolérance sévère"
            }
        }
    },

    "Groupes sanguins": {
        "rs8176719": {
            "gene": "ABO",
            "chr": "9",
            "trait": "Groupe sanguin ABO",
            "alleles": {
                "deletion": "Groupe O (homozygote pour délétion)",
                "normal": "Groupe A ou B"
            }
        }
    },

    "Maladies génétiques - Hémochromatose": {
        "rs1800562": {
            "gene": "HFE",
            "chr": "6",
            "trait": "Hémochromatose (surcharge en fer)",
            "alleles": {
                "G/G": "Normal",
                "G/A": "Porteur mutation C282Y",
                "A/A": "RISQUE ELEVE hémochromatose"
            }
        },
        "rs1799945": {
            "gene": "HFE",
            "chr": "6",
            "trait": "Hémochromatose (mutation H63D)",
            "alleles": {
                "C/C": "Normal",
                "C/G": "Porteur mutation H63D",
                "G/G": "Risque modéré hémochromatose"
            }
        }
    },

    "Maladies génétiques - Mucoviscidose": {
        "rs113993960": {
            "gene": "CFTR",
            "chr": "7",
            "trait": "Mucoviscidose (mutation F508del)",
            "alleles": {
                "normal": "Normal",
                "carrier": "Porteur F508del",
                "affected": "Mucoviscidose (homozygote)"
            }
        }
    },

    "Maladies génétiques - Thrombophilie": {
        "rs6025": {
            "gene": "F5",
            "chr": "1",
            "trait": "Facteur V Leiden (risque thrombose)",
            "alleles": {
                "C/C": "Normal",
                "C/T": "Porteur Facteur V Leiden",
                "T/T": "RISQUE ELEVE de thrombose"
            }
        }
    },

    "Maladies génétiques - Drépanocytose/Thalassémie": {
        "rs334": {
            "gene": "HBB",
            "chr": "11",
            "trait": "Anémie falciforme / Résistance au paludisme",
            "alleles": {
                "A/A": "Normal",
                "A/T": "Porteur trait drépanocytaire (résistance paludisme)",
                "T/T": "Anémie falciforme (maladie)"
            }
        }
    },

    "Maladies neurologiques - Alzheimer": {
        "rs429358": {
            "gene": "APOE",
            "chr": "19",
            "trait": "Facteur de risque Alzheimer (APOE ε4)",
            "alleles": {
                "C/C": "Pas d'allèle ε4 (risque normal)",
                "C/T": "Un allèle ε4 (risque augmenté ~3x)",
                "T/T": "Deux allèles ε4 (RISQUE TRÈS ÉLEVÉ ~12x)"
            }
        },
        "rs7412": {
            "gene": "APOE",
            "chr": "19",
            "trait": "APOE ε2 (protecteur Alzheimer)",
            "alleles": {
                "C/C": "Pas d'allèle ε2",
                "C/T": "Un allèle ε2 (protecteur contre Alzheimer)",
                "T/T": "Deux allèles ε2 (forte protection)"
            }
        }
    },

    "Maladies cardiovasculaires": {
        "rs1333049": {
            "gene": "9p21.3",
            "chr": "9",
            "trait": "Risque de maladie coronarienne",
            "alleles": {
                "C/C": "Risque normal",
                "C/G": "Risque légèrement augmenté",
                "G/G": "Risque augmenté de maladie coronarienne"
            }
        },
        "rs1801133": {
            "gene": "MTHFR",
            "chr": "1",
            "trait": "Mutation MTHFR C677T (métabolisme folates)",
            "alleles": {
                "G/G": "Normal",
                "G/A": "Hétérozygote (activité MTHFR réduite)",
                "A/A": "Homozygote (activité réduite, risque cardiovasculaire)"
            }
        }
    },

    "Métabolisme - Obésité et diabète": {
        "rs9939609": {
            "gene": "FTO",
            "chr": "16",
            "trait": "Prédisposition à l'obésité",
            "alleles": {
                "T/T": "Risque normal",
                "A/T": "Risque légèrement augmenté d'obésité",
                "A/A": "Risque augmenté d'obésité (~3kg en moyenne)"
            }
        },
        "rs7903146": {
            "gene": "TCF7L2",
            "chr": "10",
            "trait": "Diabète de type 2",
            "alleles": {
                "C/C": "Risque normal",
                "C/T": "Risque augmenté de diabète type 2",
                "T/T": "RISQUE ÉLEVÉ de diabète type 2"
            }
        },
        "rs1801282": {
            "gene": "PPARG",
            "chr": "3",
            "trait": "Sensibilité à l'insuline",
            "alleles": {
                "C/C": "Normal",
                "C/G": "Meilleure sensibilité à l'insuline (protecteur)",
                "G/G": "Protection contre diabète type 2"
            }
        }
    },

    "Système immunitaire - Maladies auto-immunes": {
        "rs2476601": {
            "gene": "PTPN22",
            "chr": "1",
            "trait": "Maladies auto-immunes (PR, diabète type 1, lupus)",
            "alleles": {
                "G/G": "Risque normal",
                "G/A": "Risque augmenté de maladies auto-immunes",
                "A/A": "Risque fortement augmenté"
            }
        },
        "rs6897932": {
            "gene": "IL7R",
            "chr": "5",
            "trait": "Sclérose en plaques",
            "alleles": {
                "C/C": "Risque normal",
                "C/T": "Risque légèrement augmenté",
                "T/T": "Risque augmenté de sclérose en plaques"
            }
        }
    },

    "Système immunitaire - Maladie coeliaque": {
        "rs2187668": {
            "gene": "HLA-DQ",
            "chr": "6",
            "trait": "Maladie coeliaque (intolérance au gluten)",
            "alleles": {
                "C/C": "Risque augmenté de maladie coeliaque",
                "C/T": "Risque intermédiaire",
                "T/T": "Risque plus faible"
            }
        }
    },

    "Détoxification et métabolisme": {
        "rs1695": {
            "gene": "GSTP1",
            "chr": "11",
            "trait": "Détoxification (glutathion S-transférase)",
            "alleles": {
                "A/A": "Activité enzymatique normale",
                "A/G": "Activité réduite",
                "G/G": "Activité fortement réduite (détox moins efficace)"
            }
        },
        "rs762551": {
            "gene": "CYP1A2",
            "chr": "15",
            "trait": "Métabolisme de la caféine",
            "alleles": {
                "A/A": "Métaboliseur RAPIDE de caféine",
                "A/C": "Métaboliseur intermédiaire",
                "C/C": "Métaboliseur LENT de caféine"
            }
        }
    },

    "Traits sensoriels": {
        "rs713598": {
            "gene": "TAS2R38",
            "chr": "7",
            "trait": "Perception du goût amer (PTC/PROP)",
            "alleles": {
                "C/C": "Super-goûteur (très sensible à l'amertume)",
                "C/G": "Goûteur moyen",
                "G/G": "Non-goûteur (insensible à certains amers)"
            }
        },
        "rs1726866": {
            "gene": "OR7D4",
            "chr": "19",
            "trait": "Perception de l'odeur (androsténone)",
            "alleles": {
                "A/A": "Perception normale",
                "A/G": "Sensibilité variable",
                "G/G": "Sensibilité modifiée aux odeurs"
            }
        }
    },

    "Morphologie et traits physiques": {
        "rs11803731": {
            "gene": "TCHH",
            "chr": "1",
            "trait": "Texture des cheveux (raides/bouclés)",
            "alleles": {
                "A/A": "Cheveux raides",
                "A/T": "Cheveux ondulés",
                "T/T": "Cheveux bouclés"
            }
        },
        "rs3827760": {
            "gene": "EDAR",
            "chr": "2",
            "trait": "Épaisseur des cheveux",
            "alleles": {
                "A/A": "Cheveux fins",
                "A/G": "Épaisseur moyenne",
                "G/G": "Cheveux épais (fréquent en Asie)"
            }
        },
        "rs6060369": {
            "gene": "GDF5",
            "chr": "20",
            "trait": "Taille adulte",
            "alleles": {
                "A/A": "Tendance à être plus petit",
                "A/G": "Taille moyenne",
                "G/G": "Tendance à être plus grand"
            }
        }
    },

    "Santé osseuse et vitamine D": {
        "rs2228570": {
            "gene": "VDR",
            "chr": "12",
            "trait": "Récepteur de la vitamine D",
            "alleles": {
                "C/C": "Activité VDR normale",
                "C/T": "Activité intermédiaire",
                "T/T": "Activité modifiée (besoin vitamine D variable)"
            }
        },
        "rs1800247": {
            "gene": "COL1A1",
            "chr": "17",
            "trait": "Risque d'ostéoporose",
            "alleles": {
                "G/G": "Risque normal",
                "G/T": "Risque légèrement augmenté",
                "T/T": "Risque augmenté d'ostéoporose"
            }
        }
    },

    "Performance et endurance": {
        "rs1815739": {
            "gene": "ACTN3",
            "chr": "11",
            "trait": "Performance musculaire (sprint vs endurance)",
            "alleles": {
                "C/C": "Muscle rapide fonctionnel (SPRINT)",
                "C/T": "Intermédiaire (mixte)",
                "T/T": "Perte protéine ACTN3 (ENDURANCE)"
            }
        },
        "rs8192678": {
            "gene": "PPARGC1A",
            "chr": "4",
            "trait": "Endurance aérobie",
            "alleles": {
                "C/C": "Réponse normale à l'entraînement",
                "C/T": "Bonne réponse à l'entraînement",
                "T/T": "Excellente réponse à l'entraînement d'endurance"
            }
        }
    },

    "Traits divers": {
        "rs17822931": {
            "gene": "ABCC11",
            "chr": "16",
            "trait": "Type de cérumen et odeur corporelle",
            "alleles": {
                "C/C": "Cérumen humide, odeur corporelle normale",
                "C/T": "Intermédiaire",
                "T/T": "Cérumen sec, odeur corporelle réduite (fréquent en Asie)"
            }
        },
        "rs4680": {
            "gene": "COMT",
            "chr": "22",
            "trait": "Tolérance à la douleur et stress",
            "alleles": {
                "G/G": "Val/Val - Métabolisme dopamine rapide (guerrier)",
                "G/A": "Val/Met - Intermédiaire",
                "A/A": "Met/Met - Métabolisme lent (inquiet, sensible douleur)"
            }
        },
        "rs1800497": {
            "gene": "ANKK1/DRD2",
            "chr": "11",
            "trait": "Système de récompense dopaminergique",
            "alleles": {
                "G/G": "Récepteurs dopamine normaux",
                "G/A": "Récepteurs légèrement réduits",
                "A/A": "Récepteurs D2 réduits (risque addiction)"
            }
        },
        "rs1801260": {
            "gene": "CLOCK",
            "chr": "4",
            "trait": "Chronotype (personne du matin/soir)",
            "alleles": {
                "T/T": "Tendance personne du matin",
                "T/C": "Intermédiaire",
                "C/C": "Tendance personne du soir"
            }
        }
    }
}


def parse_vcf_header(vcf_file):
    """Extrait les noms des échantillons du VCF"""
    with gzip.open(vcf_file, 'rt') as f:
        for line in f:
            if line.startswith('#CHROM'):
                fields = line.strip().split('\t')
                # Les échantillons commencent après la colonne FORMAT
                samples = fields[9:]
                return samples
    return []


def extract_genotypes(vcf_file, snp_dict):
    """
    Extrait les génotypes pour tous les SNPs d'intérêt
    Retourne un dictionnaire: {rs_id: {sample: genotype}}
    """
    results = defaultdict(dict)
    snps_to_find = set()

    # Créer un index rs_id -> SNP info pour recherche rapide
    snp_index = {}
    for category, snps in snp_dict.items():
        for rs_id, info in snps.items():
            snps_to_find.add(rs_id)
            snp_index[rs_id] = info

    samples = parse_vcf_header(vcf_file)
    print(f"Échantillons trouvés: {', '.join(samples)}\n")

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

            if rs_id in snps_to_find:
                # Extraire les génotypes pour chaque échantillon
                format_field = fields[8]
                sample_data = fields[9:]

                for i, sample in enumerate(samples):
                    gt = sample_data[i].split(':')[0]  # Extraire le génotype (GT)

                    # Convertir le génotype en format allèle
                    if gt == './.':
                        genotype = "Pas de données"
                    elif gt == '0/0':
                        genotype = f"{ref}/{ref}"
                    elif gt == '0/1' or gt == '1/0':
                        genotype = f"{ref}/{alt}"
                    elif gt == '1/1':
                        genotype = f"{alt}/{alt}"
                    else:
                        genotype = gt

                    results[rs_id][sample] = {
                        'genotype': genotype,
                        'ref': ref,
                        'alt': alt,
                        'chr': chrom,
                        'pos': pos
                    }

    return results, samples


def generate_report(results, samples, snp_dict):
    """Génère un rapport formaté des résultats phénotypiques"""

    print("=" * 80)
    print("RAPPORT D'ANALYSE PHÉNOTYPIQUE ET PATHOLOGIQUE")
    print("=" * 80)
    print()

    for category, snps in snp_dict.items():
        print(f"\n{'=' * 80}")
        print(f"CATÉGORIE: {category.upper()}")
        print(f"{'=' * 80}\n")

        for rs_id, snp_info in snps.items():
            if rs_id not in results:
                print(f"⚠ {rs_id} ({snp_info.get('gene', 'N/A')}) - NON TROUVÉ dans le VCF")
                print(f"  Trait: {snp_info['trait']}\n")
                continue

            print(f"SNP: {rs_id} - Gène: {snp_info.get('gene', 'N/A')}")
            print(f"Trait: {snp_info['trait']}")
            print(f"{'-' * 80}")

            for sample in samples:
                if sample in results[rs_id]:
                    data = results[rs_id][sample]
                    genotype = data['genotype']

                    print(f"\n{sample}:")
                    print(f"  Position: chr{data['chr']}:{data['pos']}")
                    print(f"  Génotype: {genotype}")

                    # Interpréter le génotype
                    if genotype == "Pas de données":
                        print(f"  Interprétation: Données manquantes")
                    elif 'alleles' in snp_info:
                        interpretation = snp_info['alleles'].get(genotype, "Interprétation non disponible")
                        print(f"  Interprétation: {interpretation}")

                        # Ajouter des alertes pour les variants pathogènes
                        if any(word in interpretation.upper() for word in ['RISQUE', 'ELEVE', 'INTOLÉRANT']):
                            print(f"  ⚠️  ATTENTION: Variant potentiellement important")

            print()

    # Résumé final
    print("\n" + "=" * 80)
    print("RÉSUMÉ DES DÉCOUVERTES IMPORTANTES")
    print("=" * 80)

    important_findings = []
    for category, snps in snp_dict.items():
        for rs_id, snp_info in snps.items():
            if rs_id in results:
                for sample in samples:
                    if sample in results[rs_id]:
                        data = results[rs_id][sample]
                        genotype = data['genotype']

                        if 'alleles' in snp_info and genotype in snp_info['alleles']:
                            interp = snp_info['alleles'][genotype]
                            if any(word in interp.upper() for word in ['RISQUE ELEVE', 'INTOLÉRANT', 'PROBABLE']):
                                important_findings.append({
                                    'sample': sample,
                                    'trait': snp_info['trait'],
                                    'gene': snp_info.get('gene', 'N/A'),
                                    'genotype': genotype,
                                    'interpretation': interp
                                })

    if important_findings:
        for finding in important_findings:
            print(f"\n{finding['sample']} - {finding['trait']} ({finding['gene']})")
            print(f"  {finding['genotype']}: {finding['interpretation']}")
    else:
        print("\nAucun variant pathogène majeur détecté.")


def main():
    if len(sys.argv) < 2:
        print("Usage: python phenotype_analysis.py <fichier.vcf.gz>")
        sys.exit(1)

    vcf_file = sys.argv[1]

    print("Extraction des génotypes en cours...\n")
    results, samples = extract_genotypes(vcf_file, PHENOTYPE_SNPS)

    # Rediriger la sortie vers un fichier
    import io
    from contextlib import redirect_stdout

    output_file = "phenotype_report.txt"

    # Sauvegarder dans un fichier ET afficher à l'écran
    with open(output_file, 'w', encoding='utf-8') as f:
        # Créer un objet qui écrit à la fois dans le fichier et stdout
        class Tee:
            def __init__(self, *files):
                self.files = files
            def write(self, data):
                for f in self.files:
                    f.write(data)
            def flush(self):
                for f in self.files:
                    f.flush()

        import sys as system
        tee = Tee(system.stdout, f)
        with redirect_stdout(tee):
            generate_report(results, samples, PHENOTYPE_SNPS)

    print(f"\n\n{'='*80}")
    print(f"Rapport sauvegardé dans: {output_file}")
    print(f"{'='*80}")


if __name__ == "__main__":
    main()
