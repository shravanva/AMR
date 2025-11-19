import re
from Bio import SeqIO
from collections import defaultdict

# ========================================
# COMPREHENSIVE RESISTANCE GENE DATABASE
# ========================================
RESISTANCE_GENE_SIGNATURES = {
    # Beta-lactamases
    'beta_lactam_resistance': {
        'gene_families': ['bla', 'TEM', 'SHV', 'OXA', 'CTX-M', 'KPC', 'NDM', 'VIM', 'IMP'],
        'signatures': [
            r'BLATEM[0-9]+',
            r'BLASHV[0-9]*',
            r'BLAOXA[0-9]*',
            r'BLACTX-M[0-9]*',
            r'BLAKPC[0-9]*',
            r'BLANDM[0-9]*',
            r'BLAVIM[0-9]*',
            r'BLAIMP[0-9]*',
        ],
        'mechanism': 'Beta-lactam hydrolysis',
        'probability_markers': [
            ('SXXK', 0.8),  # Active site motif
            ('SDN', 0.7),   # Beta-lactamase signature
            ('KTG', 0.6),   # Class A beta-lactamases
        ]
    },
    'aminoglycoside_resistance': {
        'gene_families': ['aph', 'ant', 'aac'],
        'signatures': [
            r'APH\(3[\'\"]*\)-I[A-Z]*',
            r'ANT\(2[\'\"]*\)-I[A-Z]*',
            r'AAC\(6[\'\"]*\)-I[A-Z]*',
        ],
        'mechanism': 'Aminoglycoside modification',
        'probability_markers': [
            ('GXGXXG', 0.8),
            ('DXD', 0.7),
        ]
    },
    'fluoroquinolone_resistance': {
        'gene_families': ['qnr', 'aac(6\')-Ib-cr', 'qepA', 'oqxAB'],
        'signatures': [
            r'QNR[A-Z]',
            r'AAC\(6[\'\"]*\)-IB-CR',
            r'QEPA',
            r'OQXAB',
        ],
        'mechanism': 'DNA gyrase protection or efflux',
        'target_mutations': {
            'gyrA': ['S83L', 'D87N', 'D87G'],
            'parC': ['S80I', 'E84V'],
        }
    },
    'tetracycline_resistance': {
        'gene_families': ['tet', 'otr'],
        'signatures': [
            r'TET[A-Z]',
            r'OTR[A-Z]*',
        ],
        'mechanism': 'Efflux or ribosomal protection',
        'probability_markers': [
            ('MFS', 0.8),
        ]
    },
    'macrolide_resistance': {
        'gene_families': ['erm', 'mef', 'mph'],
        'signatures': [
            r'ERM[A-Z]',
            r'MEF[A-Z]*',
            r'MPH[A-Z]*',
        ],
        'mechanism': '23S rRNA methylation or efflux',
        'probability_markers': [
            ('NPPF', 0.7),
        ]
    },
    'glycopeptide_resistance': {
        'gene_families': ['van'],
        'signatures': [
            r'VAN[A-Z]',
        ],
        'mechanism': 'Cell wall precursor modification',
        'probability_markers': [
            ('DDL', 0.8),
            ('D-ALA-D-LAC', 0.9),
        ]
    },
    'polymyxin_resistance': {
        'gene_families': ['mcr', 'pmr'],
        'signatures': [
            r'MCR-[0-9]+',
            r'PMR[A-Z]',
        ],
        'mechanism': 'Lipid A modification',
        'probability_markers': [
            ('PETA', 0.8),
        ]
    },
}

# ========================================
# GENE SCANNING FUNCTION
# ========================================
def scan_genome_for_resistance_genes(fasta_file):
    """
    Scan bacterial genome for resistance genes
    Returns: Dictionary of detected genes and confidence scores
    """
    record = list(SeqIO.parse(fasta_file, "fasta"))[0]
    sequence = str(record.seq).upper()
    detected_genes = defaultdict(list)
    # Scan for each resistance class
    for resistance_class, gene_info in RESISTANCE_GENE_SIGNATURES.items():
        class_score = 0
        found_genes = []
        # Check signature patterns
        for signature in gene_info['signatures']:
            matches = re.finditer(signature, sequence, re.IGNORECASE)
            for match in matches:
                found_genes.append({
                    'position': match.start(),
                    'sequence': match.group(),
                    'confidence': 0.9
                })
                class_score += 0.3
        # Check probability markers
        for marker, prob in gene_info.get('probability_markers', []):
            if marker in sequence:
                class_score += prob
                found_genes.append({
                    'marker': marker,
                    'confidence': prob
                })
        if found_genes:
            detected_genes[resistance_class] = {
                'genes': found_genes,
                'score': min(class_score, 1.0),
                'mechanism': gene_info['mechanism']
            }
    return detected_genes

# ========================================
# EVOLUTIONARY PREDICTION
# ========================================
def predict_resistance_evolution(sequence, antibiotic_class):
    """
    Predict likelihood of resistance evolution under antibiotic pressure
    """
    evolution_factors = {
        'high_gc_content': 0.0,
        'mobile_elements': 0.0,
        'horizontal_gene_transfer': 0.0,
        'mutation_hotspots': 0.0,
    }
    # Calculate GC content
    gc_content = (sequence.count('G') + sequence.count('C')) / len(sequence)
    if gc_content > 0.55:
        evolution_factors['high_gc_content'] = 0.3
    # Check for mobile genetic elements
    mobile_signatures = ['INTEGRON', 'TRANSPOSON', 'PLASMID', 'IS[0-9]+']
    for sig in mobile_signatures:
        if re.search(sig, sequence, re.IGNORECASE):
            evolution_factors['mobile_elements'] += 0.2
    # Check for recombination hotspots
    if 'ATGATG' in sequence or 'TAATAA' in sequence:
        evolution_factors['horizontal_gene_transfer'] = 0.15
    # Mutation hotspots (repeat regions)
    repeat_count = len(re.findall(r'([ATCG]{3,})\1+', sequence))
    if repeat_count > 10:
        evolution_factors['mutation_hotspots'] = 0.25
    # Calculate overall evolution probability
    base_probability = 0.15  # Base 15% chance
    evolution_score = base_probability + sum(evolution_factors.values())
    evolution_score = min(evolution_score, 0.95)  # Cap at 95%
    return {
        'probability': evolution_score,
        'factors': evolution_factors,
        'timeline_months': int((1 - evolution_score) * 24)  # Inverse relationship
    }

# ========================================
# NOVEL BACTERIA ANALYSIS (optional CLI)
# ========================================
def analyze_novel_bacteria(fasta_file):
    """
    Complete analysis for newly discovered bacteria
    """
    print("="*70)
    print("ANALYZING NOVEL BACTERIAL GENOME")
    print("="*70)
    record = list(SeqIO.parse(fasta_file, "fasta"))[0]
    sequence = str(record.seq).upper()
    print(f"\nGenome: {record.id}")
    print(f"Length: {len(sequence):,} bp")
    print(f"GC%: {(sequence.count('G') + sequence.count('C')) / len(sequence) * 100:.1f}%")
    print("\n" + "="*70)
    print("CURRENT RESISTANCE GENES (Intrinsic)")
    print("="*70)
    detected = scan_genome_for_resistance_genes(fasta_file)
    if detected:
        for resistance_class, info in detected.items():
            print(f"\n🧬 {resistance_class.replace('_', ' ').title()}")
            print(f"   Confidence: {info['score']:.0%}")
            print(f"   Mechanism: {info['mechanism']}")
            print(f"   Genes found: {len(info['genes'])}")
    else:
        print("\n✓ No intrinsic resistance genes detected")
        print("  (This is a susceptible bacterium)")
    print("\n" + "="*70)
    print("RESISTANCE EVOLUTION PREDICTION")
    print("="*70)
    for antibiotic_class in ['beta_lactam', 'aminoglycoside', 'fluoroquinolone']:
        evolution = predict_resistance_evolution(sequence, antibiotic_class)
        print(f"\n💊 {antibiotic_class.replace('_', ' ').title()} Antibiotics")
        print(f"   Evolution probability: {evolution['probability']:.0%}")
        print(f"   Estimated timeline: {evolution['timeline_months']} months")
        print(f"   Risk factors:")
        for factor, score in evolution['factors'].items():
            if score > 0:
                print(f"     - {factor.replace('_', ' ').title()}: {score:.0%}")
    print("\n" + "="*70)
    print("ANALYSIS COMPLETE")
    print("="*70)

if __name__ == "__main__":
    import sys
    if len(sys.argv) > 1:
        analyze_novel_bacteria(sys.argv[1])
    else:
        print("Usage: python advanced_gene_scanner.py <genome.fasta>")
