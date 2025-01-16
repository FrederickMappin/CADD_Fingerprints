import sys
from fingerprints.pubchem_fingerprint import PubchemFingerprint
from fingerprints.morgan_fingerprint import MorganFingerprint
from fingerprints.maccs_fingerprint import MACCSFingerprint
from fingerprints.atom_pairs_fingerprint import AtomPairsFingerprint
from fingerprints.avalon_fingerprint import AvalonFingerprint
from fingerprints.layered_fingerprint import LayeredFingerprint
from fingerprints.pattern_fingerprint import PatternFingerprint
from fingerprints.rdkit_fingerprint import RDKitFingerprint
from fingerprints.topological_torsions_fingerprint import TopologicalTorsionsFingerprint

def main():
    if len(sys.argv) != 4:
        print("Usage: python main.py <fingerprint_type> <input_file> <output_file>")
        sys.exit(1)
    
    fingerprint_type = sys.argv[1]
    input_file = sys.argv[2]
    output_file = sys.argv[3]
    
    if fingerprint_type == "pubchem":
        fp = PubchemFingerprint()
    elif fingerprint_type == "morgan":
        fp = MorganFingerprint()
    elif fingerprint_type == "maccs":
        fp = MACCSFingerprint()
    elif fingerprint_type == "atom_pairs":
        fp = AtomPairsFingerprint()
    elif fingerprint_type == "avalon":
        fp = AvalonFingerprint()
    elif fingerprint_type == "layered":
        fp = LayeredFingerprint()
    elif fingerprint_type == "pattern":
        fp = PatternFingerprint()
    elif fingerprint_type == "rdkit":
        fp = RDKitFingerprint()
    elif fingerprint_type == "topological_torsions":
        fp = TopologicalTorsionsFingerprint()
    else:
        print(f"Unknown fingerprint type: {fingerprint_type}")
        sys.exit(1)
    
    fp.calculate(input_file, output_file)
    print(f"Fingerprints calculated and saved to {output_file}")

if __name__ == "__main__":
    main()