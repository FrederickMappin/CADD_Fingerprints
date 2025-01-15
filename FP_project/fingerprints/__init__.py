from .base_fingerprint import Fingerprint
from .pubchem_fingerprint import PubchemFingerprint
from .other_fingerprint import OtherFingerprint
from .morgan_fingerprint import MorganFingerprint  # Add this line

__all__ = ['Fingerprint', 'PubchemFingerprint', 'OtherFingerprint', 'MorganFingerprint']