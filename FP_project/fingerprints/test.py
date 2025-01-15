from morgan_fingerprint import MorganFingerprint

# Create an instance of the class
fingerprinter = MorganFingerprint()

# Calculate fingerprints
result = fingerprinter.calculate(
    input_file='compounds.csv',
    output_file='fingerprints.csv'
)

# The result DataFrame contains the original data plus fingerprint columns
print(result.head())