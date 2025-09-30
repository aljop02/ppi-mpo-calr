from Bio import SeqIO

# Read mouse CALR fasta
record = SeqIO.read("mouse_calr.fasta", "fasta")

# Extract segments (Python is 0-based; slice end is exclusive)
seg1 = record.seq[17:238]   # residues 18–238
seg2 = record.seq[272:368]  # residues 273–368

# Save to new FASTA
with open("mouse_calr_fragments.fasta", "w") as f:
    f.write(f">mouse_CALR_18_238\n{seg1}\n")
    f.write(f">mouse_CALR_273_368\n{seg2}\n")
