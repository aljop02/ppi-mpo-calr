from pathlib import Path

ALN = Path("calr_full.aln")      # MAFFT FASTA alignment of full mouse vs full human
INP = Path("mouse_actives.txt")  # one mouse residue number per line
OUT = Path("mouse_to_human_mapping.tsv")

def read_fasta_alignment(path):
    seqs = {}
    name = None
    parts = []
    for line in path.read_text().splitlines():
        if not line.strip():  # skip blanks
            continue
        if line.startswith(">"):
            if name is not None:
                seqs[name] = "".join(parts)
            name = line[1:].strip()
            parts = []
        else:
            parts.append(line.strip().replace(" ", ""))
    if name is not None:
        seqs[name] = "".join(parts)
    return seqs

def pick_mouse_human_keys(seqs):
    keys = list(seqs.keys())
    mk = next((k for k in keys if "MOUSE" in k.upper() or "MUS" in k.upper()), None)
    hk = next((k for k in keys if "HUMAN" in k.upper() or "HOMO" in k.upper()), None)
    if mk and hk: return mk, hk
    # fallback: assume order is mouse, human
    return keys[0], keys[1]

def build_map(aln_mouse, aln_human):
    mpos = hpos = 0
    mapping = []  # per column: (mouse_pos or None, human_pos or None)
    for m, h in zip(aln_mouse, aln_human):
        mm = hh = None
        if m != "-":
            mpos += 1; mm = mpos
        if h != "-":
            hpos += 1; hh = hpos
        mapping.append((mm, hh))
    return mapping

def map_mouse_residues(mapping, mouse_list):
    out = []
    for m_req in mouse_list:
        cand = [h for (m,h) in mapping if m == m_req and h is not None]
        h_map = cand[0] if cand else None
        out.append((m_req, h_map))
    return out

def expand_pm3(nums):
    """broaden ±3 for AIRs, unique+sorted."""
    s = set()
    for n in nums:
        if n is None: continue
        for k in range(n-3, n+4):
            if k > 0: s.add(k)
    return sorted(s)

def main():
    seqs = read_fasta_alignment(ALN)
    mk, hk = pick_mouse_human_keys(seqs)
    mapping = build_map(seqs[mk], seqs[hk])

    mouse_list = [int(x.strip()) for x in INP.read_text().splitlines() if x.strip()]
    pairs = map_mouse_residues(mapping, mouse_list)

    with OUT.open("w") as f:
        f.write("mouse_res\thuman_res\n")
        for m,h in pairs:
            f.write(f"{m}\t{'' if h is None else h}\n")

    human_hits = [h for _,h in pairs if h is not None]
    broadened = expand_pm3(human_hits)

    print(f"Mapped {len(human_hits)}/{len(pairs)} residues.")
    print("Mapping written to:", OUT)
    print("Human residues (for CALR chain E) :", human_hits)
    print("Broadened CALR list (±3, for AIRs):", broadened)

if __name__ == "__main__":
    main()
