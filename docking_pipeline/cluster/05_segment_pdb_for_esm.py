from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
INPUTS = ["9WZU_apo_H.pdb", "apo_8WL6_H.pdb", "caspo_T2_H.pdb"]
CHAIN_IDS = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789"
MAX_RESIDUES_PER_CHAIN = 900

for filename in INPUTS:
    source = ROOT / "receptors" / filename
    target = source.with_name(source.stem.replace("_H", "_diffdock") + ".pdb")
    residue_order = []
    residue_index = {}
    output = []
    for line in source.read_text(encoding="ascii").splitlines(keepends=True):
        if line.startswith(("ATOM  ", "HETATM")):
            key = (line[21], line[22:26], line[26])
            if key not in residue_index:
                residue_index[key] = len(residue_order)
                residue_order.append(key)
            segment = residue_index[key] // MAX_RESIDUES_PER_CHAIN
            line = line[:21] + CHAIN_IDS[segment] + line[22:]
        output.append(line)
    target.write_text("".join(output), encoding="ascii", newline="\n")
    counts = [min(MAX_RESIDUES_PER_CHAIN, len(residue_order) - i)
              for i in range(0, len(residue_order), MAX_RESIDUES_PER_CHAIN)]
    print(f"{target.name}: {len(residue_order)} residues split as {counts}")

csv_source = ROOT / "cluster" / "diffdock_input.csv"
csv_target = ROOT / "cluster" / "diffdock_input_segmented.csv"
csv_text = csv_source.read_text(encoding="utf-8")
for filename in INPUTS:
    csv_text = csv_text.replace(filename, filename.replace("_H", "_diffdock"))
csv_target.write_text(csv_text, encoding="utf-8", newline="\n")
print(f"Wrote {csv_target.name}")
