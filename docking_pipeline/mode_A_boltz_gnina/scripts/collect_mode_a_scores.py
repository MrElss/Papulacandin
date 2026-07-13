#!/usr/bin/env python3
"""Collect isolated Mode A GNINA scores and pose-displacement RMSD."""
from __future__ import annotations

import csv
import json
import math
import statistics
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
FIELDS = ("CNNscore", "CNNaffinity", "minimizedAffinity", "minimizedRMSD_A")


def first_record(path: Path) -> str:
    return path.read_text(encoding="utf-8", errors="replace").split("$$$$", 1)[0]


def sdf_tags(record: str) -> dict[str, str]:
    lines = record.splitlines()
    values: dict[str, str] = {}
    for i, line in enumerate(lines[:-1]):
        if line.startswith(">") and "<" in line and ">" in line:
            key = line.split("<", 1)[1].split(">", 1)[0]
            values[key] = lines[i + 1].strip()
    return values


def coordinates(record: str) -> list[tuple[float, float, float]]:
    lines = record.splitlines()
    if len(lines) < 4:
        return []
    try:
        atom_count = int(lines[3][0:3])
    except ValueError:
        return []
    coords = []
    for line in lines[4 : 4 + atom_count]:
        try:
            coords.append((float(line[0:10]), float(line[10:20]), float(line[20:30])))
        except ValueError:
            return []
    return coords


def direct_rmsd(before: list[tuple[float, float, float]], after: list[tuple[float, float, float]]) -> str:
    if not before or len(before) != len(after):
        return ""
    squared = sum(
        (a[0] - b[0]) ** 2 + (a[1] - b[1]) ** 2 + (a[2] - b[2]) ** 2
        for a, b in zip(before, after)
    )
    return f"{math.sqrt(squared / len(before)):.4f}"


rows = []
for output in sorted((ROOT / "out_gnina").glob("*_min.sdf")):
    tag = output.name.removesuffix("_min.sdf")
    template, source = tag.split("__", 1)
    input_pose = ROOT / "cofold_poses" / f"{tag}.sdf"
    out_record = first_record(output)
    tags = sdf_tags(out_record)
    in_coords = coordinates(first_record(input_pose)) if input_pose.exists() else []
    out_coords = coordinates(out_record)
    rows.append(
        {
            "template": template,
            "source": source,
            "CNNscore": tags.get("CNNscore", ""),
            "CNNaffinity": tags.get("CNNaffinity", ""),
            "minimizedAffinity": tags.get("minimizedAffinity", ""),
            "minimizedRMSD_A": tags.get("minimizedRMSD", ""),
            "direct_coordinate_rmsd_A": direct_rmsd(in_coords, out_coords),
            "input_pose": str(input_pose.relative_to(ROOT)),
            "minimized_pose": str(output.relative_to(ROOT)),
        }
    )

destination = ROOT / "mode_A_scores.csv"
with destination.open("w", newline="", encoding="utf-8-sig") as handle:
    writer = csv.DictWriter(handle, fieldnames=list(rows[0]) if rows else ["template", "source", *FIELDS])
    writer.writeheader()
    writer.writerows(rows)

print(f"wrote {destination} ({len(rows)} rows)")


def write_group_summary(group_field: str, output_name: str) -> None:
    numeric = ("CNNscore", "CNNaffinity", "minimizedAffinity", "minimizedRMSD_A")
    summary_rows = []
    for group in sorted({row[group_field] for row in rows}):
        members = [row for row in rows if row[group_field] == group]
        summary = {group_field: group, "n": len(members)}
        for field in numeric:
            values = [float(row[field]) for row in members]
            summary[f"{field}_mean"] = f"{statistics.mean(values):.6f}"
            summary[f"{field}_min"] = f"{min(values):.6f}"
            summary[f"{field}_max"] = f"{max(values):.6f}"
        summary_rows.append(summary)
    output = ROOT / output_name
    with output.open("w", newline="", encoding="utf-8-sig") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(summary_rows[0]))
        writer.writeheader()
        writer.writerows(summary_rows)
    print(f"wrote {output} ({len(summary_rows)} rows)")


if rows:
    write_group_summary("template", "mode_A_scores_by_template.csv")
    write_group_summary("source", "mode_A_scores_by_sample.csv")

metrics_path = ROOT / "predictions" / "metrics.json"
if metrics_path.exists():
    payload = json.loads(metrics_path.read_text(encoding="utf-8"))
    best = payload.get("best_sample", {}).get("metrics", {})
    boltz_rows = []
    for index, result in enumerate(payload.get("all_sample_results", [])):
        metrics = result.get("metrics", {})
        boltz_rows.append(
            {
                "sample_index_by_archive_order": index,
                "structure_file": f"predictions/sample_{index}_predicted_structure.cif",
                "is_best_metrics_match": metrics == best,
                **metrics,
            }
        )
    if boltz_rows:
        output = ROOT / "boltz_metrics_summary.csv"
        with output.open("w", newline="", encoding="utf-8-sig") as handle:
            writer = csv.DictWriter(handle, fieldnames=list(boltz_rows[0]))
            writer.writeheader()
            writer.writerows(boltz_rows)
        print(f"wrote {output} ({len(boltz_rows)} rows)")
