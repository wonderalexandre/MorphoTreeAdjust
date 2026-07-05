#!/usr/bin/env python3

import argparse
import csv
import datetime as dt
import re
import statistics
import subprocess
from pathlib import Path


DEFAULT_THRESHOLDS = "10,20,30,40,49,59,69,79,89,99"
DEFAULT_DATA_DIR = Path("tests/dat/icdar_480p")


SUMMARY_PATTERNS = {
    "dynamic_subtree": re.compile(
        r"dynamic-subtree init=(?P<init>[0-9.]+) ms "
        r"iter=(?P<iter>[0-9.]+) ms "
        r"reconstruct=(?P<reconstruct>[0-9.]+) ms "
        r"total=(?P<total>[0-9.]+) ms"
    ),
    "dynamic_leaf": re.compile(
        r"dynamic-leaf init=(?P<init>[0-9.]+) ms "
        r"iter=(?P<iter>[0-9.]+) ms "
        r"reconstruct=(?P<reconstruct>[0-9.]+) ms "
        r"total=(?P<total>[0-9.]+) ms"
    ),
    "naive": re.compile(r"naive total=(?P<total>[0-9.]+) ms"),
}


MATCH_PATTERNS = {
    "dynamic_subtree_matches_naive": re.compile(r"dynamic-subtree == naive: (?P<value>true|false)"),
    "dynamic_leaf_matches_naive": re.compile(r"dynamic-leaf == naive: (?P<value>true|false)"),
    "dynamic_subtree_matches_dynamic_leaf": re.compile(
        r"dynamic-subtree == dynamic-leaf: (?P<value>true|false)"
    ),
}


IMAGE_SIZE_PATTERN = re.compile(r"Image size: (?P<rows>[0-9]+)x(?P<cols>[0-9]+)")


def repo_root() -> Path:
    return Path(__file__).resolve().parents[2]


def parse_thresholds(value: str) -> list[str]:
    thresholds = [part.strip() for part in value.split(",") if part.strip()]
    if not thresholds:
        raise ValueError("At least one threshold is required.")
    for threshold in thresholds:
        int(threshold)
    return thresholds


def default_output_path(root: Path, image_limit: int, num_thresholds: int, repeats: int) -> Path:
    stamp = dt.datetime.now().strftime("%Y%m%d")
    name = (
        f"increasing_casf_icdar480p_{image_limit}img_"
        f"{num_thresholds}t_{repeats}reps_current_{stamp}.csv"
    )
    return root / "results" / "jmiv2026" / name


def collect_images(data_dir: Path, image_limit: int) -> list[Path]:
    images = sorted(data_dir.glob("*.png"))
    if image_limit > 0:
        images = images[:image_limit]
    if not images:
        raise RuntimeError(f"No PNG images found in {data_dir}")
    return images


def parse_bool(value: str) -> bool:
    if value == "true":
        return True
    if value == "false":
        return False
    raise ValueError(f"Unexpected boolean value: {value}")


def find_line_match(pattern: re.Pattern[str], output: str, label: str) -> re.Match[str]:
    match = pattern.search(output)
    if match is None:
        raise RuntimeError(f"Could not parse {label} from dynamic_casf_apply output:\n{output}")
    return match


def parse_compare_output(output: str) -> dict[str, object]:
    row: dict[str, object] = {}

    size_match = find_line_match(IMAGE_SIZE_PATTERN, output, "image size")
    row["rows"] = int(size_match.group("rows"))
    row["cols"] = int(size_match.group("cols"))

    for prefix, pattern in SUMMARY_PATTERNS.items():
        match = find_line_match(pattern, output, prefix)
        if prefix == "naive":
            row["naive_total_ms"] = float(match.group("total"))
            continue
        row[f"{prefix}_init_ms"] = float(match.group("init"))
        row[f"{prefix}_iter_ms"] = float(match.group("iter"))
        row[f"{prefix}_reconstruct_ms"] = float(match.group("reconstruct"))
        row[f"{prefix}_total_ms"] = float(match.group("total"))

    for field, pattern in MATCH_PATTERNS.items():
        match = find_line_match(pattern, output, field)
        row[field] = parse_bool(match.group("value"))

    row["dynamic_subtree_speedup_vs_naive"] = (
        row["naive_total_ms"] / row["dynamic_subtree_total_ms"]
        if row["dynamic_subtree_total_ms"] > 0.0
        else 0.0
    )
    row["dynamic_leaf_speedup_vs_naive"] = (
        row["naive_total_ms"] / row["dynamic_leaf_total_ms"]
        if row["dynamic_leaf_total_ms"] > 0.0
        else 0.0
    )
    return row


def run_compare(binary: Path,
                image: Path,
                attribute: str,
                radio_adj: float,
                thresholds: list[str]) -> tuple[dict[str, object], str]:
    command = [
        str(binary),
        "--mode",
        "compare",
        "--attribute",
        attribute,
        "--radio-adj",
        str(radio_adj),
        "--no-output",
        str(image),
        *thresholds,
    ]
    completed = subprocess.run(command, text=True, capture_output=True, check=False)
    output = completed.stdout + completed.stderr
    if completed.returncode != 0:
        raise RuntimeError(
            "dynamic_casf_apply failed with exit code "
            f"{completed.returncode}\nCommand: {' '.join(command)}\n{output}"
        )

    row = parse_compare_output(output)
    row["command"] = " ".join(command)
    return row, output


def write_summary(csv_path: Path, rows: list[dict[str, object]]) -> Path:
    summary_path = csv_path.with_suffix(".summary.txt")
    metrics = [
        "dynamic_subtree_total_ms",
        "dynamic_subtree_iter_ms",
        "dynamic_leaf_total_ms",
        "dynamic_leaf_iter_ms",
        "naive_total_ms",
        "dynamic_subtree_speedup_vs_naive",
        "dynamic_leaf_speedup_vs_naive",
    ]

    lines = [
        f"rows={len(rows)}",
        f"images={len({row['image'] for row in rows})}",
        f"repeats={len({row['repeat'] for row in rows})}",
        f"thresholds={rows[0]['thresholds'] if rows else ''}",
        "",
    ]

    for metric in metrics:
        values = [float(row[metric]) for row in rows]
        lines.append(
            f"{metric}: sum={sum(values):.3f} "
            f"mean={statistics.fmean(values):.3f} "
            f"median={statistics.median(values):.3f}"
        )

    all_match = all(
        bool(row["dynamic_subtree_matches_naive"])
        and bool(row["dynamic_leaf_matches_naive"])
        and bool(row["dynamic_subtree_matches_dynamic_leaf"])
        for row in rows
    )
    lines.extend(["", f"all_outputs_match={str(all_match).lower()}"])
    summary_path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return summary_path


def parse_args() -> argparse.Namespace:
    root = repo_root()
    parser = argparse.ArgumentParser(
        description=(
            "Run a small increasing-CASF benchmark on ICDAR 480p images using "
            "dynamic_casf_apply compare mode."
        )
    )
    parser.add_argument("--binary", type=Path, default=root / "build-dev-tools" / "tools" / "dynamic_casf_apply")
    parser.add_argument("--data-dir", type=Path, default=DEFAULT_DATA_DIR)
    parser.add_argument("--output", type=Path, default=None)
    parser.add_argument("--image-limit", type=int, default=10)
    parser.add_argument("--repeats", type=int, default=3)
    parser.add_argument("--thresholds", default=DEFAULT_THRESHOLDS)
    parser.add_argument("--attribute", default="area")
    parser.add_argument("--radio-adj", type=float, default=1.0)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    root = repo_root()
    thresholds = parse_thresholds(args.thresholds)
    images = collect_images(args.data_dir, args.image_limit)
    output_path = args.output or default_output_path(root, len(images), len(thresholds), args.repeats)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    if not args.binary.exists():
        raise RuntimeError(f"Benchmark binary not found: {args.binary}")

    rows: list[dict[str, object]] = []
    for repeat in range(args.repeats):
        for image_index, image in enumerate(images):
            print(f"[{repeat + 1}/{args.repeats}] {image.name}", flush=True)
            parsed, _ = run_compare(args.binary, image, args.attribute, args.radio_adj, thresholds)
            parsed.update(
                {
                    "image_index": image_index,
                    "image": image.name,
                    "image_path": str(image),
                    "repeat": repeat,
                    "attribute": args.attribute,
                    "radio_adj": args.radio_adj,
                    "thresholds": ",".join(thresholds),
                    "num_thresholds": len(thresholds),
                }
            )
            rows.append(parsed)

    fieldnames = [
        "image_index",
        "image",
        "image_path",
        "repeat",
        "rows",
        "cols",
        "attribute",
        "radio_adj",
        "thresholds",
        "num_thresholds",
        "dynamic_subtree_init_ms",
        "dynamic_subtree_iter_ms",
        "dynamic_subtree_reconstruct_ms",
        "dynamic_subtree_total_ms",
        "dynamic_leaf_init_ms",
        "dynamic_leaf_iter_ms",
        "dynamic_leaf_reconstruct_ms",
        "dynamic_leaf_total_ms",
        "naive_total_ms",
        "dynamic_subtree_speedup_vs_naive",
        "dynamic_leaf_speedup_vs_naive",
        "dynamic_subtree_matches_naive",
        "dynamic_leaf_matches_naive",
        "dynamic_subtree_matches_dynamic_leaf",
        "command",
    ]

    with output_path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)

    summary_path = write_summary(output_path, rows)
    print(f"Wrote {output_path}")
    print(f"Wrote {summary_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
