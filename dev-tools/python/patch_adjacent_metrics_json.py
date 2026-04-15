#!/usr/bin/env python3
import argparse
import json
import os


METRIC_FIELDS = (
    "count_adjacent_nodes",
    "avg_area_nodes_merged",
    "max_area_nodes_merged",
)


def read_json(path):
    with open(path, "r", encoding="utf-8") as fh:
        return json.load(fh)


def find_method(data, name):
    for method in data.get("methods", []):
        if method.get("method") == name:
            return method
    return None


def build_ref_index(ref_method):
    index = {}
    base_to_name = {}
    ambiguous_bases = set()

    for image in ref_method.get("images", []):
        image_name = image.get("image")
        if not image_name:
            continue

        base = os.path.basename(image_name)
        if base:
            if base in base_to_name and base_to_name[base] != image_name:
                ambiguous_bases.add(base)
            else:
                base_to_name[base] = image_name

        by_threshold = {}
        for iteration in image.get("iterations", []):
            threshold = iteration.get("threshold")
            if threshold is None:
                continue

            phase_values = {}
            valid = True
            for phase_key in ("phase1_alg_metrics", "phase2_alg_metrics"):
                metrics = iteration.get(phase_key, {})
                if not all(field in metrics for field in METRIC_FIELDS):
                    valid = False
                    break
                phase_values[phase_key] = {
                    field: metrics[field]
                    for field in METRIC_FIELDS
                }

            if valid:
                by_threshold[threshold] = phase_values

        if by_threshold:
            index[image_name] = by_threshold

    return index, base_to_name, ambiguous_bases


def resolve_ref_iterations(image_name, ref_index, base_to_name, ambiguous_bases):
    if image_name in ref_index:
        return ref_index[image_name]

    base = os.path.basename(image_name)
    if not base or base in ambiguous_bases:
        return None

    resolved_name = base_to_name.get(base)
    if not resolved_name:
        return None
    return ref_index.get(resolved_name)


def patch_method(target_method, ref_index, base_to_name, ambiguous_bases):
    patched_iterations = 0
    missed_images = 0
    missed_thresholds = 0

    for image in target_method.get("images", []):
        image_name = image.get("image")
        if not image_name:
            continue

        ref_iterations = resolve_ref_iterations(image_name, ref_index, base_to_name, ambiguous_bases)
        if ref_iterations is None:
            missed_images += 1
            continue

        for iteration in image.get("iterations", []):
            threshold = iteration.get("threshold")
            if threshold is None:
                continue

            ref_entry = ref_iterations.get(threshold)
            if ref_entry is None:
                missed_thresholds += 1
                continue

            for phase_key in ("phase1_alg_metrics", "phase2_alg_metrics"):
                target_metrics = iteration.setdefault(phase_key, {})
                ref_metrics = ref_entry[phase_key]
                for field in METRIC_FIELDS:
                    target_metrics[field] = ref_metrics[field]

            patched_iterations += 1

    return patched_iterations, missed_images, missed_thresholds


def main():
    parser = argparse.ArgumentParser(
        description=(
            "Patch adjacent-node metrics in a benchmark JSON using a corrected reference JSON."
        )
    )
    parser.add_argument("--in", dest="input_path", required=True, help="Input JSON with broken adjacent metrics.")
    parser.add_argument("--ref", dest="ref_path", required=True, help="Reference JSON regenerated after the fix.")
    parser.add_argument("--out", dest="output_path", required=True, help="Output JSON path.")
    parser.add_argument("--target-method", default="our_subtree", help="Method name to patch in the input JSON.")
    parser.add_argument("--ref-method", default="our_subtree", help="Method name to read from the reference JSON.")
    args = parser.parse_args()

    data_in = read_json(args.input_path)
    data_ref = read_json(args.ref_path)

    target_method = find_method(data_in, args.target_method)
    if target_method is None:
        raise SystemExit(f"target method not found: {args.target_method}")

    ref_method = find_method(data_ref, args.ref_method)
    if ref_method is None:
        raise SystemExit(f"reference method not found: {args.ref_method}")

    ref_index, base_to_name, ambiguous_bases = build_ref_index(ref_method)
    if not ref_index:
        raise SystemExit("reference JSON does not contain the adjacent metrics to copy")

    patched_iterations, missed_images, missed_thresholds = patch_method(
        target_method, ref_index, base_to_name, ambiguous_bases
    )
    if patched_iterations == 0:
        raise SystemExit("no iterations patched; check image names and thresholds")

    with open(args.output_path, "w", encoding="utf-8") as fh:
        json.dump(data_in, fh, ensure_ascii=True, indent=2)
        fh.write("\n")

    print(f"patched iterations: {patched_iterations}")
    print(f"images missing in reference: {missed_images}")
    print(f"thresholds missing in reference: {missed_thresholds}")


if __name__ == "__main__":
    main()
