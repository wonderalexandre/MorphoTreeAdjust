#!/usr/bin/env python3
import argparse
import json
import os
import re
import subprocess
import sys
import tempfile


METRIC_FIELDS = (
    "count_adjacent_nodes",
    "avg_area_nodes_merged",
    "max_area_nodes_merged",
)


def read_json(path):
    with open(path, "r", encoding="utf-8") as fh:
        return json.load(fh)


def write_json(path, data):
    with open(path, "w", encoding="utf-8") as fh:
        json.dump(data, fh, ensure_ascii=True, indent=2)
        fh.write("\n")


def find_method(data, name):
    for method in data.get("methods", []):
        if method.get("method") == name:
            return method
    return None


def build_ref_index(ref_method):
    index = {}
    for image in ref_method.get("images", []):
        image_name = image.get("image")
        if not image_name:
            continue

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

    return index


def patch_method(target_method, ref_index):
    patched_iterations = 0
    missed_images = 0
    missed_thresholds = 0

    for image in target_method.get("images", []):
        image_name = image.get("image")
        ref_iterations = ref_index.get(image_name)
        if ref_iterations is None:
            missed_images += 1
            continue

        for iteration in image.get("iterations", []):
            threshold = iteration.get("threshold")
            ref_entry = ref_iterations.get(threshold)
            if ref_entry is None:
                missed_thresholds += 1
                continue

            for phase_key in ("phase1_alg_metrics", "phase2_alg_metrics"):
                target_metrics = iteration.setdefault(phase_key, {})
                for field in METRIC_FIELDS:
                    target_metrics[field] = ref_entry[phase_key][field]
            patched_iterations += 1

    return patched_iterations, missed_images, missed_thresholds


def infer_attribute_flag(method):
    attribute = method.get("attribute")
    if not attribute:
        raise ValueError("input JSON missing method.attribute")
    return str(attribute)


def infer_radio_adj(method):
    radio_adj = method.get("radio_adj")
    if radio_adj is None:
        raise ValueError("input JSON missing method.radio_adj")
    return str(radio_adj)


def infer_thresholds_arg(method):
    thresholds = method.get("thresholds")
    if not thresholds:
        raise ValueError("input JSON missing method.thresholds")
    return ",".join(str(int(v)) for v in thresholds)


def infer_images(method):
    images = []
    for image in method.get("images", []):
        image_name = image.get("image")
        if image_name:
            images.append(image_name)
    if not images:
        raise ValueError("input JSON missing method.images")
    return images


def run_benchmark(executable, method_name, method):
    thresholds_arg = infer_thresholds_arg(method)
    attribute_arg = infer_attribute_flag(method)
    radio_adj_arg = infer_radio_adj(method)
    images = infer_images(method)

    cmd = [
        executable,
        "--json",
        "--quiet",
        "--repeat", "1",
        "--warmup", "0",
        "--no-validate",
        "--method", method_name,
        "--attribute", attribute_arg,
        "--thresholds", thresholds_arg,
    ]
    for image_path in images:
        cmd.extend(["--image", image_path])

    with tempfile.NamedTemporaryFile("w", encoding="utf-8", suffix=".conf") as cfg:
        cfg.write(f"radio_adj={radio_adj_arg}\n")
        cfg.flush()
        cmd.extend(["--config", cfg.name])

        with tempfile.NamedTemporaryFile("w+", encoding="utf-8", suffix=".json") as stdout_tmp:
            proc = subprocess.Popen(
                cmd,
                stdout=stdout_tmp,
                stderr=subprocess.PIPE,
                text=True,
            )

            assert proc.stderr is not None

            stderr_chunks = []
            progress_pattern = re.compile(r"^\[\d+/\d+\] processing image: ")
            for line in proc.stderr:
                stderr_chunks.append(line)
                if progress_pattern.match(line):
                    print(line, end="", flush=True)
                else:
                    sys.stderr.write(line)
                    sys.stderr.flush()

            return_code = proc.wait()
            stdout_tmp.flush()
            stdout_tmp.seek(0)
            stdout_text = stdout_tmp.read()

    if return_code != 0:
        raise SystemExit(f"benchmark failed with exit code {return_code}")

    try:
        return json.loads(stdout_text)
    except json.JSONDecodeError as exc:
        stderr_text = "".join(stderr_chunks)
        if stderr_text:
            sys.stderr.write(stderr_text)
        raise SystemExit(f"benchmark did not emit valid JSON: {exc}") from exc


def main():
    parser = argparse.ArgumentParser(
        description=(
            "Recompute adjacent metrics by rerunning jmiv2026_benchmark from the parameters stored in an old JSON."
        )
    )
    parser.add_argument("--in", dest="input_path", required=True, help="Old JSON with zeroed adjacent metrics.")
    parser.add_argument("--out", dest="output_path", required=True, help="Output JSON path.")
    parser.add_argument(
        "--benchmark",
        dest="benchmark_path",
        required=True,
        help="Path to the corrected jmiv2026_benchmark executable.",
    )
    parser.add_argument("--target-method", default="our_subtree", help="Method name to recompute.")
    args = parser.parse_args()

    data_in = read_json(args.input_path)
    target_method = find_method(data_in, args.target_method)
    if target_method is None:
        raise SystemExit(f"target method not found: {args.target_method}")

    benchmark_data = run_benchmark(args.benchmark_path, args.target_method, target_method)
    ref_method = find_method(benchmark_data, args.target_method)
    if ref_method is None:
        raise SystemExit(f"recomputed JSON missing method: {args.target_method}")

    ref_index = build_ref_index(ref_method)
    if not ref_index:
        raise SystemExit("recomputed JSON does not contain adjacent metrics")

    patched_iterations, missed_images, missed_thresholds = patch_method(target_method, ref_index)
    if patched_iterations == 0:
        raise SystemExit("no iterations patched; image names or thresholds did not match")

    write_json(args.output_path, data_in)
    print(f"patched iterations: {patched_iterations}")
    print(f"images missing in recomputed run: {missed_images}")
    print(f"thresholds missing in recomputed run: {missed_thresholds}")


if __name__ == "__main__":
    main()
