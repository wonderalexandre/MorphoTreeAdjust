#!/usr/bin/env python3
import argparse
import json
import os
import re


# O método naive não extraí dados de flatzone. 
# Este patch adiciona os dados de flatzone para um json do método naive
#
# Exemplo:
#python3 ../../patch_flatzones_json.py \
#  --in out_1024p_our_hybrid_area_eager.json \
#  --ref out_1024p_our_area_eager.json \
#  --out out_1024p_our_hybrid_area_eagerFixed.json \
#  --target-method our_hybrid_eager --ref-method our_eager

def read_text(path):
    with open(path, "rb") as fh:
        return fh.read().decode("utf-8")


def parse_json(raw, path):
    if not raw.strip():
        raise ValueError(f"empty file: {path}")
    try:
        return json.loads(raw)
    except json.JSONDecodeError as exc:
        raise ValueError(f"invalid JSON (expected single object): {path}: {exc}") from exc



def find_method(data, name):
    methods = data.get("methods", [])
    for method in methods:
        if method.get("method") == name:
            return method
    return None


def build_ref_index(ref_method):
    index = {}
    image_metrics = {}
    base_to_name = {}
    base_ambiguous = set()
    images = ref_method.get("images", [])
    for img in images:
        image_name = img.get("image")
        if not image_name:
            continue
        base = os.path.basename(image_name)
        if base:
            if base in base_to_name and base_to_name[base] != image_name:
                base_ambiguous.add(base)
            else:
                base_to_name[base] = image_name
        iter_map = {}
        for it in img.get("iterations", []):
            threshold = it.get("threshold")
            if threshold is None:
                continue
            phase1 = it.get("phase1_metrics", {})
            phase2 = it.get("phase2_metrics", {})
            if "flatzones" not in phase1 or "flatzones" not in phase2:
                continue
            iter_map[threshold] = (phase1["flatzones"], phase2["flatzones"])
        if iter_map:
            index[image_name] = iter_map
        image_metrics[image_name] = img.get("num_flatzones_image")
    return index, image_metrics, base_to_name, base_ambiguous


def resolve_ref_entry(image_name, ref_index, base_to_name, base_ambiguous):
    entry = ref_index.get(image_name)
    if entry:
        return entry
    base = os.path.basename(image_name)
    if not base or base in base_ambiguous:
        return None
    ref_name = base_to_name.get(base)
    if not ref_name:
        return None
    return ref_index.get(ref_name)


def build_replacements(target_method, ref_index, base_to_name, base_ambiguous):
    replacements = []
    missed = 0
    for img in target_method.get("images", []):
        image_name = img.get("image")
        if not image_name:
            continue
        iter_map = resolve_ref_entry(image_name, ref_index, base_to_name, base_ambiguous)
        if not iter_map:
            missed += 1
            continue
        for it in img.get("iterations", []):
            threshold = it.get("threshold")
            if threshold is None:
                continue
            ref_vals = iter_map.get(threshold)
            if not ref_vals:
                continue
            replacements.append(ref_vals[0])
            replacements.append(ref_vals[1])
    return replacements, missed


def build_image_metrics_list(target_method, image_metrics, base_to_name, base_ambiguous):
    metrics_list = []
    missed = 0
    for img in target_method.get("images", []):
        image_name = img.get("image")
        value = image_metrics.get(image_name)
        if value is None:
            base = os.path.basename(image_name)
            if base and base not in base_ambiguous:
                ref_name = base_to_name.get(base)
                if ref_name:
                    value = image_metrics.get(ref_name)
        if value is None:
            metrics_list.append(None)
            missed += 1
        else:
            metrics_list.append(value)
    return metrics_list, missed


def find_matching_brace(text, start):
    depth = 1
    in_string = False
    escape = False
    for i in range(start + 1, len(text)):
        ch = text[i]
        if in_string:
            if escape:
                escape = False
            elif ch == "\\":
                escape = True
            elif ch == "\"":
                in_string = False
            continue
        if ch == "\"":
            in_string = True
        elif ch == "{":
            depth += 1
        elif ch == "}":
            depth -= 1
            if depth == 0:
                return i
    raise ValueError("unmatched brace while scanning object")


def find_image_objects(raw):
    pattern = re.compile(r'"images"\s*:')
    objects = []
    idx = 0
    while True:
        match = pattern.search(raw, idx)
        if not match:
            break
        pos = match.end()
        while pos < len(raw) and raw[pos] in " \t\r\n":
            pos += 1
        if pos >= len(raw) or raw[pos] != "[":
            raise ValueError("expected '[' after images key")
        pos += 1
        while pos < len(raw):
            while pos < len(raw) and raw[pos] in " \t\r\n,":
                pos += 1
            if pos >= len(raw):
                raise ValueError("unexpected end while scanning images array")
            if raw[pos] == "]":
                pos += 1
                break
            if raw[pos] != "{":
                raise ValueError("expected '{' in images array")
            end = find_matching_brace(raw, pos)
            objects.append((pos, end + 1))
            pos = end + 1
        idx = pos
    return objects


def infer_colon_spacing(obj_text):
    match = re.search(r'"[^"]+"\s*(:\s*)', obj_text)
    return match.group(1) if match else ":"


def insert_num_flatzones_image(obj_text, value):
    colon = infer_colon_spacing(obj_text)
    if "\n" not in obj_text:
        match = re.search(r'(,\s*)"total_ms"\s*:', obj_text)
        if not match:
            raise ValueError("total_ms key not found in image object")
        sep = match.group(1)
        insert = f'{sep}"num_flatzones_image"{colon}{value}{sep}'
        return obj_text[:match.start(1)] + insert + obj_text[match.end(1):]

    match = re.search(r'\n([ \t]*)"total_ms"\s*:', obj_text)
    if not match:
        raise ValueError("total_ms key not found in image object")
    indent = match.group(1)
    insert = f'\n{indent}"num_flatzones_image"{colon}{value},'
    return obj_text[:match.start()] + insert + obj_text[match.start():]


def set_num_flatzones_image(obj_text, value):
    pattern = re.compile(r'("num_flatzones_image"\s*:\s*)(-?\d+)')
    if pattern.search(obj_text):
        return pattern.sub(lambda m: f"{m.group(1)}{value}", obj_text, count=1)
    return insert_num_flatzones_image(obj_text, value)


def main():
    parser = argparse.ArgumentParser(
        description="Replace flatzones metric in a JSON benchmark summary using a reference file."
    )
    parser.add_argument("--in", dest="input_path", required=True, help="Input JSON (naive).")
    parser.add_argument("--ref", dest="ref_path", required=True, help="Reference JSON (our).")
    parser.add_argument("--out", dest="output_path", required=True, help="Output JSON path.")
    parser.add_argument("--target-method", default="naive", help="Method name to patch.")
    parser.add_argument("--ref-method", default="our", help="Reference method name.")
    args = parser.parse_args()

    raw_in = read_text(args.input_path)
    raw_ref = read_text(args.ref_path)
    data_in = parse_json(raw_in, args.input_path)
    data_ref = parse_json(raw_ref, args.ref_path)

    target_method = find_method(data_in, args.target_method)
    if target_method is None:
        raise SystemExit(f"target method not found: {args.target_method}")
    ref_method = find_method(data_ref, args.ref_method)
    if ref_method is None:
        raise SystemExit(f"reference method not found: {args.ref_method}")

    ref_index, image_metrics, base_to_name, base_ambiguous = build_ref_index(ref_method)
    if not ref_index:
        raise SystemExit("reference index empty: no flatzones metrics found")

    replacements, missed = build_replacements(target_method, ref_index, base_to_name, base_ambiguous)
    if not replacements:
        raise SystemExit("no replacements found; check images/thresholds")

    pattern = re.compile(r'("flatzones"\s*:\s*)(-?\d+)')
    total_matches = len(pattern.findall(raw_in))
    if total_matches != len(replacements):
        raise SystemExit(
            f"flatzones count mismatch: input has {total_matches}, "
            f"replacements {len(replacements)}"
        )

    it = iter(replacements)

    def repl(match):
        return f"{match.group(1)}{next(it)}"

    updated = pattern.sub(repl, raw_in)

    metrics_list, missed_images = build_image_metrics_list(
        target_method, image_metrics, base_to_name, base_ambiguous
    )
    if metrics_list:
        image_objects = find_image_objects(updated)
        if len(image_objects) != len(metrics_list):
            raise SystemExit(
                f"image object count mismatch: input has {len(image_objects)}, "
                f"expected {len(metrics_list)} (ensure one method per file)"
            )
        out = []
        last = 0
        for (start, end), value in zip(image_objects, metrics_list):
            out.append(updated[last:start])
            obj_text = updated[start:end]
            if value is not None:
                obj_text = set_num_flatzones_image(obj_text, value)
            out.append(obj_text)
            last = end
        out.append(updated[last:])
        updated = "".join(out)

    with open(args.output_path, "w", encoding="utf-8") as fh:
        fh.write(updated)

    print(
        f"patched {len(replacements)//2} iterations; images missing in reference: {missed}"
    )
    if metrics_list:
        print(f"patched num_flatzones_image; images missing metric: {missed_images}")


if __name__ == "__main__":
    main()
