#!/usr/bin/env python3
import argparse
import json
import os
import subprocess
import sys
import re

# Adiciona os dados de break-even em um json existente
#
# Exemple:
#python /Users/wonderalexandre/GitHub/MorphoTreeAdjust/tests/patch_find_break_even.py \
#  --naive /Users/wonderalexandre/GitHub/MorphoTreeAdjust/tests/build/results_news/out_480p_naive_area.json \
#  --our /Users/wonderalexandre/GitHub/MorphoTreeAdjust/tests/build/results_news/out_480p_our_subtree_area_onPixel.json \
#  --out /Users/wonderalexandre/GitHub/MorphoTreeAdjust/tests/build/results_news/out_480p_our_subtree_area_onPixel_break_even.json \
#  --graph-type on_pixel \
#  --repeat 3 --warmup 1 --radio-adj 1.5

# Exemplo extraindo somente os dados de break-even
# caffeinate -dims bash -lc '
#  for p in 480p 720p 1080p; do
#   /Users/wonderalexandre/GitHub/MorphoTreeAdjust/tests/patch_find_break_even.py \
#     --naive /Users/wonderalexandre/GitHub/MorphoTreeAdjust/tests/build/results_news/out_${p}_naive_area.json \
#     --our /Users/wonderalexandre/GitHub/MorphoTreeAdjust/tests/build/results_news/out_${p}_our_subtree_area_onPixel.json \
#     --graph-type on_pixel \
#     --repeat 3 --warmup 1 --radio-adj 1.5 \
#     --out /dev/null \
#     --cooldown-ms 20 \
#     --break-even-out /Users/wonderalexandre/GitHub/MorphoTreeAdjust/tests/build/results_news/break_even_only_${p}.json
# done
# '
#




def parse_args():
    p = argparse.ArgumentParser(description="Compute break-even index (m*) and patch our_subtree JSON.")
    p.add_argument("--naive", required=True, help="Naive JSON (with per-iteration stats)")
    p.add_argument("--our", required=True, help="our_subtree JSON to patch")
    p.add_argument("--out", help="Output JSON path (default: <our>_break_even.json)")
    p.add_argument(
        "--break-even-out",
        help="Write break-even data only to this JSON file (keeps original structure separate)",
    )
    default_runner = os.path.join(os.path.dirname(__file__), "build", "break_even_runner")
    p.add_argument("--runner", default=default_runner, help="break_even_runner binary path")
    p.add_argument("--repeat", type=int, default=3)
    p.add_argument("--warmup", type=int, default=1)
    p.add_argument("--radio-adj", type=float, default=1.5)
    p.add_argument("--graph-type", default="on_pixel", choices=["on_pixel", "on_demand", "eager", "naive", "pixel"])
    p.add_argument("--cooldown-ms", type=int, default=0, help="Cooldown between k0 evals (ms)")
    p.add_argument("--image-root", default="", help="Base directory to resolve relative image paths")
    p.add_argument("--limit", type=int, default=0, help="Limit number of images (0=all)")
    p.add_argument("--quiet", action="store_true")
    p.add_argument("--debug", action="store_true", help="Print debug stats to stderr")
    return p.parse_args()


def load_json(path):
    with open(path, "r", encoding="utf-8") as fh:
        return json.load(fh)


def build_naive_index(naive_data):
    methods = naive_data.get("methods", [])
    if not methods:
        raise ValueError("naive JSON: no methods")
    method = methods[0]
    images = method.get("images", [])
    index = {}
    for img in images:
        name = img.get("image")
        iters = img.get("iterations", [])
        if not name or not iters:
            continue
        thresholds = [it.get("threshold") for it in iters]
        totals = []
        ok = True
        for it in iters:
            total = it.get("total", {})
            if "median" not in total:
                ok = False
                break
            totals.append(float(total["median"]))
        if not ok:
            continue
        cum = []
        s = 0.0
        for v in totals:
            s += v
            cum.append(s)
        index[name] = {
            "thresholds": thresholds,
            "cum": cum,
        }
    return index


def call_runner(runner, image, thresholds, naive_cum, repeat, warmup, radio_adj, graph_type, cooldown_ms, quiet):
    thresholds_csv = ",".join(str(v) for v in thresholds)
    naive_csv = ",".join(f"{v:.6f}" for v in naive_cum)
    cmd = [
        runner,
        "--image", image,
        "--thresholds", thresholds_csv,
        "--naive-cum", naive_csv,
        "--repeat", str(repeat),
        "--warmup", str(warmup),
        "--radio-adj", str(radio_adj),
        "--graph-type", graph_type,
        "--cooldown-ms", str(cooldown_ms),
    ]
    if quiet:
        cmd.append("--quiet")
    proc = subprocess.run(cmd, capture_output=True, text=True)
    if proc.returncode != 0:
        raise RuntimeError(f"runner failed for {image}: {proc.stderr.strip()}")
    raw = proc.stdout.strip()
    if not raw:
        raise RuntimeError(f"runner returned empty output for {image}")
    return json.loads(raw)


def main():
    args = parse_args()
    naive_data = load_json(args.naive)
    our_data = load_json(args.our)
    naive_index = build_naive_index(naive_data)

    if not os.path.exists(args.runner):
        raise FileNotFoundError(f"runner not found: {args.runner}")

    image_root = args.image_root
    if not image_root:
        image_root = os.path.abspath(os.path.join(os.path.dirname(args.naive), ".."))

    methods = our_data.get("methods", [])
    if not methods:
        raise ValueError("our JSON: no methods")

    if args.out:
        out_path = args.out
    else:
        base, ext = os.path.splitext(args.our)
        if not ext:
            ext = ".json"
        out_path = f"{base}_break_even{ext}"

    # Pick the our_subtree method (any name that starts with our_subtree)
    target_method = None
    for m in methods:
        name = m.get("method", "")
        if name.startswith("our_subtree"):
            target_method = m
            break
    if target_method is None:
        raise ValueError("our JSON: no method starting with 'our_subtree'")

    images = target_method.get("images", [])
    if not images:
        raise ValueError("our JSON: no images")

    limit = args.limit if args.limit and args.limit > 0 else len(images)
    results = {}
    break_even_entries = []
    for idx, img in enumerate(images[:limit]):
        name = img.get("image")
        if not name:
            continue
        naive_info = naive_index.get(name)
        if not naive_info:
            # try basename match
            base = os.path.basename(name)
            if base:
                for key, info in naive_index.items():
                    if os.path.basename(key) == base:
                        naive_info = info
                        break
        if not naive_info:
            if not args.quiet:
                print(f"skip (naive not found): {name}", file=sys.stderr)
            continue
        if not args.quiet:
            print(f"[{idx+1}/{limit}] {name}", file=sys.stderr)
        image_path = name
        if not os.path.isabs(image_path):
            image_path = os.path.normpath(os.path.join(image_root, image_path))
        result = call_runner(
            args.runner,
            image_path,
            naive_info["thresholds"],
            naive_info["cum"],
            args.repeat,
            args.warmup,
            args.radio_adj,
            args.graph_type,
            args.cooldown_ms,
            args.quiet,
        )
        results[name] = result
        break_even_entries.append(result)
        if args.debug:
            print(f"result ok: {name}", file=sys.stderr)

    raw = None
    with open(args.our, "r", encoding="utf-8") as fh:
        raw = fh.read()
    if raw is None:
        raise ValueError("failed to read our JSON")

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

    def find_image_objects(raw_text):
        pattern = re.compile(r'"images"\s*:')
        objects = []
        idx = 0
        while True:
            match = pattern.search(raw_text, idx)
            if not match:
                break
            pos = match.end()
            while pos < len(raw_text) and raw_text[pos] in " \t\r\n":
                pos += 1
            if pos >= len(raw_text) or raw_text[pos] != "[":
                raise ValueError("expected '[' after images key")
            pos += 1
            while pos < len(raw_text):
                while pos < len(raw_text) and raw_text[pos] in " \t\r\n,":
                    pos += 1
                if pos >= len(raw_text):
                    raise ValueError("unexpected end while scanning images array")
                if raw_text[pos] == "]":
                    pos += 1
                    break
                if raw_text[pos] != "{":
                    raise ValueError("expected '{' in images array")
                end = find_matching_brace(raw_text, pos)
                objects.append((pos, end + 1))
                pos = end + 1
            idx = pos
        return objects

    def infer_indent(obj_text):
        match = re.search(r'\n([ \t]*)\"iterations\"\s*:', obj_text)
        if match:
            return match.group(1)
        match = re.search(r'\n([ \t]*)\"total_ms\"\s*:', obj_text)
        if match:
            return match.group(1)
        return "  "

    def replace_or_insert_break_even(obj_text, break_even_json):
        if "\"break_even\"" in obj_text:
            key_idx = obj_text.find("\"break_even\"")
            colon_idx = obj_text.find(":", key_idx)
            if colon_idx == -1:
                return obj_text
            i = colon_idx + 1
            while i < len(obj_text) and obj_text[i] in " \t\r\n":
                i += 1
            if i >= len(obj_text) or obj_text[i] != "{":
                return obj_text
            end = find_matching_brace(obj_text, i)
            return obj_text[:i] + break_even_json + obj_text[end + 1:]
        indent = infer_indent(obj_text)
        insert_pos = None
        m = re.search(r'\n[ \t]*\"iterations\"\s*:', obj_text)
        if m:
            insert_pos = m.start()
        else:
            insert_pos = obj_text.rfind("}")
        if insert_pos is None or insert_pos < 0:
            return obj_text
        insert = f'\\n{indent}\"break_even\":{break_even_json},'
        return obj_text[:insert_pos] + insert + obj_text[insert_pos:]

    def extract_image_name(obj_text):
        key = "\"image\""
        idx = obj_text.find(key)
        if idx == -1:
            return None
        idx = obj_text.find(":", idx + len(key))
        if idx == -1:
            return None
        idx += 1
        while idx < len(obj_text) and obj_text[idx] in " \t\r\n":
            idx += 1
        if idx >= len(obj_text) or obj_text[idx] != "\"":
            return None
        idx += 1
        start = idx
        escape = False
        while idx < len(obj_text):
            ch = obj_text[idx]
            if escape:
                escape = False
            elif ch == "\\":
                escape = True
            elif ch == "\"":
                return obj_text[start:idx]
            idx += 1
        return None

    objects = find_image_objects(raw)
    if args.debug:
        sample_key = next(iter(results), None)
        head = ""
        if objects:
            head = raw[objects[0][0]:min(objects[0][0] + 80, len(raw))]
        head_match = extract_image_name(head)
        print(
            f"debug: images={len(images)} results={len(results)} objects={len(objects)} sample={sample_key} head_match={head_match} head={head!r}",
            file=sys.stderr,
        )
    new_raw = raw
    matched = 0
    seen = 0
    for start, end in reversed(objects):
        obj_text = new_raw[start:end]
        name = extract_image_name(obj_text)
        if not name:
            continue
        if args.debug and seen < 1:
            print(f"debug: first_name={name}", file=sys.stderr)
            seen += 1
        if name not in results:
            continue
        matched += 1
        break_even_json = json.dumps(results[name], ensure_ascii=False, separators=(",", ":"))
        new_obj = replace_or_insert_break_even(obj_text, break_even_json)
        new_raw = new_raw[:start] + new_obj + new_raw[end:]
        if args.debug:
            print(f"patched: {name}", file=sys.stderr)
    if args.debug:
        print(f"debug: matched={matched}", file=sys.stderr)

    with open(out_path, "w", encoding="utf-8") as fh:
        fh.write(new_raw)
    if not args.quiet:
        print(f"Updated: {out_path}")

    if args.break_even_out:
        payload = {
            "method": target_method.get("method", ""),
            "graph_type": args.graph_type,
            "naive": args.naive,
            "our": args.our,
            "break_even": break_even_entries,
        }
        with open(args.break_even_out, "w", encoding="utf-8") as fh:
            json.dump(payload, fh, indent=2, ensure_ascii=False)
        if not args.quiet:
            print(f"Break-even only: {args.break_even_out}")


if __name__ == "__main__":
    main()
