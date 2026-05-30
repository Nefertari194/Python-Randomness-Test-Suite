from __future__ import print_function

import math
import os
import unicodedata

from bitio import iter_bit_sequences

from approximate_entropy_test import approximate_entropy_test
from binary_matrix_rank_test import binary_matrix_rank_test
from cumulative_sums_test import cumulative_sums_test
from dft_test import dft_test
from frequency_within_block_test import frequency_within_block_test
from linear_complexity_test import linear_complexity_test
from longest_run_ones_in_a_block_test import longest_run_ones_in_a_block_test
from maurers_universal_test import maurers_universal_test
from monobit_test import monobit_test
from non_overlapping_template_matching_test import non_overlapping_template_matching_test
from overlapping_template_matching_test import overlapping_template_matching_test
from random_excursion_test import random_excursion_test
from random_excursion_variant_test import random_excursion_variant_test
from runs_test import runs_test
from serial_test import serial_test


def _format_p_values(p_values):
    if p_values is None:
        return ""
    if isinstance(p_values, (list, tuple)):
        return "[" + ", ".join("%.6g" % float(x) for x in p_values) + "]"
    return "%.6g" % float(p_values)


def _normalize_test_result(success, p, plist_or_meta):
    if isinstance(plist_or_meta, dict) and plist_or_meta.get("reason"):
        return {"status": "skip", "reason": plist_or_meta.get("reason"), "p_values": None}
    if p is None and isinstance(plist_or_meta, dict) and plist_or_meta.get("p_values") is None and plist_or_meta.get("reason"):
        return {"status": "skip", "reason": plist_or_meta.get("reason"), "p_values": None}
    if isinstance(plist_or_meta, dict) and "p_values" in plist_or_meta:
        p_values = plist_or_meta.get("p_values")
    else:
        p_values = plist_or_meta
    if p_values is None and p is not None:
        p_values = p
    return {"status": "pass" if success else "fail", "reason": None, "p_values": p_values}


def _median(values):
    if not values:
        return None
    xs = sorted(values)
    n = len(xs)
    mid = n // 2
    if (n % 2) == 1:
        return xs[mid]
    return (xs[mid - 1] + xs[mid]) / 2.0


def _mean(values):
    if not values:
        return None
    return float(sum(values)) / float(len(values))


def _format_rate(x):
    if x is None:
        return "-"
    return "%.3f" % float(x)


def _format_p(x):
    if x is None:
        return "-"
    return "%.6g" % float(x)


def _disp_width(s):
    w = 0
    for ch in str(s):
        if unicodedata.combining(ch):
            continue
        if unicodedata.east_asian_width(ch) in ("F", "W"):
            w += 2
        else:
            w += 1
    return w


def _truncate_by_width(s, max_width):
    out = []
    w = 0
    for ch in str(s):
        if unicodedata.combining(ch):
            continue
        cw = 2 if unicodedata.east_asian_width(ch) in ("F", "W") else 1
        if w + cw > max_width:
            break
        out.append(ch)
        w += cw
    if _disp_width(s) > max_width and max_width >= 1:
        if max_width >= 3 and w >= 3:
            out = out[:-3] + [".", ".", "."]
        elif max_width >= 1 and out:
            out[-1] = "."
    return "".join(out)


def _render_box_table(headers, rows):
    cols = len(headers)
    widths = [_disp_width(h) for h in headers]
    for row in rows:
        for i in range(cols):
            widths[i] = max(widths[i], _disp_width(row[i]))

    def sep(ch_left, ch_mid, ch_right, ch_fill):
        parts = [ch_fill * (w + 2) for w in widths]
        return ch_left + ch_mid.join(parts) + ch_right

    def line(cells):
        out = []
        for i, c in enumerate(cells):
            s = str(c)
            out.append(" " + s + " " * (widths[i] - _disp_width(s) + 1))
        return "|" + "|".join(out) + "|"

    lines = []
    lines.append(sep("+", "+", "+", "-"))
    lines.append(line(headers))
    lines.append(sep("+", "+", "+", "-"))
    for r in rows:
        lines.append(line(r))
    lines.append(sep("+", "+", "+", "-"))
    return "\n".join(lines)


def recommend_config(sequence_length, alpha=0.01):
    n = int(sequence_length)
    params = {
        "serial_m": 4,
        "serial_max_m": 16,
        "freq_block_size": 32,
        "freq_num_blocks": None,
        "nonoverlap_group": 0,
        "nonoverlap_index": 0,
        "overlap_allow_degraded": True,
        "overlap_m": 6,
        "overlap_num_blocks": max(2, int(math.floor(n / 128))) if n > 0 else 2,
        "overlap_block_size": None,
        "random_excursion_min_J": 50,
        "random_excursion_variant_min_J": 50,
        "rank_M": 32,
        "rank_Q": 32,
        "maurers_L": None,
        "maurers_Q": None,
        "linear_M": None,
    }

    notes = [
        "random_excursion 把 J 的门槛调低了（默认 500，这里用 50），结果可能不符合标准",
        "random_excursion_variant 把 J 的门槛调低了（默认 500，这里用 50），结果可能不符合标准",
        "overlapping_template 用了降级分块方式，结果可能不符合标准",
        "non_overlapping_template 这里只测了一个固定模板，不是标准里那种全模板跑法",
    ]

    return {
        "alpha": float(alpha),
        "tests": None,
        "params": params,
        "notes": notes,
    }


def _format_param_summary(params):
    p = params or {}
    items = [
        ("serial_m", "serial 的 m（模式长度）"),
        ("freq_block_size", "block_frequency 的块大小"),
        ("freq_num_blocks", "block_frequency 的块数"),
        ("nonoverlap_group", "non_overlapping_template 的模板组"),
        ("nonoverlap_index", "non_overlapping_template 的模板下标"),
        ("overlap_m", "overlapping_template 的模板长度"),
        ("overlap_num_blocks", "overlapping_template 的块数"),
        ("random_excursion_min_J", "random_excursion 的最小 J"),
        ("random_excursion_variant_min_J", "random_excursion_variant 的最小 J"),
        ("rank_M", "matrix_rank 的 M"),
        ("rank_Q", "matrix_rank 的 Q"),
    ]
    lines = []
    for k, label in items:
        if k not in p:
            continue
        v = p.get(k)
        if v is None:
            continue
        lines.append("%s=%s（%s）" % (k, v, label))
    return lines


def run_sp800_22(config):
    key_path = config.get("key_path")
    if not key_path:
        raise ValueError("config 里得有 key_path")

    alpha = float(config.get("alpha", 0.01))
    input_format = config.get("input_format", "auto")
    bigendian = bool(config.get("bigendian", True))
    recursive = bool(config.get("recursive", False))
    verbose = bool(config.get("verbose", False))
    print_results = bool(config.get("print_results", True))
    max_items = config.get("max_items", None)

    params = config.get("params", {}) or {}
    enabled = config.get("tests", None)
    if enabled is not None:
        enabled = set(enabled)

    test_table = [
        ("monobit", lambda b: monobit_test(b, alpha=alpha, verbose=verbose)),
        ("frequency_within_block", lambda b: frequency_within_block_test(b, alpha=alpha, verbose=verbose)),
        ("runs", lambda b: runs_test(b, alpha=alpha, verbose=verbose)),
        ("dft", lambda b: dft_test(b, alpha=alpha, verbose=verbose)),
        ("cumulative_sums", lambda b: cumulative_sums_test(b, alpha=alpha, verbose=verbose)),
        ("approximate_entropy", lambda b: approximate_entropy_test(b, alpha=alpha, verbose=verbose)),
        ("serial", lambda b: serial_test(b, patternlen=params.get("serial_m"), alpha=alpha, verbose=verbose, max_m=int(params.get("serial_max_m", 16)))),
        ("binary_matrix_rank", lambda b: binary_matrix_rank_test(b, M=int(params.get("rank_M", 32)), Q=int(params.get("rank_Q", 32)), alpha=alpha, verbose=verbose)),
        ("longest_run_ones", lambda b: longest_run_ones_in_a_block_test(b, alpha=alpha, verbose=verbose)),
        ("maurers_universal", lambda b: maurers_universal_test(b, patternlen=params.get("maurers_L"), initblocks=params.get("maurers_Q"), alpha=alpha, verbose=verbose)),
        ("linear_complexity", lambda b: linear_complexity_test(b, patternlen=params.get("linear_M"), alpha=alpha, verbose=verbose)),
        ("non_overlapping_template", lambda b: non_overlapping_template_matching_test(b, alpha=alpha, verbose=verbose)),
        ("overlapping_template", lambda b: overlapping_template_matching_test(b, alpha=alpha, verbose=verbose)),
        ("random_excursion", lambda b: random_excursion_test(b, alpha=alpha, verbose=verbose)),
        ("random_excursion_variant", lambda b: random_excursion_variant_test(b, alpha=alpha, verbose=verbose)),
    ]

    results = []
    item_count = 0
    for item in iter_bit_sequences(key_path, input_format=input_format, bigendian=bigendian, recursive=recursive):
        item_count += 1
        if max_items is not None and item_count > int(max_items):
            break

        source = item["source"]
        index = item["index"]
        bits = item["bits"]

        tag = "%s#%d" % (os.path.basename(source), index)
        if not bits:
            for test_name, _ in test_table:
                if enabled is not None and test_name not in enabled:
                    continue
                row = {"tag": tag, "test": test_name, "status": "skip", "p_values": None, "reason": "这一段没读到 0/1"}
                results.append(row)
                if print_results:
                    print("%s  %-24s  skip  %s" % (tag, test_name, row["reason"]))
            continue

        for test_name, fn in test_table:
            if enabled is not None and test_name not in enabled:
                continue

            success, p, plist_or_meta = fn(bits)
            norm = _normalize_test_result(success, p, plist_or_meta)

            row = {"tag": tag, "test": test_name, "status": norm["status"], "p_values": norm["p_values"], "reason": norm["reason"]}
            results.append(row)

            if not print_results:
                continue

            if row["status"] == "skip":
                print("%s  %-24s  skip  %s" % (tag, test_name, row["reason"]))
            else:
                pv = _format_p_values(row["p_values"])
                print("%s  %-24s  %s  %s" % (tag, test_name, row["status"], pv))

    return results


def run_sp800_22_batch(config):
    key_path = config.get("key_path")
    if not key_path:
        raise ValueError("config 里得有 key_path")

    alpha = float(config.get("alpha", 0.01))
    input_format = config.get("input_format", "auto")
    bigendian = bool(config.get("bigendian", True))
    recursive = bool(config.get("recursive", False))
    verbose = bool(config.get("verbose", False))
    print_results = bool(config.get("print_results", True))
    max_items = config.get("max_items", None)

    params = config.get("params", {}) or {}
    enabled = config.get("tests", None)
    if enabled is not None:
        enabled = set(enabled)

    expected = 1.0 - alpha
    test_table = [
        ("monobit", lambda b: monobit_test(b, alpha=alpha, verbose=verbose)),
        ("frequency_within_block", lambda b: frequency_within_block_test(b, alpha=alpha, verbose=verbose, block_size=params.get("freq_block_size"), num_blocks=params.get("freq_num_blocks"))),
        ("runs", lambda b: runs_test(b, alpha=alpha, verbose=verbose)),
        ("dft", lambda b: dft_test(b, alpha=alpha, verbose=verbose)),
        ("cumulative_sums", lambda b: cumulative_sums_test(b, alpha=alpha, verbose=verbose)),
        ("approximate_entropy", lambda b: approximate_entropy_test(b, alpha=alpha, verbose=verbose)),
        ("serial", lambda b: serial_test(b, patternlen=params.get("serial_m"), alpha=alpha, verbose=verbose, max_m=int(params.get("serial_max_m", 16)))),
        ("binary_matrix_rank", lambda b: binary_matrix_rank_test(b, M=int(params.get("rank_M", 32)), Q=int(params.get("rank_Q", 32)), alpha=alpha, verbose=verbose)),
        ("longest_run_ones", lambda b: longest_run_ones_in_a_block_test(b, alpha=alpha, verbose=verbose)),
        ("maurers_universal", lambda b: maurers_universal_test(b, patternlen=params.get("maurers_L"), initblocks=params.get("maurers_Q"), alpha=alpha, verbose=verbose)),
        ("linear_complexity", lambda b: linear_complexity_test(b, patternlen=params.get("linear_M"), alpha=alpha, verbose=verbose)),
        ("non_overlapping_template", lambda b: non_overlapping_template_matching_test(b, alpha=alpha, verbose=verbose, template_group=params.get("nonoverlap_group", 0), template_index=params.get("nonoverlap_index", 0))),
        ("overlapping_template", lambda b: overlapping_template_matching_test(b, blen=params.get("overlap_m", 6), alpha=alpha, verbose=verbose, block_size=params.get("overlap_block_size"), num_blocks=params.get("overlap_num_blocks"), allow_degraded=bool(params.get("overlap_allow_degraded", False)))),
        ("random_excursion", lambda b: random_excursion_test(b, alpha=alpha, verbose=verbose, min_J=int(params.get("random_excursion_min_J", 500)))),
        ("random_excursion_variant", lambda b: random_excursion_variant_test(b, alpha=alpha, verbose=verbose, min_J=int(params.get("random_excursion_variant_min_J", 500)))),
    ]

    def ensure_stat(name):
        if name not in stats:
            stats[name] = {"pass": 0, "fail": 0, "skip": 0, "p_list": [], "skip_reasons": {}}
        return stats[name]

    stats = {}
    n_list = []
    total_items = 0

    for item in iter_bit_sequences(key_path, input_format=input_format, bigendian=bigendian, recursive=recursive):
        total_items += 1
        if max_items is not None and total_items > int(max_items):
            break

        bits = item["bits"]
        if bits:
            n_list.append(len(bits))

        if not bits:
            for test_name, _ in test_table:
                if enabled is not None and test_name not in enabled:
                    continue
                s = ensure_stat(test_name)
                s["skip"] += 1
                s["skip_reasons"]["这一段没读到 0/1"] = s["skip_reasons"].get("这一段没读到 0/1", 0) + 1
            continue

        for test_name, fn in test_table:
            if enabled is not None and test_name not in enabled:
                continue

            success, p, plist_or_meta = fn(bits)
            norm = _normalize_test_result(success, p, plist_or_meta)

            s = ensure_stat(test_name)
            if norm["status"] == "skip":
                s["skip"] += 1
                reason = norm.get("reason") or "条件不足"
                s["skip_reasons"][reason] = s["skip_reasons"].get(reason, 0) + 1
            else:
                if norm["status"] == "pass":
                    s["pass"] += 1
                else:
                    s["fail"] += 1
                p_values = norm.get("p_values")
                if isinstance(p_values, (list, tuple)):
                    for pv in p_values:
                        s["p_list"].append(float(pv))
                elif p_values is not None:
                    s["p_list"].append(float(p_values))

    note_max = int(config.get("note_max", 28))

    def finalize_row(name, s):
        p_count = s["pass"]
        f_count = s["fail"]
        sk_count = s["skip"]
        n_eff = p_count + f_count
        pass_rate = (float(p_count) / float(n_eff)) if n_eff > 0 else None
        if n_eff == 0:
            verdict = "全跳过"
        elif pass_rate is not None and pass_rate >= expected:
            verdict = "正常"
        else:
            verdict = "异常"

        p_min = min(s["p_list"]) if s["p_list"] else None
        p_avg = _mean(s["p_list"]) if s["p_list"] else None
        p_max = max(s["p_list"]) if s["p_list"] else None

        note = ""
        if sk_count > 0:
            reason = None
            if s["skip_reasons"]:
                reason = max(s["skip_reasons"].items(), key=lambda x: x[1])[0]
            if n_eff == 0:
                note = reason or "条件不足"
            else:
                note = "跳过最多: %s" % (reason or "条件不足")
        if note:
            note = _truncate_by_width(note, note_max)

        return [
            name,
            str(p_count),
            str(f_count),
            str(sk_count),
            str(n_eff),
            _format_rate(pass_rate),
            verdict,
            _format_p(p_min),
            _format_p(p_avg),
            _format_p(p_max),
            note,
        ]

    keys = list(stats.keys())
    base_order = [t[0] for t in test_table]
    def sort_key(k):
        if ":" in k:
            base, _ = k.split(":", 1)
        else:
            base = k
        try:
            idx = base_order.index(base)
        except ValueError:
            idx = 9999
        sub = 0 if ":" not in k else 1
        return (idx, sub, k)

    keys.sort(key=sort_key)

    headers = ["测试", "通过", "失败", "跳过", "有效", "通过率", "判定", "P最小", "P平均", "P最大", "备注"]
    rows = [finalize_row(k, stats[k]) for k in keys]

    summary = {
        "key_path": key_path,
        "alpha": alpha,
        "total": total_items,
        "n_min": min(n_list) if n_list else None,
        "n_med": _median(n_list) if n_list else None,
        "n_max": max(n_list) if n_list else None,
        "tests": {},
    }
    for k in keys:
        s = stats[k]
        n_eff = s["pass"] + s["fail"]
        pass_rate = (float(s["pass"]) / float(n_eff)) if n_eff > 0 else None
        if n_eff == 0:
            verdict = "全跳过"
        elif pass_rate is not None and pass_rate >= expected:
            verdict = "正常"
        else:
            verdict = "异常"
        top_reason = None
        if s["skip_reasons"]:
            top_reason = max(s["skip_reasons"].items(), key=lambda x: x[1])[0]
        summary["tests"][k] = {
            "pass": s["pass"],
            "fail": s["fail"],
            "skip": s["skip"],
            "n_eff": n_eff,
            "pass_rate": pass_rate,
            "expected": expected,
            "verdict": verdict,
            "p_min": min(s["p_list"]) if s["p_list"] else None,
            "p_avg": _mean(s["p_list"]) if s["p_list"] else None,
            "p_max": max(s["p_list"]) if s["p_list"] else None,
            "top_skip_reason": top_reason,
        }

    if print_results:
        n_info = "长度: -"
        if summary["n_min"] is not None:
            n_info = "长度: min=%s  中位=%s  max=%s" % (summary["n_min"], int(summary["n_med"]), summary["n_max"])
        head = [
            "数据: %s" % key_path,
            "序列数: %d" % total_items,
            n_info,
            "alpha=%s（p >= alpha 算通过）" % alpha,
            "通过率分母=通过+失败",
        ]
        param_lines = _format_param_summary(params)
        if param_lines:
            head.append("参数: " + "  ".join(param_lines))
        notes = config.get("notes") or []
        if notes:
            head.append("提示: " + "  ".join(str(x) for x in notes))
        print("\n".join(head))
        print(_render_box_table(headers, rows))

        detail_max = int(config.get("reason_topk", 2))
        reason_lines = []
        for k in keys:
            s = stats[k]
            if not s["skip_reasons"]:
                continue
            top = sorted(s["skip_reasons"].items(), key=lambda x: x[1], reverse=True)[:detail_max]
            chunks = ["%s(%d)" % (t[0], t[1]) for t in top]
            reason_lines.append("%s: %s" % (k, "  ".join(chunks)))
        if reason_lines:
            print("跳过原因统计（每项最多 %d 条）" % detail_max)
            for ln in reason_lines:
                print(ln)

    return summary


def run_detail(key_path, **kwargs):
    cfg = {"key_path": key_path}
    cfg.update(kwargs)
    return run_sp800_22(cfg)


def run_batch(key_path, **kwargs):
    cfg = {"key_path": key_path}
    cfg.update(kwargs)
    return run_sp800_22_batch(cfg)


def run_batch_recommended(key_path, sequence_length, **kwargs):
    base = recommend_config(sequence_length, alpha=kwargs.get("alpha", 0.01))
    base.update(kwargs)
    return run_batch(key_path, **base)

