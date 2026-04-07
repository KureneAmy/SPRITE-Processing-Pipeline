#!/usr/bin/env python3
"""
compile_report.py - SHAPE-MaP Analysis Report Generator

Parses ShapeMapper2 output files and generates professional HTML and
Markdown reports using Jinja2 templates.

Usage:
    python compile_report.py --config config.yaml --output-dir /path/to/output
"""

import argparse
import json
import logging
import os
import re
import statistics as _stats
import sys
from datetime import datetime

import yaml

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
)
logger = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _safe_float(value, default=None):
    """Convert value to float; return *default* on failure."""
    try:
        return float(value)
    except (TypeError, ValueError):
        return default


def _safe_int(value, default=None):
    """Convert value to int; return *default* on failure."""
    try:
        return int(value)
    except (TypeError, ValueError):
        return default


# ---------------------------------------------------------------------------
# Parsers
# ---------------------------------------------------------------------------

def parse_shapemapper_log(log_path):
    """
    Parse a ShapeMapper2 log file and return a dictionary of statistics.

    The log embeds output from BBmerge (paired-end merging) and Bowtie2
    (alignment), as well as ShapeMapper's own QC summary.

    Returns a dict with keys such as:
        paired_reads_total, paired_reads_merged, merge_rate,
        alignment_rate, mapped_reads_concordant_1,
        mutation_rate_modified, mutation_rate_untreated,
        mutation_rate_denatured, quality_filtered_reads,
        mean_read_depth, high_quality_profile_pct,
        shapemapper_version, bowtie2_version, bbmerge_version
    """
    stats = {}
    if not os.path.isfile(log_path):
        logger.warning("Log file not found: %s", log_path)
        return stats

    with open(log_path, "r", encoding="utf-8", errors="replace") as fh:
        content = fh.read()

    # ------------------------------------------------------------------
    # Software versions
    # ------------------------------------------------------------------
    # ShapeMapper2: "Started ShapeMapper v2.3.1 at ..."
    m = re.search(r"Started ShapeMapper\s+v?([0-9][^\s]+)", content, re.IGNORECASE)
    if m:
        stats["shapemapper_version"] = m.group(1)

    # Bowtie2 within the log: "bowtie2-align-s version 2.3.4.3"
    m = re.search(r"bowtie2[^\n]*?\bversion\s+([0-9][^\s]+)", content, re.IGNORECASE)
    if m:
        stats["bowtie2_version"] = m.group(1)

    # BBMerge: "BBMerge version 37.78" or "Bbmerge version X"
    m = re.search(r"BBMerge\s+version\s+([0-9][^\s]+)", content, re.IGNORECASE)
    if m:
        stats["bbmerge_version"] = m.group(1)

    # ------------------------------------------------------------------
    # BBmerge paired-end merging statistics
    # Format written to log by BBmerge:
    #   Pairs:           50000
    #   Joined:          42000      84.00000%
    # ------------------------------------------------------------------
    m = re.search(r"^Pairs:\s+([0-9,]+)", content, re.MULTILINE)
    if m:
        stats["paired_reads_total"] = int(m.group(1).replace(",", ""))

    m = re.search(r"^Joined:\s+([0-9,]+)", content, re.MULTILINE)
    if m:
        stats["paired_reads_merged"] = int(m.group(1).replace(",", ""))

    # Derive merge rate from the percentage on the Joined line if present
    m = re.search(r"^Joined:\s+[0-9,]+\s+([0-9.]+)%", content, re.MULTILINE)
    if m:
        stats["merge_rate"] = _safe_float(m.group(1))

    # Fallback: old-style "Paired reads merged: N / N"
    if "paired_reads_total" not in stats:
        m = re.search(r"Paired reads merged[:\s]+([0-9,]+)\s*/\s*([0-9,]+)", content)
        if m:
            stats["paired_reads_merged"] = int(m.group(1).replace(",", ""))
            stats["paired_reads_total"] = int(m.group(2).replace(",", ""))

    # ------------------------------------------------------------------
    # Bowtie2 alignment statistics
    # "97.62% overall alignment rate"  (last occurrence in file)
    # ------------------------------------------------------------------
    for m in re.finditer(r"([0-9.]+)%\s+overall alignment rate", content):
        stats["alignment_rate"] = _safe_float(m.group(1))

    # Concordantly mapped reads: "9800 (93.33%) aligned exactly 1 time"
    m = re.search(
        r"([0-9,]+)\s+\([0-9.]+%\)\s+aligned concordantly exactly 1 time",
        content,
    )
    if m:
        stats["mapped_reads_concordant_1"] = int(m.group(1).replace(",", ""))

    # ------------------------------------------------------------------
    # ShapeMapper2 QC summary — mutation rates
    # Typical formats found in the log:
    #   "Overall mutation rate in modified sample:   0.0483"
    #   "Modified channel overall mutation rate:     0.0483"
    #   "mutation rate  modified:  0.0483"
    # ------------------------------------------------------------------
    for label, key in [
        ("modified", "mutation_rate_modified"),
        ("untreated", "mutation_rate_untreated"),
        ("denatured", "mutation_rate_denatured"),
    ]:
        patterns = [
            # "Overall mutation rate in modified sample: 0.0483"
            rf"overall mutation rate in {label} sample[:\s]+([0-9.eE+\-]+)",
            # "Modified channel overall mutation rate: 0.0483"
            rf"{label} channel.*?overall mutation rate[:\s]+([0-9.eE+\-]+)",
            # "mutation rate  modified: 0.0483" or "modified: mutation rate: 0.0483"
            rf"mutation rate\s+{label}[:\s]+([0-9.eE+\-]+)",
            rf"{label}[:\s]+mutation rate[:\s]+([0-9.eE+\-]+)",
            # Generic: "modified ... 0.0483" or "modified ... 1.2e-4" inside a QC block
            rf"^\s*{label}[^\n]*?([0-9]+(?:\.[0-9]+)?(?:[eE][+\-]?[0-9]+)?)",
        ]
        for pat in patterns:
            m = re.search(pat, content, re.IGNORECASE | re.MULTILINE)
            if m:
                val = _safe_float(m.group(1))
                if val is not None:
                    stats[key] = val
                    break

    # ------------------------------------------------------------------
    # Quality-filtered reads
    # "500 reads quality-filtered" or "quality filtered: 500"
    # ------------------------------------------------------------------
    m = re.search(r"([0-9,]+)\s+reads?\s+quality[\s\-]filtered", content, re.IGNORECASE)
    if m:
        stats["quality_filtered_reads"] = int(m.group(1).replace(",", ""))
    else:
        m = re.search(r"quality[\s\-]filtered[:\s]+([0-9,]+)", content, re.IGNORECASE)
        if m:
            stats["quality_filtered_reads"] = int(m.group(1).replace(",", ""))

    # ------------------------------------------------------------------
    # Mean read depth
    # "Median nucleotide read depth: 50000" or "Mean read depth: 50000"
    # ------------------------------------------------------------------
    m = re.search(r"(?:median|mean)\s+(?:nucleotide\s+)?read\s+depth[:\s]+([0-9.]+)", content, re.IGNORECASE)
    if m:
        stats["mean_read_depth"] = _safe_float(m.group(1))

    # ------------------------------------------------------------------
    # High-quality profile percentage
    # "90.0% of nucleotides ... high-quality"
    # ------------------------------------------------------------------
    m = re.search(r"([0-9.]+)\s*%.*?nucleotides?.*?high.quality", content, re.IGNORECASE)
    if m:
        stats["high_quality_profile_pct"] = _safe_float(m.group(1))

    logger.debug("Parsed log %s: %s", log_path, stats)
    return stats


def parse_profile_txt(profile_txt_path):
    """
    Parse a ShapeMapper2 ``_profile.txt`` file and return per-channel
    read depth, mutation counts, and reactivity statistics.

    The file is tab-delimited with a header row.  Columns of interest:
      ``<Sample>_mutations``, ``<Sample>_read_depth``, ``<Sample>_rate``,
      ``Reactivity_profile``, ``Norm_profile``, ``HQ_profile``

    where <Sample> is ``Modified``, ``Untreated``, or ``Denatured``.
    """
    stats = {}
    if not os.path.isfile(profile_txt_path):
        logger.debug("Profile txt not found: %s", profile_txt_path)
        return stats

    try:
        with open(profile_txt_path, "r", encoding="utf-8") as fh:
            header_line = fh.readline()
            if not header_line.strip():
                return stats
            headers = [h.strip() for h in header_line.split("\t")]

            # Build column index map (case-insensitive)
            col_idx = {h.lower(): i for i, h in enumerate(headers)}

            rows = []
            for line in fh:
                line = line.strip()
                if not line:
                    continue
                rows.append(line.split("\t"))
    except OSError as exc:
        logger.warning("Cannot read profile txt %s: %s", profile_txt_path, exc)
        return stats

    if not rows:
        return stats

    def _col_values(name_lower):
        """Return list of valid float values for a column name."""
        idx = col_idx.get(name_lower)
        if idx is None:
            return []
        vals = []
        for row in rows:
            if idx < len(row):
                v = _safe_float(row[idx])
                if v is not None and v != -999.0:
                    vals.append(v)
        return vals

    def _col_sum(name_lower):
        vals = _col_values(name_lower)
        return int(sum(vals)) if vals else None

    # Per-channel statistics
    for label in ("modified", "untreated", "denatured"):
        mutations = _col_sum(f"{label}_mutations")
        read_depth_vals = _col_values(f"{label}_read_depth")
        read_depth = int(sum(read_depth_vals)) if read_depth_vals else None
        rates = _col_values(f"{label}_rate")

        if mutations is not None:
            stats[f"total_mutations_{label}"] = mutations
        if read_depth is not None:
            stats[f"total_read_depth_{label}"] = read_depth
        if rates:
            stats[f"mutation_rate_{label}"] = round(_stats.median(rates), 5)
        if read_depth_vals:
            stats[f"mean_read_depth_{label}"] = round(
                sum(read_depth_vals) / len(read_depth_vals), 1
            )

    # Reactivity statistics from Norm_profile column (normalised reactivities)
    norm_reactivities = _col_values("norm_profile")
    if not norm_reactivities:
        norm_reactivities = _col_values("reactivity_profile")
    if norm_reactivities:
        stats["nucleotide_count"] = len(norm_reactivities)
        stats["reactivity_mean"] = round(_stats.mean(norm_reactivities), 4)
        stats["reactivity_median"] = round(_stats.median(norm_reactivities), 4)
        stats["reactivity_stdev"] = (
            round(_stats.stdev(norm_reactivities), 4)
            if len(norm_reactivities) > 1 else 0.0
        )
        stats["reactivity_max"] = round(max(norm_reactivities), 4)
        reactive = sum(1 for r in norm_reactivities if r > 0.4)
        stats["reactive_fraction"] = round(reactive / len(norm_reactivities), 4)

    # HQ profile nucleotide count
    hq_vals = _col_values("hq_profile")
    if hq_vals:
        stats["hq_nucleotide_count"] = len(hq_vals)
        total = len(rows)
        stats["high_quality_profile_pct"] = round(100.0 * len(hq_vals) / total, 1) if total else None

    return stats


def parse_shape_file(shape_path):
    """
    Parse a ShapeMapper2 ``.shape`` file (2-column fallback).

    The .shape file has two columns: nucleotide index and normalised
    SHAPE reactivity.  Excluded positions are set to -999.

    Returns a subset of the reactivity statistics also produced by
    :func:`parse_profile_txt`.
    """
    stats = {}
    if not os.path.isfile(shape_path):
        logger.debug("Shape file not found: %s", shape_path)
        return stats

    reactivities = []
    try:
        with open(shape_path, "r", encoding="utf-8") as fh:
            for line in fh:
                line = line.strip()
                if not line or line.startswith("#"):
                    continue
                parts = line.split()
                if len(parts) >= 2:
                    val = _safe_float(parts[1])
                    if val is not None and val != -999.0:
                        reactivities.append(val)
    except OSError as exc:
        logger.warning("Cannot read shape file %s: %s", shape_path, exc)
        return stats

    if reactivities:
        stats["nucleotide_count"] = len(reactivities)
        stats["reactivity_mean"] = round(_stats.mean(reactivities), 4)
        stats["reactivity_median"] = round(_stats.median(reactivities), 4)
        stats["reactivity_stdev"] = (
            round(_stats.stdev(reactivities), 4) if len(reactivities) > 1 else 0.0
        )
        stats["reactivity_max"] = round(max(reactivities), 4)
        reactive = sum(1 for r in reactivities if r > 0.4)
        stats["reactive_fraction"] = round(reactive / len(reactivities), 4)

    return stats


def collect_all_stats(config):
    """
    Iterate over all samples and targets in *config* and collect statistics
    from log and profile files.

    Returns a dict suitable for JSON serialisation and template rendering.
    """
    output_dir = config["output_dir"]
    samples = config.get("samples", {})
    report_cfg = config.get("report_config", {})

    data = {
        "report_title": report_cfg.get("report_title", "SHAPE-MaP Analysis Report"),
        "institution": report_cfg.get("institution", ""),
        "pi_name": report_cfg.get("pi_name", ""),
        "project_id": report_cfg.get("project_id", ""),
        "analysis_date": report_cfg.get("analysis_date") or datetime.now().strftime("%Y-%m-%d"),
        "report_date": datetime.now().strftime("%Y-%m-%d"),
        "sample_count": len(samples),
        "samples": {},
        "multiqc_report": os.path.join(output_dir, "multiqc_report.html"),
        "output_dir": output_dir,
    }

    # Derive sequence/target names by scanning output directories that
    # correspond to each sample.  Fall back to an empty list if the
    # directory does not yet exist (dry-run scenario).
    for sample_name in samples:
        sample_dir = os.path.join(output_dir, sample_name)
        seq_names = []
        if os.path.isdir(sample_dir):
            seq_names = [
                d
                for d in os.listdir(sample_dir)
                if os.path.isdir(os.path.join(sample_dir, d))
            ]
        else:
            logger.warning("Sample directory not found: %s", sample_dir)

        sample_data = {"targets": {}}

        for seq_name in seq_names:
            target_dir = os.path.join(sample_dir, seq_name)
            log_path = os.path.join(
                target_dir, f"{sample_name}_{seq_name}_shapemapper_log.txt"
            )
            # Primary data source: _profile.txt (multi-column, has per-channel rates)
            profile_txt_path = os.path.join(
                target_dir, f"{sample_name}_{seq_name}_profile.txt"
            )
            # Fallback: .shape file (2-column normalised reactivities)
            shape_path = os.path.join(
                target_dir, f"{sample_name}_{seq_name}.shape"
            )
            profile_pdf = os.path.join(
                target_dir, f"{sample_name}_{seq_name}_profiles.pdf"
            )
            structure_svg = os.path.join(
                target_dir, "structure", f"{sample_name}_{seq_name}_folding.svg"
            )

            log_stats = parse_shapemapper_log(log_path)

            # Try profile.txt first; fall back to .shape for reactivity stats
            profile_stats = parse_profile_txt(profile_txt_path)
            if not profile_stats.get("nucleotide_count"):
                shape_stats = parse_shape_file(shape_path)
                # Merge: profile_stats entries take priority
                for k, v in shape_stats.items():
                    profile_stats.setdefault(k, v)

            # If mutation rates weren't found in log, promote from profile.txt
            for label in ("modified", "untreated", "denatured"):
                key = f"mutation_rate_{label}"
                if key not in log_stats and key in profile_stats:
                    log_stats[key] = profile_stats[key]

            # Similarly promote mean_read_depth from profile.txt if missing
            if "mean_read_depth" not in log_stats:
                for label in ("modified", "untreated", "denatured"):
                    depth_key = f"mean_read_depth_{label}"
                    if depth_key in profile_stats:
                        log_stats.setdefault("mean_read_depth", profile_stats[depth_key])
                        break

            sample_data["targets"][seq_name] = {
                "log_path": log_path,
                "shape_path": shape_path,
                "profile_txt_path": profile_txt_path,
                "profile_pdf": profile_pdf if os.path.isfile(profile_pdf) else None,
                "structure_svg": structure_svg if os.path.isfile(structure_svg) else None,
                "log_stats": log_stats,
                "profile_stats": profile_stats,
            }

        data["samples"][sample_name] = sample_data

    # Collect software versions from log stats (first occurrence wins)
    sw_versions = {
        "shapemapper": None,
        "bowtie2": None,
        "bbmerge": None,
    }
    for sd in data["samples"].values():
        for tgt in sd["targets"].values():
            ls = tgt["log_stats"]
            sw_versions["shapemapper"] = sw_versions["shapemapper"] or ls.get("shapemapper_version")
            sw_versions["bowtie2"] = sw_versions["bowtie2"] or ls.get("bowtie2_version")
            sw_versions["bbmerge"] = sw_versions["bbmerge"] or ls.get("bbmerge_version")
    data["software_versions"] = sw_versions

    return data


# ---------------------------------------------------------------------------
# Template rendering
# ---------------------------------------------------------------------------

def _fmt_number(value, decimals=2):
    """Format a numeric value for display in reports.

    Returns the formatted string or 'N/A' if *value* is None or invalid.
    """
    if value is None:
        return "N/A"
    try:
        return f"{float(value):,.{decimals}f}"
    except (TypeError, ValueError):
        return str(value)


def render_html_report(data, template_path, output_path):
    """Render the HTML report using Jinja2."""
    try:
        from jinja2 import Environment, FileSystemLoader, select_autoescape
    except ImportError:
        logger.error(
            "Jinja2 is required for HTML report generation.  "
            "Install it with: pip install jinja2"
        )
        sys.exit(1)

    template_dir = os.path.dirname(os.path.abspath(template_path))
    template_name = os.path.basename(template_path)

    env = Environment(
        loader=FileSystemLoader(template_dir),
        autoescape=select_autoescape(["html"]),
        keep_trailing_newline=True,
    )
    env.filters["fmt_number"] = _fmt_number

    template = env.get_template(template_name)
    html_content = template.render(**data)

    os.makedirs(os.path.dirname(os.path.abspath(output_path)), exist_ok=True)
    with open(output_path, "w", encoding="utf-8") as fh:
        fh.write(html_content)
    logger.info("HTML report written to: %s", output_path)


def render_markdown_report(data, template_path, output_path):
    """Render the Markdown report using Jinja2."""
    try:
        from jinja2 import Environment, FileSystemLoader
    except ImportError:
        logger.error(
            "Jinja2 is required for Markdown report generation.  "
            "Install it with: pip install jinja2"
        )
        sys.exit(1)

    template_dir = os.path.dirname(os.path.abspath(template_path))
    template_name = os.path.basename(template_path)

    env = Environment(
        loader=FileSystemLoader(template_dir),
        autoescape=False,
        keep_trailing_newline=True,
    )
    env.filters["fmt_number"] = _fmt_number

    template = env.get_template(template_name)
    md_content = template.render(**data)

    os.makedirs(os.path.dirname(os.path.abspath(output_path)), exist_ok=True)
    with open(output_path, "w", encoding="utf-8") as fh:
        fh.write(md_content)
    logger.info("Markdown report written to: %s", output_path)


# ---------------------------------------------------------------------------
# CLI entry point
# ---------------------------------------------------------------------------

def parse_args():
    parser = argparse.ArgumentParser(
        description="Generate SHAPE-MaP analysis reports (HTML and Markdown).",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "--config",
        required=True,
        help="Path to the Snakemake config.yaml file.",
    )
    parser.add_argument(
        "--output-dir",
        default=None,
        help=(
            "Output directory for reports.  Defaults to "
            "<config output_dir>/reports."
        ),
    )
    parser.add_argument(
        "--html-template",
        default=None,
        help="Path to the HTML Jinja2 template.",
    )
    parser.add_argument(
        "--md-template",
        default=None,
        help="Path to the Markdown Jinja2 template.",
    )
    parser.add_argument(
        "--stats-json",
        default=None,
        help=(
            "Path to write the collected statistics JSON file.  "
            "Defaults to <output-dir>/report_data.json."
        ),
    )
    parser.add_argument(
        "--no-html",
        action="store_true",
        help="Skip HTML report generation.",
    )
    parser.add_argument(
        "--no-markdown",
        action="store_true",
        help="Skip Markdown report generation.",
    )
    parser.add_argument(
        "--verbose",
        action="store_true",
        help="Enable debug logging.",
    )
    return parser.parse_args()


def main():
    args = parse_args()

    if args.verbose:
        logging.getLogger().setLevel(logging.DEBUG)

    # Load config
    config_path = os.path.abspath(args.config)
    if not os.path.isfile(config_path):
        logger.error("Config file not found: %s", config_path)
        sys.exit(1)

    with open(config_path, "r", encoding="utf-8") as fh:
        config = yaml.safe_load(fh)

    # Resolve paths relative to the config file's directory when they are
    # relative (to stay consistent with how Snakemake resolves them)
    config_dir = os.path.dirname(config_path)
    report_cfg = config.get("report_config", {})

    # Output directory
    if args.output_dir:
        reports_dir = os.path.abspath(args.output_dir)
    else:
        reports_dir = os.path.join(config["output_dir"], "reports")

    # Template paths
    script_dir = os.path.dirname(os.path.abspath(__file__))
    templates_dir = os.path.join(os.path.dirname(script_dir), "templates")

    html_template = args.html_template or os.path.join(
        templates_dir, "SHAPE-MaP_Report.html"
    )
    md_template = args.md_template or os.path.join(
        templates_dir, "SHAPE-MaP_Report.md"
    )

    # Collect statistics
    logger.info("Collecting statistics from ShapeMapper outputs …")
    stats = collect_all_stats(config)

    # Write JSON
    stats_json_path = args.stats_json or os.path.join(reports_dir, "report_data.json")
    os.makedirs(os.path.dirname(os.path.abspath(stats_json_path)), exist_ok=True)
    with open(stats_json_path, "w", encoding="utf-8") as fh:
        json.dump(stats, fh, indent=2, ensure_ascii=False)
    logger.info("Statistics JSON written to: %s", stats_json_path)

    generate_html = report_cfg.get("generate_html", True) and not args.no_html
    generate_md = report_cfg.get("generate_markdown", True) and not args.no_markdown

    # Generate HTML report
    if generate_html:
        if not os.path.isfile(html_template):
            logger.error("HTML template not found: %s", html_template)
            sys.exit(1)
        html_output = os.path.join(reports_dir, "SHAPE-MaP_Analysis_Report.html")
        render_html_report(stats, html_template, html_output)
    else:
        logger.info("HTML report generation skipped.")

    # Generate Markdown report
    if generate_md:
        if not os.path.isfile(md_template):
            logger.error("Markdown template not found: %s", md_template)
            sys.exit(1)
        md_output = os.path.join(reports_dir, "SHAPE-MaP_Analysis_Report.md")
        render_markdown_report(stats, md_template, md_output)
    else:
        logger.info("Markdown report generation skipped.")

    logger.info("Report generation complete.")


if __name__ == "__main__":
    main()
