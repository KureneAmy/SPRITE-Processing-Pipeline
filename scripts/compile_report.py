#!/usr/bin/env python3
"""
SPRITE Analysis Report Compiler
Generates HTML and Markdown reports from SPRITE pipeline outputs.

Usage:
    python compile_report.py --stats <stats_json> --html_template <html_tmpl>
                             --md_template <md_tmpl> --output_dir <dir>
                             [--samples <s1> <s2> ...]
"""

import argparse
import base64
import json
import logging
import os
import re
import sys
from datetime import datetime
from pathlib import Path

# Optional dependencies
try:
    import numpy as np
    HAS_NUMPY = True
except ImportError:
    HAS_NUMPY = False

try:
    import pandas as pd
    HAS_PANDAS = True
except ImportError:
    HAS_PANDAS = False

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s [%(levelname)s] %(message)s',
    datefmt='%Y-%m-%d %H:%M:%S'
)
logger = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Data collection helpers
# ---------------------------------------------------------------------------

def parse_ligation_efficiency(filepath):
    """Parse ligation_efficiency.txt and return a dict of stats."""
    stats = {}
    if not filepath or not os.path.isfile(filepath):
        logger.warning("Ligation efficiency file not found: %s", filepath)
        return stats
    try:
        with open(filepath) as fh:
            for line in fh:
                line = line.strip()
                if not line or line.startswith('#'):
                    continue
                # Common formats:
                #   KEY: VALUE
                #   KEY\tVALUE
                if ':' in line:
                    key, _, val = line.partition(':')
                    stats[key.strip()] = val.strip()
                elif '\t' in line:
                    parts = line.split('\t', 1)
                    stats[parts[0].strip()] = parts[1].strip()
    except OSError as exc:
        logger.error("Cannot read ligation efficiency file: %s", exc)
    return stats


def parse_cluster_file(filepath):
    """Parse a .DNA.clusters file and return basic statistics."""
    stats = {
        'total_clusters': 0,
        'size_distribution': {},
        'mean_size': 0.0,
        'median_size': 0.0,
    }
    if not filepath or not os.path.isfile(filepath):
        logger.warning("Cluster file not found: %s", filepath)
        return stats

    sizes = []
    try:
        with open(filepath) as fh:
            for line in fh:
                line = line.strip()
                if not line or line.startswith('#'):
                    continue
                # Each line is a cluster: tab-separated read IDs (or cluster descriptor)
                # cluster size = number of tab-separated fields
                parts = line.split('\t')
                size = len(parts)
                sizes.append(size)

        if sizes:
            stats['total_clusters'] = len(sizes)
            from collections import Counter
            dist = Counter(sizes)
            stats['size_distribution'] = {str(k): v for k, v in sorted(dist.items())}
            stats['mean_size'] = round(sum(sizes) / len(sizes), 2)
            sorted_sizes = sorted(sizes)
            mid = len(sorted_sizes) // 2
            if len(sorted_sizes) % 2 == 0:
                stats['median_size'] = (sorted_sizes[mid - 1] + sorted_sizes[mid]) / 2.0
            else:
                stats['median_size'] = float(sorted_sizes[mid])
    except OSError as exc:
        logger.error("Cannot read cluster file: %s", exc)
    return stats


def parse_heatmap_matrix(filepath):
    """Read a contact matrix text file and return basic stats."""
    stats = {
        'total_contacts': 0,
        'non_zero_bins': 0,
        'max_contact': 0.0,
        'mean_contact': 0.0,
    }
    if not filepath or not os.path.isfile(filepath):
        logger.warning("Heatmap matrix file not found: %s", filepath)
        return stats

    try:
        values = []
        with open(filepath) as fh:
            for line in fh:
                line = line.strip()
                if not line or line.startswith('#'):
                    continue
                for token in line.split():
                    try:
                        val = float(token)
                        if val > 0:
                            values.append(val)
                    except ValueError:
                        pass
        if values:
            stats['total_contacts'] = round(sum(values))
            stats['non_zero_bins'] = len(values)
            stats['max_contact'] = round(max(values), 2)
            stats['mean_contact'] = round(sum(values) / len(values), 4)
    except OSError as exc:
        logger.error("Cannot read heatmap matrix file: %s", exc)
    return stats


def image_to_data_uri(filepath, mime_type='image/png'):
    """Encode an image file as a base64 data URI."""
    if not filepath or not os.path.isfile(filepath):
        return ''
    try:
        with open(filepath, 'rb') as fh:
            data = base64.b64encode(fh.read()).decode('ascii')
        return f'data:{mime_type};base64,{data}'
    except OSError as exc:
        logger.warning("Cannot encode image %s: %s", filepath, exc)
        return ''


def collect_stats(args):
    """Collect all statistics from pipeline outputs and return a dict."""
    samples = args.samples or []
    out_dir = args.output_dir or ''

    logger.info("Collecting pipeline statistics for %d sample(s).", len(samples))

    report_data = {
        'report_title': args.report_title,
        'institution': args.institution,
        'pi_name': args.pi_name,
        'project_id': args.project_id,
        'analysis_date': args.analysis_date or datetime.now().strftime('%Y-%m-%d'),
        'assembly': args.assembly,
        'samples': samples,
        'num_samples': len(samples),
        'ligation_efficiency': {},
        'per_sample': {},
        'cluster_sizes_image': '',
        'heatmap_images': {},
    }

    # Ligation efficiency (single aggregated file)
    lig_file = args.ligation_efficiency
    report_data['ligation_efficiency'] = parse_ligation_efficiency(lig_file)

    # Cluster size distribution images
    cluster_sizes_png = args.cluster_sizes_png
    report_data['cluster_sizes_image'] = image_to_data_uri(cluster_sizes_png, 'image/png')

    # Per-sample data
    for sample in samples:
        sample_data = {'name': sample}

        # Cluster stats
        cluster_file = args.clusters_dir and os.path.join(
            args.clusters_dir, f'{sample}.DNA.clusters'
        )
        sample_data['cluster_stats'] = parse_cluster_file(cluster_file)

        # Heatmap matrix stats (raw and ICE-normalised)
        heatmap_dir = args.heatmap_dir or ''
        raw_file = os.path.join(heatmap_dir, f'{sample}.DNA.raw.txt') if heatmap_dir else None
        iced_file = os.path.join(heatmap_dir, f'{sample}.DNA.iced.txt') if heatmap_dir else None
        final_file = os.path.join(heatmap_dir, f'{sample}.DNA.final.txt') if heatmap_dir else None

        sample_data['raw_contacts'] = parse_heatmap_matrix(raw_file)
        sample_data['iced_contacts'] = parse_heatmap_matrix(iced_file)
        sample_data['final_contacts'] = parse_heatmap_matrix(final_file)

        # Heatmap PNG image
        heatmap_png = os.path.join(heatmap_dir, f'{sample}.DNA.final.png') if heatmap_dir else None
        sample_data['heatmap_image'] = image_to_data_uri(heatmap_png, 'image/png')

        report_data['per_sample'][sample] = sample_data
        logger.info("Collected data for sample: %s", sample)

    return report_data


# ---------------------------------------------------------------------------
# Report rendering
# ---------------------------------------------------------------------------

def _render_template(template_path, context):
    """Render a template file with the given context dict.

    Tries Jinja2 first; falls back to simple {{key}} substitution.
    """
    try:
        from jinja2 import Environment, FileSystemLoader, select_autoescape
        env = Environment(
            loader=FileSystemLoader(str(Path(template_path).parent)),
            autoescape=select_autoescape(['html']),
            keep_trailing_newline=True,
        )
        template = env.get_template(Path(template_path).name)
        return template.render(**context)
    except ImportError:
        logger.warning("Jinja2 not available; using basic template substitution.")
        with open(template_path) as fh:
            content = fh.read()
        for key, value in context.items():
            content = content.replace('{{' + key + '}}', str(value))
            content = content.replace('{{ ' + key + ' }}', str(value))
        return content


def build_html_context(report_data):
    """Build the Jinja2/substitution context for the HTML template."""
    ctx = dict(report_data)

    # Flatten ligation efficiency
    lig = report_data.get('ligation_efficiency', {})
    ctx['lig_efficiency_rows'] = ''.join(
        f'<tr><td>{k}</td><td>{v}</td></tr>' for k, v in lig.items()
    ) if lig else '<tr><td colspan="2"><em>No data available</em></td></tr>'

    # Sample summary table rows
    sample_rows = []
    for sample, data in report_data.get('per_sample', {}).items():
        cs = data.get('cluster_stats', {})
        rc = data.get('raw_contacts', {})
        row = (
            f'<tr>'
            f'<td>{sample}</td>'
            f'<td>{cs.get("total_clusters", "N/A")}</td>'
            f'<td>{cs.get("mean_size", "N/A")}</td>'
            f'<td>{cs.get("median_size", "N/A")}</td>'
            f'<td>{rc.get("total_contacts", "N/A")}</td>'
            f'</tr>'
        )
        sample_rows.append(row)
    ctx['sample_table_rows'] = ''.join(sample_rows) if sample_rows else (
        '<tr><td colspan="5"><em>No sample data available</em></td></tr>'
    )

    # Heatmap images for each sample
    heatmap_sections = []
    for sample, data in report_data.get('per_sample', {}).items():
        img_src = data.get('heatmap_image', '')
        if img_src:
            img_tag = (
                f'<figure>'
                f'<img src="{img_src}" alt="Contact heatmap for {sample}" '
                f'style="max-width:100%;height:auto;">'
                f'<figcaption>Contact heatmap: {sample}</figcaption>'
                f'</figure>'
            )
        else:
            img_tag = f'<p><em>Heatmap not available for {sample}</em></p>'
        heatmap_sections.append(img_tag)
    ctx['heatmap_sections'] = '\n'.join(heatmap_sections) if heatmap_sections else (
        '<p><em>No heatmap data available</em></p>'
    )

    # Cluster size image
    cluster_img_src = report_data.get('cluster_sizes_image', '')
    if cluster_img_src:
        ctx['cluster_sizes_figure'] = (
            f'<figure>'
            f'<img src="{cluster_img_src}" alt="Cluster size distribution" '
            f'style="max-width:100%;height:auto;">'
            f'<figcaption>Cluster size distribution across all samples</figcaption>'
            f'</figure>'
        )
    else:
        ctx['cluster_sizes_figure'] = '<p><em>Cluster size distribution image not available</em></p>'

    return ctx


def build_md_context(report_data):
    """Build the substitution context for the Markdown template."""
    ctx = dict(report_data)

    # Ligation efficiency table
    lig = report_data.get('ligation_efficiency', {})
    if lig:
        header = '| Metric | Value |\n|--------|-------|\n'
        rows = '\n'.join(f'| {k} | {v} |' for k, v in lig.items())
        ctx['ligation_efficiency_table'] = header + rows
    else:
        ctx['ligation_efficiency_table'] = '_No ligation efficiency data available._'

    # Sample summary table
    header = (
        '| Sample | Total Clusters | Mean Size | Median Size | Total Raw Contacts |\n'
        '|--------|---------------|-----------|-------------|--------------------|\n'
    )
    rows = []
    for sample, data in report_data.get('per_sample', {}).items():
        cs = data.get('cluster_stats', {})
        rc = data.get('raw_contacts', {})
        rows.append(
            f'| {sample} '
            f'| {cs.get("total_clusters", "N/A")} '
            f'| {cs.get("mean_size", "N/A")} '
            f'| {cs.get("median_size", "N/A")} '
            f'| {rc.get("total_contacts", "N/A")} |'
        )
    ctx['sample_summary_table'] = header + '\n'.join(rows) if rows else (
        '_No sample data available._'
    )

    ctx['num_samples'] = report_data.get('num_samples', 0)
    ctx['samples_list'] = ', '.join(report_data.get('samples', [])) or 'N/A'

    return ctx


def generate_report(report_data, html_template, md_template, output_dir):
    """Render and write HTML and Markdown reports."""
    Path(output_dir).mkdir(parents=True, exist_ok=True)

    html_out = os.path.join(output_dir, 'SPRITE_Analysis_Report.html')
    md_out = os.path.join(output_dir, 'SPRITE_Analysis_Report.md')
    json_out = os.path.join(output_dir, 'report_data.json')

    # Save raw stats as JSON
    # Strip embedded images from JSON to keep it readable
    json_safe = json.loads(json.dumps(report_data))
    json_safe['cluster_sizes_image'] = '<embedded>' if report_data.get('cluster_sizes_image') else ''
    for sample_data in json_safe.get('per_sample', {}).values():
        sample_data['heatmap_image'] = '<embedded>' if sample_data.get('heatmap_image') else ''
    with open(json_out, 'w') as fh:
        json.dump(json_safe, fh, indent=2)
    logger.info("Report statistics saved to %s", json_out)

    # HTML report
    if html_template and os.path.isfile(html_template):
        ctx = build_html_context(report_data)
        html_content = _render_template(html_template, ctx)
        with open(html_out, 'w', encoding='utf-8') as fh:
            fh.write(html_content)
        logger.info("HTML report written to %s", html_out)
    else:
        logger.warning("HTML template not found: %s", html_template)

    # Markdown report
    if md_template and os.path.isfile(md_template):
        ctx = build_md_context(report_data)
        md_content = _render_template(md_template, ctx)
        with open(md_out, 'w', encoding='utf-8') as fh:
            fh.write(md_content)
        logger.info("Markdown report written to %s", md_out)
    else:
        logger.warning("Markdown template not found: %s", md_template)

    return html_out, md_out, json_out


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def parse_args(argv=None):
    parser = argparse.ArgumentParser(
        description='Compile SPRITE analysis reports (HTML and Markdown).'
    )

    # Input data
    parser.add_argument('--samples', nargs='+', default=[],
                        help='Sample names to include in the report.')
    parser.add_argument('--ligation_efficiency', default=None,
                        help='Path to ligation_efficiency.txt.')
    parser.add_argument('--clusters_dir', default=None,
                        help='Directory containing .DNA.clusters files.')
    parser.add_argument('--heatmap_dir', default=None,
                        help='Directory containing heatmap matrix files.')
    parser.add_argument('--cluster_sizes_png', default=None,
                        help='Path to cluster_sizes.png.')
    parser.add_argument('--multiqc_report', default=None,
                        help='Path to multiqc_report.html (currently unused, reserved).')

    # Pre-computed stats JSON (skip collection step)
    parser.add_argument('--stats', default=None,
                        help='Path to a pre-computed report_data.json (skips data collection).')

    # Templates
    parser.add_argument('--html_template', default='templates/SPRITE_Report.html',
                        help='Path to the HTML Jinja2 template.')
    parser.add_argument('--md_template', default='templates/SPRITE_Report.md',
                        help='Path to the Markdown template.')

    # Output
    parser.add_argument('--output_dir', default='reports',
                        help='Output directory for generated reports.')

    # Report metadata
    parser.add_argument('--report_title', default='SPRITE Analysis Report',
                        help='Title for the report.')
    parser.add_argument('--institution', default='',
                        help='Institution name.')
    parser.add_argument('--pi_name', default='',
                        help='PI / investigator name.')
    parser.add_argument('--project_id', default='',
                        help='Project identifier.')
    parser.add_argument('--analysis_date', default='',
                        help='Analysis date (YYYY-MM-DD). Defaults to today.')
    parser.add_argument('--assembly', default='hg38',
                        help='Reference genome assembly (e.g. hg38, mm10).')

    return parser.parse_args(argv)


def main(argv=None):
    args = parse_args(argv)

    if args.stats and os.path.isfile(args.stats):
        logger.info("Loading pre-computed stats from %s", args.stats)
        with open(args.stats) as fh:
            report_data = json.load(fh)
    else:
        report_data = collect_stats(args)

    generate_report(
        report_data,
        html_template=args.html_template,
        md_template=args.md_template,
        output_dir=args.output_dir,
    )
    logger.info("Report generation complete.")


if __name__ == '__main__':
    main()
