"""Reference-guided assembly script

This script performs a simple iterative reference-guided polishing loop:

- maps reads to a reference with `minimap2` + `samtools`
- calls variants (using `clair3`)
- applies variants to produce a consensus with `bcftools consensus`
- repeats until the consensus stabilizes or `--max-iterations` is reached
- masks low-coverage positions in the final consensus

The environment dependencies are expected to be installed via the provided
`env/env.yaml` conda environment.
"""

from __future__ import annotations

import argparse
import logging
import shutil
import subprocess
import sys
import tempfile
import time
import gzip
from pathlib import Path

from Bio import SeqIO
import pysam


def setup_logger(log_file: Path) -> logging.Logger:
    logger = logging.getLogger("rga")
    logger.setLevel(logging.DEBUG)
    fh = logging.FileHandler(log_file)
    fh.setLevel(logging.DEBUG)
    ch = logging.StreamHandler()
    ch.setLevel(logging.INFO)
    fmt = logging.Formatter("%(asctime)s %(levelname)s: %(message)s")
    fh.setFormatter(fmt)
    ch.setFormatter(fmt)
    logger.addHandler(fh)
    logger.addHandler(ch)
    return logger


def run(cmd, logger: logging.Logger, cwd: Path | None = None, check=True):
    logger.debug("Running: %s", " ".join(cmd))
    start = time.time()
    try:
        subprocess.run(cmd, check=check, cwd=cwd)
    except subprocess.CalledProcessError as e:
        logger.error("Command failed: %s", e)
        raise
    finally:
        logger.debug("Elapsed: %.2fs", time.time() - start)


def map_reads(reference: Path, reads: Path, out_bam: Path, threads: int, w: int, k: int, logger: logging.Logger):
    tmp_bam = out_bam.with_name(out_bam.stem + ".mapped.bam")
    cmd = [
        "minimap2",
        "-ax",
        "map-ont",
        str(reference),
        str(reads),
        "-t",
        str(threads),
        "--secondary=no",
        "-w",
        str(w),
        "-k",
        str(k),
    ]
    # pipe to samtools view
    with open(tmp_bam, "wb") as fh:
        logger.info("Mapping reads with minimap2 -> %s", tmp_bam)
        p1 = subprocess.Popen(cmd, stdout=subprocess.PIPE)
        p2 = subprocess.Popen(["samtools", "view", "-F", "4", "-b", "-"], stdin=p1.stdout, stdout=fh)
        p1.stdout.close()
        p2.communicate()
        if p2.returncode != 0:
            raise subprocess.CalledProcessError(p2.returncode, cmd)

    # sort and index
    logger.info("Sorting mapped alignments -> %s", out_bam)
    run(["samtools", "sort", "--threads", str(threads), "-o", str(out_bam), str(tmp_bam)], logger)
    run(["samtools", "index", str(out_bam)], logger)
    # cleanup
    try:
        tmp_bam.unlink()
    except Exception:
        pass


def cap_coverage(in_bam: Path, out_bam: Path, reference: Path, target_coverage: int, threads: int, logger: logging.Logger):
    if target_coverage is None or target_coverage <= 0:
        logger.info("No coverage capping requested; moving input to output")
        shutil.copy2(in_bam, out_bam)
        if in_bam.with_suffix(in_bam.suffix + ".bai").exists():
            shutil.copy2(in_bam.with_suffix(in_bam.suffix + ".bai"), out_bam.with_suffix(out_bam.suffix + ".bai"))
        return

    logger.info("Capping coverage to %d using jvarkit (if available)", target_coverage)
    # Ensure samtools dict exists for jvarkit
    run(["samtools", "dict", str(reference), "-o", str(reference.with_suffix(".dict"))], logger)

    # try to run jvarkit pipeline; fallback: leave as-is
    try:
        run(
            ["jvarkit", "sortsamrefname", str(in_bam), "-o", str(in_bam.with_name("refname_sorted.bam"))],
            logger
        )
        run([
            "jvarkit",
            "biostar154220",
            "-R",
            str(reference),
            "-n",
            str(target_coverage),
            str(in_bam.with_name("refname_sorted.bam")),
            "-o",
            str(out_bam.with_name("capped.bam")),
        ], logger)
        run(["samtools", "sort", "-o", str(out_bam), str(out_bam.with_name("capped.bam"))], logger)
        run(["samtools", "index", str(out_bam)], logger)
    except Exception:
        logger.warning("jvarkit not available or failed; keeping original BAM")
        shutil.copy2(in_bam, out_bam)


def get_mapped_reads(mapped_bam: Path, capped_reads: Path, logger: logging.Logger):
    """
    Get mapped reads from BAM and write to FASTQ.
    
    :param mapped_bam: Description
    :type mapped_bam: Path
    :param capped_reads: Description
    :type capped_reads: Path
    :param logger: Description
    :type logger: logging.Logger
    """
    logger.info("Extracting mapped reads from %s -> %s using samtools fastq", mapped_bam, capped_reads)
    # Exclude unmapped (4), secondary (256) and supplementary (2048) reads -> -F 2308
    # This uses samtools C implementation which is much faster for bulk extraction.
    cmd = ["samtools", "fastq", "-F", "2308", str(mapped_bam)]
    with open(capped_reads, "wt") as out:
        subprocess.run(cmd, check=True, stdout=out)


def call_variants_with_clair3(
        reference: Path,
        bam: Path,
        vcf_out: Path,
        tmp_clair3_outdir: Path,
        model: Path,
        min_coverage: int,
        chunk_size: int,
        threads: int,
        logger: logging.Logger
):
    """Call variants with Clair3.
    """
    logger.info("Calling variants with Clair3 -> %s", vcf_out)

    # Index the reference if not already indexed
    ref_fai = reference.with_suffix(reference.suffix + ".fai")
    if not ref_fai.exists():
        logger.info("Indexing reference with samtools faidx")
        run(["samtools", "faidx", str(reference)], logger)

    cmd = [
        "run_clair3.sh",
        f"--model_path={model}",
        f"--chunk_size={chunk_size}",
        f"--min_coverage={min_coverage}",
        f"--bam={bam}",
        f"--ref_fn={reference}",
        f"--threads={threads}",
        f"--output={tmp_clair3_outdir}",
        "--platform=ont",
        "--haploid_precise",
        "--enable_variant_calling_at_sequence_head_and_tail",
        "--include_all_ctgs",
        "--enable_long_indel",
        "--no_phasing_for_fa",
    ]

    try:
        run(cmd, logger)
    except Exception as e:
        logger.exception("Clair3 failed: %s", e)
        raise

    # clair_out/merge_output.vcf.gz to vcf_out
    with gzip.open(f"{tmp_clair3_outdir}/merge_output.vcf.gz", "rt") as f, vcf_out.open("wt") as outf:
        shutil.copyfileobj(f, outf)


def filter_vcf(vcf_in: Path, vcf_out: Path, min_depth: int, min_qual: float, logger: logging.Logger):
    """Filter VCF records using pysam based on DP and QUAL.

    The function will attempt to get depth from INFO/DP or sample-level DP.
    Records that do not meet both thresholds are skipped.
    """
    logger.info("Filtering VCF %s -> %s (min_depth=%s, min_qual=%s)", vcf_in, vcf_out, min_depth, min_qual)
    invcf = pysam.VariantFile(str(vcf_in))

    # Write compressed VCF for bcftools compatibility
    if str(vcf_out).endswith('.gz'):
        outvcf = pysam.VariantFile(str(vcf_out), mode="wz", header=invcf.header)
    else:
        outvcf = pysam.VariantFile(str(vcf_out), mode="w", header=invcf.header)

    written = 0
    for rec in invcf:
        qual = rec.qual if rec.qual is not None else 0.0
        dp = None
        if "DP" in rec.info:
            dp_val = rec.info.get("DP")
            if isinstance(dp_val, (list, tuple)):
                dp = int(dp_val[0])
            else:
                dp = int(dp_val)
        else:
            # try sample-level DP
            for s in rec.samples.values():
                if "DP" in s and s["DP"] is not None:
                    dp = int(s["DP"])
                    break

        if dp is None:
            dp = 0

        if (min_qual is None or qual >= float(min_qual)) and (min_depth is None or dp >= int(min_depth)):
            outvcf.write(rec)
            written += 1

    invcf.close()
    outvcf.close()
    logger.info("Wrote %d filtered records to %s", written, vcf_out)
    return written


def make_consensus(reference: Path, vcf: Path, out_fasta: Path, logger: logging.Logger):
    logger.info("Generating consensus with bcftools consensus -> %s", out_fasta)
    # If vcf is not compressed, compress and index it using bcftools (needs .vcf.gz + .csi/ .tbi index)
    vcf_gz = vcf
    if not str(vcf).endswith(".gz"):
        vcf_gz = vcf.with_suffix(str(vcf.suffix) + ".gz")
        logger.info("Compressing VCF to %s", vcf_gz)
        run(["bcftools", "view", "-Oz", "-o", str(vcf_gz), str(vcf)], logger)

    logger.info("Indexing VCF %s", vcf_gz)
    run(["bcftools", "index", "-f", str(vcf_gz)], logger)

    run(["bcftools", "consensus", "-f", str(reference), str(vcf_gz), "-o", str(out_fasta)], logger)


def mask_low_coverage(consensus_fasta: Path, bam: Path, min_depth_final: int, logger: logging.Logger):
    if min_depth_final is None or min_depth_final <= 0:
        logger.info("No final masking requested")
        return
    logger.info("Masking positions with depth < %d", min_depth_final)
    depth_lines = subprocess.check_output(["samtools", "depth", "-a", str(bam)], text=True)
    low_pos = {}
    for line in depth_lines.splitlines():
        chrom, pos, depth = line.split()
        pos = int(pos)
        depth = int(depth)
        if depth < min_depth_final:
            low_pos.setdefault(chrom, set()).add(pos)

    records = list(SeqIO.parse(str(consensus_fasta), "fasta"))
    out_records = []
    for rec in records:
        seq = list(str(rec.seq))
        chrom = rec.id
        if chrom in low_pos:
            for p in low_pos[chrom]:
                if 1 <= p <= len(seq):
                    seq[p - 1] = "N"
        rec.seq = type(rec.seq)("".join(seq))
        out_records.append(rec)

    SeqIO.write(out_records, str(consensus_fasta), "fasta")


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Reference-guided iterative polishing")
    p.add_argument("--reference", "-r", required=True, type=Path, help="Reference FASTA")
    p.add_argument("--reads", "-i", required=True, type=Path, help="Reads FASTQ (single-end) or BAM")
    p.add_argument("--log-file", type=Path, default=None, help="Path to log file (default: tmpdir/run.log)")
    p.add_argument("--clair3-model", type=str, required=True, help="Path to clair3 model or 'auto'")
    p.add_argument("--clair-chunk-size", type=int, default=2000, help="Clair3 chunk size for parallelization")
    p.add_argument("--cov-lim", type=int, default=50, help="Coverage cap (0 or None disables)")
    p.add_argument("--minimap-window-size", type=int, default=5)
    p.add_argument("--minimap-kmer-size", type=int, default=11)
    p.add_argument("--output", "-o", type=Path, default=Path("consensus_genome.fasta"))
    p.add_argument("--min-depth", type=int, default=2)
    p.add_argument("--min-quality", type=float, default=20.0)
    p.add_argument("--min-depth-final", type=int, default=10)
    p.add_argument("--max-iterations-inner", type=int, default=20)
    p.add_argument("--max-iterations-outer", type=int, default=5)
    p.add_argument("--threads", type=int, default=4)
    return p.parse_args()


def main():
    args = parse_args()
    tmpdir = Path(tempfile.mkdtemp(prefix="rga_"))
    log_file = args.log_file if args.log_file is not None else tmpdir / "run.log"
    if args.log_file is None:
        logger_path_info = f" (temporary log: {log_file})"
    else:
        logger_path_info = ""
    logger = setup_logger(log_file)
    logger.info("Working directory: %s%s", tmpdir, logger_path_info)

    try:
        threads = args.threads
        current_ref = Path(args.reference).resolve()

        for i in range(args.max_iterations_outer):
            logger.info("Outer iteration %d/%d", i+1, args.max_iterations_outer)

            # initial mapping
            tmp_mapped_bam = tmpdir / "mapped.sorted.bam"
            map_reads(
                current_ref,
                args.reads,
                tmp_mapped_bam,
                threads,
                args.minimap_window_size,
                args.minimap_kmer_size,
                logger
            )

            # optional coverage cap
            mapped_bam = tmpdir / "mapped.capped.bam"
            cap_coverage(tmp_mapped_bam, mapped_bam, current_ref, args.cov_lim, threads, logger)

            # get reads corresponding to mapped_bam for next steps; for now we just keep the bam and re-map in each iteration
            capped_reads = tmpdir / "reads.capped.fastq"
            get_mapped_reads(mapped_bam, capped_reads, logger)

            tmp_mapped_bam.unlink()

            # iterative polishing
            for j in range(args.max_iterations_inner):
                logger.info("Iteration %d-%d", i+1, j+1)

                it = f'{i+1}_{j+1}'

                iter_bam = tmpdir / f"iter_{it}.sorted.bam"
                if j == 0:
                    logger.info("Using initial mapping for iteration %d-%d", i+1, j+1)
                    # for the first iteration, we can use the existing mapped bam and current_ref
                    # but we need to ensure it's sorted and indexed as expected by Clair3
                    shutil.move(mapped_bam, iter_bam)
                    run(["samtools", "index", str(iter_bam)], logger)
                else:
                    # for subsequent iterations, we need to re-map reads to the current reference
                    logger.info("Re-mapping reads for iteration %d-%d", i+1, j+1)
                    map_reads(
                        current_ref,
                        args.reads,
                        iter_bam,
                        threads,
                        args.minimap_window_size,
                        args.minimap_kmer_size,
                        logger
                    )

                # call variants with Clair3
                raw_vcf = tmpdir / f"iter_{it}.raw.vcf"
                call_variants_with_clair3(
                    current_ref,
                    iter_bam,
                    raw_vcf,
                    tmpdir / f"clair3_out_{it}",
                    args.clair3_model,
                    args.min_depth,
                    args.clair_chunk_size,
                    threads,
                    logger
                )

                # filter
                filt_vcf = tmpdir / f"iter_{it}.filt.vcf.gz"
                n_variants = filter_vcf(raw_vcf, filt_vcf, args.min_depth, args.min_quality, logger)

                # make consensus
                consensus = tmpdir / f"consensus_iter_{it}.fa"
                make_consensus(current_ref, filt_vcf, consensus, logger)

                # compare
                if n_variants == 0:
                    logger.info(f"Consensus stabilized at iteration {it}")
                    break
                # else set current_ref to consensus and continue
                # write consensus to file used as reference for next iteration
                next_ref = tmpdir / f"current_ref_iter_{it}.fa"
                shutil.move(consensus, next_ref)
                current_ref = next_ref

                # delete temporary files from this iteration to save space
                for f in tmpdir.glob(f"iter_{it}.*"):
                    f.unlink()

            if j == 0 and n_variants == 0:
                logger.info(f"Consensus stabilized at iteration {i+1}-{j+1} with 0 variants; stopping outer loop")
                break

        # final masking
        final_consensus = args.output
        # ensure we have a final mapping bam
        final_bam = tmpdir / "final.sorted.bam"
        map_reads(current_ref, args.reads, final_bam, threads, args.minimap_window_size, args.minimap_kmer_size, logger)
        shutil.copy2(current_ref, final_consensus)
        mask_low_coverage(final_consensus, final_bam, args.min_depth_final, logger)

        # TODO: also mask low-quality variants in the final consensus? Clair3 provides a QUAL score for each variant which could be used to mask uncertain positions as 'N' in the final consensus.

        logger.info("Final consensus written to %s", final_consensus)
        logger.info("Log file: %s", log_file)

    except Exception as e:
        logger.exception("Pipeline failed: %s", e)
        sys.exit(1)


if __name__ == "__main__":
    main()
