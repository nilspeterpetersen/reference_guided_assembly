# General stuff
# implement the fucntions (e.g. function1 and function 2 and the call of clair), that require external tools
# with subprocesses. 
#
# write outputs to log files. report running times of the different steps in the logs.

""" function 1
minimap2 -ax map-ont ref.fasta trimmed.fastq \
    -t ${task.cpus} --secondary=no \
    -w ${params.map_to_target_minimap_window_size} \
    -k ${params.map_to_target_minimap_kmer_size} \
| samtools view -F 4 -b -o mapped.bam

samtools sort --threads ${task.cpus} -o sorted.bam mapped.bam
samtools index sorted.bam
"""


""" function 2
set -euo pipefail

# Ensure jvarkit can find java where its wrapper expects it
if [[ -n "\$CONDA_PREFIX" && ! -x "\$CONDA_PREFIX/bin/java" ]]
then
    ln -sf "\$(command -v java)" "\$CONDA_PREFIX/bin/java"
fi

target_coverage=${params.sv_cap_coverage}
if [[ \$target_coverage -gt 0 ]]
then
    samtools dict ref.fasta -o ref.dict

    jvarkit sortsamrefname in.bam -o refname_sorted.bam
    jvarkit biostar154220 -R ref.fasta -n \$target_coverage refname_sorted.bam -o capped.bam

    samtools sort capped.bam -o out.bam
    samtools index out.bam
else
    mv in.bam out.bam
    mv in.bam.bai out.bam.bai
fi
"""

"""Main

---- Input
- read arguments with argparse

    - arguments
    - required:
        - reference genome (fasta) -> required
        - read
        - clair3 model.
    - optional
        - cov-lim: coverage limitation -> None if 'None', default is 50
        - minimap2 w and k (defaults are --minimap-window-size=5 and --minimap-kmer-size=11)
        - output: name of consensus sequence file (default is consensus_genome.fasta)
        - min-depth: minimum depth for a base to be called during iterations
        - min-quality: of variant call for a base to be called during iterations
        - min-depth-final: minimum depth for a base to be called in final consensus (other wise masked with N)
        - min-quality-final: of variant call for final consensus (otherwise masked with N)
        - max-iterations: default is 20

----
- create a temp directory to write all intermediate files to!

---- First mapping
- map reads to reference (like function 1)
- sample reads to coverage threshold (keep all reads if this is set to 'None')

---- Iterative variant calling
- map reads to current_reference
- run clair3 with given model
- read variants (with pysam), filter to variants due to min-depth and min-quality
- apply variants to the reference genome to build a consensus
- if consensus != current_reference: current_reference = consensus and do another iteration if max-iterations is not reached
- else we have reached the final iteration

---- Finalization
- mask all positions that are below coverage and quality threshold

"""