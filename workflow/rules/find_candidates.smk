rule find_breakends_gridss:
    input:
        alignment="resources/mason_simulated/alignment_sorted.bam",
        alignment_index="resources/mason_simulated/alignment_sorted.bam.bai",
        genome="resources/genome.fasta",
        genome_index="resources/genome.fasta.fai",
    output:
        breakends="results/candidates/breakends.vcf",
        gridss_assembly="results/candidates/gridss_assembly.bam",
    conda:
        "../envs/gridss.yaml"
    shell:
        """
        gridss --reference {input.genome} --assembly {output.gridss_assembly} --output {output.breakends} {input.alignment}
        """
