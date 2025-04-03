# TODO: Is there really no easy way to remove the temporary files created by gridss?
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
        rm -rf alignment_sorted.bam.gridss.working
        rm -rf breakends.vcf.gridss.working
        rm -rf gridss_assembly.bam.gridss.working
        rm gridss*
        rm -f libsswjni.so
        """


rule breakends_to_candidates:
    input:
        breakends="results/candidates/breakends.vcf",
        varlo=directory("resources/tools/varlociraptor"),
    output:
        candidates="results/candidates/candidates.vcf",
    conda:
        "../envs/gridss.yaml"
    shell:
        """ 
        cd {input.varlo}
        cargo run -- cnv-candidates {input.breakends} {output}
        """
