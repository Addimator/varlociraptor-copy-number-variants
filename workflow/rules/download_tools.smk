rule download_varlociraptor:
    output:
        directory("resources/tools/varlociraptor"),
    # conda:
    #     "../envs/gridss.yaml"
    shell:
        """ 
        cd resources/tools/
        git clone git@github.com:varlociraptor/varlociraptor.git
        cd varlociraptor
        git checkout copy-number-variation
        """
