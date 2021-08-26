#
# make an alignment chain if it doesn't exist
#
if "chain" not in config:

    rule make_bam:
        input:
            ref=REF,
            outgroup=OUTGROUP,
        output:
            bam=temp("results/chain/out-to-ref.bam"),
        log:
            "logs/chain/bam.log",
        conda:
            "../envs/env.yml"
        threads: 8
        shell:
            """
            minimap2 -ax asm20 -Y --eqx -t {threads} \
                {input.ref} {input.outgroup} \
                    | samtools view -u - \
                    | samtools sort -@ {threads} -m 8G - \
                > {output.bam}
            """

    rule make_psl:
        input:
            bam=rules.make_bam.output.bam,
        output:
            psl=temp("results/chain/out-to-ref.psl"),
        log:
            "logs/chain/psl.log",
        conda:
            "../envs/env.yml"
        shell:
            """
            bamToPsl {input.bam} {output.psl}  
            """

    rule make_chain:
        input:
            psl=rules.make_psl.output.psl,
        output:
            chain=CHAIN,
        log:
            "logs/chain/chain.log",
        conda:
            "../envs/env.yml"
        shell:
            """
            pslToChain {input.psl} {output.chain}
            """


#
#
#
rule prep_vcf:
    input:
        vcf=VCF,
    output:
        bcf=temp("results/vcf/input/{chrm}.bcf"),
        csi=temp("results/vcf/input/{chrm}.bcf.csi"),
    log:
        "logs/vcf/input/{chrm}.log",
    conda:
        "../envs/env.yml"
    shell:
        """
        bcftools view {input.vcf} {wildcards.chrm} \
            | bcftools sort -m 8G - \
            | bcftools +fill-tags \
            -Ob -o {output.bcf}

        bcftools index -f {output.bcf}
        """


rule prep_ref:
    input:
        reference=REF,
    output:
        ref=temp("results/ref-fasta/{chrm}.fa"),
        fai_ref=temp("results/ref-fasta/{chrm}.fa.fai"),
    log:
        "results/ref-fasta/{chrm}.log",
    conda:
        "../envs/env.yml"
    shell:
        """
        samtools faidx {input.reference} {wildcards.chrm} \
            | seqtk seq -l 60 > {output.ref}
        samtools faidx {output.ref}
        """


if "ancestor" in config:

    rule prep_ancestor:
        input:
            ancestor=ANCESTOR_FA,
        output:
            fasta=temp("results/ancestral-fasta/{chrm}.fa"),
            fai=temp("results/ancestral-fasta/{chrm}.fa.fai"),
        log:
            "results/ancestral-fasta/{chrm}.log",
        conda:
            "../envs/env.yml"
        shell:
            """
            samtools faidx {input.ancestor} {wildcards.chrm} \
                | seqtk seq -l 60 > {output.fasta}
            samtools faidx {output.fasta}
            """


else:

    rule make_ancestor_per_chr:
        input:
            bcf=rules.prep_vcf.output.bcf,
            ref=rules.prep_ref.output.ref,
            out=OUTGROUP,
            chain=CHAIN,
        output:
            fasta=temp("results/ancestral-fasta/{chrm}.fa"),
            fai=temp("results/ancestral-fasta/{chrm}.fa.fai"),
        log:
            "results/ancestral-fasta/{chrm}.log",
        conda:
            "../envs/env.yml"
        shell:
            """
            mutyper ancestor \
                {input.bcf} \
                {input.ref} \
                {input.out} \
                {input.chain} \
            {output.fasta}
            samtools faidx {output.fasta}
            """


rule annotate_vcf:
    input:
        fasta="results/ancestral-fasta/{chrm}.fa",
        bcf=rules.prep_vcf.output.bcf,
    output:
        bcf=temp("results/vcf/mutyper/{chrm}.bcf"),
        csi=temp("results/vcf/mutyper/{chrm}.bcf.csi"),
    log:
        "logs/vcf/mutyper/{chrm}.log",
    conda:
        "../envs/env.yml"
    shell:
        """
        mutyper variants {input.fasta} {input.bcf} \
            | bcftools sort -Ob -m 8G - \
            > {output.bcf}
        bcftools index -f {output.bcf}
        """


rule mutyper_vcf:
    input:
        bcf=expand(rules.annotate_vcf.output.bcf, chrm=CHRS),
        csi=expand(rules.annotate_vcf.output.csi, chrm=CHRS),
    output:
        bcf="results/vcf/mutyper.bcf",
        csi="results/vcf/mutyper.bcf.csi",
    log:
        "logs/mutyper_vcf.log",
    conda:
        "../envs/env.yml"
    shell:
        """
        bcftools concat -Ob -a \
            {input.bcf} > {output.bcf}
        bcftools index -f {output.bcf}
        """


rule mutyper_spectra:
    input:
        bcf=rules.mutyper_vcf.output.bcf,
    output:
        spectra="results/spectra/spectra.txt",
    log:
        "logs/spectra.log",
    conda:
        "../envs/env.yml"
    shell:
        """
        mutyper spectra {input.bcf} > {output.spectra}
        """


rule mutyper_spectra_stratify:
    input:
        bcf=rules.mutyper_vcf.output.bcf,
        bed=lambda wc: config["stratify"][wc.rgn],
    output:
        spectra="results/spectra/stratify/{rgn}_spectra.txt",
    log:
        "logs/spectra.{rgn}.log",
    conda:
        "../envs/env.yml"
    shell:
        """
        bcftools view --regions-file {input.bed} {input.bcf} \
            | mutyper spectra - \
            > {output.spectra}
        """


rule mutyper:
    input:
        rules.mutyper_vcf.output,
        rules.mutyper_spectra.output,
