configfile: "config.yaml"

rule all:
    input:
        expand("mapping/bowtie2/{mapping_params}/{reference}/units/{unit}.bam", \
                mapping_params = config["bowtie2_rules"]["mapping_params"], reference=config["input_files"]["references"], \
                unit=config["bowtie2_rules"]["units"])

rule extract_genomes_from_archive:
    input:
        expand("input_files/{reference}.tsv", reference=config["input_files"]["references"])
    output:
        "extracted_genomes/{reference}.tar"
    params:
        archive = config["input_files"]["archive"],
        output_dir=temp("genomes_tmp/{reference}"),
    log:
        "logs/extract_genomes_from_archive.log"
    #message: "Extract genomes from {input.archive} saving to {output.genomes}"
    shell:
        "mkdir -p {params.output_dir} ;"
        "cut -f1 {input} > {params.output_dir}/{wildcards.reference}_list ;"
        "tar -tf {params.archive} | grep -f {params.output_dir}/{wildcards.reference}_list > {params.output_dir}/{wildcards.reference}_paths ;"
        "tar -xvf {params.archive} -C {params.output_dir} -T {params.output_dir}/{wildcards.reference}_paths --strip-components=8 && "
        "tar -cvf {output} {params.output_dir}/*.gz"

rule untar_genomes:
    input:
        "extracted_genomes/{reference}.tar"
        #"extracted_genomes/{reference}.tar"
    output:
        dynamic("selected_genomes/{reference}/{n}.dna.toplevel.fa.gz")
    params:
        prefix="selected_genomes"
    shell:
        "mkdir -p {params.prefix};"
        "tar -xvf {input} -C {params.prefix} --strip-components=1 "


rule change_fasta:
    input:
        "selected_genomes/{reference}/{n}.dna.toplevel.fa.gz"
    output:
        "reference_genomes/{reference}/{n}.fasta.gz"
    message:
        "Converting fasta"
    log:
        "logs/change_fasta.log"
    run:
        import os
        import gzip
        from Bio import SeqIO
        import re
        import shutil
        import logging
        import glob

        logging.basicConfig(filename=log[0],
                            format='%(asctime)s %(levelname)-8s %(message)s',
                            level=logging.INFO,
                            datefmt='%Y-%m-%d %H:%M:%S')

        converted_count = 0
        #for input_path in input_paths:
        input_path = input[0]
        filename = os.path.basename(input_path)
        genome = re.sub(r'\.dna\.toplevel\.fa\.gz$', '', filename)
        parts = genome.split(".")
        assembly_id = parts[1]

        with gzip.open(input_path, "rt") as f_in:
            sequences = []
            for record in SeqIO.parse(f_in, "fasta"):
                record.id = "{}|{}".format(assembly_id, record.id)
                record.description = ""
                sequences.append(record)
        #output_file = os.path.join(os.path.normpath(output[0]).split(os.sep)[0],"{}.fasta".format(genome))

        output_file = os.path.join(os.path.dirname(output[0]),"{}.fasta".format(genome))
        print (output_file)
        with open(output_file,'w') as f_out:
            SeqIO.write(sequences, f_out, "fasta")
        f_out.close()

        try:
            if os.path.getsize(output_file) > 0:
                zipped_output_path = "{}.gz".format(output_file)
                with open(output_file, "rb") as f_in, gzip.open(zipped_output_path, 'wb') as f_out:
                    f_out.writelines(f_in)
                f_out.close()
                os.remove(output_file) #removing temporaty fasta
                logging.info('Finished, succesfully converted {}'.format(input_path))
            else:
                logging.warn('{} is empty file!'.format(output_file))
        except OSError as e:
            logging.error('{} does not exist, cannot convert!'.format(output_file))

rule concatenate_fasta:
    input:
        dynamic("reference_genomes/{reference}/{n}.fasta.gz")
    output:
        "references/{reference}.fasta.gz"
    params:
        prefix = "reference_genomes/{reference}"
    shell:
        #"mkdir -p {params.tmp_dir} &&"
        #"cut -f1 {input} | sed s/.dna.toplevel.fa.gz/.fasta.gz/ > {params.tmp_dir}/{wildcards.reference}.f1  &&"
        #"zcat `ls {params.input_folder}/*.fasta.gz | grep -f {params.tmp_dir}/{wildcards.reference}.f1` | gzip > {output}"
        "gunzip -c `ls {params.prefix}/*.fasta.gz` | gzip > {output}"

rule build_indexes:
    input:
        "references/{reference}.fasta.gz"
    output:
        expand("indices/{{reference}}.{index}.bt2l", index=range(1,5)),
        expand("indices/{{reference}}.rev.{index}.bt2l", index=range(1,3))
    params:
        prefix="indices/{reference}"
    # message:
    #     "Building indes for: {wildcard.reference}"
    log:
        "logs/build_indexes.log"
    threads: 20
    shell:
        "bowtie2-build --large-index --threads {threads} {input} {params.prefix}"


UNIT_TO_SAMPLE = {
     unit: sample for sample, units in config["bowtie2_rules"]["samples"].items()
     for unit in units}

#
def create_bowtie2_read_input_str(unit):
    if len(unit) == 2:
        return "-1 {unit[0]} -2 {unit[1]}".format(unit=unit)
    elif len(unit) == 1:
        return "-r {unit[0]}".format(unit=unit)
    else:
        raise(Exception("Units should either be paired library or single read library."))

from snakemake.exceptions import MissingInputException


rule bowtie2_map_large:
    input:
        lambda wildcards: config["bowtie2_rules"]["units"][wildcards.unit],
        expand("indices/{{reference}}.{index}.bt2l", index=range(1,5)),
        expand("indices/{{reference}}.rev.{index}.bt2l", index=range(1,3))
    output:
        "mapping/bowtie2/{mapping_params}/{reference}/units/{unit,\w+}.bam"
    params:
        sample=lambda wildcards: UNIT_TO_SAMPLE[wildcards.unit],
        custom=lambda wildcards: config["bowtie2_rules"]["mapping_params"][wildcards.mapping_params],
        read_input_str=lambda wildcards: create_bowtie2_read_input_str(config["bowtie2_rules"]["units"][wildcards.unit]),
        ref_idx_base="indices/{reference}"
    log:
        "mapping/bowtie2/{mapping_params}/{reference}/units/{unit}.log"
    threads: 20
    shell:
        """
        bowtie2 {params.custom} \
        --rg-id '{wildcards.unit}' \
        --rg 'SM:{params.sample}' \
        -x {params.ref_idx_base} \
        -p {threads} {params.read_input_str} \
        2> {log} | samtools view -Sbh -@ {threads} - > {output}
        """
