version 1.0

# Build Indexes for the ENCODE DCC RNA-Seq pipeline 

struct RunEnv {
  String docker
  String cpu
  String memory
  #string disks: select_first([disks,"local-disk 100 SSD"])
}

workflow build_idxs {
    input {
        File reference      # reference [fasta GZ]
        File spikeins       # spikeins [fasta GZ]
        File annotation     # annotation [GTF GZ]
        String genome       # genome (e.g 'GRCh38')
        String anno_version # annotation version (e.g 'v24')
        String docker = "ebelter/mgi:bulk-rna"
    }

    RunEnv runenv = {
      "docker": docker,
    }

    # Merged Annotation - anntoation plus spikeins (and tRNA)
    call build_merged_annotation { input :
        annotation=annotation,
        spikeins=spikeins,
        genome=genome,
        anno_version=anno_version,
        runenv=runenv,
    }

    # Transcripts - build fasta and kallisto index
    call build_transcripts_fasta { input: # docker needs gffread
        reference=reference,
        annotation=annotation,
    }
    call build_transcripts_idx { input:
        reference=build_transcripts_fasta.transcripts
    }

    # STAR
    call build_star_idx { input:
        reference=reference,
        spikeins=spikeins,
        annotation=annotation,
        anno_version=anno_version,
        genome=genome,
        runenv=runenv,
    }

    # RSEM
    call build_rsem_idx { input:
        reference=reference,
        spikeins=spikeins,
        annotation=annotation,
        anno_version=anno_version,
        genome=genome,
        runenv=runenv,
    }
}

task build_merged_annotation {
    input {
        File annotation
        File spikeins
        String genome
        String anno_version
        RunEnv runenv
    }

    String merged_annotation = "${genome + '_' + anno_version + '.merged.gtf.gz'}"
    command {
        python3 $(which merge_annotation.py) \
            ~{"--annotation " + annotation} \
            ~{"--spikeins " + spikeins} \
            ~{"--output_filename " + merged_annotation}
    }
    output {
        File merged_annotation = merged_annotation
    }

    runtime {
        docker: runenv.docker
        cpu : select_first([runenv.cpu,2])
        memory : "~{select_first([runenv.memory,'8'])} GB"
    }
}

task build_transcripts_fasta {
    input {
        File reference
        File annotation
        RunEnv runenv
    }

    command {
        $(which prep_transcripts) \
            ~{reference} \
            ~{annotation}
    }

    output {
        File transcripts = "transcripts.fasta"
    }

    runtime {
        docker: runenv.docker
        cpu: runenv.cpu
        memory: "~{select_first([runenv.memory,'8'])} GB"
    }
}

task build_transcripts_idx {
    input {
        File reference
        RunEnv runenv
    }

    command {
        $(which prep_kallisto.sh) ${reference}
    }

    output {
        File index = glob("*.idx")[0]
    }

    runtime {
        docker: runenv.docker
        cpu: runenv.cpu
        memory: "~{select_first([runenv.memory,'8'])} GB"
    }
}

task build_star_idx {
    input {
        File reference
        File spikeins
        File annotation
        String genome
        String anno_version
        RunEnv runenv
    }

    command {
        $(which prep_star.sh) \
            ~{reference} \
            ~{spikeins} \
            ~{annotation} \
            ~{anno_version} \
            ~{genome} \
            ~{runenv.cpu}
    }

    output {
        File index = glob("*.tgz")[0]
    }

    runtime {
        docker: runenv.docker
        cpu: runenv.cpu
        memory: "~{select_first([runenv.memory,'8'])} GB"
    }
}

task build_rsem_idx {
    input {
        File reference
        File spikeins
        File annotation
        String genome
        String anno_version
        RunEnv runenv
    }

    command {
        $(which prep_rsem.sh) \
            ~{reference} \
            ~{spikeins} \
            ~{annotation} \
            ~{anno_version} \
            ~{genome}
    }

    output {
        File index = glob("*.tgz")[0]
    }

    runtime {
        docker: runenv.docker
        cpu: runenv.cpu
        memory: "~{select_first([runenv.memory,'8'])} GB"
    }
}
