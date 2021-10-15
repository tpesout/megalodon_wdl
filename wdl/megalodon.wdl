version 1.0

workflow callMegalodon {

    input {
        Array[File] inputTarballOrFast5s
        File referenceFasta

        # megalodon configuration
        Array[String] megalodonOutputTypes = ["basecalls", "mod_mappings", "mods"]
        Array[String] modMotif = ["m", "CG", "0"]
        String guppyConfig = "dna_r9.4.1_450bps_modbases_5mc_hac_prom.cfg"
        Int guppyConcurrentReads = 0
        Int megalodonProcesses = 0
        Int guppyTimeout = 0
        String extraGuppyParams = ""

        # naming for final output
        String? sampleIdentifier

        # resources for megalodon
        Int memSizeGB = 128
        Int threadCount = 64
        Int gpuCount = 0
        String dockerImage = "tpesout/megalodon:latest"
    }

    scatter (inputFile in inputTarballOrFast5s) {
        call untar {
            input:
                fileToUntar=inputFile,
                dockerImage=dockerImage
        }

        Array[File] fast5s = if untar.didUntar then untar.untarredFast5s else [inputFile]

        if (gpuCount > 0) {
            scatter (fast5 in fast5s) {
                call megalodonGPU {
                    input:
                       inputFast5s = [fast5],
                       referenceFasta = referenceFasta,
                       megalodonOutputTypes = megalodonOutputTypes,
                       modMotif = modMotif,
                       guppyConfig = guppyConfig,
                       megalodonProcesses = megalodonProcesses,
                       guppyConcurrentReads = guppyConcurrentReads,
                       guppyTimeout = guppyTimeout,
                       extraGuppyParams = extraGuppyParams,
                       memSizeGB = memSizeGB,
                       threadCount = threadCount,
                       diskSizeGB = untar.fileSizeGB * 2 + 5,
                       gpuCount = gpuCount,
                       dockerImage=dockerImage
                }
            }
        }

        if (gpuCount <= 0) {
            scatter (fast5 in fast5s) {
                call megalodonCPU {
                    input:
                       inputFast5s = [fast5],
                       referenceFasta = referenceFasta,
                       megalodonOutputTypes = megalodonOutputTypes,
                       modMotif = modMotif,
                       guppyConfig = guppyConfig,
                       megalodonProcesses = megalodonProcesses,
                       guppyConcurrentReads = guppyConcurrentReads,
                       guppyTimeout = guppyTimeout,
                       extraGuppyParams = extraGuppyParams,
                       memSizeGB = memSizeGB,
                       threadCount = threadCount,
                       diskSizeGB = untar.fileSizeGB * 2 + 5,
                       dockerImage=dockerImage
                }
            }
        }

    }

    call sum {
        input:
            integers = if (gpuCount > 0) then flatten(select_all(megalodonGPU.fileSizeGB)) else flatten(select_all(megalodonCPU.fileSizeGB)),
            dockerImage=dockerImage
    }

    call mergeMegalodon {
        input:
            sampleIdentifier = sampleIdentifier,
            megalodonOutputTarballs = if (gpuCount > 0) then flatten(select_all(megalodonGPU.outputTarball)) else flatten(select_all(megalodonCPU.outputTarball)),
            megalodonOutputTypes = megalodonOutputTypes,
            diskSizeGB = sum.value * 5, #output tar, output untar, merged files, tarred merge, slop
            dockerImage=dockerImage
    }

    output {
        File mergedMegalodonResults = mergeMegalodon.mergedTarball
    }
}

task untar {
    input {
        File fileToUntar
        Int diskSizeGB = 128
        String dockerImage = "tpesout/megalodon:latest"
    }

    command <<<
        # Set the exit code of a pipeline to that of the rightmost command
        # to exit with a non-zero status, or zero if all commands of the pipeline exit
        set -o pipefail
        # cause a bash script to exit immediately when a command fails
        set -e
        # cause the bash shell to treat unset variables as an error and exit immediately
        set -u
        # echo each line of the script to stdout so we can see what is happening
        # to turn off echo do 'set +o xtrace'
        set -o xtrace

        # untar or don't
        mkdir tmp
        cd tmp
        if [[ "~{fileToUntar}" == *.tar ]] || [[ "~{fileToUntar}" == *.tar.gz ]] ; then
            tar xvf ~{fileToUntar}
            echo "true" >../untarred
        else
            echo "false" >../untarred
        fi
        cd ..

        # move everything to output
        mkdir output
        for FILE in `find tmp/ -name "*.fast5"` ; do
            mv $FILE output
        done

        # get output size
        if [[ `ls output | wc -l` == 0 ]] ; then
            OUTPUTSIZE=`du -s -BG ~{fileToUntar} | sed 's/G.*//'`
        else
            OUTPUTSIZE=`du -s -BG output/ | sed 's/G.*//'`
        fi
        echo $OUTPUTSIZE >outputsize
    >>>

    output {
        Boolean didUntar = read_boolean("untarred")
        Array[File] untarredFast5s = glob("output/*")
        Int fileSizeGB = read_int("outputsize")
    }

    runtime {
        memory: "2 GB"
        cpu: 1
        disks: "local-disk " + diskSizeGB + " SSD"
        docker: dockerImage
        preemptible: 1
    }

}

task megalodonGPU {
    input {
        # files
        Array[File] inputFast5s
        File referenceFasta
        File? customGuppyConfig

        # megalodon configuration
        Array[String] megalodonOutputTypes
        Array[String] modMotif = ["m", "CG", "0"]
        String guppyConfig = "dna_r9.4.1_450bps_modbases_5mc_hac_prom.cfg"
        Int guppyConcurrentReads = 0
        Int megalodonProcesses = 0
        Int guppyTimeout = 0
        String extraGuppyParams = ""

        # resources
        Int memSizeGB = 128
        Int threadCount = 64
        Int diskSizeGB = 128
        Int gpuCount = 1
        String gpuType = "nvidia-tesla-p100"
        String? nvidiaDriverVersion
        String dockerImage = "tpesout/megalodon:latest"
    }

	command <<<
        # Set the exit code of a pipeline to that of the rightmost command
        # to exit with a non-zero status, or zero if all commands of the pipeline exit
        set -o pipefail
        # cause a bash script to exit immediately when a command fails
        set -e
        # cause the bash shell to treat unset variables as an error and exit immediately
        set -u
        # echo each line of the script to stdout so we can see what is happening
        # to turn off echo do 'set +o xtrace'
        set -o xtrace

        # prepare input
        mkdir input_files
        for FILE in ~{ sep=' ' inputFast5s } ; do
            ln -s $FILE input_files/
        done

        # logging
        ls -lahR
        df -H

        # start constructing megalodon command
        cmd=(megalodon input_files/)
        cmd+=( --outputs ~{ sep=" " megalodonOutputTypes } )
        cmd+=( --reference ~{referenceFasta} )
        cmd+=( --mod-motif ~{ sep=" " modMotif } )
        cmd+=( --processes ~{ if megalodonProcesses > 0 then megalodonProcesses else threadCount} )
        cmd+=( --output-directory output/ )

        # cpu/gpu basecallers are different
        cmd+=( --guppy-server-path $GUPPY_GPU_DIR/bin/guppy_basecall_server )
        # user specified guppy config needs a directory
        ~{ if defined(customGuppyConfig)
            then ( "cmd+=( --guppy-config `basename " + customGuppyConfig + "` --guppy-params \"-d `dirname " + guppyConfig + "`\" )" )
            else ( "cmd+=( --guppy-config " + guppyConfig + ")" ) }

        # add GPU numbers
        cmd+=( --devices )
        G=0
        while [[ $G < ~{gpuCount} ]] ; do
            cmd+=( $G )
            G=$((G+1))
        done

        # GPU defaults are ok for GCR and GT
        ~{  if (guppyConcurrentReads > 0)
            then ("cmd+=( --guppy-concurrent-reads " + guppyConcurrentReads + " )" )
            else ("") }
        ~{  if (guppyTimeout > 0)
            then ("cmd+=( --guppy-timeout " + guppyTimeout + " )" )
            else ("") }

        # save extra agruments
        ~{  if extraGuppyParams != ""
            then "cmd+=( --guppy-params \""+extraGuppyParams+"\" )"
            else ""
        }

        # run megalodon command
        "${cmd[@]}"

        # save output
        UUID=`uuid`
        mkdir output_$UUID
        ls output/ | xargs -n1 -I{} mv output/{} output_$UUID/${UUID}_{}
        tar czvf megalodon_output_$UUID.tar.gz output_$UUID/

        # get output size
        du -s -BG output_$UUID/ | sed 's/G.*//' >outputsize

	>>>
	output {
		File outputTarball = glob("megalodon_output_*.tar.gz")[0]
        Int fileSizeGB = read_int("outputsize")
	}
    runtime {
        memory: memSizeGB + " GB"
        cpu: threadCount
        disks: "local-disk " + diskSizeGB + " SSD"
        gpuCount: gpuCount
        gpuType: gpuType
        nvidiaDriverVersion: nvidiaDriverVersion
        docker: dockerImage
    }
}

task megalodonCPU {
    input {
        # files
        Array[File] inputFast5s
        File referenceFasta
        File? customGuppyConfig

        # megalodon configuration
        Array[String] megalodonOutputTypes
        Array[String] modMotif = ["m", "CG", "0"]
        String guppyConfig = "dna_r9.4.1_450bps_modbases_5mc_hac_prom.cfg"
        Int guppyConcurrentReads = 0
        Int megalodonProcesses = 0
        Int guppyTimeout = 0
        String extraGuppyParams = ""

        # resources
        Int memSizeGB = 128
        Int threadCount = 64
        Int diskSizeGB = 128
        String dockerImage = "tpesout/megalodon:latest"
    }

    Int CPU_GUPPY_CONCURRENT_READS_DEFAULT = 4
    Int CPU_GUPPY_TIMEOUT_DEFAULT = 300

	command <<<
        # Set the exit code of a pipeline to that of the rightmost command
        # to exit with a non-zero status, or zero if all commands of the pipeline exit
        set -o pipefail
        # cause a bash script to exit immediately when a command fails
        set -e
        # cause the bash shell to treat unset variables as an error and exit immediately
        set -u
        # echo each line of the script to stdout so we can see what is happening
        # to turn off echo do 'set +o xtrace'
        set -o xtrace

        # prepare input
        mkdir input_files
        for FILE in ~{ sep=' ' inputFast5s } ; do
            ln -s $FILE input_files/
        done

        # logging
        ls -lahR
        df -H

        # start constructing megalodon command
        cmd=(megalodon input_files/)
        cmd+=( --outputs ~{ sep=" " megalodonOutputTypes } )
        cmd+=( --reference ~{referenceFasta} )
        cmd+=( --mod-motif ~{ sep=" " modMotif } )
        cmd+=( --processes ~{ if megalodonProcesses > 0 then megalodonProcesses else threadCount} )
        cmd+=( --output-directory output/ )

        # cpu/gpu basecallers are different
        cmd+=( --guppy-server-path $GUPPY_CPU_DIR/bin/guppy_basecall_server )
        # user specified guppy config needs a directory
        ~{ if defined(customGuppyConfig)
            then ( "cmd+=( --guppy-config `basename " + customGuppyConfig + "` --guppy-params \"-d `dirname " + guppyConfig + "`\" )" )
            else ( "cmd+=( --guppy-config " + guppyConfig + ")" ) }

        # CPU needs these defaults (unless set by the user)
        ~{  if (guppyConcurrentReads > 0)
            then ("cmd+=( --guppy-concurrent-reads " + guppyConcurrentReads + " )" )
            else ("cmd+=( --guppy-concurrent-reads " + CPU_GUPPY_CONCURRENT_READS_DEFAULT + " )" ) }
        ~{  if (guppyTimeout > 0)
            then ("cmd+=( --guppy-timeout " + guppyTimeout + " )" )
            else ("cmd+=( --guppy-timeout " + CPU_GUPPY_TIMEOUT_DEFAULT + " )" ) }

        # save extra agruments
        ~{  if extraGuppyParams != ""
            then "cmd+=( --guppy-params \""+extraGuppyParams+"\" )"
            else ""
        }

        # run megalodon command
        "${cmd[@]}"

        # save output
        UUID=`uuid`
        mkdir output_$UUID
        ls output/ | xargs -n1 -I{} mv output/{} output_$UUID/${UUID}_{}
        tar czvf megalodon_output_$UUID.tar.gz output_$UUID/

        # get output size
        du -s -BG output_$UUID/ | sed 's/G.*//' >outputsize

	>>>
	output {
		File outputTarball = glob("megalodon_output_*.tar.gz")[0]
        Int fileSizeGB = read_int("outputsize")
	}
    runtime {
        memory: memSizeGB + " GB"
        cpu: threadCount
        disks: "local-disk " + diskSizeGB + " SSD"
        docker: dockerImage
    }
}

task sum {
    input {
        Array[Int?] integers
        String dockerImage
    }
    Array[Int] allValidIntegers = select_all(integers)

    command <<<
        echo $((0 + ~{sep="+" allValidIntegers}))
    >>>

    output {
        Int value = read_int(stdout())
    }

    runtime {
        preemptible: 1
        docker: dockerImage
    }
}


task mergeMegalodon {
    input {
        Array[File?] megalodonOutputTarballs
        Array[String] megalodonOutputTypes
        String? sampleIdentifier
        Int threadCount = 8
        Int memSizeGB = 8
        Int diskSizeGB = 128
        String dockerImage = "tpesout/megalodon:latest"
    }
    Array[File] allValidMegalodonOutputTarballs = select_all(megalodonOutputTarballs)

    command <<<
        # Set the exit code of a pipeline to that of the rightmost command
        # to exit with a non-zero status, or zero if all commands of the pipeline exit
        set -o pipefail
        # cause a bash script to exit immediately when a command fails
        set -e
        # cause the bash shell to treat unset variables as an error and exit immediately
        set -u
        # echo each line of the script to stdout so we can see what is happening
        # to turn off echo do 'set +o xtrace'
        set -o xtrace

        # extract tarballs
        mkdir extracted
        cd extracted
        for FILE in ~{sep=' ' allValidMegalodonOutputTarballs} ; do
            tar xvf $FILE
        done
        cd ..

        # setup output
        mkdir output

        # handle each output type from megalodon
        for OUTPUT_TYPE in ~{sep=' ' megalodonOutputTypes} ; do

            if [[ $OUTPUT_TYPE == "basecalls" ]] ; then
                mkdir tmp_basecalls
                find extracted/ -name *.fastq | xargs -n1 -I{} mv {} tmp_basecalls/
                cat tmp_basecalls/*.fastq >output/merged_basecalls.fastq

            elif [[ $OUTPUT_TYPE == "mappings" ]] ; then
                mkdir tmp_mappings
                find extracted/ -name *mappings.bam | grep -v "mod_mappings" | xargs -n1 -I{} bash -c 'samtools sort -@~{threadCount} {} >tmp_mappings/$(basename {})'
                samtools merge -@~{threadCount} output/merged_mappings.bam tmp_mappings/*
                samtools index -@~{threadCount} output/merged_mappings.bam

                mkdir tmp_mapping_summary
                #TODO merge mapping summary

            elif [[ $OUTPUT_TYPE == "mod_mappings" ]] ; then
                mkdir tmp_mod_mappings
                find extracted/ -name *mod_mappings.bam | xargs -n1 -I{} bash -c 'samtools sort -@~{threadCount} {} >tmp_mod_mappings/$(basename {})'
                samtools merge -@~{threadCount} output/merged_mod_mappings.bam tmp_mod_mappings/*
                samtools index -@~{threadCount} output/merged_mod_mappings.bam

            elif [[ $OUTPUT_TYPE == "mods" ]] ; then
                mkdir tmp_mods
                find extracted/ -name *db | xargs -n1 -I{} bash -c 'UUID=$(uuid) ; mkdir tmp_mods/$UUID ; mv {} tmp_mods/$UUID/per_read_modified_base_calls.db'
                megalodon_extras merge modified_bases tmp_mods/* --max-processes ~{threadCount}
                mv megalodon_merge_mods_results/*.db output/merged_per_read_modified_base_calls.db

                mkdir tmp_bed_methyl/
                find extracted/ -name *bed | xargs -n1 -I{} bash -c 'cat {} | sort -k1,1V -k2,2n >tmp_bed_methyl/$(basename {})'
                megalodon_extras merge aggregated_modified_bases --sorted-inputs --output-bed-methyl-file output/merged_bed_methyl.bed tmp_bed_methyl/*

            else
                echo "Unrecognized Megalodon output type: $OUTPUT_TYPE"
            fi
        done


        # finalize output
        cd output/
        tar czvf ~{if defined(sampleIdentifier) then sampleIdentifier + "." else ""}merged_megalodon.tar.gz *
        mv *merged_megalodon.tar.gz ..
    >>>

    output {
        File mergedTarball = glob("*merged_megalodon.tar.gz")[0]
    }


    runtime {
        memory: memSizeGB + "GB"
        cpu: threadCount
        disks: "local-disk " + diskSizeGB + " SSD"
        docker: dockerImage
        preemptible: 1
    }


}

