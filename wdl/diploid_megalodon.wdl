version 1.0

import "megalodon.wdl" as megalodon

workflow diploidMegalodon {

    input {
        Array[File] inputTarballOrFast5s
        File referenceFastaH1
        File referenceFastaH2
        File readIdsH1
        File readIdsH2

        # megalodon configuration
        Array[String] megalodonOutputTypes = ["basecalls", "mod_mappings", "mods", "per_read_mods"]
        Array[String] modMotif = ["m", "CG", "0"]
        String guppyConfig = "dna_r9.4.1_450bps_modbases_5mc_hac_prom.cfg"
        Int guppyConcurrentReads = 6
        Int megalodonProcesses = 0
        Int guppyTimeout = 500
        String extraGuppyParams = ""

        # naming for final output
        String? sampleIdentifier

        # resources for megalodon
        Int memSizeGB = 64
        Int threadCount = 12
        Int gpuCount = 1
        String gpuType = "nvidia-tesla-v100"
        String nvidiaDriverVersion = "418.87.00"
        Int maxRetries = 4 # workaround for Terra failure to initilize drivers
        Array[String] zones = 	[ "us-west1-b" ]
        String dockerImage = "tpesout/megalodon:latest"
    }

    scatter (inputFile in inputTarballOrFast5s) {
        call megalodon.untar {
            input:
                fileToUntar=inputFile,
                zones=zones,
                dockerImage=dockerImage
        }

        Array[File] fast5s = if untar.didUntar then untar.untarredFast5s else [inputFile]

        call diploidMegalodonGPU {
            input:
               inputFast5s = fast5s,
               referenceFastaH1 = referenceFastaH1,
               referenceFastaH2 = referenceFastaH2,
               readIdsFileH1 = readIdsH1,
               readIdsFileH2 = readIdsH2,
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
               gpuType = gpuType,
               nvidiaDriverVersion = nvidiaDriverVersion,
               zones=zones,
               maxRetries = maxRetries,
               dockerImage=dockerImage
        }

    }

    call megalodon.sum as sizeH1 {
        input:
            integers = diploidMegalodonGPU.fileSizeGBH1,
            dockerImage=dockerImage
    }

    call megalodon.sum as sizeH2 {
        input:
            integers = diploidMegalodonGPU.fileSizeGBH2,
            dockerImage=dockerImage
    }

    call megalodon.mergeMegalodon as mergeH1 {
        input:
            sampleIdentifier = if defined(sampleIdentifier) then sampleIdentifier + ".H1" else "H1",
            megalodonOutputTarballs = diploidMegalodonGPU.outputTarballH1,
            megalodonOutputTypes = megalodonOutputTypes,
            diskSizeGB = sizeH1.value * 5, #output tar, output untar, merged files, tarred merge, slop
            zones=zones,
            dockerImage=dockerImage
    }

    call megalodon.mergeMegalodon as mergeH2 {
        input:
            sampleIdentifier = if defined(sampleIdentifier) then sampleIdentifier + ".H2" else "H2",
            megalodonOutputTarballs = diploidMegalodonGPU.outputTarballH2,
            megalodonOutputTypes = megalodonOutputTypes,
            diskSizeGB = sizeH2.value * 5, #output tar, output untar, merged files, tarred merge, slop
            zones=zones,
            dockerImage=dockerImage
    }

    output {
        File mergedMegalodonResultsH1 = mergeH1.mergedTarball
        File mergedMegalodonResultsH2 = mergeH2.mergedTarball
    }
}



task diploidMegalodonGPU {
    input {
        # files
        Array[File] inputFast5s
        File referenceFastaH1
        File referenceFastaH2
        File readIdsFileH1
        File readIdsFileH2
        File? customGuppyConfig

        # megalodon configuration
        Array[String] megalodonOutputTypes
        Array[String] modMotif = ["m", "CG", "0"]
        String guppyConfig = "dna_r9.4.1_450bps_modbases_5mc_hac_prom.cfg"
        Int guppyConcurrentReads = 4
        Int megalodonProcesses = 0
        Int guppyTimeout = 500
        String extraGuppyParams = ""

        # resources
        Int memSizeGB = 64
        Int threadCount = 12
        Int diskSizeGB = 128
        Int gpuCount = 1
        String gpuType = "nvidia-tesla-v100"
        String nvidiaDriverVersion = "418.87.00"
        Int maxRetries = 4 # workaround for Terra failure to initilize drivers
        Array[String] zones = 	[ "us-west1-b" ]
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

        # prep
        UUID=`uuid`
        mkdir input_files
        for FILE in ~{ sep=' ' inputFast5s } ; do
            ln -s $FILE input_files/
        done

        ################
        ### PATERNAL ###
        ################

        # start constructing megalodon command
        cmdH1=(megalodon input_files/)
        cmdH1+=( --outputs ~{ sep=" " megalodonOutputTypes } )
        cmdH1+=( --reference ~{referenceFastaH1} )
        cmdH1+=( --mod-motif ~{ sep=" " modMotif } )
        cmdH1+=( --processes ~{ if megalodonProcesses > 0 then megalodonProcesses else threadCount} )
        cmdH1+=( --output-directory output/ )
        cmdH1+=( --read-ids-filename ~{readIdsFileH1} )
        # cpu/gpu basecallers are different
        cmdH1+=( --guppy-server-path $GUPPY_GPU_DIR/bin/guppy_basecall_server )
        # user specified guppy config needs a directory
        ~{ if defined(customGuppyConfig)
            then ( "cmdH1+=( --guppy-config `basename " + customGuppyConfig + "` --guppy-params \"-d `dirname " + guppyConfig + "`\" )" )
            else ( "cmdH1+=( --guppy-config " + guppyConfig + ")" ) }

        # add GPU numbers
        cmdH1+=( --devices )
        G=0
        while [[ $G < ~{gpuCount} ]] ; do
            cmdH1+=( $G )
            G=$((G+1))
        done

        # GPU defaults are ok for GCR and GT
        ~{  if (guppyConcurrentReads > 0)
            then ("cmdH1+=( --guppy-concurrent-reads " + guppyConcurrentReads + " )" )
            else ("") }
        ~{  if (guppyTimeout > 0)
            then ("cmdH1+=( --guppy-timeout " + guppyTimeout + " )" )
            else ("") }

        # save extra agruments
        ~{  if extraGuppyParams != ""
            then "cmdH1+=( --guppy-params \""+extraGuppyParams+"\" )"
            else ""
        }

        # run megalodon command
        "${cmdH1[@]}"

        # save output
        mkdir output_${UUID}_H1
        ls output/ | xargs -n1 -I{} mv output/{} output_${UUID}_H1/${UUID}_H1_{}
        tar czvf megalodon_output_${UUID}_H1.tar.gz output_${UUID}_H1/

        # get output size
        du -s -BG output_${UUID}_H1/ | sed 's/G.*//' >outputsize_H1



        ################
        ### MATERNAL ###
        ################

        # start constructing megalodon command
        cmdH2=(megalodon input_files/)
        cmdH2+=( --outputs ~{ sep=" " megalodonOutputTypes } )
        cmdH2+=( --reference ~{referenceFastaH2} )
        cmdH2+=( --mod-motif ~{ sep=" " modMotif } )
        cmdH2+=( --processes ~{ if megalodonProcesses > 0 then megalodonProcesses else threadCount} )
        cmdH2+=( --output-directory output/ )
        cmdH2+=( --read-ids-filename ~{readIdsFileH2} )
        # cpu/gpu basecallers are different
        cmdH2+=( --guppy-server-path $GUPPY_GPU_DIR/bin/guppy_basecall_server )
        # user specified guppy config needs a directory
        ~{ if defined(customGuppyConfig)
            then ( "cmdH2+=( --guppy-config `basename " + customGuppyConfig + "` --guppy-params \"-d `dirname " + guppyConfig + "`\" )" )
            else ( "cmdH2+=( --guppy-config " + guppyConfig + ")" ) }

        # add GPU numbers
        cmdH2+=( --devices )
        G=0
        while [[ $G < ~{gpuCount} ]] ; do
            cmdH2+=( $G )
            G=$((G+1))
        done

        # GPU defaults are ok for GCR and GT
        ~{  if (guppyConcurrentReads > 0)
            then ("cmdH2+=( --guppy-concurrent-reads " + guppyConcurrentReads + " )" )
            else ("") }
        ~{  if (guppyTimeout > 0)
            then ("cmdH2+=( --guppy-timeout " + guppyTimeout + " )" )
            else ("") }

        # save extra agruments
        ~{  if extraGuppyParams != ""
            then "cmdH2+=( --guppy-params \""+extraGuppyParams+"\" )"
            else ""
        }

        # run megalodon command
        "${cmdH2[@]}"

        # save output
        mkdir output_${UUID}_H2
        ls output/ | xargs -n1 -I{} mv output/{} output_${UUID}_H2/${UUID}_H2_{}
        tar czvf megalodon_output_${UUID}_H2.tar.gz output_${UUID}_H2/

        # get output size
        du -s -BG output_${UUID}_H2/ | sed 's/G.*//' >outputsize_H2

	>>>
	output {
		File outputTarballH1 = glob("megalodon_output_*_H1.tar.gz")[0]
		File outputTarballH2 = glob("megalodon_output_*_H2.tar.gz")[0]
        Int fileSizeGBH1 = read_int("outputsize_H1")
        Int fileSizeGBH2 = read_int("outputsize_H2")
	}
    runtime {
        memory: memSizeGB + " GB"
        cpu: threadCount
        disks: "local-disk " + diskSizeGB + " SSD"
        gpuCount: gpuCount
        gpuType: gpuType
        maxRetries: maxRetries
        nvidiaDriverVersion: nvidiaDriverVersion
        docker: dockerImage
        zones: zones
    }
}
