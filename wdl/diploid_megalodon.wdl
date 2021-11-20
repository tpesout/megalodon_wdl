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
        Int megalodonProcesses = 6
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

        scatter (fast5 in fast5s) {
            call megalodon.megalodonGPU as megalodonH1 {
                input:
                   inputFast5s = [fast5],
                   referenceFasta = referenceFastaH1,
                   readIdsFilename = readIdsH1,
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
                   nvidiaDriverVersion = nvidiaDriverVersion,
                   zones=zones,
                   maxRetries = maxRetries,
                   dockerImage=dockerImage
            }
            call megalodon.megalodonGPU as megalodonH2 {
                input:
                   inputFast5s = [fast5],
                   referenceFasta = referenceFastaH2,
                   readIdsFilename = readIdsH2,
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
                   nvidiaDriverVersion = nvidiaDriverVersion,
                   zones=zones,
                   maxRetries = maxRetries,
                   dockerImage=dockerImage
            }
        }

    }

    call megalodon.sum as sumH1 {
        input:
            integers = flatten(select_all(megalodonH1.fileSizeGB)),
            dockerImage=dockerImage
    }

    call megalodon.sum as sumH2 {
        input:
            integers = flatten(select_all(megalodonH2.fileSizeGB)),
            dockerImage=dockerImage
    }

    call megalodon.mergeMegalodon as mergeH1 {
        input:
            sampleIdentifier = sampleIdentifier,
            megalodonOutputTarballs = flatten(select_all(megalodonH1.outputTarball)),
            megalodonOutputTypes = megalodonOutputTypes,
            diskSizeGB = sumH1.value * 5, #output tar, output untar, merged files, tarred merge, slop
            zones=zones,
            dockerImage=dockerImage
    }

    call megalodon.mergeMegalodon as mergeH2 {
        input:
            sampleIdentifier = sampleIdentifier,
            megalodonOutputTarballs = flatten(select_all(megalodonH2.outputTarball)),
            megalodonOutputTypes = megalodonOutputTypes,
            diskSizeGB = sumH2.value * 5, #output tar, output untar, merged files, tarred merge, slop
            zones=zones,
            dockerImage=dockerImage
    }

    output {
        File mergedMegalodonResultsH1 = mergeH1.mergedTarball
        File mergedMegalodonResultsH2 = mergeH2.mergedTarball
    }
}
