version 1.0

import "megalodon.wdl" as megalodon

workflow diploidMegalodon {

    input {
        Array[File] inputTarballOrFast5s
        File referenceFastaPat
        File referenceFastaMat
        File readIdsPat
        File readIdsMat

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

        scatter (fast5 in fast5s) {
            call megalodon.megalodonGPU as megalodonPat {
                input:
                   inputFast5s = [fast5],
                   referenceFasta = referenceFastaPat,
                   readIdsFilename = readIdsPat,
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
            call megalodon.megalodonGPU as megalodonMat {
                input:
                   inputFast5s = [fast5],
                   referenceFasta = referenceFastaMat,
                   readIdsFilename = readIdsMat,
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

    }

    call megalodon.sum as sumPat {
        input:
            integers = flatten(select_all(megalodonPat.fileSizeGB)),
            dockerImage=dockerImage
    }

    call megalodon.sum as sumMat {
        input:
            integers = flatten(select_all(megalodonMat.fileSizeGB)),
            dockerImage=dockerImage
    }

    call megalodon.mergeMegalodon as mergePat {
        input:
            sampleIdentifier = if defined(sampleIdentifier) then sampleIdentifier + ".paternal" else "paternal",
            megalodonOutputTarballs = flatten(select_all(megalodonPat.outputTarball)),
            megalodonOutputTypes = megalodonOutputTypes,
            diskSizeGB = sumPat.value * 5, #output tar, output untar, merged files, tarred merge, slop
            zones=zones,
            dockerImage=dockerImage
    }

    call megalodon.mergeMegalodon as mergeMat {
        input:
            sampleIdentifier = if defined(sampleIdentifier) then sampleIdentifier + ".maternal" else "maternal",
            megalodonOutputTarballs = flatten(select_all(megalodonMat.outputTarball)),
            megalodonOutputTypes = megalodonOutputTypes,
            diskSizeGB = sumMat.value * 5, #output tar, output untar, merged files, tarred merge, slop
            zones=zones,
            dockerImage=dockerImage
    }

    output {
        File mergedMegalodonResultsPat = mergePat.mergedTarball
        File mergedMegalodonResultsMat = mergeMat.mergedTarball
    }
}
