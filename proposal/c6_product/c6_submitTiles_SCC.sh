

parameters="/usr3/graduate/mkmoon/GitHub/TAscience/c6_Parameters.json"
tileList="h12v04"

workDir=$( jq --raw-output .SCC.workDir $parameters )
dataDir=$( jq --raw-output .SCC.dataDir $parameters )

jobTime="$(date +%Y_%m_%d_%H_%M_%S)"

logDir=$( jq --raw-output .SCC.logDir $parameters ) 
mkdir -p $logDir


#Download NBAR first
if [ $(jq .setup.downloadImagery $parameters ) == true ] 
then
    while read -r tile
    do
        nameArg="-N DL_${tile}"
        logArg_download="-o ${logDir}Download_${tile}_${jobTime}.txt"
        downloadArg="-l download"
        imgStartYr=$(( $( jq .setup.imgStartYr $parameters ) )) 
        imgEndYr=$(( $( jq .setup.imgEndYr $parameters ) ))     
        qsub $nameArg $logArg_download $downloadArg runDownloadNBAR.sh $tile $baseDir $imgStartYr $imgEndYr
    done < $tileList

    while read -r tile
    do
        tileDir="${workDir}${tile}/"
        paramName="${tileDir}parameters_${jobTime}.json"
        nameArg="-N R_${tile}"
        logArg="-o ${logDir}Run_${tile}_${jobTime}.txt"
        holdArg="-hold_jid DL_${tile}"
        qsub $nameArg $logArg $nodeArgs $holdArg MSLSP_runTile_SCC.sh $tile $paramName $jobTime
    done < $tileList
else
    while read -r tile
    do
        tileDir="${workDir}${tile}/"
        paramName="${tileDir}parameters_${jobTime}.json"
        nameArg="-N R_${tile}"
        logArg="-o ${logDir}Run_${tile}_${jobTime}.txt"
        qsub $nameArg $logArg $nodeArgs MSLSP_runTile_SCC.sh $tile $paramName $jobTime
    done < $tileList
fi 

