#!/bin/bash
# Run commands in batch mode. 

## Basic information. 
args=("${@}")
myBasename=${0##*/}
path2Me=${0}

# Configuration
nProcF="Nproc"
if [ -n "$nProcF" ] && [ -f "$nProcF" ]
then
	while read proc_1time
	do
		break
	done <"$nProcF"
fi

## input : user defined. 
proc_1time=1
inLis=${args[0]} # input command list file. one command name per line. 
ss=1  # For the first line. 
ee=9999  # For the last line. 
[ -n "${args[1]}" ] && ss="${args[1]}"
[ -n "${args[2]}" ] && ee="${args[2]}"


## check if input is OK. 
if [[ -n "$inLis" && -n "$ss" && -n "$ee" && -f "$inLis" ]]
then
	echo "[Rec][$(date)] Input list=[$inLis], Start-End liens are [${ss}-${ee}]"; 
else
	echo "[Err]"
	echo "[Err]/bin/bash $myBasename input_cmd_list BeginLineNum EndLineNum"
	echo "[Err]----------------------------------------------------------"
	echo "[Err]Path to script: $path2Me"
	echo "[Err]"
	echo "[Err]inLis=$inLis,ss=$ss,ee=$ee"
	echo "[Err]Use file [${nProcF}] to control process number executed one time."
	echo "[Err]"
	exit 1
fi
## carry out iterations. 
files=()
while read line
do
	files+=("$line")
done <"$inLis"

used_pid=()    # storing the process IDs currently running. 
ii=$((ss-1))   # start index of iteration. 
if [[ "$ee" -gt "${#files[@]}" ]]
then
	ee="${#files[@]}"
	echo "[Rec][$(date)] Start-End lines restricted to [${ss}-${ee}]"
fi

while [ "$ii" -lt "$ee" ]
do
	new_used_pid=()
	for cur_pid in ${used_pid[@]}
	do
		if ps -p $cur_pid > /dev/null
		then
			# If this process ID is running, record it
			new_used_pid+=( "$cur_pid" )
		fi
	done
	used_pid=("${new_used_pid[@]}")

	if [ -n "$nProcF" ] && [ -f "$nProcF" ]
	then
		while read proc_1time
		do
			break
		done <"$nProcF"
	fi

	if [ "${#used_pid[@]}" -lt "$proc_1time" ]
	then
		## start Commands to be processed . 
		if [ -z "$files[$ii]" ] 
		then
			# If the line is null, just skip this line. 
			echo "[Wrn][$(date)]Line [$realII] is null. Skipped."
			((ii++))
			continue
		fi

		realII=$(( ii+1 ))
		echo "[Msg][$(date)]Process [$realII] command line."
		echo "[Msg][CMD]Cmd line: ${files[$ii]}"
		${files[$ii]} & 
		used_pid+=($!)
		echo "[Msg][PID]${used_pid[@]}" 
		## end   Commands to be processed . 

		((ii++))
	fi
	# Wait some time for the next throwing. 
	sleep 1
done

echo "[Rec][$(date)]All commands finished."; 




