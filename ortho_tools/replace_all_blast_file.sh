pid_toKill=$1
file_toRM=$2
file_toUse=$3
if [ -n "$file_toUse" ]
then
	rm $file_toRM
	chmod a-w $file_toUse
	cp -p $file_toUse $file_toRM
	kill -9 $pid_toKill
else
	echo ""
	echo "bash $0   PID_toKill   File_toRM   File_toUse" 
	echo ""
	exit; 
fi

