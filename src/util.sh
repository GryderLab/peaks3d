
km(){
	echo $1 | sed s/[mM]/000k/ | sed s/[kK]/000/
}

nz(){ ## non-zero file
	if [[ -f $1 && `stat -c %s $1` -gt 0 ]];then echo 1; else echo ""; fi
}

progress-bar(){
    let _progress=(${1}*100/${2}*100)/100
    let _done=(${_progress}*4)/10
    let _left=40-$_done
    _fill=$(printf "%${_done}s")
    _empty=$(printf "%${_left}s")
printf "\rProgress : [${_fill// /#}${_empty// /-}] ${_progress}%%"
        if [ $_progress -eq  100 ];then echo "..done";fi
}


