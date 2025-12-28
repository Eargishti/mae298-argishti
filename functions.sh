echoStep(){
	echo "["`date`"]:" "$1"
}

displayVars(){
	for var in "$@"; do
		echoStep "	The value of $var = ${!var}"
	done
}
