#!/bin/bash 
BASEDIR=$(dirname $(realpath "$0"));
for f in $BASEDIR/src/*.sh;do
	. $f
done

