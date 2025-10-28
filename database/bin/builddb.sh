#!/bin/bash
FILE=`readlink -e $0`
bindir=`dirname $FILE`
rootdir=`dirname $bindir`
cd $rootdir

chmod +x bin/*

$bindir/download.sh
$bindir/curate.sh
