#!/bin/bash
############################################################################
#	This script uses cdd2go.pl from https://github.com/aleimba/bac-genomics-scripts/tree/master/cdd2cog
#	to annotate COG categories after rps-blast
#	Author Arturo Vera
#	Apr 2020
#	avera@ccg.unam.mx
############################################################################

print_usage(){
	echo "Usage: $0 COG_database path_to_cdd2cog"
	echo "example: bash ~/bin/COG_differential_analysis/scripts/ccd2cog.sh ~/database ~/bin/bac-genomics-scripts/cdd2cog/"
}

if [ $# -le 0 ]
	then
	print_usage
	exit 1
fi

database=$1
cdd2go=$2
wd=$(pwd)
for i in *.out;
        do
        mkdir $i.dir
        cd $i.dir
        ln -s ../$i .
        perl $cdd2go/cdd2go.pl -r $i -c $database/cddid.tbl -f $database/fun.txt -w $database/whog
        cd $wd
done
