#!/usr/bin/bash
#################################################################################
#	This script pars the COG IDs and values of an scpecific COG category
#	Author Arturo Vera
#	April 2020
#	Dependencies:
#		-Directroies created by cdd2go.pl https://github.com/aleimba/bac-genomics-scripts/tree/master/cdd2cog#cdd2cog
#		-Perl
##################################################################################

print_usage(){
	echo "Usage: $0 COG_Category"
	echo "Example: $0 L"
}

if [ $# -le 0 ]
	then
	print_usage
	exit 1;
	fi
	

wd=$(pwd);
catego=$1;
#Obtain the COG tables for specific COG category
for i in *.dir ;
	do 
	name=$(echo $i|cut -d . -f1);
	cat $i/results/protein-id_cog.txt |\
	perl -e '($file,$catego)=@ARGV;open(FILE,$file);while(<FILE>){chomp;@col=split(/\t/);if($col[2]=~$catego){print "$_\n";}}' - $catego|\
 	cut -f 2|fgrep -f - $i/results/cog_stats.txt|\
       	perl -e '($file,$name)=@ARGV;open(FILE,$file);print "COG\tAnnot\tNumber.Of.Genes.$name\n";while(<FILE>){chomp;print "$_\n";}' - $name \
	> $wd/$name.$catego.catego.tab;
done

#Obtain the COGS annotation

for i in *$catego.catego.tab ;
	do cat $i|awk '{if (NR!=1) {print}}'|cut -f 1,2; 
done |sort|uniq > COGS.$catego.tsv

