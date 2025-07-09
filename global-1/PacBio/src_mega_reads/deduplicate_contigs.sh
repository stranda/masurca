#!/bin/bash
#this script aligns the assembly to itself and de-duplicates contigs, assumes masurca on the path

ASM_DIR=$1;
ASM_PREFIX=$2;
NUM_THREADS=$3;
PLOIDY=$4;
if [ "$5" = "" ];then 
ASM_SUBDIR="9-terminator";
else
ASM_SUBDIR=$5;
fi

if [ $PLOIDY -gt 1 ];then
MERGE_LEN=20000
HAP_SIM_RATE=90
else
MERGE_LEN=10000
HAP_SIM_RATE=94
fi

MYPATH="`dirname \"$0\"`"
MYPATH="`( cd \"$MYPATH\" && pwd )`"
export PATH=$MYPATH:$PATH;

set -e

#here we map the contigs against themselves to figure out which ones are redundant
if [ ! -e $ASM_DIR/self_map.contigs.success ];then
rm -f $ASM_DIR/filter_map.contigs.success
perl -ane 'BEGIN{$seq=""}{if($F[0] =~ /^\>/){if(length($seq)>500){print length($seq)," $rn $seq\n";}$seq="";$rn=$F[0];}else{$seq.=$F[0];}}END{print length($seq)," $rn $seq\n";}' $ASM_DIR/$ASM_SUBDIR/$ASM_PREFIX.scf.fasta | sort -S 50% -nrk1 | awk '{print $2"\n"$3}' >  $ASM_DIR/scaffolds.ref.fa && \
nucmer -t $NUM_THREADS --batch $(($(stat -c%s "$ASM_DIR/scaffolds.ref.fa")/25)) -l 31 -c 200 -b 100 -p $ASM_DIR/sasm_to_sasm  $ASM_DIR/scaffolds.ref.fa  $ASM_DIR/$ASM_SUBDIR/$ASM_PREFIX.scf.fasta && \
touch $ASM_DIR/self_map.contigs.success || exit
fi

if [ ! -e $ASM_DIR/filter_map.contigs.success ];then
awk 'BEGIN{p=1;}{if($1 ~/^>/){if(substr($1,2)==$2) p=0; else p=1;} if(p==1) print $0;}' $ASM_DIR/sasm_to_sasm.delta > $ASM_DIR/sasm_to_sasm.noself.delta && \
delta-filter -i $HAP_SIM_RATE -q $ASM_DIR/sasm_to_sasm.noself.delta > $ASM_DIR/sasm_to_sasm.noself.fdelta && \
show-coords -lcHr $ASM_DIR/sasm_to_sasm.noself.fdelta | awk '{if($12>$13) print $0}' |merge_matches_and_tile_coords_file_new.pl $MERGE_LEN | perl -ane '{$cov{$F[18]}+=$F[15];}END{foreach $k(keys %cov){print $k,"\n" if($cov{$k}>60);}}' > $ASM_DIR/sduplicates.txt && \
awk 'BEGIN{p=1;}{if($1 ~/^>/){if(substr($1,2)==$2) p=0; else p=1;} if(p==1) print $0;}' $ASM_DIR/sasm_to_sasm.delta| show-coords -lcH /dev/stdin | awk '{if($12>$13 && $10>int("'$HAP_SIM_RATE'") && $16>90) print $NF}' >> $ASM_DIR/sduplicates.txt && \
ufasta extract -v -f $ASM_DIR/sduplicates.txt $ASM_DIR/$ASM_SUBDIR/$ASM_PREFIX.scf.fasta | ufasta format > $ASM_DIR/primary.$ASM_PREFIX.scf.fasta && \
ufasta extract -f $ASM_DIR/sduplicates.txt $ASM_DIR/$ASM_SUBDIR/$ASM_PREFIX.scf.fasta | ufasta format > $ASM_DIR/alternative.$ASM_PREFIX.scf.fasta && \
touch $ASM_DIR/filter_map.contigs.success || exit
fi

