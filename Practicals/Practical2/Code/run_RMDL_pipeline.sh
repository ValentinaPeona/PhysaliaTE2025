# need to have blast, mafft and t-coffee in your $PATH or change how to call them

CONS=/path/to/consensus.fasta
REF=/path/to/reference.fasta
DIR=/path/to/output/dir

cd $DIR
makeblastdb -in $REF -out $REF -dbtype nucl -parse_seqids
awk '{print$1;}' $CONS > consensi.fa.classified.clean
CONS=consensi.fa.classified.clean
sed -i 's/\//_/g' $CONS

perl RMDL_curation_pipeline.pl $REF $REF $CONS

rm -f *emp.out
mkdir final
cd aligned; for i in $(ls *.fa); do name=`ls $i | cut -f1 -d "."`; cat $i | perl -ne 'chomp;s/>\s+/>/;if(/>(\S+)/){$id{$1}++;$id2=$1;}if($id{$id2}==1){print "$_\n"}' >../final/$name.fa; done; cd ../
cd final; for i in $(ls *.fa); do name=`ls $i | cut -f1 -d "."`; t_coffee -other_pg seq_reformat -in $i -action +rm_gap 95 >$name.gaps95.fa; done; cd ../