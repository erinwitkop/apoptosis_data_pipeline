#isolate MSTRG IDs in the cluster into a smaller file that can be used by R locally 
#06_MSTRG_isolate_Sig_nonSig.sh


#Find line with Matching stringtie geneID found in the StringTie --merge output file
#Make file for all the MSTRGID_Gene_Sig.txt with all the significantly differentially expressed MSTRG transcripts

#cut -f2 MSTRGID_Gene_Sig_oshv1 > MSTRGID_Gene_Sig_oshv1.txt

args=($(cat MSTRGID_Gene_Sig_oshv1.txt))
while read -r line
do
    for i in ${args[@]}
    do
        case "$line" in
            *"$i"*) echo "$line";;
        esac
    done
done <"Bac_viral_stringtie_merged.gtf"