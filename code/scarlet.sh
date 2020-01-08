

if [ "$#" -le 3 ]; then
    echo "USAGE: scarlet.sh [read count file] [cn tree file] [output prefix] [plotting style (OPTIONAL)]"
    exit 1
fi

readcount_file=$1
cntree_file=$2
output_prefix=$3
plotting_style=$4
if [ -z "$plotting_style" ];
then
    plotting_style='ALL'
fi

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"

python $DIR/scarlet.py $readcount_file $cntree_file $output_prefix

python $DIR/plot_tree.py ${output_prefix}.B_ancestor $cntree_file $plotting_style $output_prefix 

dot -Tpdf ${output_prefix}.dot -o ${output_prefix}.pdf -v
