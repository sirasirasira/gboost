#./src/lpboost_correct -m 1 -x 8 -n 0.5 data/cpdbtrain
#./eval/evaluator model data/cpdbtest

./src/lpboost_correct -x 3 -m 1 -n 0.4 ~/desktop/experiment/data/cpdb_buckets/train0.gsp >a
./eval/evaluator model ~/desktop/experiment/data/cpdb_buckets/test0.gsp >a
#echo "acc: `python ~/desktop/experiment/bin/grid_search/auc_eval.py temp.txt 0`"

