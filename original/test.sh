#./src/lpboost_correct -m 1 -x 8 -n 0.5 data/cpdbtrain
#./eval/evaluator model data/cpdbtest

./src/lpboost -m 1 -x 6 -n 0.5 ~/desktop/experiment/data/NCI47_comp_buckets/train0.gsp
./eval/evaluator model ~/desktop/experiment/data/NCI47_comp_buckets/test0.gsp

